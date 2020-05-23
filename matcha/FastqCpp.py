import collections
import concurrent.futures
import itertools
from pathlib import Path
from string import Formatter

import numpy as np
import pandas as pd

import _matcha

class FastqCpp:
    """
    Fastq reading, barcode matching, and optional export of barcoded fastqs

    Attributes:
        matches (Dict[str, MatchQuality]): After calling read_chunk, holds the quality of matches for each barcode_name.

    """
    MatcherConfig = collections.namedtuple("MatcherConfig", ["sequence_name", "barcode_name", "matcher", "match_start"])
    

    def __init__(self, threads=None):
        self._started_reading = False
        self._inputs = {} # sequence_name -> input path
        self._outputs = {} # sequence_name -> output path
        self._barcodes = [] # list of barcode configs
        self._fastq_files = {} # sequence_name -> c++ FastqFile object
        
        # Note: Output names are calculated by alternating a literal from _pattern_literals with a 
        # looked up barcode match / read name field from _pattern_fields (see self.set_output_names)
        self._name_literals = [] # List of strings with any literal text to be included in output names
        self._name_field_indexes = [] # List of field indices in output pattern. -1 == read_name, -i = (i-2)th field when splitting the read name by ":"
        self._barcode_name_to_index = {} # Dict from barcode_name to field_index
        self.set_output_names("{read_name}")

        # Dictionary of barcode_name -> match results for most recent chunk
        self.matches = {}

        if threads:
            self._executor = concurrent.futures.ThreadPoolExecutor(threads)
            self._map = lambda *args: list(self._executor.map(*args))
        else:
            self._map = lambda *args: list(map(*args))
        
    
    def add_sequence(self, sequence_name, input_path, output_path=""):
        """
        Add fastq input and/or output files for barcode matching 

        Args:
            sequence_name (str): Name of sequence (typically R1, R2, I1, or I2)
            input_path (str): Path of fastq file for sequence
            output_path (str): Path of output fastq containing matching reads (optional)
        """
        if self._started_reading:
            raise Exception("Can't modify FastqReader settings after calling read_chunk")

        if isinstance(input_path, Path):
            input_path = str(input_path)
        if isinstance(output_path, Path):
            output_path = str(output_path)
        
        self._inputs[sequence_name] = input_path
        self._outputs[sequence_name] = output_path
        

    def add_barcode(self, barcode_name, matcher, sequence_name, match_start=0):
        """
        Add barcode matcher on a sequence

        Args:
            barcode_name (str): Name of barcode (used to access match data from this barcode)
            matcher (Matcher): matcha.Matcher object holding the valid barcodes
            sequence_name (str): Name of sequence to match on (typically R1, R2, I1, or I2)
            match_start (int): 0-based index to start matching from in the sequence
        """
        if self._started_reading:
            raise Exception("Can't modify FastqReader settings after calling read_chunk")

        if barcode_name in self._parsed_attributes or barcode_name == "read_name":
            raise ValueError("Can't add barcode with reserved name read_name, lane, tile, x, or y")

        config = self.MatcherConfig(sequence_name, barcode_name, matcher, match_start)
        self._barcodes.append(config)
    
    _parsed_attributes = {"lane": 3, "tile": 4, "x": 5, "y": 6} #0-based indices of attributes in bcl2fastq2 name when split by ':'

    def set_output_names(self, pattern):
        """
        Set pattern for read names in the output fastq.

        Args:
            pattern (str): Python format string. Available variables for substition are: 
                - read_name: The original read name from the input fastq
                - Any barcodes that have been added via inputs.add_barcodes: outputs the label for best match barcode
                - lane, tile, x, y: derived from Illumina read name position info
        """
        if self._started_reading:
            raise Exception("Can't modify FastqReader settings after calling read_chunk")

        literals = []
        fields = []
        
        carryover_string = ""
        for literal, field, _, _ in Formatter().parse(pattern):
            if field is not None:
                literals.append(carryover_string + literal)
                fields.append(field)
                carryover_string = ""
            else:
                caryover_string += literal
        literals.append(carryover_string)

        self._name_literals = literals # List of strings with any literal text to be included in output names

        self._name_field_indexes = [] # List of field indices in output pattern. -1 == read_name, -i = (i-2)th field when splitting the read name by ":"
        self._barcode_name_to_index = {} # Dict from barcode_name to field_index

        for f in fields:
            if f == "read_name":
                self._name_field_indexes.append(-1)
            elif f in self._parsed_attributes:
                self._name_field_indexes.append(-self._parsed_attributes[f] - 2)
            else:
                if f not in self._barcode_name_to_index:
                    new_idx = max(self._barcode_name_to_index.values(), default = -1) + 1
                    self._barcode_name_to_index[f] = new_idx
                self._name_field_indexes.append(self._barcode_name_to_index[f])
        
    def _validate_config(self):
        """
        Validate that the configuration is consistent: 
            - Barcodes rely only on valid sequence names
            - Output pattern relies only on valid barcode names

        Call when starting reading (later calls are no-ops)
        """
        # Only validate if we have not yet started reading
        if self._started_reading:
            return
    
        # Matchers rely only on valid sequence names
        for m in self._barcodes:
            assert m.sequence_name in self._inputs

        # Output pattern relies only on valid barcode names
        barcode_names = [b.barcode_name for b in self._barcodes]
        for var in self._barcode_name_to_index:
            assert var in barcode_names or var in self._parsed_attributes or var == "read_name"

        self._init_cpp_objects()

        self._started_reading = True

    def _init_cpp_objects(self):
        for read_name in self._inputs:
            fastq_file = _matcha.FastqFile(
                self._inputs[read_name], 
                self._name_literals, 
                self._name_field_indexes, 
                self._outputs[read_name]
            )
            self._fastq_files[read_name] = fastq_file
    
    def _read_fastq(self, sequence_name, max_chunk_size):
        fastq_file = self._fastq_files[sequence_name]
        return fastq_file.read_chunk(max_chunk_size)

    def _match_barcode(self, barcode):
        b = barcode
        fastq_file = self._fastq_files[b.sequence_name]
        match_results = fastq_file.match(b.matcher._matcher, b.match_start, b.match_start + b.matcher.sequence_length)
        self.matches[b.barcode_name] =  b.matcher.process_matches(match_results)

    def read_chunk(self, max_chunk_size):
        """
        Read and match barcodes on a chunk of data 
        
        Args:
            max_chunk_size (int): Maximum number of reads to process in the chunk

        Returns:
            Number of records read in chunk
        """
        self._validate_config()

        all_records_read = self._map(
            self._read_fastq, 
            self._fastq_files, 
            itertools.repeat(max_chunk_size))
        all_records_read = list(all_records_read)

        records_read = all_records_read[0]
        if any(r != records_read for r in all_records_read):
            raise Exception("Unequal number of records read from input fastq files") 
        
        if records_read == 0:
            return 0

        self._map(self._match_barcode, self._barcodes)
        
        return records_read

    def _write_fastq(self, read_name, filter, match_indexes, matchers):
        if not self._outputs[read_name]:
            return # Don't write output unless requested
        self._fastq_files[read_name].write_chunk(filter, match_indexes, matchers)
        
    def write_chunk(self, filter):
        """
        Write the current chunk of data in fastq format to output sequence files.

        Args:
            filter (array_like[bool]): Boolean array of same length as last chunk_size. 
                For each element set to True, the corresponding read from the last-read chunk will be output.
        """
        # If all reads are masked, return instantly
        if not np.any(filter):
            return

        # Lookup barcode names if needed
        field_count = max(self._barcode_name_to_index.values()) + 1
        match_indexes = [None] * field_count
        matchers = [None] * field_count

        for b in self._barcodes:
            if b.barcode_name not in self._barcode_name_to_index:
                continue
            index = self._barcode_name_to_index[b.barcode_name]
            matchers[index] = b.matcher._matcher
            match_indexes[index] = self.matches[b.barcode_name].match
        
        self._map(
            self._write_fastq,
            self._fastq_files,
            itertools.repeat(filter),
            itertools.repeat(match_indexes),
            itertools.repeat(matchers)
        )
    
    def close(self):
        """Close output files"""
        for f in self._fastq_files.values():
            f.close()
