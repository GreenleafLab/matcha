import collections
import gzip
from pathlib import Path
import re
from string import Formatter
import subprocess

import numpy as np
import pandas as pd

class FastqReader:
    """
    Fastq reading, barcode matching, and optional export of barcoded fastqs

    Attributes:
        matches (Dict[str, MatchQuality]): After calling read_chunk, holds the quality of matches for each barcode_name.

    """
    MatcherConfig = collections.namedtuple("MatcherConfig", ["sequence_name", "barcode_name", "matcher", "match_start"])
    

    def __init__(self):
        self._started_reading = False
        self._inputs = {}
        self._outputs = {}
        self._barcodes = []
        self._output_pattern = None
        self.matches = {}
        self._sequences = {}
        self._parse_read_attributes = False
    
    def add_sequence(self, sequence_name, input_path, output_path=None):
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
        
        self._inputs[sequence_name] = read_gzip(input_path)
        if output_path:
            self._outputs[sequence_name] = write_gzip(output_path)

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

        if barcode_name in self._parsed_attributes + ["read_name"]:
            raise ValueError("Can't add barcode with reserved name lane, tile, x, or y")

        config = self.MatcherConfig(sequence_name, barcode_name, matcher, match_start)
        self._barcodes.append(config)
    
    _parsed_attributes = ["lane", "tile", "x", "y"]

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

        self._output_pattern = pattern
        self._pattern_literals = literals
        self._pattern_fields = fields
        
        # Check if any attributes that require special parsing are in the pattern
        if any(attr in fields for attr in self._pattern_fields):
            self._parse_read_attributes = True
        
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
        if self._output_pattern:
            barcode_names = [b.barcode_name for b in self._barcodes]
            for var in self._pattern_fields:
                assert var in barcode_names or var in self._parsed_attributes or var == "read_name"

    def read_chunk(self, max_chunk_size):
        """
        Read and match barcodes on a chunk of data 
        
        Args:
            max_chunk_size (int): Maximum number of reads to process in the chunk

        Returns:
            Number of records read in chunk
        """
        self._validate_config()

        sequences = {
            name: read_fastq(f, max_reads=max_chunk_size) for name, f in self._inputs.items()
        }
        
        num_reads = [len(s.seq) for s in sequences.values()]
        assert num_reads.count(num_reads[0]) == len(num_reads) # Make sure the read count was equal from all files

        if num_reads[0] == 0:
            return 0

        for b in self._barcodes:
            self.matches[b.barcode_name] = b.matcher.match_all(sequences[b.sequence_name].seq, b.match_start)

        self._sequences = sequences
        return num_reads[0]
        


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
        if self._output_pattern:
            match_names = {}
            for b in self._barcodes:
                if b.barcode_name not in self._pattern_fields:
                    continue
                match_names[b.barcode_name] = b.matcher.labels[self.matches[b.barcode_name].match]
                            
        lane, tile, x, y = None, None, None, None

        substitution_vars = {key: pd.Series(val)[filter] for key, val in match_names.items()}

        for name, f in self._outputs.items():
            # Filter the reads
            seq = pd.Series(self._sequences[name].seq)[filter].str.decode("ascii")
            qual = pd.Series(self._sequences[name].qual)[filter].str.decode("ascii")
            name = pd.Series(self._sequences[name].name)[filter].str.decode("ascii")
            # Create the names
            if self._output_pattern:
                substitution_vars["read_name"] = name
                if self._parse_read_attributes and lane is None:
                    parse_name_frame = name.str.split(":", expand=True)
                    substitution_vars["lane"] = parse_name_frame[3]
                    substitution_vars["tile"] = parse_name_frame[4]
                    substitution_vars["x"]    = parse_name_frame[5]
                    substitution_vars["y"]    = parse_name_frame[6]
                name = self._pattern_literals[0]
                for literal, field in zip(self._pattern_literals[1:], self._pattern_fields):
                    name += substitution_vars[field] + literal
            
            # Make the output 
            output = "@" + name + "\n" + seq + "\n+\n" + qual + "\n"
            f.write(output.str.cat().encode())
    
    def close(self):
        """Close output files"""
        for f in self._outputs.values():
            f.close()

FastqData = collections.namedtuple("FastqData", ["name", "seq", "qual"])
def read_fastq(file, max_reads = -1):
    lines = [file.readline() for line in range(max_reads*4)]
    # Handle hitting the end of the input
    if not lines[-1]:
        if not lines[0]:
            return FastqData([], [], [])
        last_read = -1
        while not lines[last_read]:
            last_read -= 1
        lines = lines[:last_read+1]
        assert len(lines) % 4 == 0
    
    name = [l.strip()[1:] for l in lines[::4]]
    seq = [l.strip() for l in lines[1::4]]
    qual = [l.strip() for l in lines[3::4]]
    return FastqData(name, seq, qual)
    
def read_gzip(filename):
    if filename.endswith(".gz"):
        # Use subprocess so that the decompression can happen in a parallel
        ps = subprocess.Popen(('gzip', '-d', '-c', filename), stdout=subprocess.PIPE)
        return ps.stdout
    else:
        return open(filename, 'rb')

def write_gzip(filename):
    if filename.endswith(".gz"):
        ps = subprocess.Popen(('gzip', '-c'), stdout=open(filename, 'wb'), stdin=subprocess.PIPE)
        return ps.stdin
    else:
        return open(filename, 'wb')