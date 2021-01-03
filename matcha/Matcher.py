import collections
import itertools
import re
import subprocess

import numpy as np

import _matcha

class Matcher:
    """
    Barcode matcher python wrapper

    Args:
        sequences (List[str]): Barcode DNA sequences
        labels (List[str]): Labels for barcode sequences (optional)
        _matcher (_matcha.Matcher): C++ matcher object to use
    """
    def __init__(self, sequences, _matcher, labels=None):     
        self.sequences = np.array(sequences)
        self.binary_sequences = np.array([_matcha.stringToBinary(seq) for seq in self.sequences])
        if labels is None:
            labels = self.sequences
        self.labels = np.array(labels)
        self.sequence_length = len(self.sequences[0])

        self._matcher = _matcher
        self._matcher.add_sequences(self.sequences)
        self._matcher.add_labels(labels)

    def match_all(self, sequences, start=0):
        """Match all sequences in a list
        
        Args:
            sequences (List[str]): Query sequences
            start (int): 0-based position of first base to use in barcode match
        
        Returns:
            collections.namedtuple: Named tuple. Field seq is binary encoding of best matching sequence.
                Other fields give match quality and may vary based on algorithm.
                Fields by algorithm:
                    list: dist, second_best_dist -- hamming distance to best and second best matches respectively
        """
        end = start + self.sequence_length
        result = self._matcher.match_all(sequences, start, end)
        return self.process_matches(result)

    def process_matches(self, match_result):
        """Process a match result quality based on the type of algorithm used"""
        best_match, raw_quality = match_result[0], match_result[1]
        return MatchResult(best_match, raw_quality & 63, np.right_shift(raw_quality, 6), self.labels)
   
class ListMatcher(Matcher):
    """
    List matcher iterates through possible matches in a list. Slow for many valid barcodes,
    but has no limits on the maximum number of mismatches. As a rule of thumb, this is the best choice if you
    have under 20 valid sequences, but you should probably not use it to match over 100 sequences

    Args:
        sequences (List[str]): Barcode DNA sequences
        labels (List[str]): Labels for barcode sequences (optional)
    """
    def __init__(self, sequences, labels=None):   
        _matcher = _matcha.ListMatcher()
        super().__init__(sequences, _matcher, labels)  

class MatchResult:
    """
    Container type for match results.

    Attributes:
        label (numpy.ndarray): String array of labels for the best match. Only valid if labels were provided while creating the matcher object
        dist (numpy.ndarray): Integer array of distances (number of mismatches) to the best match
        second_best_dist (numpy.ndarray): Integer array of distances (number of mismatches) to the second-best match
        match (numpy.ndarray): Integer array of indexes of the best match in the Matcher object's list of valid sequences
    """
    def __init__(self, match, dist, second_best_dist, labels):
        self.match = match
        self.dist = dist
        self.second_best_dist = second_best_dist
        self._labels = labels
    
    @property
    def label(self):
        return self._labels[self.match]

class HashMatcher(Matcher):
    """
    Hash matcher uses hash tables of subsequences for barcode search, using the algorithm of Norouzi et al. https://arxiv.org/pdf/1307.2982.pdf. 
    Fast for many valid barcodes and few possible mismatches. Increasing the max_mismatches will decrease performance. For
    10x-style barcodes (16bp, ~1M valid barcodes), recommended settings are max_mismatches=1, subsequence_count=2

    Args:
        sequences (List[str]): Barcode DNA sequences
        max_mismatches (int): Maximum mismatches to match against
        subsequence_count (int): Number of subsequence indexes to use. 
            In general, use lower subsequence_count for larger number of valid labels,
            and higher subsequence_count for searching against a larger number of mismatches.
        labels (List[str]): Labels for barcode sequences (optional)

    """
    def __init__(self, sequences, max_mismatches, subsequence_count, labels=None):   
        # Get masks for extracting subsequences
        # Stripe base indexes for each subsequence, in case there are contiguous chunks of similarity in valid barcodes
        valid_bases = list(range(len(sequences[0])))
        subsequence_indexes = [valid_bases[i::subsequence_count] for i in range(subsequence_count)]
        subsequence_indexes.sort(key = lambda l: len(l))
        subsequence_masks = [self.get_mask(indexes) for indexes in subsequence_indexes]

        # Calculate masks for finding neigbors during lookup
        # Terminology as from equation 3 in: https://arxiv.org/pdf/1307.2982.pdf 
        mismatch_masks = []
        r_prime = max_mismatches // subsequence_count 
        a = max_mismatches % subsequence_count
        for i in range(subsequence_count):
            mismatch_range = r_prime if i <= a else r_prime - 1
            mismatch_masks.append(self.get_mismatch_masks(subsequence_indexes[i], mismatch_range))

        _matcher = _matcha.HashMatcher(subsequence_masks, mismatch_masks, max_mismatches)
        super().__init__(sequences, _matcher, labels)  

    @staticmethod
    def get_mask(indexes):
        mask = 0
        for i in indexes:
            mask |= 3 << (i*2)
        return mask

    @staticmethod
    def get_mismatch_masks(indexes, mismatch_range):
        """Return a (xor) mask for every possible sequence neighbor within a given mismatch range"""
        mismatch_masks = []

        mismatch_positions = itertools.chain(
            *(itertools.combinations(indexes, m) for m in range(mismatch_range+1))
        )
        for pos in mismatch_positions:
            for vals in itertools.product([1,2,3], repeat=len(pos)):
                mask = 0
                for p, v in zip(pos, vals):
                    mask |= v << (p*2)
                mismatch_masks.append(mask)
        return mismatch_masks