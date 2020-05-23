import random

import pytest
import numpy as np

import matcha
import _matcha

from .utils import hamming_dist, random_sequence, random_mismatches


def get_neighbors(seqs, distance):
    if type(seqs) is str:
        seqs = set([seqs])
    if distance == 0:
        return seqs

    ret = set()
    for s in seqs:
        for i in range(len(s)):
            for c in ["A", "T", "G", "C"]:
                ret.add(s[:i] + c + s[i+1:])
    ret.update(seqs)
    return get_neighbors(ret, distance - 1)

def test_neighbors():
    sequence = "ATGCAGTAC"
    distance = 3
    neighbors1 = get_neighbors(sequence, distance)
    neighbor_masks = matcha.HashMatcher.get_mismatch_masks(list(range(len(sequence))), distance)
    assert len(neighbor_masks) == len(neighbors1)
    
    neighbors2 = set()
    binary_seq, flag = _matcha.stringToBinary(sequence)
    for n in neighbor_masks:
        neighbors2.add(_matcha.binaryToString(binary_seq ^ n, len(sequence), flag))
    
    assert neighbors1 == neighbors2

def test_neighbors2():
    sequence = "ATGCAGTAC"
    distance = 1
    indexes = [0, 4]
    correct = [
        "ATGCAGTAC",
        "TTGCAGTAC", "GTGCAGTAC", "CTGCAGTAC",
        "ATGCTGTAC", "ATGCGGTAC", "ATGCCGTAC",
    ]

    binary_seq, flag = _matcha.stringToBinary(sequence)
    neighbor_masks = matcha.HashMatcher.get_mismatch_masks(indexes, distance)
    neighbors = [_matcha.binaryToString(binary_seq ^ n, len(sequence), flag) for n in neighbor_masks]
    assert set(neighbors) == set(correct)


def test_basic_matching():
    barcode_sequences = ["ATGC", "TGAC", "ACAA", "CGAT"]
    query_sequences   = ["ATGC", "TCAC", "ACAA", "CAAG"]
    dists             = [0, 1, 0, 2]
    labels = ["one", "two", "three", "four"]    

    m = matcha.HashMatcher(barcode_sequences, 2, 3, labels)
    results = m.match_all(query_sequences, 0)

    assert list(m.labels[results.match]) == labels
    assert list(results.dist) == dists


def assert_match_results_equal(reference, r, max_mismatches, sequence_len):
    within_best = reference.dist <= max_mismatches
    within_second_best = reference.second_best_dist <= max_mismatches

    # r must either have the correct value, or for cases where the distance is over max_mismatches it can have the given default
    assert all((r.match == reference.match) | ((r.match == (2**64-1)) & ~within_best))
    assert all((r.dist == reference.dist) | ((r.dist == (2**6-1)) & ~within_best))
    assert all((r.second_best_dist == reference.second_best_dist) | ((r.second_best_dist == (2**6-1)) & ~within_second_best))

def test_full_matching():
    # 1. Generate barcode sequences (at least 3)
    # 2. Generate list to match against (Should contain a range of mismatch counts)
    # 3. Create hashmatchers with different max mismatches and subsequence counts (up to 3 subsequences)
    sequence_len = 8
    barcode_count = 10
    sequence_count = 100

    barcode_sequences =  [random_sequence(sequence_len, "ATGC") for i in range(barcode_count)]

    mismatch_counts = [random.randint(0, sequence_len) for i in range(sequence_count)]
    mismatch_against = random.choices(barcode_sequences, k=sequence_count)
    sequences = [random_mismatches(b, m) for b, m in zip(mismatch_against, mismatch_counts)]

    ref_matcher = matcha.ListMatcher(barcode_sequences)
    ref_results = ref_matcher.match_all(sequences)

    hash_matchers = [
        matcha.HashMatcher(barcode_sequences, mismatch, subseqs) 
        for mismatch in range(sequence_len)
        for subseqs in range(1,4)
    ]

    for max_mismatches in range(sequence_len):
        for subseqs in range(1,4):
            r = matcha.HashMatcher(barcode_sequences, max_mismatches, subseqs).match_all(sequences)
            assert_match_results_equal(ref_results, r, max_mismatches, sequence_len)

def test_all_different():
    sequence_len = 10
    barcode_sequences = [c * sequence_len for c in "ATG"]

    for subseqs in range(1,4):
        r = matcha.HashMatcher(barcode_sequences, sequence_len, subseqs).match_all(["C" * sequence_len])
        assert list(r.match) == [0]
        assert list(r.dist) == [sequence_len]
        assert list(r.second_best_dist) == [sequence_len]