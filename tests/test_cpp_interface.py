import random

import numpy as np
import pytest

import _matcha

from .utils import hamming_dist, random_sequence


def test_matcher_exception_N():
    m = _matcha.ListMatcher()
    with pytest.raises(RuntimeError, match="Sequence .* has N's"):
        m.add_sequences(["ANA"])

def test_matcher_exception_length():
    m = _matcha.ListMatcher()
    m.add_sequences(["A"])
    with pytest.raises(RuntimeError, match="Sequence .* does not match size"):
        m.add_sequences(["AA"])

def test_binary_conversion():
    sequence_length = 10
    n_sequences = 30
    sequences = [random_sequence(sequence_length) for i in range(n_sequences)]
    
    binary = [_matcha.stringToBinary(seq) for seq in sequences]

    ret_sequences = [_matcha.binaryToString(seq[0], sequence_length, seq[1]) for seq in binary]

    assert sequences == ret_sequences

def test_binary_conversion_vectorized():
    sequence_length = 10
    n_sequences = 30
    start = 1
    end = 6
    sequences = [random_sequence(sequence_length) for i in range(n_sequences)]

    no_vector = [_matcha.stringToBinary(seq[start:end]) for seq in sequences]
    vector = _matcha.stringsToBinary(sequences, start, end)

    assert np.all(np.array(no_vector).T == vector)

def mismatch_compare(sequence_length):
    # Test matching against 3 sequences, make sure match, top dist, and 2nd best dist are right
    seqs = [random_sequence(sequence_length, "ATGC") for i in range(3)]
            
    query = random_sequence(sequence_length)

    m = _matcha.ListMatcher()
    m.add_sequences(seqs)
    res = m.match_all([query], 0, sequence_length)
    matcha_best_match = seqs[res[0,0]]

    distances = [hamming_dist(seq, query) for seq in seqs]
    best_match = seqs[np.argmin(distances)]

    assert matcha_best_match in seqs
    assert hamming_dist(matcha_best_match, query) == min(distances)
    assert res[1] & 63 == min(distances)
    # Check that 2nd-closest distance is right also
    assert res[1] >> 6 == np.sort(distances)[-2]


def test_mismatches_random_short():
    random.seed("happyseed")
    for i in range(3000):
        mismatch_compare(3)

def test_mismatches_random_long():
    random.seed("happyseed2")
    for i in range(3000):
        mismatch_compare(10)
