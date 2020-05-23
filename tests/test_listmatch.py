import random

import pytest
import numpy as np

import matcha
import _matcha

def test_basic_matching():
    barcode_sequences = ["ATGC", "TGAC", "ACAA", "CGAT"]
    query_sequences   = ["ATGC", "TCAC", "ACAA", "CAAG"]
    dists             = [0, 1, 0, 2]
    labels = ["one", "two", "three", "four"]    

    m = matcha.ListMatcher(barcode_sequences, labels)
    results = m.match_all(query_sequences, 0)

    assert list(m.labels[results.match]) == labels
    assert list(results.dist) == dists
