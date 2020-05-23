import random
import timeit

import numpy as np

import matcha
import _matcha

from utils import hamming_dist, random_sequence, random_mismatches, results_equal


def construct_benchmark_seqs(sequence_len = 8, barcode_count = 10, sequence_count = 1000, seed="benchmarkme"):
    random.seed(seed)
    barcode_sequences =  [random_sequence(sequence_len, "ATGC") for i in range(barcode_count)]

    mismatch_counts = [random.randint(0, sequence_len) for i in range(sequence_count)]
    mismatch_against = random.choices(barcode_sequences, k=sequence_count)
    sequences = [random_mismatches(b, m) for b, m in zip(mismatch_against, mismatch_counts)]

    binary_sequences = _matcha.stringsToBinary(sequences, 0, sequence_len)
    return barcode_sequences, binary_sequences, sequences

def binary_results_equal(ref_results, binary_results, max_mismatches, sequence_len):
    best_match, raw_quality = binary_results[0], binary_results[1]
    other_result = matcha.MatchResult(best_match, raw_quality & 63, np.right_shift(raw_quality, 6), [])
    return results_equal(ref_results, other_result, max_mismatches, sequence_len)

def benchmark(max_mismatches, subsequence_counts, **kwargs):
    """Benchmark a list of """
    sequence_len = kwargs["sequence_len"]
    barcode_seqs, binary_seqs, seqs = construct_benchmark_seqs(**kwargs)
    list_matcher = matcha.ListMatcher(barcode_seqs)

    ref_results = list_matcher.match_all(seqs)
    
    print("### Running timing experiments ###")
    print("ListMatcher: ", end="")
    time_results = np.zeros_like(binary_seqs)
    t = timeit.timeit(stmt = "list_matcher._matcher.match_raw(binary_seqs, time_results)", number=1, globals=locals())
    if not binary_results_equal(ref_results, time_results, sequence_len, sequence_len):
        print("(failed correctness check)")
    else:
        print(f"{t:.3f}")
    
    for mismatch, subseq in zip(max_mismatches, subsequence_counts):
        print(f"HashMatcher(max_mismatches = {mismatch}, subsequences = {subseq}) ", end="")
        m = matcha.HashMatcher(barcode_seqs, mismatch, subseq)
        time_results = np.zeros_like(binary_seqs)
        t = timeit.timeit(stmt = "m._matcher.match_raw(binary_seqs, time_results)", number=1, globals=locals())
        if not binary_results_equal(ref_results, time_results, mismatch, sequence_len):
            print("(failed correctness check)")
        else:
            print(f"{t:.3f}")

if __name__ == "__main__":
    print("Benchmarking sequence_len 8, 72 barcodes, 1M sequences")
    benchmark([0,1,1,2,2,3,3],[1,1,2,2,3,2,3], sequence_len=8, barcode_count=1000, sequence_count=1000000)
    print("Benchmarking sequence_len 8, 1000 barcodes, 1M sequences")
    benchmark([0,1,2],[1,2,2], sequence_len=8, barcode_count=1000, sequence_count=1000000)
    print("Benchmarking sequence_len 16, 1000 barcodes, 1M sequences")
    benchmark([0,1,1,2,2,3,3],[1,1,2,2,3,2,3], sequence_len=8, barcode_count=1000, sequence_count=1000000)
