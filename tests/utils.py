import random

def hamming_dist(a, b):
    return sum(c1 != c2 for c1, c2 in zip(a,b))

def random_sequence(k, alphabet="ATGCN"):
    return "".join(random.choices(alphabet, k=k))

def random_mismatches(reference, mismatches, alphabet="ATGCN"):
    mismatch_indices = sorted(random.sample(range(len(reference)), mismatches))
    mismatch_values = [random.choice(alphabet) for i in range(mismatches)]
    a = list(reference)
    for i, val in zip(mismatch_indices, mismatch_values):
        a[i] = val
    return "".join(a)

def results_equal(reference, r, max_mismatches, sequence_len):
    within_best = reference.dist <= max_mismatches
    within_second_best = reference.second_best_dist <= max_mismatches

    # r must either have the correct value, or for cases where the distance is over max_mismatches it can have the given default
    if not all((r.match == reference.match) | ((r.match == (2**64-1)) & ~within_best)):
        return False
    if not all((r.dist == reference.dist) | ((r.dist == (2**6-1)) & ~within_best)):
        return False
    if not all((r.second_best_dist == reference.second_best_dist) | ((r.second_best_dist == (2**6-1)) & ~within_second_best)):
        return False
    return True