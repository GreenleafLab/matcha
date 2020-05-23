def choose(n, k):
    """n choose k"""
    prod = 1
    for i in range(1, k+1):
        prod *= n - i + 1
        prod /= i
    return prod

def cost(n, k, b, r):
    """n = # index elements, k = bp length, b = # of bins to use, r = Max # mismatches to tolerate
    Note: assumes cost of mismatch check = cost of hash table lookup
    """
    rprime = r // b
    a = r % b
    #
    s = k // b # Length of sequence to use (some will be s+1)
    long_b = k - b*s
    short_b = b - long_b
    #
    assert short_b * s + long_b * (s+1) == k
    #
    bin_len = [s] * short_b + [s+1] * long_b
    local_r = [r] * (a + 1) + [r-1] * (b - a - 1)
    #
    def cost_per_bin(b, r):
        """Cost for lookups in bin b, radius r"""
        return (1 + n/4**b) * sum(3 * choose(b, i) for i in range(r+1))
    #
    return sum(cost_per_bin(b, r) for b, r in zip(bin_len, local_r))

def optimal_bins(n, k, r):
    """Estimate the optimal number of bins to use for a barcode matching algorithm
    n = # of valid barcodes to match against
    k = bp length of barcode
    r = Maximum number of mismatches to tolerate
    """
    best_cost = cost(n, k, 1, r)
    b = 1
    while True:
        b += 1
        if cost(n, k, b, r) < best_cost:
            best_cost = cost(n, k, b, r)
        else:
            return b - 1