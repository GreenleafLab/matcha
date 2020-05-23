#include "ListMatcher.h"

void ListMatcher::add_sequence(uint64_t seq) {
    sequences.push_back(seq);
}

uint64_t ListMatcher::match(uint64_t seq, uint64_t flag, uint64_t &qual) {
    uint64_t best_match;
    uint64_t best_dist = max_dist;
    uint64_t next_dist = max_dist;

    for (size_t i = 0; i < sequences.size(); i ++) {
        uint64_t mismatches = hammingDistance(seq, flag, sequences[i]);
        if (mismatches < best_dist) {
            best_match = (uint64_t) i;
            next_dist = best_dist;
            best_dist = mismatches;
        } else if (mismatches < next_dist) {
            next_dist = mismatches;
        }
    }

    //qual format is: bottom N bits = # mismatches to best match, next N bits = # mismatches to 2nd best match
    qual = next_dist << dist_bits | best_dist;
    return best_match;
}