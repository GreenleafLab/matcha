#ifndef MATCHA_HASH_MATCHER_H
#define MATCHA_HASH_MATCHER_H

#include <algorithm>
#include <unordered_map>

#include "Matcher.h"

using std::unordered_multimap;

// Worker backend of HashMatcher algorithm as described here: https://arxiv.org/pdf/1307.2982.pdf
class HashMatcher: public Matcher {
private:
    uint max_mismatches; // Only used to limit what matches are returned
    vector<uint64_t> chunk_masks;
    vector<vector<uint64_t>> mismatch_masks; // Masks to xor with lookup chunk to get neighboring sequences
    vector<unordered_multimap<uint64_t, uint32_t>> chunk_indexes;
public:
    // chunk_masks -- List of masks to be bitwise-anded to extract chunks of input sequences
    // mismatch_masks -- List lists of masks to be xor-ed with with chunks to get neighboring mismatches
    HashMatcher(vector<uint64_t> chunk_masks, vector<vector<uint64_t>> mismatch_masks, uint max_mismatches);
    void add_sequence(uint64_t seq) override;
    uint64_t match(uint64_t seq, uint64_t flag, uint64_t &qual) override; //qual format is: bottom 6 bits = # mismatches to best match, next 6 bits = # mismatches to 2nd best match
};

#endif // MATCHA_HASH_MATCHER_H