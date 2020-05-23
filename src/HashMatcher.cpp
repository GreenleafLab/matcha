#include "HashMatcher.h"

// #include <iostream>
// #include <bitset>
// using std::cerr;
// using std::endl;
// using std::bitset;

HashMatcher::HashMatcher(vector<uint64_t> chunk_masks, vector<vector<uint64_t>> mismatch_masks, uint max_mismatches) {
    this->chunk_masks = chunk_masks;
    this->mismatch_masks = mismatch_masks;
    this->max_mismatches = max_mismatches;
    if (chunk_masks.size() != mismatch_masks.size()) {
        runtime_error("chunk_masks and mismatch_masks have different lengths");
    }
    // cerr << endl;
    for (size_t i = 0; i < chunk_masks.size(); i++) {
        chunk_indexes.push_back(unordered_multimap<uint64_t, uint32_t>());
        // cerr << "Chunk (" << i << ") = " << binaryToString(chunk_masks[i], 10, chunk_masks[i]) << endl;
        // cerr << "Mismatch (" << i << ") = " << endl;
        //for (auto m : mismatch_masks[i]) {
            // cerr << "\t" << binaryToString(m, 10, ~chunk_masks[i]) << endl;
        //}
    }
}

void HashMatcher::add_sequence(uint64_t seq)  {
    size_t seq_index = sequences.size();
    sequences.push_back(seq);

    for (size_t i = 0; i < chunk_masks.size(); i++) {
        // cerr << "Inserting seq_index " << seq_index << " into chunk " << i << " under " << binaryToString(seq, k, ~chunk_masks[i]) << endl;
        chunk_indexes[i].insert({seq & chunk_masks[i], seq_index});
    }
}

//qual format is: 
//  - bottom 6 bits = # mismatches to best match, 
//  - next 6 bits = # mismatches to 2nd best match

uint64_t HashMatcher::match(uint64_t seq, uint64_t flag, uint64_t &qual)  {
    uint64_t best_match = -1;
    uint64_t best_dist = max_dist;
    uint64_t next_dist = max_dist;
    // cerr << "Matching seq = " << binaryToString(seq, k, flag) << endl;
    for (size_t i = 0; i < chunk_masks.size(); i++) {
        uint64_t chunk_mask = chunk_masks[i];
        for (uint64_t mismatch_mask : mismatch_masks[i]) {
            uint64_t query_seq = (seq ^ mismatch_mask) & chunk_mask;
                // cerr << "\tmismatch_mask = " << binaryToString(mismatch_mask, k, ~chunk_mask) << 
                // " query_seq = " << binaryToString(query_seq, k, ~chunk_mask) <<
                // " count = " << chunk_indexes[i].count(query_seq) << endl;
            auto ret = chunk_indexes[i].equal_range(query_seq);
            for (auto it = ret.first; it != ret.second; it++) {
                uint32_t candidate_idx = it->second;
                if (candidate_idx == best_match) continue;
                uint64_t mismatches = hammingDistance(seq, flag, sequences[candidate_idx]);
                // cerr << "\t\tcandidate_idx = " << candidate_idx << " distance = " << mismatches << endl;
                if (mismatches > max_mismatches) {
                    continue;
                } else if (mismatches == best_dist) {
                    best_match = std::min(best_match, (uint64_t) candidate_idx);
                    next_dist = best_dist;
                } else if (mismatches < best_dist) {
                    best_match = (uint64_t) candidate_idx;
                    next_dist = best_dist;
                    best_dist = mismatches;
                } else if (mismatches < next_dist) {
                    next_dist = mismatches;
                }       
            }
        }
    }
    
    //qual format is: bottom N bits = # mismatches to best match, next N bits = # mismatches to 2nd best match
    qual = next_dist << dist_bits | best_dist;
    return best_match;
}
