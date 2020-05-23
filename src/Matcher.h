#ifndef MATCHA_MATCHER_H
#define MATCHA_MATCHER_H

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "BinaryConverter.h"



namespace py = pybind11;

using std::uint64_t;
using std::string;
using std::vector;
using std::unordered_map;
using std::runtime_error;

class Matcher {
protected: 
    size_t k = 0; // Length of barcode
    vector<uint64_t> sequences; // List of barcode sequences
    vector<string> labels; // (optional) List of names for sequences (same order as sequences vector)
public:
    virtual ~Matcher() {}
    void add_sequences(vector<string> sequences); // Add all sequences to matcher
    vector<string> get_sequences(); // Get list of sequences in matcher
    py::array_t<uint64_t> matchAll(vector<string> strings, const size_t start, const size_t end); // Match all sequences in a list
    void matchRaw(py::array_t<uint64_t> seqs, py::array_t<uint64_t> output); // Used for benchmarking

    bool has_labels();
    void add_label(string label);
    void add_labels(vector<string> label);
    string get_label(uint64_t index);
    vector<string> get_labels(vector<uint64_t> indexes);

    virtual void add_sequence(uint64_t seq) {throw runtime_error("Not Implemented");}; // Add barcode sequence to match against
    virtual uint64_t match(uint64_t seq, uint64_t flag, uint64_t &qual) {throw runtime_error("Not Implemented");}; // Return the index of closest matching barcode to seq + quality
private:
    void _matchAll(const vector<string> &strings, const size_t start, const size_t end, uint64_t *out); //Inner worker for matchAll, safe without holding GIL
};


// Compute hamming distance between seq+flag and barcode
inline uint64_t hammingDistance(uint64_t seq, uint64_t flag, uint64_t barcode) {
    uint64_t diff = barcode ^ seq; // at least 1 bit set for each mismatched group of 2
    diff = (diff | diff >> 1 | flag) & 0x5555555555555555; // Lower bit set in each group of 2 with a mismatch or N
    
    uint64_t mismatches = __builtin_popcountll(diff);
    
    return mismatches;
}

const uint dist_bits = 6;
const uint max_dist = (1 << dist_bits) - 1;

#endif // MATCHA_MATCHER_H