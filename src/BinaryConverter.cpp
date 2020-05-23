#include "BinaryConverter.h"

// Adapted from kallisto
// https://github.com/pachterlab/kallisto/blob/master/src/BUSData.cpp However,
// flag is a 2-bit encoded binary mask: 00 if base != N, 01 if base == N
// A = 00, T = 11, G = 10, C = 01
// Lowest bits are at the beginnig of the string (e.g. 1100 = "AT", 1001="CG")
uint64_t stringToBinary(const std::string &s, uint64_t &flag) {
    return stringToBinary(s.c_str(), s.size(), flag);
}

uint64_t stringToBinary(const char *s, const size_t len, uint64_t &flag) {
    uint64_t r = 0;
    
    flag = 0;

    size_t k = len;
    if (k > 32) {
        k = 32;
    }

    for (size_t i = 0; i < k; i++) {
        uint64_t x = ((*s) & 4) >> 1;
        if (((*s) & 3) == 2) {
            flag |= 1 << (2 * i);
        }
        r |= (x + ((x ^ (*s & 2)) >> 1)) << (2 * i);
        s++;
    }
    return r;
}

py::array_t<uint64_t> stringsToBinary(vector<string> strings, const size_t start, const size_t end) {
    py::array_t<uint64_t, py::array::c_style> result({(size_t) 2,  strings.size()});
    auto r = result.mutable_unchecked<2>();

    size_t len = end - start;
    for (size_t i = 0; i < strings.size(); i++) {
        uint64_t flag = 0;
        uint64_t seq = stringToBinary(strings[i].c_str() + start, len, flag);
        r(0, i) = seq;
        r(1, i) = flag;
    } 
    return result;
}

char binary_decoder[4] = {'A', 'C', 'G', 'T'};
string binaryToString(const uint64_t seq, const size_t len, const uint64_t flag) {
    size_t k = len;
    if (k > 32) {
        k = 32;
    }

    string s(k,'N');

    for (size_t i = 0; i < k; i++) {
        if (flag & 1 << (2*i)) continue; // We have an N flagged, so leave as N
        else s[i] = binary_decoder[seq >> (2*i) & 3];
    }

    return s;
}