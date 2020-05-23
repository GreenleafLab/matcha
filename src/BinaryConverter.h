#ifndef MATCHA_BINARY_CONVERTER_H
#define MATCHA_BINARY_CONVERTER_H

#include <cstdint>
#include <string>

#include <pybind11/numpy.h>
namespace py = pybind11;

using std::uint64_t;
using std::string;
using std::vector;

// Adapted from kallisto
// https://github.com/pachterlab/kallisto/blob/master/src/BUSData.cpp However,
// flag is a binary mask of whether a base == N
uint64_t stringToBinary(const string &s, uint64_t &flag);

uint64_t stringToBinary(const char *s, const size_t len, uint64_t &flag);

py::array_t<uint64_t> stringsToBinary(vector<string> strings, const size_t start, const size_t end);

string binaryToString(const uint64_t seq, const size_t len, const uint64_t flag);

#endif // MATCHA_BINARY_CONVERTER_H