#ifndef MATCHA_FASTQ_FILE_H
#define MATCHA_FASTQ_FILE_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "gzstream/gzstream.h"
#include "Matcher.h"

namespace py = pybind11;

using std::string;
using std::vector;
using std::tuple;
using std::istream;
using std::ostream;
using std::ofstream;

namespace py = pybind11;

class FastqFile {
protected: 
    vector<string> name;
    vector<string> seq;
    vector<string> qual;
    igzstream in;
    bool use_out_gz;
    ofstream out_txt;
    ogzstream out_gz;
    vector<string> pattern_literals;
    vector<int> pattern_fields;
    vector<int> name_fields;
    vector<int> name_fields_lookup;
public:
    FastqFile(string in_path, vector<string> literals, vector<int> fields, string out_path = "");
    size_t read_chunk(size_t max_records);
    py::array_t<uint64_t> match(Matcher &m, const size_t start, const size_t end); // Match all sequences from last chunk read

    tuple<vector<string>, vector<string>, vector<string> > inspect_reads(); // Returns a tuple of the (name, seq, qual) vectors
    void write_chunk(py::array_t<bool> mask, vector<py::array_t<uint64_t>> sequence_matches, vector<Matcher> matchers);
    void close();
};

#endif // MATCHA_FASTQ_FILE_H