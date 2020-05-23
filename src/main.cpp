#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include "Matcher.h"
#include "ListMatcher.h"
#include "HashMatcher.h"
#include "BinaryConverter.h"
#include "FastqFile.h"

namespace py = pybind11;

PYBIND11_MODULE(_matcha, m) {
    m.doc() = R"pbdoc(
        Matcha C++ backend
        -----------------------

        .. currentmodule:: matcha

        .. autosummary::
           :toctree: _generate

    )pbdoc";

    m.def("stringToBinary", [](string s) {
        uint64_t flag = 0; 
        uint64_t seq = stringToBinary(s, flag); 
        return std::make_tuple(seq, flag);
    });

    m.def("stringsToBinary", &stringsToBinary);

    m.def("binaryToString", &binaryToString);

    //bindings to Matcher class
    py::class_<Matcher>matcher(m, "Matcher");
    matcher.def("add_sequences", &Matcher::add_sequences)
        .def("get_sequences", &Matcher::get_sequences)
        .def("match_all", &Matcher::matchAll)
        .def("match_raw", &Matcher::matchRaw)
        .def("has_labels", &Matcher::has_labels)
        .def("add_label", &Matcher::add_label)
        .def("add_labels", &Matcher::add_labels)
        .def("get_label", &Matcher::get_label)
        .def("get_labels", &Matcher::get_labels);

    py::class_<ListMatcher>(m, "ListMatcher", matcher)
        .def(py::init<>()) 
        .def("match", &ListMatcher::match);

    py::class_<HashMatcher>(m, "HashMatcher", matcher)
        .def(py::init<vector<uint64_t>, vector<vector<uint64_t>>, uint>()) 
        .def("match", &HashMatcher::match);

    py::class_<FastqFile>(m, "FastqFile")
        .def(py::init<string, vector<string>, vector<int>, string >())
        .def("read_chunk", &FastqFile::read_chunk, py::call_guard<py::gil_scoped_release>())
        .def("match", &FastqFile::match)
        .def("inspect_reads", &FastqFile::inspect_reads)
        .def("write_chunk", &FastqFile::write_chunk)
        .def("close", &FastqFile::close);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
