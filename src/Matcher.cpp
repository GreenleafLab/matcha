#include "Matcher.h"

void Matcher::add_sequences(vector<string> new_sequences) {
    uint64_t flag = 0;

    for (string s : new_sequences) {
        if (k == 0)
            k = s.size();
        else if (k != s.size())
            throw runtime_error("Sequence " + s + " does not match size");

        uint64_t seq = stringToBinary(s, flag);
        if (flag) throw runtime_error("Sequence " + s + " has N's");

        add_sequence(seq);
    }
}

vector<string> Matcher::get_sequences() {
    vector<string> ret;
    ret.reserve(sequences.size());
    for (uint64_t seq : sequences) {
        ret.push_back(binaryToString(seq, k, 0));
    }
    return ret;
}


void Matcher::_matchAll(const vector<string> &strings, const size_t start, const size_t end, uint64_t *out) {
    size_t n = strings.size();

    size_t len = end - start;
    for (size_t i = 0; i < strings.size(); i++) {
        uint64_t flag = 0;
        uint64_t seq = stringToBinary(strings[i].c_str() + start, len, flag);
        uint64_t qual = 0;
        uint64_t match_idx = match(seq, flag, qual);
        out[i] = match_idx;
        out[n + i] = qual;
    } 
}

py::array_t<uint64_t> Matcher::matchAll(vector<string> strings, const size_t start, const size_t end) {
    size_t n = strings.size();
    
    uint64_t *out = new uint64_t[n*2];

    py::capsule free_when_done(out, [](void *f) {
      auto mem = reinterpret_cast<uint64_t *>(f);
      delete mem;
    });

    // Allow the work to run in parallel
    py::gil_scoped_release release;
    Matcher::_matchAll(strings, start, end, out);
    py::gil_scoped_acquire acquire;    
    
    // Create numpy array that takes ownership of the data buffer
    py::array_t<uint64_t, py::array::c_style> result({(size_t) 2,  strings.size()}, out, free_when_done);
    return result;
}

void Matcher::matchRaw(py::array_t<uint64_t> seqs, py::array_t<uint64_t> output) {
    if (seqs.ndim() != 2 || output.ndim() != 2) throw runtime_error("Seqs and output must have ndim == 2");
    if (seqs.shape(1) != output.shape(1)) throw runtime_error("Seqs and output must have same number of columns");
    if (seqs.shape(0) != 2 || output.shape(0) != 2) throw runtime_error("Seqs and output must have 2 rows each");

    auto seq = seqs.unchecked<2>();
    auto res = output.mutable_unchecked<2>();
    for (size_t i = 0; i < seq.shape(1); i++) {
        uint64_t qual = 0;
        res(0,i) = match(seq(0,i), seq(1,i), qual);
        res(1,i) = qual;
    }
}

bool Matcher::has_labels() {
    return labels.size() == sequences.size();
}

void Matcher::add_label(string label) {
    labels.push_back(label);
}

void Matcher::add_labels(vector<string> new_labels) {
    for (string l : new_labels) {
        add_label(l);
    }
}

string Matcher::get_label(uint64_t index) {
    return labels[index];
}

vector<string> Matcher::get_labels(vector<uint64_t> indexes) {
    vector<string> ret;
    for (auto i : indexes) {
        ret.push_back(get_label(i));
    }
    return ret;
}