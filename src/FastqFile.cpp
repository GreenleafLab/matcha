#include "FastqFile.h"

using namespace std;

// From Stack Overflow 
static bool endsWith(const std::string& str, const std::string& suffix)
{
    return str.size() >= suffix.size() && 0 == str.compare(str.size()-suffix.size(), suffix.size(), suffix);
}

FastqFile::FastqFile(string in_path, vector<string> literals, vector<int> fields, string out_path) {
    in.open(in_path.c_str());
    if (in.fail()) throw invalid_argument("Could not open file: " + in_path);
    if (out_path.size() > 0) {
        if (endsWith(out_path, ".gz")) {
            use_out_gz = true;
            out_gz.open(out_path.c_str());
            if (out_gz.fail()) throw invalid_argument("Could not open file: " + out_path);
        } else {
            use_out_gz = false;
            out_txt.open(out_path.c_str());
            if (out_txt.fail()) throw invalid_argument("Could not open file: " + out_path);
        }
    }

    pattern_literals = literals;
    pattern_fields = fields;
    for (auto f : fields) {
        if (f < -1) {
            size_t field_idx = -f - 2;
            //Double-check this name field hasn't been seen already
            if(std::find(name_fields.begin(), name_fields.end(), field_idx) != name_fields.end()) continue;
            name_fields.push_back(field_idx);
        }
    }
    if (name_fields.size() == 0) return;
    
    // Ensure name_fields[name_fields_lookup[x]] = x for all used fields
    sort(name_fields.begin(), name_fields.end());
    name_fields_lookup.resize(name_fields[name_fields.size() - 1]);
    
    size_t current_field = 0;
    for (size_t i = 0; i < name_fields.size(); i++) {
        if (i == (size_t) name_fields[current_field]) {
            name_fields_lookup[i] = current_field;
            current_field += 1;
        }
    }
}

size_t FastqFile::read_chunk(size_t max_records) {
    if (name.size() != max_records) {
        name.resize(max_records);
        seq.resize(max_records);
        qual.resize(max_records);
    }
    size_t lines_read = 0;
    string temp_name, temp_seq, temp_qual;

    for (;lines_read < max_records; lines_read++) {
        getline(in, temp_name);
        getline(in, temp_seq);
        getline(in, temp_qual); // Dummy to read the + symbol
        getline(in, temp_qual);

        if(in.eof()) {
            break;
        }

        name[lines_read] = temp_name.substr(1);
        seq[lines_read] = temp_seq;
        qual[lines_read] = temp_qual;
    }

    if (lines_read != max_records) {
        name.resize(lines_read);
        seq.resize(lines_read);
        qual.resize(lines_read);
    }
    return lines_read;
}

py::array_t<uint64_t> FastqFile::match(Matcher &m, const size_t start, const size_t end) {
    return m.matchAll(seq, start, end);
}

tuple<vector<string>, vector<string>, vector<string> > FastqFile::inspect_reads() {
    return make_tuple(name, seq, qual);
}

void FastqFile::write_chunk(py::array_t<bool> mask, vector<py::array_t<uint64_t>> raw_matches, vector<Matcher> matchers) {
    //Make a buffer to read in parts of the read name
    vector<string> parsed_name_fields(name_fields.size());
    
    vector<py::detail::unchecked_reference<uint64_t,1>> matches;
    matches.reserve(raw_matches.size());
    for (size_t i = 0; i < raw_matches.size(); i++) {
        matches.push_back(raw_matches[i].unchecked<1>());
    }

    ostream *out;
    if (use_out_gz) out = &out_gz;
    else out = &out_txt;

    auto m = mask.unchecked<1>();

    py::gil_scoped_release release;
    for (int i = 0; i < m.shape(0); i++) {
        if (!m[i]) continue;
        // Make the output name
        *out << "@" << pattern_literals[0];

        // Parse fields from the read name
        if (name_fields.size()) {
            int current_field = 0;
            size_t pos = 0;
            for (int j = 0; j < name_fields[name_fields.size() - 1]; j++) {
                size_t nextpos = name[i].find(":", pos);
                if (name_fields[current_field] == j) {
                    parsed_name_fields[current_field] = name[i].substr(pos, nextpos);
                }
                pos = nextpos;
            }
        }
        for (size_t j = 0; j < pattern_fields.size(); j++) {
            //output field
            int f = pattern_fields[j];
            if (f == -1) {
                *out << name[i];
            } else if (f < -1) {
                *out << parsed_name_fields[name_fields_lookup[-f - 2]];
            } else {
                *out << matchers[f].get_label(matches[f][i]);
            }
            //output final literal
            *out << pattern_literals[j+1];
        }

        *out << "\n" << seq[i] << "\n+\n" << qual[i] << "\n";
    }
    *out << flush;
    py::gil_scoped_acquire acquire;   
}

void FastqFile::close() {
    in.close();
    if (use_out_gz)
        out_gz.close();
    else
        out_txt.close();
}