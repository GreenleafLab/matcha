#ifndef MATCHA_LIST_MATCHER_H
#define MATCHA_LIST_MATCHER_H

#include "Matcher.h"

class ListMatcher: public Matcher {
public:
    void add_sequence(uint64_t seq) override;
    uint64_t match(uint64_t seq, uint64_t flag, uint64_t &qual) override; //qual format is: bottom 6 bits = # mismatches to best match, next 6 bits = # mismatches to 2nd best match
};

#endif // MATCHA_LIST_MATCHER_H