//
// Created by t8 on 2/18/21.
//

#ifndef COLLAPSER_COLLAPSER_H
#define COLLAPSER_COLLAPSER_H

int sequence_counter(char*);

struct stat statbuf;
struct c_string_comparator {
    bool operator()(const char* a, const char* b) const {
        return strcmp(a, b) < 0;
    }
};

#endif //COLLAPSER_COLLAPSER_H
