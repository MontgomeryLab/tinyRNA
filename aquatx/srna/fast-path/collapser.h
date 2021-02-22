//
// Created by t8 on 2/18/21.
//

#ifndef COLLAPSER_COLLAPSER_H
#define COLLAPSER_COLLAPSER_H

#include <unordered_map>

int sequence_counter(char*);

struct stat statbuf;

struct cstr{
    char *s = nullptr;
    long h = -1;

    void heap(){
        int len = strlen(s);
        s = new char[len+1];
        memcpy(s, s, len);
    }
};

#endif //COLLAPSER_COLLAPSER_H
