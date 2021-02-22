#include <map>
#include <deque>
#include <cstring>
#include <iostream>

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include "collapser.h"

using namespace std;


struct cstr_equal{
    bool operator()(const char *a, const char *b) const {
        return strcmp(a, b) == 0;
    }
};

struct cstr_hash{
    // This is a simplification of cPython's __hash__() for str type
    long operator()(const char* str) const {
        int len;
        const char *p;
        long x;

        len = strlen(str);
        p = str;
        x = *p << 7;
        while (--len >= 0)
            x = (1000003*x) ^ *p++;
        x ^= len;

        return x;
    }
};

int main(int argv, char* args[]) {
    sequence_counter(args[1]);
    return 0;
}

int sequence_counter(char *fastq_file) {
    int fd = open (fastq_file, ios::in | ios::binary);
    if(fd == -1) {
        cerr << "Error on file open." << endl;
        return 0;
    } else {
        fstat (fd,&statbuf);
    }

    /* Advise the kernel of our access pattern.  */
    posix_fadvise (fd, 0, 0, POSIX_FADV_SEQUENTIAL);
    unordered_map<char*, size_t, cstr_hash, cstr_equal> counter(statbuf.st_size/(200));
    deque<__detail::_Node_iterator<pair<char *const, unsigned long>, false, true>> order;

    size_t BUFFER_SIZE = 16*1024;
    char* buf = new char[BUFFER_SIZE + 1];

    char* linestart;
    size_t lines = 1;
    size_t partial_read = 0;
    size_t bytes_read;

    while((bytes_read = read(fd, buf, BUFFER_SIZE)) > 0) {
        char *q = buf;
        char *end = q + bytes_read;
        if (!partial_read) linestart = q;

        /* memchr is more efficient than ptr search with lines >15.  */
        while ((q = (char *) memchr(q, '\n', end - q))) {

            if (lines % 2 == 0 && lines % 4 != 0) {

                buf[q - buf] = '\0';

                char *seq = new char[q - linestart + 1];
                memcpy(seq, linestart, q - linestart + 1);

                pair<__detail::_Node_iterator<pair<char *const, unsigned long>, false, true>, bool>
                where = counter.insert({seq, 1});

                if (!where.second) {
                    ++where.first->second;
                    free(seq);
                } else {
                    order.push_back(where.first);
                }
            }

            ++q;
            ++lines;
            linestart = q;
        }

        // Handle lines split across chunks here
        size_t last_partial = partial_read;
        partial_read = end - linestart;
        if (buf[BUFFER_SIZE] != '\n'){
            // "Unshrink" the buffer since it no longer needs to accommodate the last partial
            buf -= last_partial;
            BUFFER_SIZE += last_partial;
            // Copy the partial line to the beginning of the buffer
            // Correct linestart for the new location of the partial line
            // Advance the buffer pointer so that read() doesn't overwrite it
            // And shrink BUFFER_SIZE to accommodate the new partial line
            strncpy(buf, linestart, partial_read);
            linestart = buf;
            buf += partial_read;
            BUFFER_SIZE -= partial_read;
        } else {
            // There was no partial line this time
            // Unshrink the buffer if we had to accommodate a partial last time
            buf -= last_partial;
            BUFFER_SIZE += last_partial;
            partial_read = 0;
        }
    }


    /*size_t rec_len = strlen(counter.begin()->first) + strlen(": ") + 5 + 1;
    char *report = new char[counter.size() * rec_len];
    char *a = report;

    for (const auto& rec: order){
        size_t true_length = sprintf(a, "%s: %lu\n", rec->first, rec->second);
        a += true_length;
        free(rec->first);
    }

    cout << report;*/
    close(fd);
    return lines;
}
