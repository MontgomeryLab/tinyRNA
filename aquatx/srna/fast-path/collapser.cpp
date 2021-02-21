#include <map>
#include <cstring>
#include <iostream>

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/mman.h>

#include "collapser.h"

using namespace std;

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
    map<char*, size_t, c_string_comparator> counter;

    size_t BUFFER_SIZE = 16*1024;
    char* buf = new char[BUFFER_SIZE + 1];
    size_t lines = 1;
    size_t last_partial = 0;
    size_t partial_read = 0;

    char* linestart;
    size_t bytes_read;

    while((bytes_read = read(fd, buf, BUFFER_SIZE)) > 0) {

        char *q = buf;
        char *end = q + bytes_read;
        if (!partial_read) linestart = q;

        /* memchr is more efficient with longer lines.  */
        while ((q = (char *) memchr(q, '\n', end - q))) {

            if (lines % 2 == 0 && lines % 4 != 0) {
                // Rather than allocating and copying twice, change this sequence's newline to a null character
                // Then provide the map with a pointer to sequence start; mem-mapped file happily supports
                // Lower_bound and emplace_hint help avoid a double log(n) search and unnecessary copy
                buf[q - buf] = '\0';
                char *seq = new char[q - linestart];
                strcpy(seq, linestart);

                auto where = counter.lower_bound(linestart);
                if (where == counter.end() || strcmp(where->first, linestart) != 0) {
                    counter.emplace_hint(
                            where,
                            piecewise_construct,
                            forward_as_tuple(seq),
                            forward_as_tuple(1)
                    );
                } else {
                    ++where->second;
                }
            }

            ++q;
            ++lines;
            linestart = q;
        }

        // Handle sequences split across chunks here
        last_partial = partial_read;
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
            last_partial = 0;
        }
    }

    for (auto rec = counter.rbegin(); rec != counter.rend(); ++rec){
        cout << rec->first << ": " << rec->second << endl;
    }

    close(fd);
    return lines;
}
