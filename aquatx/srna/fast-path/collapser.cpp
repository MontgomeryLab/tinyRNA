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

uintmax_t sequence_counter(char *fastq_file) {
    int fd = open (fastq_file, ios::in | ios::binary);
    if(fd == -1) {
        cerr << "Error on file open." << endl;
        return 0;
    } else {
        fstat (fd,&statbuf);
    }

    char *fastq = (char*) mmap(nullptr, statbuf.st_size, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);
    if (fastq == MAP_FAILED){
        cerr << "Memory mapping the input fastq file failed." << endl;
        close(fd);
        return 0;
    }

    /* Advise the kernel of our access pattern.  */
    madvise(fastq, statbuf.st_size, MADV_SEQUENTIAL);
    map<char*, size_t, c_string_comparator> counter;

    char *q = fastq;
    char *seqstart = q;
    char *end = fastq + statbuf.st_size;
    uintmax_t lines = 1;

    /* memchr is more efficient with longer lines.  */
    while ((q = (char*) memchr (q, '\n', end - q)))
    {
        if (lines % 2 == 0 && lines % 4 != 0){
            // Rather than allocating and copying twice, change this sequence's newline to a null character
            // Then provide the map with a pointer to sequence start; mem-mapped file happily supports
            // Lower_bound and emplace_hint help avoid a double log(n) search and unnecessary copy
            fastq[q - fastq] = '\0';

            auto where = counter.lower_bound(seqstart);
            if (where == counter.end() || strcmp(where->first, seqstart) != 0){
                counter.emplace_hint(
                        where,
                        piecewise_construct,
                        forward_as_tuple(seqstart),
                        forward_as_tuple(1)
                );
            } else {
                ++where->second;
            }
        }

        ++q;
        ++lines;
        seqstart = q;
    }

    munmap(fastq, statbuf.st_size);
    close(fd);
    return lines;
}
