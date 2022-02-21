

cdef inline void read_fastq(char* fname, vector[PyObject *] &out):
    cdef:
        FILE *cfile = fopen(fname, "rb")
        char *line = NULL
        size_t count = 1
        size_t l = 0

    if cfile == NULL:
        raise Exception("Blamo!")

    while getline(&line, &l, cfile) != -1:
        if count % 2 == 0 and count % 4 != 0:
            seq = ShortSeqFactory.from_chars(line)
            Py_XINCREF(<PyObject *>seq)
            out.push_back(<PyObject *>seq)

        count += 1


cdef inline void read_fastq_raw(char * fname, vector[char *] &out) nogil:
    cdef:
        FILE *cfile = fopen(fname, "rb")
        char *linecpy = NULL
        char *line = NULL
        size_t count = 1
        size_t l = 0

    if cfile == NULL:
        raise Exception("Blamo!")

    while getline(&line, &l, cfile) != -1:
        if count % 2 == 0 and count % 4 != 0:
            linecpy = strdup(line)
            out.push_back(linecpy)

        count += 1