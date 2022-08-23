from ._stepvector import StepVector
from HTSeq._HTSeq import GenomicInterval

cdef class ChromVector(object):
    """Counting vector covering a chromosome.

    This class supports only the tinyRNA Cython implementation of the
    HTSeq StepVector
    """

    cdef public object array
    cdef public GenomicInterval iv
    cdef public int offset
    cdef public bint is_vector_of_sets
    cdef public str _storage
    cdef public str typecode
    cdef public str memmap_dir

    @classmethod
    def create(cls, GenomicInterval iv, str typecode, str storage, str memmap_dir=""):
        """Create ChromVector from GenomicInterval

        Args:
            iv (GenomicInterval): A GenomicInterval describing the chromosome
              vector.
            typecode ('d', 'i', 'l', 'b', or 'O'): What kind of data will be
              stored inside this chromosome vector. 'd' for double, 'i' for int,
              'l' for long int, 'b' for boolean, 'O' for arbitrary objects
              (e.g. sets).
            storage ('step', 'ndarray', or 'memmap'): What kind of storage to
              use. 'ndarray' is appropriate for short chromosomes and stores
              each position in the genome into memory. 'memmap' stores all
              positions, but maps the memory onto disk for larger chromosomes.
              'step' is a sparse representation similar to CSR matrices whereby
              only the boundaries between genomic stretches with differing
              data content are stored - see HTSeq.StepVector.
            memmap_dir (str): If using 'memmap' storage, what folder to store
              the memory maps. These can get quite big.

        Returns:
            An instance of ChromVector with the requested options.

        """
        ncv = cls()
        ncv.iv = iv
        ncv.array = StepVector.create()

        ncv._storage = storage
        ncv.typecode = typecode
        # NOTE: As long as autochromosomes in GenomicArray are infinite length
        # this has pretty limited use, but that might change
        ncv.offset = iv.start
        ncv.is_vector_of_sets = False
        ncv.memmap_dir = memmap_dir
        return ncv

    @classmethod
    def _create_view(cls, ChromVector vec, GenomicInterval iv):
        if iv.length == 0:
            raise IndexError("Cannot subset to zero-length interval.")
        v = cls()
        v.iv = iv
        v.array = vec.array
        v.offset = vec.offset
        v.is_vector_of_sets = vec.is_vector_of_sets
        v._storage = vec._storage
        return v

    def extend_to_include(self, iv):
        if iv.strand != self.iv.strand:
            raise ValueError(
                'The new interval must match the current strandedness',
            )

        # Step 1: extend the interval
        length = self.iv.length
        startdiff = max(self.iv.start - iv.start, 0)
        self.iv.extend_to_include(iv)
        self.offset = self.iv.start

        # Step 2: extend the array if needed, and shift-copy the old values
        if self._storage == 'ndarray':
            if self.typecode != 'O':
                array = numpy.zeros(shape=(self.iv.length,), dtype=self.typecode)
            else:
                array = numpy.empty(shape=(self.iv.length,), dtype=self.typecode)
                array[:] = None
            array[startdiff: startdiff + length] = self.array[:]
        elif self._storage == 'memmap':
            array = numpy.memmap(
                shape=(self.iv.length,), dtype=self.typecode,
                filename=os.path.join(
                    self.memmap_dir,
                    self.iv.chrom + self.iv.strand + str(self.iv.start) + '_' \
                        + str(self.iv.length) + ".nmm"),
                mode='w+',
            )
            array[startdiff: startdiff + length] = self.array[:]
        else:
            # The StepVector is created in ChromVector.create without explicit
            # boundaries, so it's already bound by 0, +inf. So we do not need
            # to extend it here, but rather just set the slice to the right
            # value
            array = self.array
        self.array = array

    def __getitem__(self, index):
        """Index or slice the chromosome.

        The index can be a few things:
        - an integer: get the value of the vector at that chromosome coordinate
        - a 1-step slice e.g "4:7": get a view of the chromosome region
          between those coordinates. The array data are not copied.
        - a GenomicInterval: similar to slices, with the additional choice of
          strandedness. If this argument is stranded but the chromosome itself
          is not stranded, a nonstranded view of the chromosome region is
          returned.

        """
        cdef slice index_slice
        cdef long int index_int
        cdef long int start, stop
        cdef GenomicInterval iv

        if isinstance(index, int):
            index_int = index
            if index_int < self.iv.start or index_int >= self.iv.end:
                raise IndexError
            return self.array[index_int - self.offset]

        elif isinstance(index, slice):
            index_slice = index
            if index_slice.start is None:
                start = self.iv.start
            else:
                start = index_slice.start
                if start < self.iv.start:
                    raise IndexError("start too small")

            if index_slice.stop is None:
                stop = self.iv.end
            else:
                stop = index_slice.stop
                if stop > self.iv.end:
                    raise IndexError("stop too large")

            iv = GenomicInterval(self.iv.chrom, start, stop, self.iv.strand)

            if not self.iv.contains(iv):
                raise IndexError
            return ChromVector._create_view(self, iv)

        elif isinstance(index, GenomicInterval):
            if not self.iv.contains(index):
                raise IndexError

            if self.iv.strand is strand_nostrand and \
                    index.strand is not strand_nostrand:
                iv = index.copy()   # Is this correct now?
                iv.strand = strand_nostrand
            else:
                iv = index

            return ChromVector._create_view(self, iv)

        else:
            raise TypeError("Illegal index type")

    def __setitem__(self, index, value):
        cdef slice index_slice
        cdef long int start, stop

        if isinstance(value, ChromVector):
            if self.array is value.array and value.iv.start == index.start and \
                    value.iv.end == index.stop and (index.step is None or index.step == 1):
                return
            else:
                raise NotImplementedError(
                    "Required assignment signature not yet implemented.")

        if isinstance(index, int):
            self.array[index - self.iv.start] = value

        elif isinstance(index, slice):
            index_slice = index
            if index_slice.start is not None:
                start = index_slice.start
                if start < self.iv.start:
                    raise IndexError("start too small")
            else:
                start = self.iv.start
            if index_slice.stop is not None:
                stop = index_slice.stop
                if stop > self.iv.end:
                    raise IndexError("stop too large")
            else:
                stop = self.iv.end
            if start > stop:
                raise IndexError("Start of interval is after its end.")
            if start == stop:
                raise IndexError("Cannot assign to zero-length interval.")
            self.array[start - self.offset: stop -
                       self.iv.start: index.step] = value

        elif isinstance(index, GenomicInterval):
            if index.chrom != self.iv.chrom:
                raise KeyError("Chromosome name mismatch.")
            if self.iv.strand is not strand_nostrand and \
                    self.iv.strand is not self.index.strand:
                raise KeyError("Strand mismatch.")
            self.array[index.iv.start - self.iv.start,
                       index.iv.end - self.iv.start] = value
        else:
            raise TypeError("Illegal index type")

    def __iadd__(self, value):
        if not self.is_vector_of_sets:
            self.array[self.iv.start - self.offset: self.iv.end -
                       self.offset].__iadd__(value)
        else:
            def addval(x):
                y = x.copy()
                y.add(value)
                return y

            self.apply(addval)
        return self

    def __iter__(self):
        return self.values()

    def values(self):
        return iter(self.array[self.iv.start - self.offset: self.iv.end - self.offset])

    def steps(self):
        return _HTSeq_internal.ChromVector_steps(self)

    def apply(self, fun):
        for iv, value in self.steps():
            self.array[iv.start - self.offset: iv.end -
                       self.offset] = fun(value)

    def __repr__(self):
        return "<%s object, %s, %s>" % (self.__class__.__name__, str(self.iv), self._storage)

    def __reduce__(self):
        assert self.__class__ is ChromVector
        return(_ChromVector_unpickle,
               (self.array, self.iv, self.offset, self.is_vector_of_sets, self._storage))


def _ChromVector_unpickle(array, iv, offset, is_vector_of_sets, _storage):
    cv = ChromVector()
    cv.array = array
    cv.iv = iv
    cv.offset = offset
    cv.is_vector_of_sets = is_vector_of_sets
    cv._storage = _storage
    return cv