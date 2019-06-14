'''
Module that generates quasi-random Sobol sequences. They are needed for better space coverage.
'''

"""C: Can we have tests for expected values (e.g. from sampler)?"""

"""C: Better to load/save directions.py as numpy array."""


import numpy as np
import math
import sys

#local files

"""C: from .directions import directions"""

from . import directions

"""C: Delete, no need for py2 compatibility"""

if sys.version_info[0] > 2:
    long = int

#Class that generates Sobol sequences
class SobolSample:

    global scale
    scale = 31

    def __init__(self,N,D):

        """C: if len(directions.directions) < D:"""

        if D > len(directions.directions) + 1:
            raise ValueError("Error in Sobol sequence: not enough dimensions")

        L = int(math.ceil(math.log(N) / math.log(2)))

        if L > scale:
            raise ValueError("Error in Sobol sequence: not enough bits")

        """Can do multiple assignment at once, e.g.

        self.N, self.D, self.L = N, D, L

        """

        self.N = N
        self.D = D
        self.L = L
        self.current = 0
        self.Y = np.array([long(0)]*self.D)
        self.V = self.generate_V()

    def index_of_least_significant_zero_bit(self,value):
        index = 1
        while((value & 1) != 0):
            value >>= 1
            index += 1
        return index

    def generate_V(self):

        N = self.N
        D = self.D
        L = self.L

        """C: dtype = int"""

        V = np.zeros([L+1,D], dtype=long)
        """C: First element already zero.

        V[1:,0] = [1 << (scale - j) for j in range(1, L + 1)]

        """
        V[:,0] = [0]+[1 << (scale - j) for j in range(1, L + 1)]

        """C: ???

        for i in range(D):
            m = np.array(directions.directions[i], dtype=int)

        """

        """C: Faster to already have as array:

        In [1]: from directions import directions
           ...:

        In [2]: import numpy as np
           ...:

        In [3]: a = np.zeros((len(directions), max([len(row) for row in directions])))
           ...:

        In [4]: for i, row in enumerate(directions):
           ...:     a[i, :len(row)] = row
           ...:

        In [5]: def get_with_mask(array):
           ...:     for row in array:
           ...:         return row[row > 0]
           ...:

        In [6]: def get_with_count(array):
           ...:     for row in array:
           ...:         return row[:np.count_nonzero(row)]
           ...:

        In [7]: def convert(lst):
           ...:     for row in lst:
           ...:         np.array(row, dtype=int)
           ...:

        In [8]: %timeit get_with_mask(a)
           ...:
        1.95 µs ± 51.6 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)

        In [9]: %timeit convert(directions)
           ...:
        49 ms ± 3.04 ms per loop (mean ± std. dev. of 7 runs, 10 loops each)

        In [10]: %timeit get_with_count(a)
            ...:
        792 ns ± 18.6 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)

        """

        for i in range(1,D):
            m = np.array(directions.directions[i - 1], dtype=int)

            """C: Maybe easier to use these value directly instead of creating
            new variables. If you really want new variables, give them names
            people can understand."""

            a = m[0]
            s = len(m) - 1

            # The following code discards the first row of the ``m`` array
            # Because it has floating point errors, e.g. values of 2.24e-314
            if L <= s:
                V[:,i] = [0]+[1 << (scale-j) for j in range(1, L+1)]
            else:
                V[1:s+1,i] = [m[j] << (scale-j) for j in range(1,s+1)]
                for j in range(s + 1, L + 1):
                    V[j,i] = V[j-s,i] ^ (V[j-s,i] >> s)
                    for k in range(1, s):
                        V[j,i] ^= ((a >> (s - 1 - k)) & 1) * V[j-k][i]
        return V

    def generate_sample(self):
        sample_one = np.zeros(self.D)

        """C: Huge gains here from using numpy functions at once instead of a Python loop"""

        for i in range(self.D):
            self.Y[i] ^= self.V[self.index_of_least_significant_zero_bit(self.current - 1),i]
            sample_one[i] = float(self.Y[i] / math.pow(2, scale))
        self.current += 1
        return sample_one

    def __iter__(self):
        return self

    def __next__(self):
        if self.current > self.N-1:
            raise StopIteration
        elif self.current == 0:
            self.current = 1
            return self.Y
        else:
            return self.generate_sample()

    def generate_all_samples(self):

        N = self.N
        D = self.D
        V = self.V

        sample_all = np.zeros([N,D])

        X = long(0)
        for j in range(1, N):
            X ^= V[self.index_of_least_significant_zero_bit(j - 1)]
            sample_all[j][:] = [float(x / math.pow(2, scale)) for x in X]
        return sample_all








