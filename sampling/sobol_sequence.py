"""
Module that generates quasi-random Sobol sequences. They are needed for a more uniform space coverage.
"""

"""C: Can we have tests for expected values (e.g. from sampler)?"""
"""S: In which format? I already ran stress tests."""

import numpy as np
import math
import sys

#local files
from .directions import directions


class SobolSample:
    """
    Generate Sobol quasi random sequences of size n_dimension one at a time using method __next__ or 
    all samples simultaneously as a matrix of size n_samples x n_dimensions.

    Literature
    ----------
    [2003] Joe and Kuo
        Remark on Algorithm 659: Implementing Sobol's quasirandom sequence generator
        https://doi.org/10.1145/641876.641879
    [2008] Joe and Kuo
        Constructing Sobol sequences with better two-dimensional projections
        https://doi.org/10.1137/070709359

    """

    def __init__(self,n_samples,n_dimensions,scale=31):

        if n_dimensions > len(directions) + 1:
            raise ValueError("Error in Sobol sequence: not enough dimensions")

        L = int(math.ceil(math.log(n_samples) / math.log(2)))

        if L > scale:
            raise ValueError("Error in Sobol sequence: not enough bits")

        self.n_samples, self.n_dimensions, self.scale, self.L = n_samples, n_dimensions, scale, L
        self.current = 0
        self.Y = np.array([int(0)]*self.n_dimensions)
        self.V = self.generate_V()


    def index_of_least_significant_zero_bit(self,value):
        """
        Generate index of the least significant zero bit of a value.

        Parameters
        ----------
        Value : int

        Returns
        -------
        Index : int
            Index of the least significant zero bit

        """

        index = 1
        while((value & 1) != 0):
            value >>= 1
            index += 1
        return index


    def generate_V(self):
        """
        Compute matrix V that is needed for Sobol quasi random sequences generation.

        Returns
        -------
        V : ndarray

        """

        n_samples, n_dimensions, L = self.n_samples, self.n_dimensions, self.L

        V = np.zeros([L+1,n_dimensions], dtype=int)
        V[1:,0] = [1 << (self.scale - j) for j in range(1, L + 1)]

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
           ...:         row[row > 0]
           ...:

        In [6]: def get_with_count(array):
           ...:     for row in array:
           ...:         row[:np.count_nonzero(row)]
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
        """S: this part I didn't understand :("""

        for i in range(n_dimensions-1):

            m = np.array(directions[i], dtype=int)
            s = len(m) - 1

            # The following code discards the first row of the ``m`` array
            # Because it has floating point errors, e.g. values of 2.24e-314
            if L <= s:
                V[1:,i+1] = [m[j] << (self.scale-j) for j in range(1, L+1)]
            else:
                V[1:s+1,i+1] = [m[j] << (self.scale-j) for j in range(1,s+1)]
                for j in range(s + 1, L + 1):
                    V[j,i+1] = V[j-s,i+1] ^ (V[j-s,i+1] >> s)
                    for k in range(1, s):
                        V[j,i+1] ^= ((m[0] >> (s - 1 - k)) & 1) * V[j-k][i+1]

        return V


    def generate_sample(self):
        """
        Generate an array of size n_dimensions that contains one samples of Sobol quasi random sequence.

        Returns
        -------
        sample_one : array
            Array of size n_dimensions that contains one Sobol quasi random sequence

        """

        sample_one = np.zeros(self.n_dimensions)

        """C: Huge gains here from using numpy functions at once instead of a Python loop"""
        """S: How..?"""

        for i in range(self.n_dimensions):
            self.Y[i] ^= self.V[self.index_of_least_significant_zero_bit(self.current - 1),i]
            sample_one[i] = float(self.Y[i] / math.pow(2, self.scale))
        self.current += 1
        return sample_one


    def __iter__(self):
        return self


    def __next__(self):
        if self.current > self.n_samples-1:
            raise StopIteration
        elif self.current == 0:
            self.current = 1
            return self.Y
        else:
            return self.generate_sample()


    def generate_all_samples(self):
        """
        Generate a matrix that contains n_samples number of samples for n_dimensions number of dimensions.

        Returns
        -------
        sample_all : ndarray
            Matrix of size n_samples x n_dimensions that contains Sobol quasi random sequences

        """

        n_samples, n_dimensions, V = self.n_samples, self.n_dimensions, self.V
        sample_all = np.zeros([n_samples,n_dimensions])

        X = int(0)
        for j in range(1, n_samples):
            X ^= V[self.index_of_least_significant_zero_bit(j - 1)]
            sample_all[j][:] = [float(x / math.pow(2, self.scale)) for x in X]
        return sample_all








