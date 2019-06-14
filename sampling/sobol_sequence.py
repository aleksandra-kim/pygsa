'''
Module that generates quasi-random Sobol sequences. They are needed for better space coverage.
'''

import numpy as np
import math
import sys

#local files
from . import directions

if sys.version_info[0] > 2:
    long = int

#Class that generates Sobol sequences
class SobolSample:
    
    def __init__(self,n_samples,n_dimensions,scale=31):

        if n_dimensions > len(directions.directions) + 1:
            raise ValueError("Error in Sobol sequence: not enough dimensions")

        L = int(math.ceil(math.log(n_samples) / math.log(2)))

        if L > scale:
            raise ValueError("Error in Sobol sequence: not enough bits")

        self.n_samples.   = n_samples
        self.n_dimensions = n_dimensions
        self.scale = scale
        self.L = L
        self.current = 0
        self.Y = np.array([long(0)]*self.n_dimensions)
        self.V = self.generate_V()
           
    def index_of_least_significant_zero_bit(self,value):
        index = 1
        while((value & 1) != 0):
            value >>= 1
            index += 1
        return index
        
    def generate_V(self):

        n_samples    = self.n_samples
        n_dimensions = self.n_dimensions    
        L = self.L    

        V = np.zeros([L+1,n_dimensions], dtype=long)
        V[:,0] = [0]+[1 << (scale - j) for j in range(1, L + 1)]
        
        for i in range(1,n_dimensions):
            m = np.array(directions.directions[i - 1], dtype=int)
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
        sample_one = np.zeros(self.n_dimensions)
        for i in range(self.n_dimensions):
            self.Y[i] ^= self.V[self.index_of_least_significant_zero_bit(self.current - 1),i]
            sample_one[i] = float(self.Y[i] / math.pow(2, scale))
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

        n_samples    = self.n_samples
        n_dimensions = self.n_dimensions    
        V = self.V 

        sample_all = np.zeros([n_samples,n_dimensions])

        X = long(0)
        for j in range(1, n_samples):
            X ^= V[self.index_of_least_significant_zero_bit(j - 1)]
            sample_all[j][:] = [float(x / math.pow(2, scale)) for x in X]
        return sample_all








