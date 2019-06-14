'''
Functions that run Monte Carlo (MC) simulations and compute Sobol indices
Input: Sampler class that can generate samples (inputs) for the model,
        samples can be generated one at a time using 'next' and all simultaneously as a matrix of size NxD,
        where N is the # of MC runs and D - # of inputs
       Model, such that output = Model(sample)
'''

from sampling.sobol_sequence import SobolSample
import numpy as np
import time


#Compute Sobol indices while generating samples separately
def sobol_indices_one(Samp,Model):
    """C: Every function should have a docstring (like this) which tells what it does (if not clear from the name), the inputs arguments, and what the function returns."""

    """C: Can we say Sampler instead of Samp?"""

	#we generate N number of samples, actual number of MC runs is N//2

    """C: Can't we just use Samp.N?"""
	N = Samp.N
	D = Samp.D

    """C: Note that this will round down. May not be desired behaviour."""

    nruns = N//2

	next(Samp) #TODO change later, some bug

    """C: Consider using numpy arrays here, e.g. np.zeros(nruns). Much better
    performance for numerical stuff plus broadcasting."""

    """C: Ideally these would be defined in the docstring"""

	y_A = [0]*nruns
	y_B = [0]*nruns
	y_j = np.zeros([nruns,D])

	total_index = [0]*D
	first_index = [0]*D

    """C: Could also be:

    for i, sample_A, sample_B in zip(range(nruns), Samp, Samp):

    """

	for i in range(nruns): #compute total indices for all parameters

		sample_A = next(Samp)
		sample_B = next(Samp)

        """C: Need more documentation on what is returned here (vector? float?)"""

		y_A[i]   = Model(sample_A)
		y_B[i]   = Model(sample_B)

		for j in range(D):

            """C: Maybe need to make a copy?"""

			sample_j = sample_A
			sample_j[j] = sample_B[j]
			y_j[i][j] = Model(sample_j)

	for j in range(D):

        """C: This would be much faster if they were numpy arrays"""

        total_index[j] = [(y_A[i]-y_j[i,j])**2 for i in range(nruns)]

        """C: Normally would have some space here, e.g. x / y * z"""

		total_index[j] = sum(total_index[j])/2/nruns

		first_index[j] = [y_B[i]*(y_j[i,j]-y_A[i]) for i in range(nruns)]
		first_index[j] = sum(first_index[j])/nruns

	return first_index,total_index


#Compute Sobol indices while generating all samples simultaneously
def sobol_indices_all(Samp,Model):

	#we generate N number of samples, actual number of MC runs is N//2
    N = Samp.N
    D = Samp.D
    nruns = N//2

    samples = Samp.generate_all_samples() #change later, there's a bug

    y_A = [0]*nruns
    y_B = [0]*nruns
    y_j = [0]*nruns

    total_index = [0]*D
    first_index = [0]*D

    A = samples[:nruns]
    B = samples[nruns:]

    for i in range(nruns):
        y_A[i] = Model(A[i])
        y_B[i] = Model(B[i])

    for j in range(D): #compute total indices for all parameters

        J = A
        J[:,j] = B[:,j]

        for i in range(nruns):
            y_j[i] = Model(J[i])

        total_index[j] = [(y_A[k]-y_j[k])**2 for k in range(nruns)]
        total_index[j] = sum(total_index[j])/2/nruns

        first_index[j] = [y_B[k]*(y_j[k]-y_A[k]) for k in range(nruns)]
        first_index[j] = sum(first_index[j])/nruns

    return first_index,total_index


#wrapper function that chooses how to calculate sobol indices depending on the problem at stake
#for now, if problem has <1000 dimensions, generate all samples simultaneously. Otherwise, one sample at a time
def sobol_indices(N,D,Model,Samp=None):

	if D < 5000:

        """C: if Samp is None"""

        if Samp==None:
			Samp = SobolSample(N*2,D)
		t0 = time.time()
		first_index,total_index = sobol_indices_all(Samp,Model)
		t = time.time()-t0
		return t
	else:
		if Samp==None:
			Samp = SobolSample(N*2+1,D) #TODO +1 occurs due to a bug
		t0 = time.time()
		first_index,total_index = sobol_indices_one(Samp,Model)
		t = time.time()-t0
		return t














