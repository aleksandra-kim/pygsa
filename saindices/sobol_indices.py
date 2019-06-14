'''
Functions that run Monte Carlo (MC) simulations and compute Sobol indices
Input: Sampler class that can generate samples (inputs) for the model, 
        samples can be generated one at a time using 'next' and all simultaneously as a matrix of size n_samples x n_dimensions, 
        where n_samples is the # of MC runs and n_dimensions - # of dimensions
       Model, such that output = Model(sample)
'''

from sampling.sobol_sequence import SobolSample
import numpy as np
import time


#Compute Sobol indices while generating samples separately 
def sobol_indices_one(Sampler,Model):

	n_samples    = Sampler.n_samples
	n_dimensions = Sampler.n_dimensions
	n_samples_2  = n_samples//2

	next(Sampler) #TODO change later, needed because of a bug

	y_A = [0]*n_samples_2
	y_B = [0]*n_samples_2
	y_j = np.zeros([n_samples_2,n_dimensions])

	total_index = [0]*n_dimensions
	first_index = [0]*n_dimensions

	for i in range(n_samples_2): #compute total indices for all parameters

		sample_A = next(Sampler)
		sample_B = next(Sampler)

		y_A[i]   = Model(sample_A)
		y_B[i]   = Model(sample_B)

		for j in range(n_dimensions):

			sample_j = sample_A
			sample_j[j] = sample_B[j]
			y_j[i][j] = Model(sample_j)

	for j in range(n_dimensions):
		total_index[j] = [(y_A[i]-y_j[i,j])**2 for i in range(n_samples_2)]
		total_index[j] = sum(total_index[j])/2/n_samples_2

		first_index[j] = [y_B[i]*(y_j[i,j]-y_A[i]) for i in range(n_samples_2)]
		first_index[j] = sum(first_index[j])/n_samples_2

	return first_index,total_index


#Compute Sobol indices while generating all samples simultaneously
def sobol_indices_all(Sampler,Model):

    n_samples    = Sampler.n_samples
    n_dimensions = Sampler.n_dimensions
    n_samples_2  = n_samples//2

    samples = Sampler.generate_all_samples()

    y_A = [0]*n_samples_2
    y_B = [0]*n_samples_2
    y_j = [0]*n_samples_2

    total_index = [0]*n_dimensions
    first_index = [0]*n_dimensions
    
    A = samples[:n_samples_2]
    B = samples[n_samples_2:]

    for i in range(n_samples_2):
        y_A[i] = Model(A[i])
        y_B[i] = Model(B[i])

    for j in range(n_dimensions): #compute total indices for all parameters
        
        J = A
        J[:,j] = B[:,j]

        for i in range(n_samples_2):
            y_j[i] = Model(J[i])

        total_index[j] = [(y_A[k]-y_j[k])**2 for k in range(n_samples_2)]
        total_index[j] = sum(total_index[j])/2/n_samples_2

        first_index[j] = [y_B[k]*(y_j[k]-y_A[k]) for k in range(n_samples_2)]
        first_index[j] = sum(first_index[j])/n_samples_2 
            
    return first_index,total_index


#wrapper function that chooses how to calculate sobol indices depending on the problem
#for now, if problem has <1000 dimensions, generate all samples simultaneously. Otherwise, one sample at a time
def sobol_indices(n_samples,n_dimensions,Model,Sampler=None):

	if n_dimensions < 1000: 
		if Sampler==None:
			Sampler = SobolSample(n_samples*2,n_dimensions)
		t0 = time.time()
		first_index,total_index = sobol_indices_all(Sampler,Model)
		t = time.time()-t0
		return t
	else: 
		if Sampler==None:
			Sampler = SobolSample(n_samples*2+1,n_dimensions) #TODO +1 occurs due to a bug
		t0 = time.time()
		first_index,total_index = sobol_indices_one(Sampler,Model)
		t = time.time()-t0
		return t














