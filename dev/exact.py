from SALib.sample import sobol_sequence
import numpy as np
from copy import copy
import pandas as pd

def exact_first_total(problem, N, Model):

    num_vars = problem['num_vars']

    # sample
    skip_values = 1000
    base_sequence = sobol_sequence.sample(N+skip_values, 2*num_vars)
    A = base_sequence[skip_values:, :num_vars]
    B = base_sequence[skip_values:, num_vars:]
    
    # variables to store lca scores
    y_A = np.zeros(N)
    y_B = np.zeros(N)
    y_j = np.zeros(N)

    # Monte Carlo simulations without resampling
    y_A = [Model(A[i,:]) for i in range(N)]
    y_B = [Model(B[i,:]) for i in range(N)]

    var_y = np.var( np.array([y_A,y_B]).flatten() )
    
    first_index = np.zeros(num_vars)
    total_index = np.zeros(num_vars)

    # First order
    for j in range(num_vars):
        J1 = copy(B) # all factors change but j-th
        expectation_1 = np.zeros(N)
        for i in range(N):
            J1[:, j] = copy(A[i, j]) # j-th factor stays the same
            y_J1 = [Model(J1[k,:]) for k in range(N)]
            expectation_1[i] = np.mean(y_J1)
        first_index[j] = np.var(expectation_1) / var_y
    
    # Total order
    for j in range(num_vars):
        expectation_T = np.zeros(N)
        for i in range(N):
            JT_row = copy(A[i,:]) # all factors stay the same but j-th
            JT = np.kron(np.ones(N), JT_row).reshape(N,num_vars)
            JT[:,j] = copy(B[:, j])
            y_JT = [Model(JT[k,:]) for k in range(N)]
            expectation_T[i] = np.mean(y_JT)
        total_index[j] = 1-np.var(expectation_T) / var_y

    df = pd.DataFrame([first_index, total_index], index = ['S1 exact', 'ST exact'])
    df = df.transpose()

    return df