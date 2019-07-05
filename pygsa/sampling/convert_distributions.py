from stats_arrays import LognormalUncertainty, NormalUncertainty, \
    DiscreteUniform, UniformUncertainty, TriangularUncertainty
from ..constants import Q_LOW, Q_HIGH    
import numpy as np
from itertools import compress

# TODO implement other distributions in stats_arrays

# TODO do we need this anymore?
def get_distr_indices_params(params,id_distr):
    list_ = params['uncertainty_type']==id_distr
    indices = list(compress(range(len(list_)), list_))
    params = params[indices]
    indices = np.array(indices,dtype=int)
    return indices,params



def convert_sample_to_lognor(params,sample):
    """
    Convert uniform in [0,1] to LOGNORMAL distribution
    """

    #Make sure all parameters have proper distributions
    assert np.all(params['uncertainty_type'] == LognormalUncertainty.id)
    
    if params["minimum"]:
        Q_MIN = LognormalUncertainty.cdf(params, params["minimum"])
    else:
        Q_MIN = Q_LOW
    if params["maximum"]:
        Q_MAX = LognormalUncertainty.cdf(params, params["maximum"])
    else:
        Q_MIN = Q_HIGH
        
    q = (Q_MAX - Q_MIN) * sample + Q_MIN
    
    params_converted = LognormalUncertainty.ppf(params, q)
    return params_converted
    



def convert_sample_to_normal(params,sample):
    """
    Convert uniform in [0,1] to NORMAL distribution
    """

    #Make sure all parameters have proper distributions
    assert np.all(params['uncertainty_type'] == NormalUncertainty.id)
    
    
    if params["minimum"]:
        Q_MIN = NormalUncertainty.cdf(params, params["minimum"])
    else:
        Q_MIN = Q_LOW
    if params["maximum"]:
        Q_MAX = NormalUncertainty.cdf(params, params["maximum"])
    else:
        Q_MIN = Q_HIGH
        
    q = (Q_MAX - Q_MIN) * sample + Q_MIN

    params_converted  = NormalUncertainty.ppf(params, q)
    return params_converted



def convert_sample_to_triang(params,sample):
    """
    Convert uniform in [0,1] to TRIANGULAR distribution
    """

    #Make sure all parameters have proper distributions
    assert np.all(params['uncertainty_type'] == TriangularUncertainty.id)

    params_converted = TriangularUncertainty.ppf(params, sample)
    return params_converted



def convert_sample_to_unifor(params,sample):
    """
    Convert uniform in [0,1] to UNIFORM distribution
    """

    #Make sure all parameters have proper distributions
    assert np.all(params['uncertainty_type'] == UniformUncertainty.id)

    params_converted = UniformUncertainty.ppf(params, sample)
    return params_converted



def convert_sample_to_dscr_u(params,sample):
    """
    Convert uniform in [0,1] to DISCRETE UNIFORM distribution
    """

    #Make sure all parameters have proper distributions
    assert np.all(params['uncertainty_type'] == DiscreteUniform.id)

    params_converted = DiscreteUniform.ppf(params, sample) 
    return params_converted







