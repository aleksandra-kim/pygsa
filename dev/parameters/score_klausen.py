# Imports. ARE THERE OTHERS?
import numpy as np
import stats_arrays as sa

id_distr = {
            "lognormal":        sa.LognormalUncertainty.id,
            "normal":           sa.NormalUncertainty.id,
            "uniform":          sa.UniformUncertainty.id,
            "triangular":       sa.TriangularUncertainty.id,
            "discrete uniform": sa.DiscreteUniform.id,
           }

# Functions to convert from unifr [0,1] to other distributions
def convert_to_normal (params, sample):
    q = (q_high-q_low)*sample + q_low
    return norm.ppf(q,loc=params['loc'],scale=params['scale'])
         
def convert_to_lognormal (params, sample):
    q = (q_high-q_low)*sample + q_low
    return lognorm.ppf(q, s=params['scale'], scale=np.exp(params['loc']))
    
def convert_to_triang (params, sample):
    q     = sample[:len(indices_type["triang"])]
    loc   = params_type["triang"]['minimum']
    scale = params_type["triang"]['maximum']-params_type["triang"]['minimum']
    c     = (params_type["triang"]['loc']-loc)/scale
    params_triang_new = triang.ppf(q,c=c,loc=loc,scale=scale)

def convert_to_uniform (params, sample):
    return (sample * (params["maximum"] - params["minimum"])) + params["minimum"]

def convert_to_discrete_uniform (params, 8sample):   
    return np.rint((sample * (params["maximum"]-1 - params["minimum"])) #-1 because stats array excludes the upper bound
                           +  params["minimum"])

def convert_to_distr (samples, params, model, distributions = id_distr):
    # This function converts samples into the distributions 
    # of each parameter defined in the NamedParameter object params.
    # And then return the new parameters according to the function model.
    
    # Model needs to return new paramaters and the linked activities in the following format:
    # [
    #    (input, output, values)
    # ....
    # ]
    
    # For now id_distr is not a parameter
    new_params = {}
    
    
    for i, key in enumerate(params):
        if params[key]["ucertainty_type"] == id_distr["normal"]:
            new_params[key] = convert_to_normal(params, samples[i])
        if params[key]["ucertainty_type"] == id_distr["lognormal"]:
            new_params[key] = convert_to_lognormal (params, samples[i])
        if params[key]["ucertainty_type"] == id_distr["triangular"]:
            new_params[key] = convert_to_triang (params, samples[i])
        if params[key]["ucertainty_type"] == id_distr["uniform"]:
            new_params[key] = convert_to_uniform (params, samples[i])
        if params[key]["ucertainty_type"] == id_distr["discrete uniform"]:
            new_params[key] = convert_to_discrete_uniform (params, samples[i])
        
        #Apply model and return new_parameters
        return model(new_params)

def find_bio_tech (new_params_from_model):
    new_params_tech, new_params_bio = []
    
    for p in new_params_from_model:
        if p[0][0] != "biosphere3":
            new_params_tech.append()
        else:
            new_params_bio.append(p)
    
    return new_params_tech, new_params_bio

def get_indices_from_model(array, new_params_type):
    
    # new_params_type is either new_params_tech or new_params_bio
    # array is either tech_params or bio_params
    # This functions works with model that returns value as described above      
    
    index = []
    
    for p in new_params_type:
            mask_output = array["col"] == lca.activity_dict[p[1]]
            mask_input  = array["row"] == lca.activity_dict[p[0]]
            index_d = list(compress(range(len(mask_output)), mask_output))
            index_i = list(compress(range(len(mask_input)), mask_input))
            index.append(set(index_d) & set(index_i))
    
    return index
    
def replace_params(array, new_params_type, indices_type):
   
   # array is either tech_params or bio_params
   
   vector_amount = array["amount"]
    
   for i, p in enumerate(new_params_type):
       vector_new = np.put(vector_amount, indices[i], p[1])
         
   return vector_new
           

def score_from_sample(lca, bio_vector_new, tech_vector_new):
    
    # This function calculates new scores 
                                                 
    lca.rebuild_biosphere_matrix(bio_vector_new):
    lca.rebuild_technosphere_matrix(tech_vector_new)
    
    c = sum(lca.characterization_matrix)
    score = ((c*lca.biosphere_matrix))*spsolve(lca.technosphere_matrix,d)[0]
    
    return score
    
             
           