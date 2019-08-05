#TODO we need not consider exchanges in db that are not being used
#TODO remove repetitive exchanges


import numpy as np
import pandas as pd
from pypardiso import spsolve
from pygsa.sampling.convert_distributions import *
from klausen.named_parameters import NamedParameters
import brightway2 as bw


class GSAinLCA:
    """
    Perform Global Sensivitiy Analysis (GSA) in Life Cycle Assessment (LCA).
    """

    def __init__(self, lca, inputs, parameters=None, ParametersModel=None):
        """
        Initialize class by defining all databases, running initial lcia, creating inputs and parameters dictionary.

        Attributes
        ----------
        lca : instance of bw.LCA class
            Contains information on the demand and LCIA method
        inputs : list of strings
            Contains strings that indicate user's choice of input factors for Sensitivity Analysis. 
            Can contain the following: 'technosphere', 'demand_exc', specific databases names, 'biosphere'.
            In all cases, all exchanges that correspond to the exchanges from the 'inputs' list will be included in the analysis.
            If 'biosphere' option is chosen, all biosphere flows that correspond to the 'inputs' activities will be included.
        parameters : dictionary or NamedParameters object
            Contains parameters for the parameterization of exchanges and flows.
            Should contain the following fields: 
                'name / code / ID' or other parameter identifier
                'uncertainty_type': ID of the distribution consistent with the stats_arrays package
                relevant characteristics of the distribution, such as 'loc', 'scale', 'minimum', etc
        ParametersModel : object
            Object that returns exchanges when running ParametersModel(parameters).
            Exchanges should be in the np.array format with the fields similar to np.dtype([ ('input_db', '<U20'),
                                                                                             ('input_code','<U32'),
                                                                                             ('output_db', '<U20'),
                                                                                             ('output_code', '<U32'),
                                                                                             ('amount', '<f4') ])

        Returns
        -------
        A new GSAinLCA object that contains self.: inputs, parameters, ParametersModel, databases, lca, 
                                                   inputs_dict, parameters_dict


        """

        self.inputs, self.parameters, self.ParametersModel = inputs, parameters, ParametersModel

        # 1 Run lci
        lca.lci()
        lca.lcia()
        lca.build_demand_array()
        self.lca = lca
        print('Initial LCA score: ' + str(lca.score))
        self.scores = np.array([])

        # 2 Extract databases names, TODO there should be a simpler way 
        databases = set()
        databases.update( set([key[0] for key in lca.activity_dict.keys()]) )
        databases.update( set([key[0] for key in lca.biosphere_dict.keys()]) )
        self.databases = databases

        # 3 Generate parameters dictionary
        if parameters != None and ParametersModel != None:
            if type(self.parameters) is not NamedParameters:
                #Initiate NamedParameters object
                self.parameters = NamedParameters(self.parameters)
            self.parameters.static()
            self.convert_named_parameters_to_array()
            self.obtain_parameterized_exchanges()

        # 4 Generate inputs dictionary
        self.split_inputs()

        #Generate pandas dataframe
        self.sa_pandas_init()


    def split_inputs(self):
        """
        Split information for each inputs option. Needed to identify which tech_params and bio_params are used in GSA.

        Returns
        -------
        A GSAinLCA object that contains self.inputs_dict, where inputs_dict is of the form: 
        { input_1: {
            tech_params:        np.array in the same format as lca.tech_params
            tech_params_where:  np.array(dtype=int),
            tech_n_params:      len(tech_params) TODO delete if not used
            bio_params:         np.array in the same format as lca.bio_params
            bio_params_where:   np.array(dtype=int),
            bio_n_params:       len(bio_params) TODO delete if not used
            },
          input_2: {...},
          input_3: {...},
          ...
        }
        """

        lca    = self.lca
        inputs = self.inputs

        inputs_dict = {}  # Only store exchanges with uncertainty

        # Keep track of which tech_params and bio_params are already included to the analysis
        # Needed to avoid running sa indices computations twice for the same tech or bio params. 
        # Initialize with parameterized exchanges
        if self.parameters != None and self.ParametersModel != None:
            indices_tech_all = self.parameters_dict['tech_params_where']
            indices_bio_all  = self.parameters_dict['bio_params_where']
        else:
            indices_tech_all = np.array([], dtype=int)
            indices_bio_all  = np.array([], dtype=int)

        for input_ in inputs:

            inputs_dict[input_] = {}

            indices_tech = np.array([], dtype=int)
            indices_bio  = np.array([], dtype=int)

            if input_ == 'technosphere':
                indices_tech = np.where(lca.tech_params['uncertainty_type']!=0)[0]
                if 'biosphere' in inputs:
                    indices_bio = np.where(lca.bio_params['uncertainty_type']!=0)[0]

            elif input_ == 'demand_exc':
                # Select all products that pertain to activities in the given demand vector
                for act_index in np.nonzero(lca.demand_array)[0]:
                    mask_tech = np.all([lca.tech_params['uncertainty_type']!=0, lca.tech_params['col']==act_index], axis=0)
                    indices_tech = np.concatenate([indices_tech, np.where(mask_tech)[0]])
                    if 'biosphere' in inputs:
                        mask_bio = np.all([lca.bio_params['uncertainty_type']!=0, lca.bio_params['col']==act_index], axis=0)
                        indices_bio = np.concatenate([indices_bio, np.where(mask_bio)[0]])

            elif input_ in self.databases:
                # Select all products and flows that are linked to the given database
                # Indices corresponding to exchanges in the tech_params depending on the given database
                db_act_indices_tech = [val for key,val in lca.activity_dict.items()  if key[0]==input_]
                if len(db_act_indices_tech) > 0:
                    db_act_index_min_tech = db_act_indices_tech[0]
                    db_act_index_max_tech = db_act_indices_tech[-1]
                    mask = lambda i : np.all( [lca.tech_params['uncertainty_type']!=0, 
                                               lca.tech_params['col']==i,
                                               lca.tech_params['amount']!=0], axis=0 )
                    indices_tech = [ np.where( mask(i) ) [0] for i in range(db_act_index_min_tech, db_act_index_max_tech+1) ]
                    indices_tech = np.concatenate(indices_tech)

                # Indices corresponding to flows in the biosphere params depending on the given database
                if 'biosphere' in inputs:
                    mask = lambda j : np.all( [lca.bio_params['uncertainty_type']!=0, lca.bio_params['col']==j], axis=0 )
                    indices_bio = [ np.where(mask(j))[0] for j in range(db_act_index_min_tech, db_act_index_max_tech+1) ]
                    indices_bio = np.concatenate(indices_bio)

            indices_tech = np.sort(indices_tech)
            indices_bio  = np.sort(indices_bio)

            # Do not add indices_tech that are already in the indices_tech_all
            indices_tech_same = np.intersect1d(indices_tech, indices_tech_all)
            pos_tech = np.array([ np.where(indices_tech==s)[0] for s in indices_tech_same ]).flatten()
            indices_tech = np.delete(indices_tech, pos_tech)
            np.append(indices_tech_all, indices_tech)

            # Do not add indices_bio that are already in the indices_bio_all
            indices_bio_same = np.intersect1d(indices_bio, indices_bio_all)
            pos_bio = np.array([ np.where(indices_bio==s)[0] for s in indices_bio_same ]).flatten()
            indices_bio = np.delete(indices_bio, pos_bio)
            np.append(indices_bio_all, indices_bio)
            
            inputs_dict[input_]['tech_params']       = lca.tech_params[indices_tech] #TODO maybe remove later, indices should be sufficient
            inputs_dict[input_]['tech_params_where'] = indices_tech
            inputs_dict[input_]['tech_n_params']     = len(indices_tech) #TODO remove later

            inputs_dict[input_]['bio_params']       = lca.bio_params[indices_bio] #TODO maybe remove later
            inputs_dict[input_]['bio_params_where'] = indices_bio
            inputs_dict[input_]['bio_n_params']     = len(indices_bio)


        self.indices_tech_all = indices_tech_all #TODO remove later
        self.indices_bio_all  = indices_bio_all
        self.inputs_dict = inputs_dict



    def convert_sample_to_proper_distribution(self,params,sample):
        """
        Map uniform samples on [0,1] to certain params and convert this sample to the correct distribution
        that is specified in the params np.array.

        Attributes
        ----------
        params : np.array 
            params dtype should contain 'uncertainty_type' and uncertainty/distribution information consistent with stats_arrays.
            Can be in the same format as lca.tech_params.
        sample : np.array
            Array that contains uniform samples on [0,1] with the same length as params.

        Returns
        -------
        converted_sample : np.array
            Sample with the correct distribution as specified in the params.

        """

        # Make sure that sample length is the same as the number of parameters #TODO change for group sampling
        assert len(params) == len(sample)

        uncertainties_dict = dict([ (uncert_choice, params[u'uncertainty_type'] == uncert_choice.id) \
                                    for uncert_choice in sa.uncertainty_choices \
                                    if any(params[u'uncertainty_type'] == uncert_choice.id) ])

        converted_sample = np.empty(len(params))

        for key in uncertainties_dict:
            mask = uncertainties_dict[key]
            converted_sample[mask] = convert_sample(params[mask], sample[mask]).flatten()

        return converted_sample



    def replace_non_parameterized_exchanges(self,sample):
        """
        Replace non parameterized exchanges, namely replace values in self.amounts_tech and self.amounts_bio 
        with the new sample values for all self.inputs. 
        self.i_sample iterates over a sample to select subsamples of the correct length for each option in inputs.

        Attributes
        ----------
        sample : np.array
            Array that contains uniform samples on [0,1] with the values for all params in all inputs.

        Returns
        -------
        A GSAinLCA object that contains resampled values of self.amounts_tech and self.amounts_bio.

        """

        for input_ in self.inputs:          

            # 1 Technosphere
            tech_params    = self.inputs_dict[input_]['tech_params']
            tech_n_params  = self.inputs_dict[input_]['tech_n_params']
            tech_subsample = sample[self.i_sample : self.i_sample+tech_n_params]

            self.i_sample += tech_n_params

            converted_tech_params = self.convert_sample_to_proper_distribution(tech_params,tech_subsample)
            np.put(self.amount_tech, self.inputs_dict[input_]['tech_params_where'], converted_tech_params)

            # 2 Biosphere
            bio_params    = self.inputs_dict[input_]['bio_params']
            bio_n_params  = self.inputs_dict[input_]['bio_n_params']
            bio_subsample = sample[self.i_sample : self.i_sample+bio_n_params]

            self.i_sample += bio_n_params

            converted_bio_params = self.convert_sample_to_proper_distribution(bio_params,bio_subsample)
            np.put(self.amount_bio, self.inputs_dict[input_]['bio_params_where'], converted_bio_params)



    def convert_named_parameters_to_array(self):
        """
        Convert parameters that are used in the parameterized exchanges to an np.array that contains uncertainty information
        of the parameters.

        Returns
        -------
        A GSAinLCA object that contains self.parameters_array with the dtype as specified in dtype_parameters.

        """

        dtype_parameters = np.dtype([ ('name', '<U40'), #TODO change hardcoded 40 here
                                      ('uncertainty_type', 'u1'), 
                                      ('amount', '<f4'),
                                      ('loc', '<f4'), 
                                      ('scale', '<f4'), 
                                      ('shape', '<f4'), 
                                      ('minimum', '<f4'), 
                                      ('maximum', '<f4'),
                                      ('negative', '?')  ])

        parameters_array = np.zeros(len(self.parameters), dtype_parameters)
        parameters_array[:] = np.nan
 
        for i, name in enumerate(self.parameters):
            parameters_array[i]['name']     = name
            parameters_array[i]['negative'] = False
            for k,v in self.parameters.data[name].items():
                parameters_array[i][k] = v

        self.parameters_array = parameters_array



    def obtain_parameterized_exchanges(self):
        """
        Get information about parameterized exchanges, which are obtained by running ParametersModel(parameters)
        in the form of a dictionary.

        Returns
        -------
        A GSAinLCA object that contains self.parameters_dict in the same format as self.inputs_dict.

        """

        lca = self.lca

        exchanges = self.ParametersModel.run(self.parameters)

        indices_tech = np.array([], dtype=int)
        indices_bio  = np.array([], dtype=int)

        get_input  = lambda exc: (exc['input_db'],  exc['input_code'])
        get_output = lambda exc: (exc['output_db'], exc['output_code'])

        exc_tech = np.array([exc for exc in exchanges if get_input(exc) in lca.activity_dict])
        if exc_tech.shape[0] != 0:
            mask_tech    = lambda i,j : np.where( np.all([lca.tech_params['row']==i, lca.tech_params['col']==j], axis=0) )
            indices_tech = np.hstack([ mask_tech( lca.activity_dict[get_input(exc)],lca.activity_dict[get_output(exc)] ) \
                                              for exc in exc_tech]) [0]

        exc_bio = np.array([exc for exc in exchanges if get_input(exc) in lca.biosphere_dict])
        if exc_bio.shape[0] != 0:
            mask_bio    = lambda i,j : np.where( np.all([lca.bio_params['row']==i, lca.bio_params['col']==j], axis=0) )
            indices_bio = np.hstack([ mask_bio( lca.biosphere_dict[get_input(exc)],lca.activity_dict[get_output(exc)] ) \
                                  for exc in exc_bio]) [0]
        parameters_dict = {}

        parameters_dict['tech_params_where']   = indices_tech
        parameters_dict['tech_params_amounts'] = np.array([ exc['amount'] for exc in exc_tech ])
        parameters_dict['tech_n_params'] = len(indices_tech)

        parameters_dict['bio_params_where']   = indices_bio
        parameters_dict['bio_params_amounts'] = np.array([ exc['amount'] for exc in exc_bio ])
        parameters_dict['bio_n_params'] = len(indices_bio)

        # TODO remove this check later on maybe
        assert indices_tech.shape[0] == parameters_dict['tech_n_params']
        assert indices_bio.shape[0]  == parameters_dict['bio_n_params']

        self.parameters_dict = parameters_dict



    def update_parameterized_exchanges(self, parameters):

        """
        Update parameterized exchanges by running ParametersModel(parameters) with the new converted parameters. 

        Attributes
        ----------
        parameters : NamedParameters object
            Contains parameters and their uncertainty information for the parameterization of exchanges. 

        Returns
        -------
        A GSAinLCA object that contains updated self.parameters_dict.

        """

        exchanges = self.ParametersModel.run(parameters)

        get_input  = lambda exc: (exc['input_db'],  exc['input_code'])

        exc_tech = np.array([exc for exc in exchanges if get_input(exc) in self.lca.activity_dict])
        exc_bio  = np.array([exc for exc in exchanges if get_input(exc) in self.lca.biosphere_dict])

        self.parameters_dict['tech_params_amounts'] = np.array([exc['amount'] for exc in exc_tech ])
        self.parameters_dict['bio_params_amounts']  = np.array([exc['amount'] for exc in exc_bio  ])




    def replace_parameterized_exchanges(self,sample):
        """
        Replace parameterized exchanges, namely replace self.amounts_tech and self.amounts_bio 
        after running the ParametersModel(parameters) with the new sample values for parameters. 

        Attributes
        ----------
        sample : np.array
            Array that contains uniform samples on [0,1] with the same length as parameters.

        Returns
        -------
        A GSAinLCA object that contains resampled values of self.amounts_tech and self.amounts_bio.

        """

        parameters_subsample = sample[self.i_sample : self.i_sample+len(self.parameters_array)]
        self.i_sample += len(self.parameters_array)

        # Convert uniform [0,1] sample to proper parameters distributions
        converted_parameters = self.convert_sample_to_proper_distribution(self.parameters_array,parameters_subsample)

        new_parameters = {}

        # Put converted values to parameters class, order of converted_parameters is the same as in parameters_array
        for i in range(len(self.parameters_array)):
            name = self.parameters_array[i]['name']
            new_parameters[name] = converted_parameters[i]

        # Update parameterized exchanges with the new converted values of parameters
        self.update_parameterized_exchanges(new_parameters)

        # Replace values in self.amounts_tech and self.amounts_bio
        np.put(self.amount_tech, self.parameters_dict['tech_params_where'], self.parameters_dict['tech_params_amounts'])
        np.put(self.amount_bio,  self.parameters_dict['bio_params_where'],  self.parameters_dict['bio_params_amounts'])
    


    def model(self,sample):

        """
        Rerun LCIA with a new sample and return score.

        Attributes
        ----------
        sample : np.array
            Array that contains uniform samples on [0,1] with the values for all parameters and all params from all inputs.

        Returns
        -------
        score : np.array
            Contains LCIA score for all LCIA methods. TODO probably can only do scalars atm

        """

        lca = self.lca
        
        self.amount_tech = lca.tech_params['amount']
        self.amount_bio  = lca.bio_params['amount']

        self.i_sample = 0
        self.replace_non_parameterized_exchanges(sample)
        self.replace_parameterized_exchanges(sample)

        lca.rebuild_technosphere_matrix(self.amount_tech)
        lca.rebuild_biosphere_matrix(self.amount_bio)

        score = (sum(lca.characterization_matrix)*lca.biosphere_matrix) * \
                spsolve(lca.technosphere_matrix,lca.demand_array)

        np.append(self.scores, score)

        return score

    

    def sa_pandas_init(self):
        """
        Initialize a dataframe to store sensitivity indices later on. 

        Returns
        -------
        A GSAinLCA object that contains self.sensitivity_indices_df dataframe with 
          columns: 'Products or flows' and 'Activities' corresponding to inputs and outputs of exchanges resp.
                   For parameters these values coincide.
          index:   consecutive numbers of the varied exchanges/parameters.

        """

        lca = self.lca

        ind_activity  = 0
        ind_product   = 1
        ind_biosphere = 2

        cols = []
        rows = []
        inputs = []

        #All exchanges in inputs
        for input_ in self.inputs:
            for i in self.inputs_dict[input_]['tech_params']:
                act  = lca.reverse_dict() [ind_activity] [i['col']]
                prod = lca.reverse_dict() [ind_product]  [i['row']]
                cols += [ bw.get_activity(act) ['name'] ]
                rows += [ bw.get_activity(prod)['name'] ]
                inputs += [input_]
            for j in self.inputs_dict[input_]['bio_params']:
                act = lca.reverse_dict() [ind_activity]  [j['col']]
                bio = lca.reverse_dict() [ind_biosphere] [j['row']]
                cols += [ bw.get_activity(act) ['name'] ]
                rows += [ bw.get_activity(prod)['name'] ]
                inputs += [input_]

        if self.parameters != None:
            # All parameters
            parameters_names_list = [name for name in self.parameters_array['name']]
            cols += parameters_names_list
            rows += parameters_names_list
            inputs += ['Parameters']*len(parameters_names_list)

        df = pd.DataFrame([inputs, rows, cols], index = ['Inputs', 'Products or flows', 'Activities'])
        df = df.transpose()

        self.sensitivity_indices_df = df



    def sa_pandas_append(self, sa_dict):
        """
        Update a dataframe with the new sensitivity indices taken from sa_dict dictionary.

        Attributes
        ----------
        sa_dict : dictionary
            Dictionary that contains sensitivity indices from some method.

        Returns
        -------
        A GSAinLCA object with the updated self.sensitivity_indices_df dataframe, where new columns include 
        computed sensitivity indices for all exchanges/parameters of interest.

        """

        df  = self.sensitivity_indices_df
        df_sa = pd.DataFrame(sa_dict, columns = sa_dict.keys(), index = df.index)
        self.sensitivity_indices_df = df.join(df_sa)










