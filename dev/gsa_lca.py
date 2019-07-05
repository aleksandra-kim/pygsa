import numpy as np
import pandas as pd
from stats_arrays import LognormalUncertainty, NormalUncertainty, \
    DiscreteUniform, UniformUncertainty, TriangularUncertainty, \
    UncertaintyBase, uncertainty_choices
from pypardiso import spsolve
from pygsa.sampling.convert_distributions import *
from klausen.named_parameters import NamedParameters


#TODO for non_parameterized activities: need to check that order of activities, 
#order of conversion from sampling and order when constructing database is the same.

class GSAinLCA:
    """
    Perform GSA in LCA
    """

    def __init__(self,lca,inputs,parameters=None,parameters_model=None,options=None):

        self.inputs, self.parameters, self.parameters_model = inputs, parameters, parameters_model

        #Extract databases names, TODO implement smth simpler 
        db = set()
        db.update( set([key[0] for key in lca.activity_dict.keys()]) )
        db.update( set([key[0] for key in lca.biosphere_dict.keys()]) )
        self.db = db

        #Run lci
        lca.lci()
        lca.build_demand_array()
        self.lca = lca

        #Generate inputs dictionary
        self.split_inputs()
        
        if self.parameters != None:
            if type(self.parameters) is not NamedParameters:
                
                #Initiate NamedParameters object
                self.parameters=NamedParameters(self.parameters)
                
            self.parameters.static()  
            
            #Generate parameters dictionary
            self.obtain_parameterized_activities()
            
        
        

    def split_inputs(self):
        
        #Split information of each input group
        
        lca = self.lca
        inputs = self.inputs

        inputs_dict = {} #only store inputs with uncertainty

        for input_ in inputs:

            inputs_dict[input_] = {}
            indices_tech = []
            indices_bio  = []

            if input_ == 'technosphere':
                indices_tech = np.where(lca.tech_params['uncertainty_type']!=0)[0]
                if 'biosphere' in inputs:
                    indices_bio = np.where(lca.bio_params['uncertainty_type']!=0)[0]

            elif input_ == 'demand_acts':
                #select all products that pertain to activities in the given demand
                # TODO these are being initialized two times.
                indices_tech = np.array([],dtype=int)
                indices_bio  = np.array([],dtype=int)
                for act_index in np.nonzero(lca.demand_array)[0]:
                    mask_tech = np.all([lca.tech_params['uncertainty_type']!=0, lca.tech_params['col']==act_index], axis=0)
                    indices_tech = np.concatenate([indices_tech,np.where(mask_tech)[0]])
                    if 'biosphere' in inputs:
                        mask_bio = np.all([lca.bio_params['uncertainty_type']!=0, lca.bio_params['col']==act_index], axis=0)
                        indices_bio = np.concatenate([indices_bio,np.where(mask_bio)[0]])
                indices_tech = np.sort(indices_tech)
                indices_bio = np.sort(indices_bio)
            
            #TODO we need not consider activities in db that are not being used
            elif input_ in self.db:
                #select all products that are linked to the given database
                #indices corresponding to the given database in the activity dictionary
                db_act_indices_tech = [val for key,val in lca.activity_dict.items()  if key[0]==input_]
                if len(db_act_indices_tech) > 0:
                    db_act_index_min_tech = db_act_indices_tech[0]
                    db_act_index_max_tech = db_act_indices_tech[-1]
                    mask = lambda i : np.all( [lca.tech_params['uncertainty_type']!=0, lca.tech_params['col']==i] , axis=0 )
                    indices_tech = [np.where( mask(i) ) [0] for i in range(db_act_index_min_tech,db_act_index_max_tech+1)]
                    indices_tech = np.sort(np.concatenate(indices_tech))

                #indices corresponding to flows in the biosphere params depending on the technosphere activities
                if 'biosphere' in inputs:
                    mask = lambda j : np.all( [lca.bio_params['uncertainty_type']!=0, lca.bio_params['col']==j] , axis=0 )
                    indices_bio = [np.where(mask(j))[0] for j in range(db_act_index_min_tech,db_act_index_max_tech+1)]
                    indices_bio = np.sort(np.concatenate(indices_bio))

                
            inputs_dict[input_]['tech_params']       = lca.tech_params[indices_tech]
            inputs_dict[input_]['tech_params_where'] = indices_tech
            inputs_dict[input_]['tech_n_params']     = len(indices_tech)

            inputs_dict[input_]['bio_params']       = lca.bio_params[indices_bio]
            inputs_dict[input_]['bio_params_where'] = indices_bio
            inputs_dict[input_]['bio_n_params']     = len(indices_bio)


        self.inputs_dict = inputs_dict



    def convert_sample_to_proper_distribution(self,params,sample):
              
        #Make sure that sample length is the same as the number of parameters #TODO change for group sampling
        assert len(params) == len(sample)
        
        # TODO When all distributions in stats_arrays are implemented, use:
        # uncertainties_dict = dict([(choices, params[u'uncertainty_type'] == choices.id) for choices in uncertainty_choices if any(params[u'uncertainty_type'] == choices.id)])
        
        # For now, we only include some uncertainty types
        uncertainty_types = LognormalUncertainty, NormalUncertainty, TriangularUncertainty, \
            UniformUncertainty, DiscreteUniform
        uncertainties_dict = dict([(types, params[u'uncertainty_type'] == types.id) for types in uncertainty_types if any(params[u'uncertainty_type'] == types.id)])
        
        converted_sample = np.zeros(len(params))
        
        # TODO is there a quicker and more elegant way to do this?
        for key in uncertainties_dict :
            if key ==  LognormalUncertainty:
                mask = uncertainties_dict[key]
                converted_sample[mask] = convert_sample_to_lognor(params[mask], sample[mask]).flatten()
            if key ==  NormalUncertainty:
                mask = uncertainties_dict[key]
                converted_sample[mask] = convert_sample_to_normal(params[mask], sample[mask]).flatten()
            if key ==  TriangularUncertainty:
                mask = uncertainties_dict[key]
                converted_sample[mask] = convert_sample_to_triang(params[mask], sample[mask]).flatten()
            if key ==  UniformUncertainty:
                mask = uncertainties_dict[key]
                converted_sample[mask] = convert_sample_to_unifor(params[mask], sample[mask]).flatten()
            if key ==  DiscreteUniform:
                mask = uncertainties_dict[key]
                converted_sample[mask]  = convert_sample_to_dscr_u(params[mask], sample[mask]).flatten()
                
        return converted_sample



    def replace_non_parameterized(self,sample):

        for input_ in self.inputs:            

            #1. Technosphere
            tech_params    = self.inputs_dict[input_]['tech_params']
            tech_n_params  = self.inputs_dict[input_]['tech_n_params']
            tech_subsample = sample[self.i_sample : self.i_sample+tech_n_params]

            self.i_sample += tech_n_params

            converted_tech_params = self.convert_sample_to_proper_distribution(tech_params,tech_subsample)
            np.put(self.amount_tech, self.inputs_dict[input_]['tech_params_where'], converted_tech_params)

            #2. Biosphere
            bio_params    = self.inputs_dict[input_]['bio_params']
            bio_n_params  = self.inputs_dict[input_]['bio_n_params']
            bio_subsample = sample[self.i_sample : self.i_sample+bio_n_params]

            self.i_sample += bio_n_params

            converted_bio_params = self.convert_sample_to_proper_distribution(bio_params,bio_subsample)
            np.put(self.amount_bio, self.inputs_dict[input_]['bio_params_where'], converted_bio_params)



    def convert_named_parameters_to_array(self):
             
        parameters_names = sorted([key for key in self.parameters])
        parameters_array = UncertaintyBase.from_dicts(*[self.parameters.data[key] for key in parameters_names])
        
        self.parameters_names = parameters_names
        self.parameters_array = parameters_array



    def obtain_parameterized_activities(self):

        lca = self.lca

        activities = self.parameters_model(self.parameters)

        ind_row = 0
        ind_col = 1
        ind_amt = 2

        acts_tech    = [act for act in activities if act[ind_row] in lca.activity_dict]
        mask_tech    = lambda i,j : np.where( np.all([lca.tech_params['row']==i, lca.tech_params['col']==j], axis=0) )
        indices_tech = np.hstack([ mask_tech( lca.activity_dict[act[ind_row]],lca.activity_dict[act[ind_col]] ) \
                                          for act in acts_tech]) [0]

        acts_bio    = [act for act in activities if act[ind_row] in lca.biosphere_dict]
        mask_bio    = lambda i,j : np.where( np.all([lca.bio_params['row']==i, lca.bio_params['col']==j], axis=0) )
        indices_bio = np.hstack([ mask_bio( lca.biosphere_dict[act[ind_row]],lca.activity_dict[act[ind_col]] ) \
                                  for act in acts_bio]) [0]

        parameters_dict = {}

        parameters_dict['ind_row'] = ind_row
        parameters_dict['ind_col'] = ind_col
        parameters_dict['ind_amt'] = ind_amt

        parameters_dict['tech_params_where']   = indices_tech
        parameters_dict['tech_params_amounts'] = np.array([ act[ind_amt] for act in acts_tech ])
        parameters_dict['tech_n_params'] = len(indices_tech)

        parameters_dict['bio_params_where']   = indices_bio
        parameters_dict['bio_params_amounts'] = np.array([ act[ind_amt] for act in acts_bio ])
        parameters_dict['bio_n_params'] = len(indices_bio)

        #TODO remove this check later on maybe
        assert indices_tech.shape[0] == parameters_dict['tech_n_params']
        assert indices_bio.shape[0]  == parameters_dict['bio_n_params']

        self.parameters_dict = parameters_dict



    def update_parameterized_activities(self):

        activities = self.parameters_model(self.new_parameters)

        ind_row = self.parameters_dict['ind_row']
        ind_amt = self.parameters_dict['ind_amt']

        acts_tech = [act for act in activities if act[ind_row] in self.lca.activity_dict]
        acts_bio  = [act for act in activities if act[ind_row] in self.lca.biosphere_dict]

        self.parameters_dict['tech_params_amounts'] = np.array([ act[ind_amt] for act in acts_tech ])
        self.parameters_dict['bio_params_amounts']  = np.array([ act[ind_amt] for act in acts_bio ])




    def replace_parameterized(self,sample):
        
        self.convert_named_parameters_to_array()

        parameters_subsample = sample[self.i_sample : self.i_sample+len(self.parameters_array)]
        self.i_sample += len(self.parameters_array)

        #Convert uniform [0,1] sample to proper parameters distributions
        converted_parameters = self.convert_sample_to_proper_distribution(self.parameters_array,parameters_subsample)
        
        #Make dictionary of new_parameters to replicate output of klausen 
        new_parameters = {}
   
        #Put converted values to parameters class, order of converted_parameters is the same as in parameters_array
        for i in range(len(self.parameters_array)):
                 
            name = self.parameters_names[i]
            new_parameters[name] = converted_parameters[i]
            self.new_parameters = new_parameters
                 
        #Obtain dictionary of parameterized tech_params and bio_params given the parameters_model
        self.update_parameterized_activities()

        np.put(self.amount_tech, self.parameters_dict['tech_params_where'], self.parameters_dict['tech_params_amounts'])
        np.put(self.amount_bio,  self.parameters_dict['bio_params_where'],  self.parameters_dict['bio_params_amounts'])
        


    def model(self,sample):

        lca = self.lca
        
        self.amount_tech = lca.tech_params['amount']
        self.amount_bio  = lca.bio_params['amount']

        self.i_sample = 0
        self.replace_non_parameterized(sample)
        self.replace_parameterized(sample)

        lca.rebuild_technosphere_matrix(self.amount_tech)
        lca.rebuild_biosphere_matrix(self.amount_bio)

        score = (sum(lca.characterization_matrix)*lca.biosphere_matrix) * \
                spsolve(lca.technosphere_matrix,lca.demand_array)

        return score











