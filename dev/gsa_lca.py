import numpy as np
import pandas as pd
from pypardiso import spsolve
from pygsa.sampling.convert_distributions import *
from klausen.named_parameters import NamedParameters

#Local files
from pygsa.constants import *


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
				print(indices_tech)

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
		#params should have uncertainty information in the same format as stats_array

		#Make sure that sample length is the same as the number of parameters #TODO change for group sampling
		assert len(params) == len(sample)

		#Info on distributions
		indices_lognor, params_lognor = get_distr_indices_params(params,ID_LOGNOR)
		indices_normal, params_normal = get_distr_indices_params(params,ID_NORMAL)
		indices_triang, params_triang = get_distr_indices_params(params,ID_TRIANG)
		indices_unifor, params_unifor = get_distr_indices_params(params,ID_UNIFOR)
		indices_dscr_u, params_dscr_u = get_distr_indices_params(params,ID_DSCR_U)
		n_lognor = len(params_lognor)
		n_normal = len(params_normal)
		n_triang = len(params_triang)
		n_unifor = len(params_unifor)
		n_dscr_u = len(params_dscr_u)

		#Lognormal
		i_start, i_end = 0, n_lognor
		params_lognor_converted = convert_sample_to_lognor(params_lognor,sample[i_start:i_end])

		#Normal
		i_start, i_end = i_end, i_end+n_normal
		params_normal_converted = convert_sample_to_normal(params_normal,sample[i_start:i_end])

		#Triangular
		i_start, i_end = i_end, i_end+n_triang
		params_triang_converted = convert_sample_to_triang(params_triang,sample[i_start:i_end])		

		#Uniform
		i_start, i_end = i_end, i_end+n_unifor
		params_unifor_converted = convert_sample_to_unifor(params_unifor,sample[i_start:i_end])		

		#Discrete uniform
		i_start, i_end = i_end, i_end+n_dscr_u
		params_dscr_u_converted = convert_sample_to_unifor(params_dscr_u,sample[i_start:i_end])	

		all_indices = np.concatenate([ indices_lognor,
									   indices_normal,
									   indices_triang,
									   indices_unifor,
									   indices_dscr_u  ])

		all_params  = np.concatenate([ params_lognor_converted,
									   params_normal_converted,
									   params_triang_converted,
									   params_unifor_converted,
									   params_dscr_u_converted  ])

		#Construct converted sample
		converted_sample = np.zeros(len(params))
		np.put(converted_sample,all_indices,all_params)

		return converted_sample



	def replace_non_parameterized(self,sample):

		#TODO remove repetitive params

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

		dtype_parameters = np.dtype([ ('name', '<U40'), #TODO change hardcoded 40 here
  									  ('uncertainty_type', 'u1'), 
  									  ('amount', '<f4'),
									  ('loc', '<f4'), 
									  ('scale', '<f4'), 
									  ('shape', '<f4'), 
									  ('minimum', '<f4'), 
									  ('maximum', '<f4')  ])

		parameters_array = np.empty(len(self.parameters),dtype_parameters)
		parameters_array[:] = np.nan

		for i,name in enumerate(self.parameters):
			parameters_array[i]['name'] = name
			for k,v in self.parameters.data[name].items():
				#to avoid confusion when parameters amounts were generated before, disregard amounts completely
				if not k == 'amount':
					parameters_array[i][k] = v

		self.parameters_array = parameters_array



	def obtain_parameterized_activities(self):

		lca = self.lca

		self.parameters.static() #need to run because this is how klausen works
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

		parameters_dict['tech_params_where']   = indices_tech
		parameters_dict['tech_params_amounts'] = np.array([ act[ind_amt] for act in acts_tech ])
		parameters_dict['tech_n_params'] = len(indices_tech)

		parameters_dict['bio_params_where']    = indices_bio
		parameters_dict['bio_params_amounts']  = np.array([ act[ind_amt] for act in acts_bio ])
		parameters_dict['bio_n_params'] = len(indices_bio)

		#TODO remove this check later on maybe
		assert indices_tech.shape[0] == parameters_dict['tech_n_params']
		assert indices_bio.shape[0]  == parameters_dict['bio_n_params']

		self.parameters_dict = parameters_dict



	def replace_parameterized(self,sample):

		if type(self.parameters) is NamedParameters:

			self.convert_named_parameters_to_array()

			parameters_subsample = sample[self.i_sample : self.i_sample+len(self.parameters_array)]
			self.i_sample += len(self.parameters_array)

			#Convert uniform [0,1] sample to proper parameters distributions
			converted_parameters = self.convert_sample_to_proper_distribution(self.parameters_array,parameters_subsample)

			#Put converted values to parameters class, order of converted_parameters is the same as in parameters_array
			for i in range(len(self.parameters_array)):
				name = self.parameters_array[i]['name']
				self.parameters.data[name]['amount'] = converted_parameters[i]

			#Obtain dictionary of parameterized tech_params and bio_params given the parameters_model
			# self.obtain_parameterized_activities() #TODO decide where to put this, maybe in init is better

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

		return (sum(lca.characterization_matrix)*lca.biosphere_matrix) * \
				spsolve(lca.technosphere_matrix,lca.demand_array)











