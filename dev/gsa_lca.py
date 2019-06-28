#TODO ask people what should vary in foreground and background, wrt databases


import numpy as np
from pypardiso import spsolve
from pygsa.sampling.convert_distributions import *
from klausen.named_parameters import NamedParameters

#Local files
from pygsa.constants import *


class GSAinLCA:
	"""
	Perform GSA in LCA
	"""

	def __init__(self,lca,inputs,options=None):

		self.inputs = inputs

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



	def split_inputs(self):
		#Split information of each input group

		lca = self.lca
		inputs = self.inputs

		dtype_tech = lca.tech_params.dtype
		dtype_bio  = lca.bio_params.dtype

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



	def replace_parameterized(self,parameters,sample):

		if type(params) is NamedParameters:
			

		np.put(self.amount_tech, self.inputs_dict[input_]['tech_params_where'], converted_tech_params)

		return
		


	def model(self,sample):

		lca = self.lca
		
		self.amount_tech = lca.tech_params['amount']
		self.amount_bio  = lca.bio_params['amount']

		self.i_sample = 0
		self.replace_non_parameterized(sample)
		self.replace_parameterized(sample)


		lca.rebuild_technosphere_matrix(self.amount_tech)
		lca.rebuild_biosphere_matrix(self.amount_bio)

		A = lca.technosphere_matrix
		B = lca.biosphere_matrix
		c = sum(lca.characterization_matrix)
		d = lca.demand_array

		score = (c*B)*spsolve(A,d) #run it before MC to factorize matrix A

		return score















