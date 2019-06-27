#TODO ask people what should vary in foreground and background, wrt databases


import numpy as np
from scipy.stats import lognorm, norm, triang
from itertools import compress
from pypardiso import spsolve


#Local files
from pygsa.constants import *


class GSAinLCA:
	"""
	Perform GSA in LCA
	"""

	def __init__(self,lca,inputs,options=None):

		self.inputs = inputs

		#Extract databases names
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

		self.q_low  = (1-THREE_SIGMA_Q)/2
		self.q_high = (1+THREE_SIGMA_Q)/2



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



	def get_distr_indices_params(self,params,id_distr):
		list_ = params['uncertainty_type']==id_distr
		indices = list(compress(range(len(list_)), list_))
		params = params[indices]
		return indices,params



	def convert_sample_to_proper_distribution(self,params,sample):
		#params should have uncertainty information in the same format as stats_array

		assert len(params) == len(sample)

		#Info on distributions
		indices_lognor,params_lognor = self.get_distr_indices_params(params,ID_LOGNOR)
		indices_normal,params_normal = self.get_distr_indices_params(params,ID_NORMAL)
		indices_triang,params_triang = self.get_distr_indices_params(params,ID_TRIANG)
		n_lognor = len(params_lognor)
		n_normal = len(params_normal)
		n_triang = len(params_triang)

		converted_sample = np.zeros(len(params))

		#Convert uniform to lognormal distribution
		q = (self.q_high-self.q_low)*sample[ : n_lognor] + self.q_low
		params_lognor_converted = lognorm.ppf(q,s=params_lognor['scale'],scale=np.exp(params_lognor['loc']))
		del q

		#Convert uniform to normal distribution
		q = (self.q_high-self.q_low)*sample[n_lognor : n_lognor+n_normal] + self.q_low
		params_normal_converted  = norm.ppf(q,loc=params_normal['loc'],scale=params_normal['scale'])
		del q

		#Convert uniform to triangular distribution
		q = sample[n_lognor+n_normal : n_lognor+n_normal+n_triang]
		loc   = params_triang['minimum']
		scale = params_triang['maximum']-params_triang['minimum']
		c     = (params_triang['loc']-loc)/scale
		params_triang_converted = triang.ppf(q,c=c,loc=loc,scale=scale)

		#Construct converted sample
		np.put(converted_sample,indices_lognor,params_lognor_converted)
		np.put(converted_sample,indices_normal,params_normal_converted)
		np.put(converted_sample,indices_triang,params_triang_converted)

		return converted_sample



	def replace_non_parameterized():
		return amount_tech, amount bio



	def replace_parameterized():
		return amount_tech, amount bio



	def model(self,sample):

		lca = self.lca
		inputs_dict = self.inputs_dict

		amount_tech = lca.tech_params['amount']
		amount_bio  = lca.bio_params['amount']

		i_sample = 0

		for input_ in self.inputs:

			#Technosphere
			tech_params    = inputs_dict[input_]['tech_params']
			tech_n_params  = inputs_dict[input_]['tech_n_params']
			tech_subsample = sample[i_sample : i_sample+tech_n_params]

			i_sample += tech_n_params

			converted_tech_params = self.convert_sample_to_proper_distribution(tech_params,tech_subsample)
			np.put(amount_tech, inputs_dict[input_]['tech_params_where'], converted_tech_params)

			#Biosphere
			bio_params    = inputs_dict[input_]['bio_params']
			bio_n_params  = inputs_dict[input_]['bio_n_params']
			bio_subsample = sample[i_sample : i_sample+bio_n_params]

			i_sample += bio_n_params

			converted_bio_params = self.convert_sample_to_proper_distribution(bio_params,bio_subsample)
			np.put(amount_bio, inputs_dict[input_]['bio_params_where'], converted_bio_params)


		lca.rebuild_technosphere_matrix(amount_tech)
		lca.rebuild_biosphere_matrix(amount_bio)

		A = lca.technosphere_matrix
		B = lca.biosphere_matrix
		c = sum(lca.characterization_matrix)
		d = lca.demand_array

		score = (c*B)*spsolve(A,d) #run it before MC to factorize matrix A

		return score















