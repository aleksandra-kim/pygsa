import brightway2 as bw
from gsa_lca import *
from pygsa.saindices.sobol_indices import sobol_indices

#Local files
from parameters.cge_func import geothermal_conventional_model
from parameters.cge_klausen import parameters
from parameters.lookup_func import *

#Set current project
bw.projects.set_current('Geothermal')

from parameters.replace import *

#Choose demand
_, _, _, _, _, _, _, _, _, _, _, _, electricity_prod_conv, _ = lookup_geothermal()
demand = {electricity_prod_conv: 1}

#Choose LCIA method
ILCD_CC = [method for method in bw.methods if "ILCD 2.0 2018 midpoint no LT" in str(method) and "climate change total" in str(method)]
ILCD_HH = [method for method in bw.methods if "ILCD 2.0 2018 midpoint no LT" in str(method) and "human health" in str(method)]
ILCD_EQ = [method for method in bw.methods if "ILCD 2.0 2018 midpoint no LT" in str(method) and "ecosystem quality" in str(method)]
ILCD = ILCD_CC + ILCD_HH + ILCD_EQ
method = ILCD[0]

#Run LCIA
lca = bw.LCA(demand,method)
lca.lci()
lca.lcia()
print('Initial LCA score: ' + str(lca.score))

#Run GSA
# inputs = ['demand_acts','geothermal energy']
inputs = []
gsa_in_lca = GSAinLCA(lca,inputs,parameters,geothermal_conventional_model)

n_dimensions = len(gsa_in_lca.parameters.data) #Number of parameters
print( '# of parameters: ' + str(len(gsa_in_lca.parameters.data)))

print( 'parameterized activities: ' +
	   '# of tech and bio params -> ' 
		+ str([ gsa_in_lca.parameters_dict['tech_n_params'], gsa_in_lca.parameters_dict['bio_n_params'] ]) )
for input_ in inputs:
	print( input_ + ': ' +
		   '# of tech and bio params -> ' 
			+ str([ gsa_in_lca.inputs_dict[input_]['tech_n_params'], gsa_in_lca.inputs_dict[input_]['bio_n_params'] ]) )

	n_dimensions += gsa_in_lca.inputs_dict[input_]['tech_n_params'] + gsa_in_lca.inputs_dict[input_]['bio_n_params']


def sobol_pandas(gsa_in_lca,first,total):

	lca = gsa_in_lca.lca
	ind_reverse_act  = 0
	ind_reverse_prod = 1
	ind_reverse_bio  = 2
	ind_db_name  = 0
	ind_act_code = 1

	sum_first = np.sum(first)
	sum_total = np.sum(total)
	normalized_first = np.array([ first[i]/sum_first for i in range(len(first)) ])
	normalized_total = np.array([ total[i]/sum_total for i in range(len(total)) ])

	# get_act_prod_bio = lambda ind_reverse, p:
	activities, products_flows = [], []
	for input_ in gsa_in_lca.inputs:

		#Technosphere activities, columns
		act_mask_tech = gsa_in_lca.inputs_dict[input_]['tech_params_where']
		act_pos_tech  = lca.tech_params[act_mask_tech]['col']
		act_tech      = [ bw.Database(lca.reverse_dict()[ind_reverse_act][p][ind_db_name]).get(lca.reverse_dict()[ind_reverse_act][p][ind_act_code])['name'] \
					      for p in act_pos_tech]

		#Products, rows
		prod_mask = gsa_in_lca.inputs_dict[input_]['tech_params_where']
		prod_pos  = lca.tech_params[prod_mask]['row']
		products  = [ bw.Database(lca.reverse_dict()[ind_reverse_prod][p][ind_db_name]).get(lca.reverse_dict()[ind_reverse_prod][p][ind_act_code])['name'] \
					  for p in prod_pos ]

		#Technosphere activities wrt biosphere flows, columns
		act_mask_bio = gsa_in_lca.inputs_dict[input_]['bio_params_where']
		act_pos_bio  = lca.bio_params[act_mask_bio]['col']
		act_bio      = [ bw.Database(lca.reverse_dict()[ind_reverse_act][p][ind_db_name]).get(lca.reverse_dict()[ind_reverse_act][p][ind_act_code])['name'] \
					      for p in act_pos_bio]

		#Biosphere flows, rows
		bio_mask  = gsa_in_lca.inputs_dict[input_]['bio_params_where']
		bio_pos   = lca.bio_params[bio_mask]['row']
		bio_flows = [ bw.Database(lca.reverse_dict()[ind_reverse_bio][p][ind_db_name]).get(lca.reverse_dict()[ind_reverse_bio][p][ind_act_code])['name'] \
					   for p in bio_pos]

		activities     += act_tech + act_bio
		products_flows += products + bio_flows
		
	data = { 'Products or flows': products_flows + list(gsa_in_lca.parameters),
			 'Activities': activities + list(gsa_in_lca.parameters),
			 'First': first, 
			 'Total': total, 
			 'Normalized first': normalized_first,
			 'Normalized total': normalized_total  }

	# Creates pandas DataFrame
	df = pd.DataFrame(data) 

	gsa_in_lca.sobol_df = df


print("Number of dimensions: "+str(n_dimensions))
n_runs = 2

first, total = sobol_indices(n_runs, n_dimensions, gsa_in_lca.model, Sampler=None)
sobol_pandas(gsa_in_lca,first,total)























  






