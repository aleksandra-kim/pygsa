import pandas as pd
import numpy as np
import brightway2 as bw

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
		
	data = { 'Products or flows': products_flows + sorted(list(gsa_in_lca.parameters)),
			 'Activities': activities + sorted(list(gsa_in_lca.parameters)),
			 'First': first, 
			 'Total': total, 
			 'Normalized first': normalized_first,
			 'Normalized total': normalized_total  }

	# Creates pandas DataFrame
	df = pd.DataFrame(data) 

	gsa_in_lca.sobol_df = df