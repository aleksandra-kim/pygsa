# Set working directory
# Note that because pygsa is not installed, we need set pygsa as working directory and add dev. for scripts in dev folder.
import os
os.chdir(r'C:\Users\andre\Documents\Github\pygsa')

import brightway2 as bw
from dev.gsa_lca import GSAinLCA
from pygsa.saindices.sobol_indices import sobol_indices

#Local files
from dev.parameters.lookup_func import lookup_geothermal
from dev.parameters.cge_klausen import parameters
from dev.sobol_pandas import sobol_pandas

# Need to set current project before importing the model
bw.projects.set_current('Geothermal')

# On GitIgnore
from dev.parameters.cge_model import *

#Choose demand
_, _, _, _, _, _, _, _, _, _, _, _, _, _, electricity_prod_conv, _ = lookup_geothermal()
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
#inputs = ['demand_acts']
inputs = []
gsa_in_lca = GSAinLCA(lca,inputs,parameters,geothermal_conventional_model)

#Print some info
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


print("Number of dimensions: "+ str(n_dimensions))

#Choose number of runs
n_runs = 2

#Perform GSA
first, total = sobol_indices(n_runs, n_dimensions, gsa_in_lca.model, Sampler=None)
sobol_pandas(gsa_in_lca,first,total)

#Export to excel
gsa_in_lca.sobol_df.to_excel('dev\sobol.xlsx')






















  






