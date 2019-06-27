import brightway2 as bw
from gsa_lca import *
from pygsa.saindices.sobol_indices import sobol_indices

###1. LCA model
# bw.projects.set_current('Sobol indices')
bw.projects.set_current('Thrive')

method = ('IPCC 2013', 'climate change', 'GWP 100a')
ei = bw.Database('cutoff35')
bs = bw.Database('biosphere3')
demand_act = [act for act in ei if 'market for electricity, high voltage'==act['name'] and 'CH' in act['location']][0]
demand_amt = 1
demand = {demand_act: demand_amt}

lca = bw.LCA(demand,method)
lca.lci()
lca.lcia()

inputs = ['demand_acts']
gsa_in_lca = GSAinLCA(lca,inputs)

n_dimensions = 0
for input_ in inputs:
	print( [ gsa_in_lca.inputs_dict[input_]['tech_n_params'], gsa_in_lca.inputs_dict[input_]['bio_n_params'] ])
	n_dimensions += gsa_in_lca.inputs_dict[input_]['tech_n_params'] + gsa_in_lca.inputs_dict[input_]['bio_n_params']

n_runs = 2

first, total = sobol_indices(n_runs, n_dimensions, gsa_in_lca.model, Sampler=None)
print([first, total])




