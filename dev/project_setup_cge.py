import brightway2 as bw
from gsa_lca import *
from pygsa.saindices.sobol_indices import sobol_indices
from SALib.sample import saltelli, fast_sampler
from SALib.analyze import sobol, fast
import time

#Set current project
bw.projects.set_current('Geothermal')

#Local files
from parameters.cge_model import GeothermalConventionalModel
from parameters.cge_klausen import parameters
from parameters.lookup_func import lookup_geothermal
from exact import exact_first_total

#Set up project parameters, model, etc
gt_model = GeothermalConventionalModel()

#Choose demand
_, _, _, _, _, _, _, _, _, _, _, _, _, _, electricity_prod_conv, _ = lookup_geothermal()
demand = {electricity_prod_conv: 1}

#Choose LCIA method
ILCD_CC = [method for method in bw.methods if "ILCD 2.0 2018 midpoint no LT" in str(method) and "climate change total" in str(method)]
ILCD_HH = [method for method in bw.methods if "ILCD 2.0 2018 midpoint no LT" in str(method) and "human health" in str(method)]
ILCD_EQ = [method for method in bw.methods if "ILCD 2.0 2018 midpoint no LT" in str(method) and "ecosystem quality" in str(method)]
ILCD_RE = [method for method in bw.methods if "ILCD 2.0 2018 midpoint no LT" in str(method) and "resources" in str(method)]
ILCD = ILCD_CC + ILCD_HH + ILCD_EQ + ILCD_RE
method = ILCD[0]

# Run LCIA
lca = bw.LCA(demand,method)

# Run GSA
inputs = []
gsa_in_lca = GSAinLCA(lca,inputs,parameters,gt_model)

# Activities in the gsa_in_lca.inputs_dict
n_dimensions = 0
for input_ in inputs:

    if input_ == 'biosphere':
        continue

    print( input_ + ': ' +
           '# of tech and bio params -> ' 
            + str([ gsa_in_lca.inputs_dict[input_]['tech_n_params'], gsa_in_lca.inputs_dict[input_]['bio_n_params'] ]) )
    n_dimensions += gsa_in_lca.inputs_dict[input_]['tech_n_params'] + gsa_in_lca.inputs_dict[input_]['bio_n_params']

# Parameters
n_dimensions += len(gsa_in_lca.parameters.data) # Number of parameters
print( '# of parameters: ' + str(len(gsa_in_lca.parameters.data)))
print( 'parameterized activities: ' +
       '# of tech and bio params -> ' 
        + str([ gsa_in_lca.parameters_dict['tech_n_params'], gsa_in_lca.parameters_dict['bio_n_params'] ]) )

################################
##### SALib implementation #####
################################
N = 50
calc_second_order = False
problem = {
  'num_vars': n_dimensions,
  'names': np.arange(n_dimensions),
  'bounds': np.array([[0,1]]*n_dimensions)
}

# Sobol indices
t0 = time.time()
X  = saltelli.sample(problem, N, calc_second_order=calc_second_order)
t1 = time.time()
print(str(t1-t0) + ' sec for Sobol sampling')

t0 = time.time()
Y  = np.zeros(X.shape[0])
for i in range(X.shape[0]):
  sample = X[i,:]
  Y[i]   = gsa_in_lca.model(sample)
t1 = time.time()
print(str(t1-t0) + ' sec for running the LCA model with Sobol sampling')

t0 = time.time()
sa_sobol = sobol.analyze(problem, Y, print_to_console=False, calc_second_order=calc_second_order)
t1 = time.time()
print(str(t1-t0) + ' sec for Sobol analyze')

gsa_in_lca.sa_pandas_append(sa_sobol)



###########################
###### Exact formulas #####
###########################
t0 = time.time()
sa_exact_formula = exact_first_total(problem, N, gsa_in_lca.model)
t1 = time.time()

gsa_in_lca.sa_pandas_append(sa_exact_formula)
print(str(t1-t0) + ' sec for exact formulas')

## Sort results and save dataframe as .xlsx
#gsa_in_lca.sensitivity_indices_df = gsa_in_lca.sensitivity_indices_df.sort_values(['S1'], ascending = False)
#column_order = ['Inputs', 'Activities', 'Products or flows', 'S1', 'S1 exact', 'ST', 'ST exact', 'S1_conf', 'ST_conf']
#gsa_in_lca.sensitivity_indices_df = gsa_in_lca.sensitivity_indices_df.reindex(column_order, axis=1)


name= "cge_sa_indices_N" + str(N) + ".xlsx"
gsa_in_lca.sensitivity_indices_df.to_excel(name)




# # 2. FAST - minimum number of runs is 70 so far, at 70 runs very slowly!
# N = 70

# t0 = time.time()
# X = fast_sampler.sample(problem, N)
# t1 = time.time()
# print(str(t1-t0) + ' sec for FAST sampling')

# Y = np.zeros(X.shape[0])
# for i in range(X.shape[0]):
#   sample = X[i,:]
#   Y[i] = gsa_in_lca.model(sample)

# t0 = time.time()
# sa_fast = fast.analyze(problem, Y, print_to_console=False)
# t1 = time.time()
# print(str(t1-t0) + ' sec for FAST analyze')

# gsa_in_lca.sa_pandas_append(sa_fast)


# gsa_in_lca.sensitivity_indices_df.to_excel('sa_indices.xlsx')




# def sobol_pandas(gsa_in_lca,first,total):

#     lca = gsa_in_lca.lca
#     ind_reverse_act  = 0
#     ind_reverse_prod = 1
#     ind_reverse_bio  = 2
#     ind_db_name  = 0
#     ind_act_code = 1

#     sum_first = np.sum(first)
#     sum_total = np.sum(total)
#     normalized_first = np.array([ first[i]/sum_first for i in range(len(first)) ])
#     normalized_total = np.array([ total[i]/sum_total for i in range(len(total)) ])

#     # get_act_prod_bio = lambda ind_reverse, p:
#     activities, products_flows = [], []
#     for input_ in gsa_in_lca.inputs:

#         #Technosphere activities, columns
#         act_mask_tech = gsa_in_lca.inputs_dict[input_]['tech_params_where']
#         act_pos_tech  = lca.tech_params[act_mask_tech]['col']
#         act_tech      = [ bw.Database(lca.reverse_dict()[ind_reverse_act][p][ind_db_name]).get(lca.reverse_dict()[ind_reverse_act][p][ind_act_code])['name'] \
#                           for p in act_pos_tech]

#         #Products, rows
#         prod_mask = gsa_in_lca.inputs_dict[input_]['tech_params_where']
#         prod_pos  = lca.tech_params[prod_mask]['row']
#         products  = [ bw.Database(lca.reverse_dict()[ind_reverse_prod][p][ind_db_name]).get(lca.reverse_dict()[ind_reverse_prod][p][ind_act_code])['name'] \
#                       for p in prod_pos ]

#         #Technosphere activities wrt biosphere flows, columns
#         act_mask_bio = gsa_in_lca.inputs_dict[input_]['bio_params_where']
#         act_pos_bio  = lca.bio_params[act_mask_bio]['col']
#         act_bio      = [ bw.Database(lca.reverse_dict()[ind_reverse_act][p][ind_db_name]).get(lca.reverse_dict()[ind_reverse_act][p][ind_act_code])['name'] \
#                           for p in act_pos_bio]

#         #Biosphere flows, rows
#         bio_mask  = gsa_in_lca.inputs_dict[input_]['bio_params_where']
#         bio_pos   = lca.bio_params[bio_mask]['row']
#         bio_flows = [ bw.Database(lca.reverse_dict()[ind_reverse_bio][p][ind_db_name]).get(lca.reverse_dict()[ind_reverse_bio][p][ind_act_code])['name'] \
#                        for p in bio_pos]

#         activities     += act_tech + act_bio
#         products_flows += products + bio_flows

#     parameters_list = gsa_in_lca.parameters_array['name'].tolist()
        
#     data = { 'Products or flows': products_flows + parameters_list,
#              'Activities': activities + parameters_list,
#              'First': first, 
#              'Total': total, 
#              'Normalized first': normalized_first,
#              'Normalized total': normalized_total  }

#     # Creates pandas DataFrame
#     df = pd.DataFrame(data) 

#     gsa_in_lca.sobol_df = df


# print("Number of dimensions: "+str(n_dimensions))
# n_runs = 2

# first, total = sobol_indices(n_runs, n_dimensions, gsa_in_lca.model, Sampler=None)
# sobol_pandas(gsa_in_lca,first,total)


# gsa_in_lca.sobol_df.to_excel('sobol.xlsx')






















  






