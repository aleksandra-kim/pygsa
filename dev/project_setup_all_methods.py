import brightway2 as bw
from gsa_lca import *
from SALib.sample import saltelli, fast_sampler
from SALib.analyze import sobol, fast
import time

# Set current project
bw.projects.set_current('Geothermal')

# Local files
from parameters.cge_model import GeothermalConventionalModel
from parameters.cge_klausen import parameters
from parameters.lookup_func import lookup_geothermal

# Set up project parameters, model, etc
gt_model = GeothermalConventionalModel()

# Choose demand
_, _, _, _, _, _, _, _, _, _, _, _, _, _, electricity_prod_conv, _ = lookup_geothermal()
demand = {electricity_prod_conv: 1}

# Choose LCIA method
ILCD_CC = [method for method in bw.methods if "ILCD 2.0 2018 midpoint no LT" in str(method) and "climate change total" in str(method)]
ILCD_HH = [method for method in bw.methods if "ILCD 2.0 2018 midpoint no LT" in str(method) and "human health" in str(method)]
ILCD_EQ = [method for method in bw.methods if "ILCD 2.0 2018 midpoint no LT" in str(method) and "ecosystem quality" in str(method)]
ILCD = ILCD_CC + ILCD_HH + ILCD_EQ


# Number of MC runs
N = 2
inputs = []
calc_second_order = False


##############################
##### SA for all methods #####
##############################

sa_sobol = {}

for method in ILCD:

    lca = bw.LCA(demand,method)
    gsa_in_lca = GSAinLCA(lca,inputs,parameters,gt_model)

    # Count number of dimensions
    ## a) Exchanges
    n_dimensions = 0
    for input_ in inputs:
        if input_ == 'biosphere':
            continue
        n_dimensions += gsa_in_lca.inputs_dict[input_]['tech_n_params'] + gsa_in_lca.inputs_dict[input_]['bio_n_params']
    ## b) Parameters
    n_dimensions += len(gsa_in_lca.parameters.data) # Number of parameters
    
    problem = {
      'num_vars': n_dimensions,
      'names': np.arange(n_dimensions),
      'bounds': np.array([[0,1]]*n_dimensions)
    }

    # Sobol indices
    X  = saltelli.sample(problem, N, calc_second_order=calc_second_order)
    Y  = np.zeros(X.shape[0])
    for i in range(X.shape[0]):
      sample = X[i,:]
      Y[i]   = gsa_in_lca.model(sample)
    sa_sobol[method] = sobol.analyze(problem, Y, print_to_console=False, calc_second_order=calc_second_order)

#########################
##### Visualization #####
#########################
# 1. Change full names of methods to smth shorter for the plots 
methods_mapping = {
    ('ILCD 2.0 2018 midpoint no LT', 'climate change', 'climate change total'): 
        'climate change total',
    ('ILCD 2.0 2018 midpoint no LT', 'human health', 'carcinogenic effects'):
        'HH carcinogenic effects',
    ('ILCD 2.0 2018 midpoint no LT', 'human health', 'ionising radiation'):
        'HH ionising radiation',
    ('ILCD 2.0 2018 midpoint no LT', 'human health', 'non-carcinogenic effects'):
        'HH non-carcinogenic effects',
    ('ILCD 2.0 2018 midpoint no LT', 'human health', 'ozone layer depletion'):
        'HH ozone layer depletion',
    ('ILCD 2.0 2018 midpoint no LT', 'human health', 'photochemical ozone creation'):
        'HH photochemical ozone creation',
    ('ILCD 2.0 2018 midpoint no LT', 'human health', 'respiratory effects, inorganics'):
        'HH respiratory effects, inorganics',
    ('ILCD 2.0 2018 midpoint no LT', 'ecosystem quality', 'freshwater and terrestrial acidification'):
        'EQ freshwater and terrestrial acidification',
    ('ILCD 2.0 2018 midpoint no LT', 'ecosystem quality', 'freshwater ecotoxicity'):
        'EQ freshwater ecotoxicity',
    ('ILCD 2.0 2018 midpoint no LT', 'ecosystem quality', 'freshwater eutrophication'):
        'EQ freshwater eutrophication',
    ('ILCD 2.0 2018 midpoint no LT', 'ecosystem quality', 'marine eutrophication'):
        'EQ marine eutrophication',
    ('ILCD 2.0 2018 midpoint no LT', 'ecosystem quality', 'terrestrial eutrophication'):
        'EQ terrestrial eutrophication'
}

names = gsa_in_lca.sensitivity_indices_df['Activities'].values #TODO this is specific to our setup

# 2. Extract first and total order indices
def extract_indices(sa_dict, names, sa_index, methods_mapping):
    df = pd.DataFrame([], columns=list(methods_mapping.values()), index=names)
    for k, v in sa_dict.items():
        method_short = methods_mapping[k]
        df[method_short] = v[sa_index]
    return df

df_first = extract_indices(sa_sobol, names, 'S1', methods_mapping)
df_total = extract_indices(sa_sobol, names, 'ST', methods_mapping)

# 2. Convert indices to rankings
def sa_to_rankings(sa_df):
    df = pd.DataFrame([], index=sa_df.index, columns=sa_df.columns)
    for col in sa_df.columns:
        df[col] = np.argsort(np.argsort(sa_df[col]))
    return df

df_first_rankings = sa_to_rankings(df_first)
df_total_rankings = sa_to_rankings(df_total)


# 3. Plot results for all methods
import plotly.graph_objs as go

df = df_first_rankings

dimensions = [0]*df.shape[1]
i = 0
for col in df.columns:
    dimensions[i] = dict( label = col,
                          values = df[col].values)
    i += 1

data = [
    go.Parcoords(
        line = dict(colorscale = 'Jet',
                    showscale = True,
                    reversescale = True,
                    cmin = -4000,
                    cmax = -100),
        dimensions = dimensions,
        labelangle = 90              
    )
]

fig = go.Figure(data)
fig.write_image('par_coor.png')












