from parameters.cge_klausen import parameters
from parameters.cge_func import geothermal_conventional_model
from parameters.lookup_func import *


# CONVENTIONAL GEOTHERMAL
parameters.static()
params_sta_conv=geothermal_conventional_model(parameters)

#Lookup activities
_, _, _, _, _, _, _, _, collection_pipelines, _, _, _, electricity_prod_conventional, _, = lookup_geothermal()

act = bw.get_activity(electricity_prod_conventional)

if not bw.Database("geothermal energy").search(act["name"] + " zeros"):
    act.copy(name = act["name"] + " (zeros)")

# Delete all exchanges
for exc in act.exchanges():
    exc.delete()

# Insert new exchanges      
for inp in params_sta_conv:
    if inp[0][0] != "biosphere3":
        act.new_exchange(input = inp[0], amount = float(inp[2]), type= "technosphere").save()
    else:
        act.new_exchange(input = inp[0], amount = float(inp[2]), type= "biosphere").save()  