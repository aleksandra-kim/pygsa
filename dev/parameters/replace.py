import brightway2 as bw
from .lookup_func import *

def replace(parameters,gt_model):

	# CONVENTIONAL GEOTHERMAL
	parameters.static()
	gt_model.run(parameters)
	params_sta_conv = gt_model.array_io

	#Lookup activities
	_, _, _, _, _, _, _, _, _, _, _, _, _, _, electricity_prod_conventional, _, = lookup_geothermal()

	act = bw.get_activity(electricity_prod_conventional)

	if not bw.Database("geothermal energy").search(act["name"] + " zeros"):
		act.copy(name = act["name"] + " (zeros)")

	# Delete all exchanges
	for exc in act.exchanges():
	    exc.delete()

	# Insert new exchanges      
	for inp in params_sta_conv:
		if inp['input_db'] != "biosphere3":
			print(inp)
			# act.new_exchange(input = (inp['input_db'],inp['input_code']), amount = float(inp['amount']), type= "technosphere").save()
		else:
			print( type( tuple ( (str(inp['input_db']),str(inp['input_code'])) ) ))
			print(float(inp['amount']))
			# act.new_exchange(input = tuple(inp['input_db'],inp['input_code']), amount = float(inp['amount']), type= "biosphere").save() 