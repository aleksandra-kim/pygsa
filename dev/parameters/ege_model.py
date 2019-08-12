import numpy as np

# REMEMBER TO ADD . IF USING PYGSA AND REMOVE IF NOT
from .lookup_func import lookup_geothermal


class GeothermalEnhancedModel:
    
    def __init__(self, exploration = True):
        
        # Init constants
        self.drilling_waste_per_metre = 450 # kilogram (for open hole diameter of 8.5 in and assume factor 3 production line to total volume drilled)

        
        self.wellhead, self.diesel, self.steel, self.cement, self.water, \
        self.drilling_mud, self.drill_wst, self.wells_closr, self.coll_pipe, \
        self.plant, self.ORC_fluid, self.ORC_fluid_wst, self.diesel_stim, self.co2, \
        _, self.electricity_prod = lookup_geothermal()

        # Make-up wells 
        
        # We are assuming that:
        # there is a need for 3 exploratory wells
        # these are equivalent to 0.3 normal wells.
        # exploratory wells have same average depth and same need of materials as normal wells 
        # THE IMPACT OF THIS NEEDS TO BE INVESTIGATED

        if exploration:
            number_of_expl_wells = 3
        else:
            number_of_expl_wells = 0

        equivalency_expl_to_normal_wells = 0.3
        self.number_of_eq_expl_wells = number_of_expl_wells * equivalency_expl_to_normal_wells

        self.input_output()

    def input_output(self):

        #TODO careful with the hardcoded values in the string length '<U20'
        dtype_io = np.dtype([ ('input_db', '<U20'),
                              ('input_code','<U40'),
                              ('output_db', '<U20'),
                              ('output_code', '<U40'),
                              ('amount', '<f4') ])

        input_  = np.array([self.wellhead, self.diesel, self.steel, self.cement, self.water, self.drilling_mud,
                  self.drill_wst, self.wells_closr, self.coll_pipe, self.plant, self.ORC_fluid, self.ORC_fluid_wst,
                  self.diesel_stim])
        output_ = self.electricity_prod
        amounts = np.zeros(len(input_),dtype=float)

        array_io = np.empty(len(input_),dtype=dtype_io)

        array_io['input_db']    = input_[:,0]
        array_io['input_code']  = input_[:,1]
        array_io['output_db']   = output_[0]
        array_io['output_code'] = output_[1]
        array_io['amount']      = amounts

        self.array_io = array_io
        
    def run(self, params):
           
        lifetime_electricity_generated = (params["installed_capacity"] *
                                          params["capacity_factor"] *
                                          (1 - params["auxiliary_power"]) *
                                          params["lifetime"] * 8760000) # kilowatt hour
        
        total_metres_drilled = (params["number_of_wells"] * params["average_depth_of_wells"]) # metres
        
        diesel_consumption = (total_metres_drilled * 
                              params["specific_diesel_consumption"]) # MJ of thermal energy as diesel burned
        
        steel_consumption = (total_metres_drilled * 
                             params["specific_steel_consumption"])  # kilogram
        
        cement_consumption = (total_metres_drilled * 
                              params["specific_cement_consumption"]) # kilogram
        
        water_cem_consumption = (cement_consumption * (1/0.65))  # kilogram
        
        drilling_mud_consumption = (total_metres_drilled *
                                    params["specific_drilling_mud_consumption"])  # cubic meter
        
        drilling_waste = - (total_metres_drilled * self.drilling_waste_per_metre) # kilogram (minus because waste)
        
        total_collection_pipelines = (params["number_of_wells"] * params["collection_pipelines"])  # metres
        
        # This parameter generates a value that is between 1 and max of number_of_wells. 
        # However, NOTE that if random number 0to1 is 0, then this results in number of wells stimulated equal to 0.
        number_of_wells_stim = np.round ( 
                0.5 + # So that at leat one well is stimulatd
                params["number_of_wells_stimulated_0to1"] * params["number_of_wells"])  
                                                                                              
        water_stim_consumption = number_of_wells_stim * params["water_stimulation"]
        
        total_water_consumption = water_cem_consumption + water_stim_consumption
        
        diesel_for_stim = ((water_stim_consumption / 1000) * params["specific_electricity_stimulation"]) / 0.3 # To convert from electricity to thermal energy
        
        ORC_fluid_consumption = 300 * params["installed_capacity"]
       
        # Amounts per kwh generated
        amounts = np.array([ params["number_of_wells"], 
                             diesel_consumption, 
                             steel_consumption, 
                             cement_consumption,
                             total_water_consumption, 
                             drilling_mud_consumption, 
                             drilling_waste, 
                             total_metres_drilled,
                             total_collection_pipelines, 
                             params["installed_capacity"],
                             ORC_fluid_consumption,
                             ORC_fluid_consumption,
                             diesel_for_stim]) / lifetime_electricity_generated 
        
        self.array_io['amount'] = amounts

        return self.array_io
    
    
    def run_ps(self, params):
        
        lifetime_electricity_generated = (params["installed_capacity"] *
                                          params["capacity_factor"] *
                                          (1 - params["auxiliary_power"]) *
                                          params["lifetime"] * 8760000) # kilowatt hour
        
        total_metres_drilled = (params["number_of_wells"] * params["average_depth_of_wells"]) # metres
        
        diesel_consumption = (total_metres_drilled * 
                              params["specific_diesel_consumption"]) # MJ of thermal energy as diesel burned
        
        steel_consumption = (total_metres_drilled * 
                             params["specific_steel_consumption"])  # kilogram
        
        cement_consumption = (total_metres_drilled * 
                              params["specific_cement_consumption"]) # kilogram
        
        water_cem_consumption = (cement_consumption * (1/0.65))  # kilogram
        
        drilling_mud_consumption = (total_metres_drilled *
                                    params["specific_drilling_mud_consumption"])  # cubic meter
        
        drilling_waste = - (total_metres_drilled * drilling_waste_per_metre) # kilogram (minus because waste)
        
        total_collection_pipelines = (params["number_of_wells"] * params["collection_pipelines"])  # metres
        
        # This parameter generates a value that is between 1 and max of number_of_wells. 
        # However, NOTE that if random number 0to1 is 0, then this results in number of wells stimulated equal to 0.
        number_of_wells_stim = np.round ( 
                0.5 + # So that at leat one well is stimulatd
                params["number_of_wells_stimulated_0to1"] * params["number_of_wells"])  
                                                                                              
        water_stim_consumption = number_of_wells_stim * params["water_stimulation"]
        
        total_water_consumption = water_cem_consumption + water_stim_consumption
        
        diesel_for_stim = ((water_stim / 1000) * params["specific_electricity_stimulation"]) / 0.3 # To convert from electricity to thermal energy
        
        ORC_fluid_consumption = 300 * params["installed_capacity"]
              
        #TODO Need to find a way to include array in column "amount" in the structured array "array_io"
        return [
            (self.wellhead, self.electricity_prod, np.array(number_of_wells/lifetime_electricity_generated)),
            (self.diesel, self.electricity_prod, np.array(diesel_consumption/lifetime_electricity_generated)),
            (self.steel, self.electricity_prod, np.array(steel_consumption/lifetime_electricity_generated)),
            (self.cement, self.electricity_prod, np.array(cement_consumption/lifetime_electricity_generated)),
            (self.water, self.electricity_prod, np.array(total_water_consumption/lifetime_electricity_generated)),
            (self.drilling_mud, self.electricity_prod, np.array(drilling_mud/lifetime_electricity_generated)),
            (self.drill_wst, self.electricity_prod , np.array(drilling_waste/lifetime_electricity_generated)),
            (self.wells_closr, self.electricity_prod, np.array(total_metres_drilled/lifetime_electricity_generated)), # well closure
            (self.coll_pipe, self.electricity_prod, np.array(total_collection_pipelines/lifetime_electricity_generated)),
            (self.plant, self.electricity_prod, np.array(params["installed_capacity"]/lifetime_electricity_generated)),
            (self.ORC_fluid, self.electricity_prod, np.array(ORC_fluid_consumption/lifetime_electricity_generated)),
            (self.ORC_fluid_wst, self.electricity_prod, np.array(ORC_fluid_consumption/lifetime_electricity_generated)),
            (self.diesel_stim, self.electricity_prod, np.array(diesel_stim/lifetime_electricity_generated))
    ] 