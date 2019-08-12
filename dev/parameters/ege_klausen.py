import stats_arrays as sa
import numpy as np
import klausen

# Insert parameters distribution and generate klausen instance
# NOTE: need to amend uncertainty distribution of capacity factor

parameters = {
        "number_of_wells": {
                "minimum": 2,
                "maximum": 4,  # DiscreteUniform does not include maximum
                "uncertainty_type": sa.DiscreteUniform.id
        },
        "average_depth_of_wells": {
                "minimum": 2500,
                "maximum": 6000,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "collection_pipelines": {
                "minimum": 50,
                "maximum": 200,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "installed_capacity": {
                "minimum": 0.4,
                "maximum": 11,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "lifetime": {
                "minimum": 20,
                "maximum": 40,
                "uncertainty_type": sa.NormalUncertainty.id,
                "loc": 30,
                "scale": 5
        },
        "capacity_factor": {
                "minimum": 0.85,
                "maximum": 0.95,
                "uncertainty_type": sa.UniformUncertainty.id,
        },
        "auxiliary_power": {
                "minimum": 0.12,
                "maximum": 0.28,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "specific_diesel_consumption": {
                "minimum": 3000,
                "maximum": 14000,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "specific_steel_consumption": {
                "minimum": 75,
                "maximum": 150,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "specific_cement_consumption": {
                "minimum": 16,
                "maximum": 100,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "specific_drilling_mud_consumption": {
                "minimum": 0.5,
                "maximum": 0.8,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "number_of_wells_stimulated_0to1": {
                "minimum": 0,
                "maximum": 1,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "water_stimulation": {
                "minimum": 10000 * 1000,
                "maximum": 60000 * 1000,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "specific_electricity_stimulation": {
                "minimum": 10,
                "maximum": 140,
                "uncertainty_type": sa.UniformUncertainty.id
        }
                
}

params = klausen.NamedParameters(parameters)
