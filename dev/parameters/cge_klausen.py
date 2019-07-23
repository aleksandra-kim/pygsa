import stats_arrays as sa
import numpy as np
import klausen

# Insert parameters distribution and generate klausen instance
# NOTE: need to amend uncertainty distribution of capacity factor

parameters = {
        "gross_power_per_well": {
                "minimum": 4,
                "maximum": 16,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "average_depth_of_wells": {
                "minimum": 660,
                "maximum": 4000,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "collection_pipelines": {
                "minimum": 250,
                "maximum": 750,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "installed_capacity": {
                "minimum": 10,
                "maximum": 130,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "co2_emissions": {
                "minimum": 100.004,
                "maximum": 100.740,
                "uncertainty_type": sa.LognormalUncertainty.id,
                "loc": np.log(100.07718487920206497),
                "scale": 100.985026884879192
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
                "minimum": 0.032,
                "maximum": 0.048,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "specific_diesel_consumption": {
                "minimum": 1600,
                "maximum": 2800,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "specific_steel_consumption": {
                "minimum": 80,
                "maximum": 130,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "specific_cement_consumption": {
                "minimum": 30,
                "maximum": 50,
                "uncertainty_type": sa.UniformUncertainty.id
        },
        "specific_drilling_mud_consumption": {
                "minimum": 0.5,
                "maximum": 0.8,
                "uncertainty_type": sa.UniformUncertainty.id
        }
}

parameters = klausen.NamedParameters(parameters)

