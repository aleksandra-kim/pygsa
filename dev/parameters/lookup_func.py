# Lookup activities function
import brightway2 as bw

# Conventional
def lookup_geothermal():

    db_geothe = bw.Database("geothermal energy")
    db_ecoinv = bw.Database("ecoinvent 3.5 cutoff")
    db_biosph = bw.Database("biosphere3")
    
    wellhead      = db_geothe.search("geothermal wellhead")[0].key
    diesel        = db_ecoinv.search("market diesel, burned diesel-electric generating set")[1].key 
    steel         = db_ecoinv.search("market steel low alloyed")[0].key
    cement        = db_ecoinv.search("market cement portland", filter={"location": "ROW"})[0].key
    water         = db_ecoinv.search("market tap water", filter={"location": "ROW"})[0].key
    drilling_mud  = db_geothe.search("drilling mud")[0].key
    drill_wst     = db_ecoinv.search("market drilling waste")[0].key
    wells_closr   = db_ecoinv.search("deep well closure")[0].key
    coll_pipe     = db_geothe.search("collection pipelines")[0].key
    plant         = db_geothe.search("geothermal plant, double flash (electricity)")[0].key
    ORC_fluid     = db_ecoinv.search("market perfluoropentane")[0].key
    ORC_fluid_wst = db_ecoinv.search("treatment perfluoropentane")[0].key
    diesel_stim   = db_ecoinv.search("market diesel, burned diesel-electric generating set")[0].key
    co2           = db_biosph.search("Carbon dioxide, fossil")[0].key
    electricity_prod_conventional = db_geothe.search("electricity production, geothermal, conventional")[0].key
    electricity_prod_enhanced     = db_geothe.search("electricity production, geothermal, enhanced")[0].key
    
    return wellhead,     \
           diesel,       \
           steel,        \
           cement,       \
           water,        \
           drilling_mud, \
           drill_wst,    \
           wells_closr,  \
           coll_pipe,    \
           plant,        \
           diesel_stim,  \
           co2,          \
           electricity_prod_conventional, \
           electricity_prod_enhanced


"""
S: does search always return activities in the same order?
"""