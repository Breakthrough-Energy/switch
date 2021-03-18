# -*- coding: utf-8 -*-

"""
Created on Wed Mar 16 12:27:29 2021

@author: YifanLi @ Breakthrough Energy Science
"""

# This code transfers generation and transmission data from PowerSimData scenario grid to Switch2.0 input format.
# Generate data inputs needed for all basic modules of Switch2 expansion optimization model in .csv format.
# No geographical reduction will be performed, all data maintains original nodal model resolution.
# Refer to TemporalSlicing.py for time series related input data development.
# Investment cost data from NREL ATB.


import os, sys, csv
import pandas as pd
import numpy as np
import math as m

# user_defined inputs.

base_year = input("Please enter base study year (normally PowerSimData scenario year): ")
inv_period = input("Please enter investment period year, separate by space: ")

inv_period = inv_period.split()
inv_stage = len(inv_period)

while inv_stage < 1:
    print("Invalid investment period, please re-enter.")
    inv_period = input("Please enter investment period year, separate by space: ")
    inv_period = inv_period.split()
    inv_stage = len(inv_period)

if inv_stage == 1:
    print("Single stage expansion identified.")
else:
    print("Multi stage expansion identified.")
    
period_start = input("Please enter start year for each period, separate by space: ")
period_end = input("Please enter end year for each period, separate by space: ")

period_start = period_start.split()
period_end= period_end.split()
    
    
# Pre-defined financial parameters.

discount_rate = 0.079
interest_rate = 0.029

# load in the PowerSimData scenario user specified.

scenarioID = input("Please enter the scenario ID: ")

from powersimdata.scenario.scenario import Scenario
scenario = Scenario(scenarioID)
grid = scenario.state.get_grid()

# Retired unit list, units with 0 Pmax.

ret_unit_list = (grid.plant['Pmax'][grid.plant['Pmax'] == 0]).index

# Fuel (suggest to separate type and fuel in PowerSimData later.):

type_psd = grid.plant.type.drop_duplicates().tolist()
fuel = ['Coal', 'NaturalGas', 'Uranium']
non_fuel = ['Wind', 'Solar', 'Water', 'Geothermal']
gen_invcost = [5983000, 1956000, 1616000, 1600000, 7186000, 5255000, 7696000, 5342000]
dic_fuel = {}

for i in range(len(type_psd)):
    if type_psd[i] == 'wind' or type_psd[i] == 'wind_offshore':
        dic_fuel[type_psd[i]] = non_fuel[0]
    elif type_psd[i] == 'solar':
        dic_fuel[type_psd[i]] = non_fuel[1]
    elif type_psd[i] == 'hydro':
        dic_fuel[type_psd[i]] = non_fuel[2]
    elif type_psd[i] == 'geothermal':
        dic_fuel[type_psd[i]] = non_fuel[3]
    elif type_psd[i] == 'coal':
        dic_fuel[type_psd[i]] = fuel[0]
    elif type_psd[i] == 'ng':
        dic_fuel[type_psd[i]] = fuel[1]
    elif type_psd[i] == 'nuclear':
        dic_fuel[type_psd[i]] = fuel[2]

dic_gen_cap_lim = {}

for i in range(len(type_psd)):
    if type_psd[i] == 'coal':
        dic_gen_cap_lim[type_psd[i]] = float(0)
    else:
        dic_gen_cap_lim[type_psd[i]] = float(5000)

dic_gen_age = {}

for i in range(len(type_psd)):
    if type_psd[i] == 'hydro':
        dic_gen_age[type_psd[i]] = 60
    elif type_psd[i] == 'coal' or type_psd[i] == 'nuclear':
        dic_gen_age[type_psd[i]] = 40
    else:
        dic_gen_age[type_psd[i]] = 20

dic_variable = {}

for i in range(len(type_psd)):
    if type_psd[i] == 'hydro' or 'wind' or 'solar' or 'wind_offshore':
        dic_variable[type_psd[i]] = 1
    else:
        dic_variable[type_psd[i]] = 0

dic_baseload = {}
        
for i in range(len(type_psd)):
    if type_psd[i] == 'coal' or 'nuclear':
        dic_baseload[type_psd[i]] = 1
    else:
        dic_baseload[type_psd[i]] = 0
        
dic_geninvcost = {}

for i in range(len(type_psd)):
    dic_geninvcost[type_psd[i]] = gen_invcost[i]
        
        
# Breakdown generation cost from grid data into FVOM and fuel costs.
# Need to update if generation cost modeling has been improved later in PowerSimData.
# Assume 70% of variable cost is from fuel cost.
# If there are different cost for same fuel type at the same bus (not common), get average. 

bus_type_fuelcost = grid.plant[["bus_id", "type", "GenFuelCost"]]
bus_type_fuelcost['fuel'] = bus_type_fuelcost['type'].map(dic_fuel)

df_ave_fuel_cost = pd.pivot_table(bus_type_fuelcost, values = 'GenFuelCost', 
                   index=['bus_id'], columns = ['fuel'], aggfunc = np.mean, fill_value = 0.0) 

gencost_splitter = 0.7

df_gen_cost_ratio = (grid.gencost['after']['f2'].sub(grid.gencost['after']['f1'])).divide(
    grid.gencost['after']['p2'].sub(grid.gencost['after']['p1']), fill_value = 0.0)
df_gen_cost_ratio = df_gen_cost_ratio.fillna(value = 0.0)

df_gen_fix_cost = (grid.gencost['after']['f1']).divide(grid.plant['Pmax'], fill_value = 0.0)
df_gen_fix_cost = df_gen_fix_cost.multiply(8760, fill_value = 0.0)
df_gen_fix_cost = df_gen_fix_cost.fillna(value = 0.0)

df_gen_fuel_related = df_gen_cost_ratio.multiply(0.7)
df_gen_vom = df_gen_cost_ratio.multiply(0.3)

bus_type_fuelcost['AveFuelRelated'] = df_gen_fuel_related
bus_type_fuelcost['BusAveFuelCost'] = ""

for i in range(len(bus_type_fuelcost.index)):
    temp_bus_id = bus_type_fuelcost.loc[bus_type_fuelcost.index[i]].at['bus_id']
    temp_fuel = bus_type_fuelcost.loc[bus_type_fuelcost.index[i]].at['fuel']
    bus_type_fuelcost.at[bus_type_fuelcost.index[i], 
                         'BusAveFuelCost'] = df_ave_fuel_cost.loc[temp_bus_id].at[temp_fuel]
    
bus_type_fuelcost['Calc_HeatRate'] = bus_type_fuelcost['AveFuelRelated'].divide(
    bus_type_fuelcost['BusAveFuelCost'], fill_value = 0.0)
bus_type_fuelcost = bus_type_fuelcost.fillna(value = 0.0)


def AveFuelCost(bus_id, fuel, year):
    
    """Calculate average fuel cost by type on each bus.
       :param int64 bus_id: bus number ID.
       :param str fuel: fuel type.
       :param str year: investment year.
       :return: (*float*) -- average fuel cost for a certain fuel type on each bus.
    """
    
    ave_fuel_cost = df_ave_fuel_cost.at[bus_id, fuel]
    ave_fuel_cost = ave_fuel_cost * (interest_rate**(int(year) - int(base_year)))
    ave_fuel_cost = round(ave_fuel_cost,2)
    return ave_fuel_cost

def CalcHeatRate(plant_id):
    
    """Calculate heat rate for a certain plant.
       :param int64 plant_id: plant ID.
       :return: (*float*) -- heat rate for a plant.
    """    
    
    calc_heat_rate = bus_type_fuelcost.at[plant_id, 'Calc_HeatRate']
    calc_heat_rate = round(calc_heat_rate,2)
    if plant_id in ret_unit_list:
        calc_heat_rate = 999
    if calc_heat_rate == 0.0:
        calc_heat_rate = '.'
    return calc_heat_rate

def GIS_distance(branch_id):
    
    """Calculate length for a given branch ID.
       :param int64 branch_id: branch ID.
       :return: (*float*) -- straight line distance for a branch.
    """
    
    pi = float (3.14159265354)
    r = float (3958.756)
    e = float (1.60934)
    lat1 = grid.branch.at[branch_id, 'from_lat']
    long1 = grid.branch.at[branch_id, 'from_lon']
    lat2 = grid.branch.at[branch_id, 'to_lat']
    long2 = grid.branch.at[branch_id, 'to_lon']
    if lat1 == lat2 and long1 == long2:
        distance = 0.001
    else:
        distance = m.acos (m.sin (lat1 / 180 * pi) * m.sin (lat2 / 180 * pi) + 
                           m.cos (lat1 / 180 * pi) * m.cos (lat2 / 180 * pi) * 
                           m.cos (long1 / 180 * pi - long2 / 180 * pi) ) * r * e
    distance = round(distance, 2)
    return distance

def GIS_distance2(from_bus_id, to_bus_id):
    
    """Calculate distance based on two given busses.
       :param int64 from_bus_id: from bus number ID.
       :param int64 to_bus_id: to bus number ID.
       :return: (*float*) -- straight line distance connecting from and to busses.
    """
    
    pi = float (3.14159265354)
    r = float (3958.756)
    e = float (1.60934)
    lat1 = grid.bus.at[from_bus_id, 'lat']
    long1 = grid.bus.at[from_bus_id, 'lon']
    lat2 = grid.bus.at[to_bus_id, 'lat']
    long2 = grid.bus.at[to_bus_id, 'lon']
    if lat1 == lat2 and long1 == long2:
        distance = 0.001
    else:
        distance = m.acos (m.sin (lat1 / 180 * pi) * m.sin (lat2 / 180 * pi) + 
                           m.cos (lat1 / 180 * pi) * m.cos (lat2 / 180 * pi) * 
                           m.cos (long1 / 180 * pi - long2 / 180 * pi) ) * r * e
    distance = round(distance, 2)
    return distance

def Branch_efficiency(branch_id):
    
    """Identify power efficiency for a given branch ID.
       :param int64 branch_id: branch ID.
       :return: (*float*) -- efficiency rate for a branch.
    """
    
    from_kv = grid.bus.at[(grid.branch.at[branch_id, 'from_bus_id']), 'baseKV']
    to_kv = grid.bus.at[(grid.branch.at[branch_id, 'to_bus_id']), 'baseKV']
    if from_kv == to_kv:
        if from_kv == 115:
            branch_efficiency = 0.9
        elif from_kv == 138:
            branch_efficiency = 0.94
        elif from_kv == 161:
            branch_efficiency = 0.96
        elif from_kv == 230:
            branch_efficiency = 0.97
        elif from_kv == 345:
            branch_efficiency = 0.98
        elif from_kv == 500 or 765:
            branch_efficiency = 0.99
    else:
        branch_efficiency = 0.99
    return branch_efficiency


# write SWITCH2.0 input files.
# write financials.csv file.

with open('financials.csv', 'w', newline='') as financials_file:
    writer = csv.writer(financials_file)
    writer.writerow(["base_financial_year", "discount_rate", "interest_rate"])
    writer.writerow([base_year, discount_rate, interest_rate])

# write fuels.csv file.

with open('fuels.csv', 'w', newline='') as fuels_file:
    writer = csv.writer(fuels_file)
    writer.writerow(["fuel", "co2_intensity", "upstream_co2_intensity"])
    
    for i in range(len(fuel)):
        writer.writerow([fuel[i], ".", "."])

# write non_fuel_energy_source.csv file.

with open('non_fuel_energy_source.csv', 'w', newline='') as non_fuel_energy_source_file:
    writer = csv.writer(non_fuel_energy_source_file)
    writer.writerow(["energy_source"])
    
    for i in range(len(non_fuel)):
        writer.writerow([non_fuel[i]])

# write fuel_cost.csv file. 
   
with open('fuel_cost.csv', 'w', newline='') as fuel_cost_file:
    writer = csv.writer(fuel_cost_file)
    writer.writerow(["load_zone", "fuel", "period", "fuel_cost"])
    
    for i in range(inv_stage):
        for j in range(len(df_ave_fuel_cost.index)):
            for k in range(len(fuel)):
                
# Make sure no zero fuel costs.
                
                if AveFuelCost(df_ave_fuel_cost.index[j], fuel[k], inv_period[i]) == 0:
                    pass
                else:
                    writer.writerow([df_ave_fuel_cost.index[j], fuel[k], inv_period[i], 
                                    AveFuelCost(df_ave_fuel_cost.index[j], 
                                    fuel[k], inv_period[i])])                                   
                                    
# write load_zones.csv file.
# assume no local transmission and distributed transmission consideration.

with open('load_zones.csv', 'w', newline='') as load_zone_file:
    writer = csv.writer(load_zone_file)
    writer.writerow(["LOAD_ZONE", "dbid", "existing_local_td", "local_td_annual_cost_per_mw"])
    
    for i in range(len(grid.bus.index)):
        writer.writerow([grid.bus.index[i], i+1, 99999, 0])

# write generation_projects_info.csv file.
# assume no generation interconnection local network upgrades.

with open('generation_projects_info.csv', 'w', newline='') as gen_proj_info_file:
    writer = csv.writer(gen_proj_info_file)
    writer.writerow(["GENERATION_PROJECT", "gen_tech", "gen_load_zone", 
                     "gen_connect_cost_per_mw", "gen_capacity_limit_mw", 
                     "gen_full_load_heat_rate", "gen_variable_om", "gen_max_age", 
                     "gen_min_build_capacity", "gen_scheduled_outage_rate", 
                     "gen_forced_outage_rate", "gen_is_variable", 
                     "gen_is_baseload", "gen_is_cogen", "gen_energy_source", 
                     "gen_unit_size", "gen_ccs_capture_efficiency", 
                     "gen_ccs_energy_load", "gen_storage_efficiency", 
                     "gen_store_to_release_ratio"])
    
    for i in range(len(grid.plant.index)):
        writer.writerow([''.join(["g", str(grid.plant.index[i]), "i"]), 
                         grid.plant["type"].iloc[i], grid.plant["bus_id"].iloc[i], 0, 
                         dic_gen_cap_lim[grid.plant["type"].iloc[i]], 
                         CalcHeatRate(grid.plant.index[i]), df_gen_vom.iloc[i], 
                         dic_gen_age[grid.plant["type"].iloc[i]], 0, 0, 0, 
                         dic_variable[grid.plant["type"].iloc[i]], 
                         dic_baseload[grid.plant["type"].iloc[i]], 
                         0, dic_fuel[grid.plant["type"].iloc[i]], 
                         '.', '.', '.', '.', '.'])
    
        writer.writerow([''.join(["g", str(grid.plant.index[i])]), 
                         grid.plant["type"].iloc[i], grid.plant["bus_id"].iloc[i], 0, 
                         '.', 
                         CalcHeatRate(grid.plant.index[i]), df_gen_vom.iloc[i], 
                         dic_gen_age[grid.plant["type"].iloc[i]], 0, 0, 0, 
                         dic_variable[grid.plant["type"].iloc[i]], 
                         dic_baseload[grid.plant["type"].iloc[i]], 
                         0, dic_fuel[grid.plant["type"].iloc[i]], 
                         '.', '.', '.', '.', '.'])

# write gen_build_predetermined.csv file.

with open('gen_build_predetermined.csv', 'w', newline='') as gen_exist_file:
    writer = csv.writer(gen_exist_file)
    writer.writerow(["GENERATION_PROJECT", "build_year", "gen_predetermined_cap"])
    
    for i in range(len(grid.plant.index)):
        writer.writerow([''.join(["g", str(grid.plant.index[i])]), 2019, 
                         grid.plant['Pmax'].iloc[i]])

# write gen_build_costs.csv file.

with open('gen_build_costs.csv', 'w', newline='') as gen_buildcost_file:
    writer = csv.writer(gen_buildcost_file)
    writer.writerow(["GENERATION_PROJECT", "build_year", "gen_overnight_cost", "gen_fixed_om"])
    
    for i in range(len(grid.plant.index)):
        writer.writerow([''.join(["g", str(grid.plant.index[i])]), 2019, 
                         dic_geninvcost[grid.plant["type"].iloc[i]], 
                         df_gen_fix_cost.iloc[i]])
    
    for i in range(inv_stage):
        for j in range(len(grid.plant.index)):
            writer.writerow([''.join(["g", str(grid.plant.index[j]), "i"]), inv_period[i], 
                             dic_geninvcost[grid.plant["type"].iloc[j]], 
                             df_gen_fix_cost.iloc[j]])
        
# write transmission_lines.csv file.
# no circuit ID in PowerSimData.

with open('transmission_lines.csv', 'w', newline='') as transmission_file:
    writer = csv.writer(transmission_file)
    writer.writerow(["TRANSMISSION_LINE", "trans_lz1", "trans_lz2", "trans_length_km", 
                     "trans_efficiency", "existing_trans_cap"])

    for i in range(len(grid.dcline.index)):
        writer.writerow([''.join(["l", str(grid.dcline.index[i])]), 
                         grid.dcline["from_bus_id"].iloc[i], grid.dcline["to_bus_id"].iloc[i], 
                         GIS_distance2(grid.dcline["from_bus_id"].iloc[i], 
                                      grid.dcline["to_bus_id"].iloc[i]), 0.99, 
                         grid.dcline["Pmax"].iloc[i]])
    
    for i in range(len(grid.branch.index)):
        writer.writerow([''.join(["l", str(grid.branch.index[i])]), 
                         grid.branch["from_bus_id"].iloc[i], grid.branch["to_bus_id"].iloc[i], 
                         GIS_distance(grid.branch.index[i]), 
                         Branch_efficiency(grid.branch.index[i]), 
                         grid.branch["rateA"].iloc[i]])
    

# write trans_params.csv file.

with open('trans_params.csv', 'w', newline='') as trans_params_file:
    writer = csv.writer(trans_params_file)
    writer.writerow(["trans_capital_cost_per_mw_km", "trans_lifetime_yrs", 
                     "trans_fixed_om_fraction", "distribution_loss_rate", ])
    writer.writerow([621, 40, 0, 0])


# write periods.csv file.

with open('periods.csv', 'w', newline='') as periods_file:
    writer = csv.writer(periods_file)
    writer.writerow(["INVESTMENT_PERIOD", "period_start", "period_end"])
    
    for i in range(inv_stage):
        writer.writerow([inv_period[i], period_start[i], period_end[i]])


# end.
                         


