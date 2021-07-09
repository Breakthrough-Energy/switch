# Copyright (c) 2015-2019 The Switch Authors. All rights reserved.
# Licensed under the Apache License, Version 2.0, which is in the LICENSE file.
"""
Add a soft transmission cap to the model. Specifying transmission_softcap will impose 
a system-wide cap on the total MW-km of transmission in the system. When any additional
transmission is built beyond this cap, a financial cost equal to transmission_overbuild_penalty
will be applied to every additional increment of transmission capacity. 


"""
from __future__ import division
import os
from pyomo.environ import Set, Param, Expression, Constraint, Var, Suffix
from pyomo.environ import *
import pandas as pd
import switch_model.reporting as reporting

def define_components(model):
    model.TxPenalty = Var(model.PERIODS, within=NonNegativeReals)
    model.transmission_softcap = Param(model.PERIODS, default=1e10,
        doc="A soft cap above which a penalty is applied to additional tranmsission buildout.")
    model.transmission_overbuild_penalty= Param(model.PERIODS, default=0,
        doc="Penalty applied to transmission buildout per MW-km above the cap.")
    model.TransmissionOverCosts = Expression(model.PERIODS,
        rule=lambda model, period: model.transmission_overbuild_penalty[period] * (sum(
            model.TxCapacityNameplate[tx, period] * model.trans_length_km[tx]
            for tx in model.TRANSMISSION_LINES
            ) - model.transmission_softcap[period]),
        doc=("Calculates the penalty for transmission overbuild emissions."))
    model.ApplyTxPenalty = Constraint(model.PERIODS,
        rule=lambda model, period:
            model.TxPenalty[period] >= model.TransmissionOverCosts[period])
    model.Cost_Components_Per_Period.append('TxPenalty')


def load_inputs(model, switch_data, inputs_dir):
     """
     A single file can be used to specify both transmission_softcap and transmission_overbuild_penalty.

     Expected input files:
     transmission_cap.csv
         PERIOD, transmission_softcap, transmission_overbuild_penalty

     """
     switch_data.load_aug(
         filename=os.path.join(inputs_dir, 'transmission_cap.csv'),
         autoselect=True,
         index=model.PERIODS,
         param=(model.transmission_softcap, model.transmission_overbuild_penalty,))

