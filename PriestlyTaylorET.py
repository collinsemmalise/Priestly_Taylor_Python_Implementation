#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 10:32:31 2018

@author: Emma Collins
Last Updated: August 16, 2018
Priestly Taylor Modeling for Evapotranspiration
Based on the Equations presented by Priestly and Taylor (1972)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates
import datetime

#input: albedo constant, solar radiation ((MJ/m**2)/d), cloudiness factor, net
#       emissivity, Temperature (Kelvin)
#output: net radiation ((MJ/m**2)/day)
def calculate_net_radiation(albedo, Rs, C, epsilon, T_k):
    o = 4.89*10**-9
    Rn = (1-albedo)*Rs-C*epsilon*o*((T_k)**4)
    return Rn

#input: solar radiation (MJ/m**2/d), extraterrestrial radiation (MJ/m**2/d), ac constant, bc constant
#output: cloudiness factor
def calculate_cloudiness(Rs, Ra,ac=0.72, bc=0.28):
    Ra = np.transpose(Ra)[0]
    R = np.divide(Rs, Ra, out=np.zeros_like(Rs), where=Ra!=0)
    return ac*R+bc

#input: Temperature (˚C)
#output: slope of vaporization curve at temp T (kPa/˚C)
def calculate_Delta(T):
    Delta = 4098*(0.6108*np.exp((17.27*T)/(T+237.3)))/(T+237.3)**2
    return Delta

#input: psychrometric constant, temperature of top and bottom, vapor
#       pressure of top and bottom
#output: Bowen's ratio (K/kPa) or (˚C/kPa)
def calculate_Bowen(gamma, T2,T1,e2,e1):
    return gamma*((T2-T1)/(e2-e1))

#input: slope of vaporization curve, psychrometric constant, Bowen's Ratio
#output: saturation deficit factor
def calculate_alpha(Delta, gamma, beta):
    return (Delta+gamma)/(Delta*(1+beta))

#input: Pressure (kPa), latent heat of vaporization (lambda)
#output: psychrometric constant (kPa/˚C)
def psychrometric_constant(P, lam):
    return 0.00163*P/lam

#input: temperature (deg Celsius)
#output: vapor pressure (kPa)
def calculate_vapor_pressure(T):
    e = 0.6108*np.exp(17.27*T/(237.3+T))
    return e

#input: air temperature (deg Celsius), relative humidity (%)
#output: actual vapor pressure (kPa)
def calculate_actual_vapor_pressure(T, RH):
    es = calculate_vapor_pressure(T)
    e = RH*es/100
    return e

#input: saturation deficit constant, slope of vaporization curve,
#       psychrometric constant, net solar radiation
#output: evapotranspiration (MJ/m**2/day)
def Priestley_Taylor(alpha, Delta, gamma, Rn):
    return alpha*((Delta)/(Delta+gamma))*Rn

# Calculate Net emissivity
# input: temperature
# output: net emissivity
def net_Emissivity(T):
    return 0.261*np.exp(-7.77*10**-4*T**2)-0.02


#input: air temp(Celsius), Fuel temp (Celsius), elevation(m), RH(%)
#       fuel moisture (%), net solar radiation, extraterrestrial solar radiation
#       albedo constant
#output: ET
def PET(air_T, fuel_T, elevation,RH, fuel_moist, Rs, Ra, albedo):
    P = 101.3 - 0.01055*elevation
    lam = 2.501 - 0.002361*air_T
    Delta = calculate_Delta(air_T)
    gamma = psychrometric_constant(P, lam)
    air_vape = calculate_actual_vapor_pressure(air_T, RH)
    fuel_vape = calculate_actual_vapor_pressure(fuel_T, fuel_moist)
    beta = calculate_Bowen(gamma, air_T,fuel_T,air_vape,fuel_vape)
    alpha = calculate_alpha(Delta, gamma, beta)
    C = calculate_cloudiness(Rs, Ra,ac=0.72, bc=0.28)
    T_k = air_T-273.15
    epsilon = net_Emissivity(air_T)
    Rn = calculate_net_radiation(albedo, Rs, C, epsilon, T_k)
    ET = Priestley_Taylor(alpha, Delta, gamma, Rn)
    return ET

#Graphs the results of PET estimation
# input: dates for the x values (as yyyy-mm-dd hh:mm), PET estimates for y values, optional title, y label, and xlabel
def graph_ET_results(date_time, evapotranspiration, title = "Evapotranspiration", ylabel = 'mm', xlabel = 'Date'):
    dates = []
    for i in range(len(date_time)):
        h = datetime.datetime.strptime(date_time[i], "%Y-%m-%d %H:%M")
        dates.append(datetime.datetime.date(h))
    dates = sorted(list(set(dates)))
    dates_plot = matplotlib.dates.date2num(dates)
    plt.figure(figsize = (15,10))
    plt.plot_date(dates_plot, np.asarray(evapotranspiration),'-', ydate = False)
    plt.title(title)
    plt.xticks(rotation = 'vertical')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.show()
    return dates
