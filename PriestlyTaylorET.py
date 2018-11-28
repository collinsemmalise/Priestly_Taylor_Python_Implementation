#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 10:32:31 2018

@author: Emma Collins
Last Updated: November 5, 2018
Priestly Taylor Modeling for Potential Evapotranspiration
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates
import datetime


def calculate_net_radiation(albedo, Rs, C, epsilon, T_k):
    """
    Calculates Net Solar Radiation in MJ/m^2/day and returns value or list
    
    keyword arguments:
        albedo = albedo of area
        Rs = solar radiation in MJ/m^2/day
        C = cloudiness factor
        epsilon = Net Emissivity
        T_k = Temperature (Kelvin)
    """
    o = 4.89*10**-9
    Rn = (1-albedo)*Rs-C*epsilon*o*((T_k)**4)
    return Rn


def calculate_cloudiness(Rs, Ra,ac=0.72, bc=0.28):
    """
    Calculates Cloudiness Factor and returns value or numpy array
    
    keyword arguments:
        Rs = Solar radiation in MJ/m^2/day
        Ra = Extraterrestrial Solar Radiation in MJ/m^2/day
        ac = constant (default is 0.72)
        bc = constant (default is 0.28)
    """
    Ra = np.transpose(Ra)[0]
    R = np.divide(Rs, Ra, out=np.zeros_like(Rs), where=Ra!=0)
    return ac*R+bc


def calculate_Delta(T):
    """ 
    Calculates Slope of Vaporization Curve (in kPa/˚C) 
        at Temperature T (from Tetens, 1930) and returns value or numpy array
    
    key arguments:
        T = Temperature in Celsius 
    """
    Delta = 4098*(0.6108*np.exp((17.27*T)/(T+237.3)))/(T+237.3)**2
    return Delta


def calculate_Bowen(gamma, T2,T1,e2,e1):
    """
    Calculates Approximation of Bowen's Ratio (Drexler, 2004) 
    and returns value or list
    
    key arguments:
        gamma = psychrometric constant
        T2 = temperature of higher elevation (Celsius or Kelvin)
        T1 = temperature of lower elevation (Celsius or Kelvin)
        e2 = vapor pressure of higher elevation (kPa)
        e1 = vapor pressure of lower elevation (kPa)
    """
    return gamma*((T2-T1)/(e2-e1))


def calculate_alpha(Delta, gamma, beta):
    """ 
    Calculates Saturation Deficit Factor and returns value or list
    
    key arguments:
        Delta = Slope of Vaporization Curve
        gamma = psychrometric constant
        beta = Bowen's Ratio
    """
    return (Delta+gamma)/(Delta*(1+beta))


def psychrometric_constant(P, lam):
    """
    Calculates Psychrometric Constant and returns value or list (kPa/˚C)
    
    key arguments:
        P = pressure (kPa)
        lam = latent heat of vaporization
    """
    return 0.00163*P/lam


def calculate_vapor_pressure(T):
    """
    Calculates Saturation Vapor Pressure (kPa) and returns value or numpy array
    
    key arguments:
        T = Temperature (Celsius) 
    """
    e = 0.6108*np.exp(17.27*T/(237.3+T))
    return e


def calculate_actual_vapor_pressure(T, RH):
    """ 
    Calculates Actual Vapor (kPa) and returns value or numpy array
    
    key arguments:
        T = Tempertaure (Celsius)
        RH = Relative Humidity  (%)
    """
    es = calculate_vapor_pressure(T)
    e = RH*es/100
    return e


def Priestley_Taylor(alpha, Delta, gamma, Rn):
    """ 
    Calculates Evapotranspiration (MJ/m^2/day) and returns value or list
    
    key arguments:
        alpha = saturation deficit constant
        Delta = slope of vaporization curve
        gamma = psychrometric constant
        Rn = net solar radiation
    """
    return alpha*((Delta)/(Delta+gamma))*Rn

def net_Emissivity(T):
    """
    Calculates Net Emissivity and returns value or numpy array
    
    key arguments:
        T = Temperature (Celsius)
    """
    return 0.261*np.exp(-7.77*10**-4*T**2)-0.02


def PET(air_T, fuel_T, elevation,RH, fuel_moist, Rs, Ra, albedo):
    """
    Calculates Potential ET using Priestley-Taylor Equations 
    (Priestley and Taylor 1972) and returns value or list
    
    key arguments:
        air_T = air temperature (Celsius)
        fuel_T = fuel temperature or temperature of ground (Celsius)
        elevation = elevation of ground surface (meters)
        RH = relative humidity (%)
        fuel_moist = fuel moisture or moisture of ground (%)
        Rs = net solar radiation (MJ/m^2/day)
        Ra = extraterrestrial solar radiation (MJ/m^2/day)
        albedo = albedo constant
    """
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

def graph_ET_results(date_time, evapotranspiration, title = "Evapotranspiration", ylabel = 'mm', xlabel = 'Date', verbose = True):
    """
    Graphs PET results and returns matplotlib dates object from x axis
    
    key arguments:
        date_time = pandas series of Timestamp objects 
        evapotranspiration = PET estimates corresponding with Timestamps
        title = title for graph (default is "Evapotranspiration")
        ylabel = y axis label for graph (default is 'mm')
        xlabel = x axis label for graph (default is 'Date')
        verbose(bool) if True, displays graph (default is true)
    """
    
    dates = []
    for i in range(len(date_time)):
        dates.append(datetime.datetime.date(date_time.iloc[i]))
    dates = sorted(list(set(dates)))
    dates_plot = matplotlib.dates.date2num(dates)
    if(verbose):
        plt.figure(figsize = (15,10))
        plt.plot_date(dates_plot, np.asarray(evapotranspiration),'-', ydate = False)
        plt.title(title)
        plt.xticks(rotation = 'vertical')
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.show()
    return dates
    
