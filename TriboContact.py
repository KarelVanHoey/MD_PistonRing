# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 13:18:44 2020
@ProjectTitle : Thermal Mixed-Hydrodynamic Lubrication Solver for simulation of a compression ring - cylinder liner contact in an IC Engine.
@author: Dieter Fauconnier
@version : 1.0
"""

import numpy as np
import scipy.special as special # Sparse Matrix Definitions and operations
import scipy.integrate as integral # Sparse Matrix Definitions and operations

class TriboContact:
    
    def __init__(self,Engine):
    
        self.Engine=Engine
        
        """ Equivalent Young's modulus of Hertzian contact"""
        self.YoungsModulus=1.0/((1.0-Engine.Cylinder.Material.PoissonModulus**2.0)/Engine.Cylinder.Material.YoungsModulus + (1.0-Engine.CompressionRing.Material.PoissonModulus**2)/Engine.CompressionRing.Material.YoungsModulus)
        self.Domain=np.array([-Engine.CompressionRing.Thickness/2,Engine.CompressionRing.Thickness/2])
        
        """ Roughness parameters """
        self.Roughness=np.sqrt(Engine.Cylinder.Roughness**2.0 + Engine.CompressionRing.Roughness**2.0)
        self.Zeta=97.0e9
        self.Kappa=1.56e-6
        self.Tau0=2.0e6
        self.f_b=0.3
        self.RoughnessParameter=self.Zeta*self.Kappa*self.Roughness
        self.RoughnessSlope=self.Roughness/self.Kappa
        self.Lambda_c = 2.2239

        """Geometry of the cylinder and ring"""
        self.L = 2 * np.pi * self.Engine.Cylinder.Radius
        self.b = self.Engine.CompressionRing.Thickness
        self.delta = self.Engine.CompressionRing.CrownHeight
        self.R_cylinder = Engine.Cylinder.Radius
        self.R_piston = Engine.Piston.Radius
    
        
        """Wear Coefficients"""
        self.WearCoefficient_Cylinder=2.5e-10
        self.WearCoefficient_CompressionRing=1.25e-10

    def I2(self,l): 
        I2 = (0.5*(l**2+1)*special.erfc(l/np.sqrt(2.0)) - (l/np.sqrt(2.0*np.pi))*np.exp(-l**2.0/2.0))/np.sqrt(l)
        return I2
    
    def I52(self,l):
        I52 = ((1.0/(8.0*np.sqrt(np.pi)))*np.exp(-l**2.0/4.0)*(l**(3.0/2.0))*((2.0*l**2.0+3.0)*special.kv(3.0/4.0,l**2.0/4.0)-(2.0*l**2.0+5.0)*special.kv(1.0/4.0,l**2.0/4.0)))/np.sqrt(l)
        return I52
    
    def I2_lambda(self, l):
        return self.I2(l) * (l ** (-0.5))
    
    def I52_lambda(self, l):
        return self.I52(l) * (l ** (-0.5))


#################
##### TO DO #####
#################
    def AsperityContact(self,StateVector,time):

        Lambda=StateVector[time].Lambda

        AsperityArea = (np.pi**2) * ((self.RoughnessParameter) ** 2) * self.L * np.sqrt(self.Roughness * (self.b ** 2) * 0.25 / self.delta) * integral.quad(self.I2_lambda, Lambda, self.Lambda_c)[0]
        AsperityLoad = 1.06666667 * np.sqrt(2) * np.pi * ((self.RoughnessParameter) ** 2) * np.sqrt(self.Roughness / self.Kappa) * self.YoungsModulus * np.sqrt(self.Roughness * (self.b ** 2) * 0.25 / self.delta) * integral.quad(self.I52_lambda, Lambda, self.Lambda_c)[0]
        AsperityFriction = self.Tau0 * AsperityArea / self.L + self.f_b * AsperityLoad

        R_eq = 1 / ((1 / self.R_cylinder) + (1 / self.R_piston))

        StateVector[time].AsperityArea= AsperityArea
        StateVector[time].AsperityLoad= AsperityLoad
        StateVector[time].AsperityFriction= AsperityFriction
        StateVector[time].AsperityContactPressure= AsperityLoad / AsperityArea
        StateVector[time].HertzianContactPressure= (np.pi / 4) * np.sqrt((AsperityLoad * self.YoungsModulus) / (np.pi * R_eq))
        
        
#################
##### TO DO #####
#################       
    def Wear(self,Ops,Time,StateVector,time):
        
        # Calculate Wear Depth on the Piston Ring  
        StateVector[time].WearDepthRing= # accumulated wear depth on the ring         
        # Calculate The Wear Depth on the Cylinder wall
        StateVector[time].WearLocationsCylinder= # array of unique Positions where the pistion passes by  
        StateVector[time].WearDepthCylinder= #incremental wear depth on the positions in the array above