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
        # self.R_cylinder = Engine.Cylinder.Radius
        # self.R_piston = Engine.Piston.Radius
        
        """Wear Coefficients"""
        self.WearCoefficient_Cylinder=2.5e-10
        # self.WearCoefficient_Cylinder=29e3 * 3.2e-6 / 1000 # COATING
        self.WearCoefficient_CompressionRing=1.25e-10
        # self.WearCoefficient_CompressionRing = 29e3 * 3.2e-6 / 1000 # COATING
        

    def I2(self,l): 
        I2 = (0.5*(l**2+1)*special.erfc(l/np.sqrt(2.0)) - (l/np.sqrt(2.0*np.pi))*np.exp(-l**2.0/2.0))/np.sqrt(l)
        return I2
    
    def I52(self,l):
        I52 = ((1.0/(8.0*np.sqrt(np.pi)))*np.exp(-l**2.0/4.0)*(l**(3.0/2.0))*((2.0*l**2.0+3.0)*special.kv(3.0/4.0,l**2.0/4.0)-(2.0*l**2.0+5.0)*special.kv(1.0/4.0,l**2.0/4.0)))/np.sqrt(l)
        return I52


#################
##### TO DO #####
#################
    def AsperityContact(self,StateVector,time):

        # Lambda = min(StateVector[time].Lambda, self.Lambda_c)
        Lambda = StateVector[time].Lambda

        if Lambda < self.Lambda_c:      # Asperities make contact with cylinder --> Asperity contact area, load and pressure > 0

            AsperityContactArea = ((np.pi * self.RoughnessParameter) ** 2) * self.L * np.sqrt(self.Roughness * (self.b ** 2) * 0.25 / self.delta) * integral.quad(self.I2, Lambda, self.Lambda_c,limit=100)[0]
            AsperityLoad = 1.06666667 * np.sqrt(2) * np.pi * (self.RoughnessParameter ** 2) * np.sqrt(self.Roughness / self.Kappa) * self.YoungsModulus * np.sqrt(self.Roughness * (self.b ** 2) * 0.25 / self.delta) * integral.quad(self.I52, Lambda, self.Lambda_c,limit=100)[0]
            AsperityFriction = self.Tau0 * AsperityContactArea / self.L + self.f_b * AsperityLoad
            StateVector[time].AsperityContactPressure= AsperityLoad / AsperityContactArea       # Nowhere used, but can be plotted
        
        else:                           # Chance of asperity contact is zero --> no asperity load

            AsperityContactArea = 0.0
            AsperityLoad = 0.0
            AsperityFriction = 0.0
            StateVector[time].AsperityContactPressure= 0.0       # Nowhere used, but can be plotted

        # R_eq = self.Engine.CompressionRing.Curvature   # Uit paper van Gore et al.: vergeleken met curvature gedefinieerd in Engine, beetje verwarrend want hier: curvature == R_eq uit Engine (normaal gezien: curvature = R^(-1))

        StateVector[time].AsperityContactArea= AsperityContactArea
        StateVector[time].AsperityLoad= AsperityLoad
        StateVector[time].AsperityFriction= AsperityFriction
        StateVector[time].HertzianContactPressure= (np.pi / 4) * np.sqrt((AsperityLoad * self.YoungsModulus) / (np.pi * self.Engine.CompressionRing.Curvature))
        
        
#################
##### TO DO #####
#################       
    def Wear(self,Ops,Time,StateVector,time):
        
        # p_tmin1 = StateVector[time-1].HertzianContactPressure
        p_t = StateVector[time].HertzianContactPressure
        s_t = Ops.SlidingDistance[time-1]

        if time > 1:  
            s_tmin1 = Ops.SlidingDistance[time-2]

        elif time == 1:
            s_tmin1 = 0.0

        # Calculate Wear Depth on the Piston Ring
        # StateVector[time].WearDepthRing= StateVector[time-1].WearDepthRing + self.WearCoefficient_CompressionRing / self.Engine.CompressionRing.Material.Hardness * (p_tmin1 + p_t) * 0.5 * (s_t - s_tmin1) # accumulated wear depth on the ring         
        StateVector[time].WearDepthRing = StateVector[time-1].WearDepthRing + self.WearCoefficient_CompressionRing / self.Engine.CompressionRing.Material.Hardness * p_t * (s_t - s_tmin1) # accumulated wear depth on the ring   
        # Calculate The Wear Depth on the Cylinder wall

            # assumptions: From asperity contact size, the height of the contact zone is determined. Uniform wear due to p_Hertz is assumed in that zone. The position of the piston in the cylinder is known, so the wear zone is known.
        # WearLength = StateVector[time].AsperityContactArea / self.L
        # dL = WearLength / 20
        # StateVector[time].WearLocationsCylinder= np.arange(Ops.PistonPosition[time] - WearLength, Ops.PistonPosition[time] - WearLength + dL, step=dL) # array of unique positions where the piston passes by
        # StateVector[time].WearDepthCylinder= self.WearCoefficient_Cylinder / self.Engine.Cylinder.Material.Hardness * p_t * np.abs(Ops.PistonVelocity[time]) * Time.dt * np.ones(np.size(StateVector[time].WearLocationsCylinder)) # incremental wear depth on the positions in the array above
        
        StateVector[time].WearLocationsCylinder = np.unique(np.round(Ops.PistonPosition,8))
        WearIncrement = np.zeros(np.size(StateVector[time].WearLocationsCylinder))
        # WearIncrement[StateVector[time].WearLocationsCylinder.index(np.round(Ops.PistonPosition[time],8))] += self.WearCoefficient_Cylinder / self.Engine.Cylinder.Material.Hardness * p_t * np.abs(Ops.PistonVelocity[time]) * Time.dt
        WearIncrement[np.where(StateVector[time].WearLocationsCylinder == np.round(Ops.PistonPosition[time],8))[0][0]] += self.WearCoefficient_Cylinder / self.Engine.Cylinder.Material.Hardness * p_t * np.abs(Ops.PistonVelocity[time]) * Time.dt
        StateVector[time].WearDepthCylinder = StateVector[time - 1].WearDepthCylinder + WearIncrement