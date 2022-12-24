# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 17:21:10 2021
@ProjectTitle : Thermal Mixed-Hydrodynamic Lubrication Solver for simulation of a compression ring - cylinder liner contact in an IC Engine.
@author: Dieter Fauconnier
@version : 1.0
"""


import numpy as np 
import matplotlib.pyplot as plt


    
def Report_PT(Grid,State,time=None): # initiatlization

    f1, ax1 = plt.subplots()
    color = 'tab:blue'
    ax1.set_xlabel('$x [m]$')
    if time is not None:
        ax1.set_ylabel('$P [MPa]$ at time =' + str(time*5/100) + 'ms',color=color)
    else:
        ax1.set_ylabel('$P [MPa]$',color=color)
    ax1.plot(Grid.x,State.Pressure/1e6,'x-', linewidth=1,color=color)
    ax1.tick_params(axis='y')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('$T [^\circ C]$',color=color)  # we already handled the x-label with ax1
    ax2.plot(Grid.x,State.Temperature-273.15,'x-', linewidth=1,color=color)
    ax2.tick_params(axis='y')
    # f1.tight_layout()  # otherwise the right y-label is slightly clipped
    # plt.show()
    return f1
    

    


def Report_Ops(Time,Ops,time):     
    
    f2, ax1 = plt.subplots()
    color = 'tab:blue'
    # ax1.set_xlabel('$t [s]$')
    ax1.set_xlabel('$Crank angle [rad]$')
    ax1.set_ylabel('$U [m/s]$',color=color)
    # ax1.plot(Time.t,Ops.SlidingVelocity,'-',Time.t[time],Ops.SlidingVelocity[time],'ko', linewidth=1,color=color)
    ax1.plot(Ops.CranckAngle,Ops.SlidingVelocity,'-',Ops.CranckAngle[time],Ops.SlidingVelocity[time],'ko', linewidth=1,color=color)
    plt.xlabel('Crank angle $\psi$ [rad]')
    ax1.tick_params(axis='y')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('$F [N/m]$',color=color)  # we already handled the x-label with ax1
    # ax2.plot(Time.t,Ops.CompressionRingLoad,'-',Time.t[time],Ops.CompressionRingLoad[time],'ko',linewidth=1,color=color)
    ax2.plot(Ops.CranckAngle,Ops.CompressionRingLoad,'-',Ops.CranckAngle[time],Ops.CompressionRingLoad[time],'ko',linewidth=1,color=color)
    ax2.tick_params(axis='y')
    f2.tight_layout()  # otherwise the right y-label is slightly clipped
    ax1.vlines(Ops.CranckAngle[time],-15,15,'k','--', linewidth=.6)
    pi = np.pi
    psi = np.arange(0, 4 * pi + pi/2, step=(pi/2))
    plt.xticks(psi,['0','π/2', 'π', '3π/2', '2π','5π/2', '3π', '7π/2', '4π'])
    # ax2.set_xticks(psi,['0','π/2', 'π', '3π/2', '2π','5π/2', '3π', '7π/2', '4π'])
    ax1.set_ylim([-15,15])
    # plt.show()
    return f2

def Report_Ops_PT(Time,Ops,time, Grid,State):
    fig, axs = plt.subplots(2)
    color = 'tab:blue'
    axs[0].set_xlabel('$t [s]$')
    axs[0].set_ylabel('$U [m/s]$',color=color)
    axs[0].plot(Time.t,Ops.SlidingVelocity,'-',Time.t[time],Ops.SlidingVelocity[time],'ko', linewidth=1,color=color)
    axs[0].tick_params(axis='y')
    ax2 = axs[0].twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax2.set_ylabel('$F [N/m]$',color=color)  # we already handled the x-label with ax1
    ax2.plot(Time.t,Ops.CompressionRingLoad,'-',Time.t[time],Ops.CompressionRingLoad[time],'ko',linewidth=1,color=color)
    ax2.tick_params(axis='y')
    

    color = 'tab:blue'
    axs[1].set_xlabel('$x [m]$')
    axs[1].set_ylabel('$P [MPa]$',color=color)
    axs[1].plot(Grid.x,State.Pressure/1e6,'x-', linewidth=1,color=color)
    axs[1].tick_params(axis='y')
    ax3 = axs[1].twinx()  # instantiate a second axes that shares the same x-axis
    color = 'tab:red'
    ax3.set_ylabel('$T [^\circ C]$',color=color)  # we already handled the x-label with ax1
    ax3.plot(Grid.x,State.Temperature-273.15,'x-', linewidth=1,color=color)
    ax3.tick_params(axis='y')
    # ax3.set_ylim([0,5.5])
    axs[1].set_ylim([0,3.5])
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()
    return fig
