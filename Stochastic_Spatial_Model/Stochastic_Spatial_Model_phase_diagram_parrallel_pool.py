# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a stochastic-spatial model of intraspecific competition
"""
from __future__ import division 
#import division for division?
import glob

#import maptlotlib for plotting
import matplotlib.pyplot as plt
#import numpy for arrays
import numpy as np
#import pandas for saving data
import pandas as pd
#import colors from matplotlib for custom color scheme
from matplotlib import colors
from multiprocessing import Pool
import time
#import os to make directories for saving simulation outputs
import os
#%%
#set font for graphs
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times New Roman"
#names for folders
experiment_title = 'phase_diagram'

if os.path.isdir(experiment_title)==False:
    os.system('mkdir '+experiment_title)
#%%
#set colors for different species
cmap = colors.ListedColormap(['black','#9813DE', '#E2E604', '#008F8F'])
bounds=[0,1,2,3,4]
norm = colors.BoundaryNorm(bounds, cmap.N)

###global variables for simulation run
##frame counter to determine frame numbers
frame_counter=0
L=50
N=L**2

ratio_list = list(np.flip(np.arange(0,1.1,.1)))
growth_list = list(np.arange(.03,.06003,.003))

x1=0.05
num_iter=20000
dt=1
tmax=dt*num_iter
moves = np.array([[-1,0],[1,0],[0,-1],[0,1]])
diffusion=1
ranger=4
#init experiment 

#function
def pbc_reset(old_pos,L):
    new_pos = old_pos
    for i in range(2):
        if old_pos[i] >= L:
            new_pos[i]=0
        if old_pos[i] < 0:
            new_pos[i]= L-1
    return new_pos

simulation_iter = 0

gpa = 0.042
gec = 0.036
gef = 0.109

kpa = (0.7198-1)/10
kec = (0.6078-1)/10
kef = (0.1942-1)/10

file_list = glob.glob(experiment_title+'/*.xlsx', recursive=False)

for ratio in ratio_list:
    for growth in growth_list:
        tester_str = str(ratio)+'_'+str(growth)
        for file in file_list:
            if tester_str in file:
                if ratio in ratio_list:
                    ratio_list.remove(ratio)
                if growth in growth_list:
                    growth_list.remove(growth)
                continue
#%%


def run_simulation(ratio, growth, rep):

    A=np.array([[0,0,0,0],[gec,kec,kec,kec],[growth,kpa,kpa,kpa],[gef,kef,kef,kef]])

    #this runs the simulation
    dispersal=0
    spc_3=[]
    spc_2=[]
    spc_1=[]
    dispersal_counter = []
    simulation_time = []
    
    x = np.random.choice([0,1,2,3],[L,L], p=[0.7,0.1,0.1,0.1])
    
    for i in range(num_iter):
        species=[0,0,0]
        for j in range(N):
            foc_pos = np.random.choice(range(L),[1,2])[0]
            foc_id =x[foc_pos[0],foc_pos[1]]
            payoff_list = []
            for move in moves:
                opp_pos = foc_pos + move
                opp_pos = pbc_reset(opp_pos,L)
                opp_id= x[opp_pos[0], opp_pos[1]]
                payoff_list.append(A[foc_id,opp_id])
            payoff = np.mean(payoff_list)
            if foc_id > 0:
                species[foc_id-1]+=1
            if payoff > 0 :
                birth = np.random.uniform() < payoff*dt
                if birth:
                    rep_pos=foc_pos+moves[np.random.choice(range(4))]
                    rep_pos=pbc_reset(rep_pos,L)
                    x[rep_pos[0] , rep_pos[1]]=foc_id
            elif payoff < 0:
                death = np.random.uniform()<np.abs(payoff)*dt
                if death:
                   rep_pos=foc_pos+moves[np.random.choice(range(4))]
                   rep_pos=pbc_reset(rep_pos,L)
                   x[foc_pos[0],foc_pos[1]]=x[rep_pos[0],rep_pos[1]]

        if species[1]/N >= 0.7198*0.9:
            dispersal=1
            # print('pa14 dispersal event')
            if ratio == 0.0:
                A=np.array([[0,0,0,0],[gec,kec,kec,kec],[-1,-1,-1,-1],[gef,kef,kef,kef]])
            else:
                A=np.array([[0,0,0,0],[-ratio,-ratio,-ratio,-ratio],[-1,-1,-1,-1],[-ratio,-ratio,-ratio,-ratio]])
            # A=np.array([[0,0,0,0],[-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1]])
        if species[1]/N <= 0.7198*0.1:
            dispersal=0
            # print('pa14 dispersal event over')
            A=np.array([[0,0,0,0],[gec,kec,kec,kec],[growth,kpa,kpa,kpa],[gef,kef,kef,kef]])
        #save data for later graphing
        if i%20 == 0:
            spc_1.append(species[0]/(species[0]+species[1]+species[2]))
            spc_2.append(species[1]/(species[0]+species[1]+species[2]))
            spc_3.append(species[2]/(species[0]+species[1]+species[2]))
            dispersal_counter.append(dispersal)
            simulation_time.append(0)  
    #save data to csv for later 
    volume_data = dict()
    volume_data['Species1']=spc_1

    volume_data['Species2']=spc_2
    volume_data['Species3']=spc_3
    volume_data['Dispersal'] = dispersal
    volume_dataframe = pd.DataFrame(data=volume_data)
    volume_dataframe.to_excel(experiment_title+'/'+experiment_title+'_'+str(ratio)+'_'+str(growth)+'_'+str(rep)+'_'+str(num_iter)+'.xlsx', index=False)

#%%
if __name__ == "__main__":
    # first way, using multiprocessing
    arguments = [(ratio, growth, rep_number) for ratio in ratio_list for growth in growth_list for rep_number in np.arange(0, 20, 1)]
    start_time = time.perf_counter()
    with Pool(processes=30) as pool:
        results = pool.starmap(run_simulation, arguments)
    finish_time = time.perf_counter()
    print("Program finished in {} seconds - using multiprocessing".format(finish_time - start_time))