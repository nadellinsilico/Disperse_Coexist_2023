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
from itertools import product
#import os to make directories for saving simulation outputs
import os
#%%
#set font for graphs
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times New Roman"
#names for folders
experiment_title = 'ternary_plot_f'

if os.path.isdir(experiment_title)==False:
    os.system('mkdir '+experiment_title)

#%%
"""
This sections gets starting points
"""
# Define the range and increments
start = 0
stop = 1
increment = 0.05
reps = 3

# Calculate the number of steps
num_steps = int((stop - start) / increment) + 1
# Generate all combinations
combinations = list(product(range(num_steps), repeat=3))
# Convert combinations to actual values
X = [[start + c[0] * increment, start + c[1] * increment, start + c[2] * increment] for c in combinations]

to_test = list(filter(lambda num: np.sum(num) != 0, X))

print(len(to_test))
ec_start_list = []
pa_start_list = []
ef_start_list = []
identifier_start_list = []

identity = 0
for value in to_test:
    if round(np.sum(value), 2) == 1:
        for c in np.arange(0,reps):
            identity+=1
            ec_start_list.append(value[0])
            pa_start_list.append(value[1])
            ef_start_list.append(value[2])
            identifier_start_list.append(identity)
        
print(len(pa_start_list))
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

x1=0.05
num_iter=100
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
#%%


def run_simulation(ec_start, pa_start, ef_start, identifier, gpa=gpa, gec=gec, gef=gef, kpa=kpa, kec=kec, kef=kef):

    A=np.array([[0,0,0,0],[gec,kec,kec,kec],[gpa,kpa,kpa,kpa],[gef,kef,kef,kef]])

    #this runs the simulation
    dispersal=0
    spc_3=[]
    spc_2=[]
    spc_1=[]
    dispersal_counter = []
    simulation_time = []
    
    x = np.random.choice([0,1,2,3],[L,L], p=[0.7,ec_start,pa_start,ef_start])
    
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
            A=np.array([[0,0,0,0],[gec,kec,kec,kec],[-1,-1,-1,-1],[gef,kef,kef,kef]])
            # A=np.array([[0,0,0,0],[-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1]])
        if species[1]/N <= 0.7198*0.1:
            dispersal=0
            # print('pa14 dispersal event over')
            A=np.array([[0,0,0,0],[gec,kec,kec,kec],[gpa,kpa,kpa,kpa],[gef,kef,kef,kef]])
        #save data for later graphing
        if i%5 == 0:
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
    volume_dataframe.to_excel(experiment_title+'/'+experiment_title+'_'+str(identifier)+'.xlsx', index=False)

#%%
if __name__ == "__main__":
    # first way, using multiprocessing
    arguments = [(ec_start, pa_start, ef_start, identifier) 
                 for ec_start in ec_start_list for pa_start in pa_start_list for ef_start in ef_start_list for identifier in identifier_start_list]
    start_time = time.perf_counter()
    with Pool(processes=30) as pool:
        results = pool.starmap(run_simulation, arguments)
    finish_time = time.perf_counter()
    print("Program finished in {} seconds - using multiprocessing".format(finish_time - start_time))