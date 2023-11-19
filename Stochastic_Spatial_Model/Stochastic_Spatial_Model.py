# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a stochastic-spatial model of intraspecific competition
"""
#import division for division?
from __future__ import division 
#import maptlotlib for plotting
import matplotlib.pyplot as plt
#import numpy for arrays
import numpy as np
#import pandas for saving data
import pandas as pd
#import os to make directories for saving simulation outputs
import os
#import colors from matplotlib for custom color scheme
from matplotlib import colors
#%%
#set font for graphs
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Times New Roman"

#set seed
np.random.seed(seed=1769)

#names for folders
experiment_title = 'm_biased_dispersal_New'
folder = 'm_images_biased_dispersal_New'
graph_folder = 'm_graphs_biased_dispersal_New'

#make folder for saving output
if os.path.isdir(folder)==False:
    os.system('mkdir '+folder)
if os.path.isdir(graph_folder)==False:
    os.system('mkdir '+graph_folder)

#lists to save data to
spc_3=[]
spc_2=[]
spc_1=[]
time_array=[]
spc_3_b=[]
spc_2_b=[]
spc_1_b=[]
spc_3_d=[]
spc_2_d=[]
spc_1_d=[]
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

A=np.array([[0,0,0,0],[0.036,(0.6078-1)/10,(0.6078-1)/10,(0.6078-1)/10],[0.042,(0.7198-1)/10,(0.7198-1)/10,(0.7198-1)/10],[0.109,(0.1942-1)/10,(0.1942-1)/10,(0.1942-1)/10]])


x1=0.5
num_iter=3500
dt=1
tmax=dt*num_iter
moves = np.array([[-1,0],[1,0],[0,-1],[0,1]])
diffusion=1
ranger=4
#init experiment 
x = np.random.choice([0,1,2,3],[L,L], p=[0.7,.1,0.1,0.1])

#function
def pbc_reset(old_pos,L):
    new_pos = old_pos
    for i in range(2):
        if old_pos[i] >= L:
            new_pos[i]=0
        if old_pos[i] < 0:
            new_pos[i]= L-1
    return new_pos

#this runs the simulation
for i in range(num_iter):
    species=[0,0,0]
    species_births=[0,0,0]
    species_neighbors=[[],[],[]]
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
    #save data for later graphing
    spc_1.append(species[0])
    spc_2.append(species[1])
    spc_3.append(species[2])
    time_array.append(i)
    spc_1_b.append(species_births[0])
    spc_2_b.append(species_births[1])
    spc_3_b.append(species_births[2])
    spc_1_d.append(species_neighbors[0])
    spc_2_d.append(species_neighbors[1])
    spc_3_d.append(species_neighbors[2])
    
    if species[1]/N >= 0.7198*0.9:
        print('pa14 dispersal event')
        A=np.array([[0,0,0,0],[0.036,(0.6078-1)/10,(0.6078-1)/10,(0.6078-1)/10],[-1,-1,-1,-1],[0.109,(0.1942-1)/10,(0.1942-1)/10,(0.1942-1)/10]])
        # A=np.array([[0,0,0,0],[-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1]])
    if species[1]/N <= 0.7198*0.1:
        print('pa14 dispersal event over')
        A=np.array([[0,0,0,0],[0.036,(0.6078-1)/10,(0.6078-1)/10,(0.6078-1)/10],[0.042,(0.7198-1)/10,(0.7198-1)/10,(0.7198-1)/10],[0.109,(0.1942-1)/10,(0.1942-1)/10,(0.1942-1)/10]])
    
    # x =np.take(x,np.random.permutation(x.shape[0]),axis=0,out=x)
    if i%10==0:
        rng = np.random.default_rng()
        x = rng.permuted(x)
        #display images

    if i%10==0 or i==799 or i==0:
        frame_counter+=1
        plt.figure()
        plt.axis('off')
        plt.title('gen %d'%i)
        plt.imshow(x,interpolation='nearest', origin='lower', cmap=cmap, norm=norm)
        if len(str(frame_counter))==1:
            plt.savefig(folder+'/gen_00%d'%frame_counter, dpi=300)
        elif len(str(frame_counter))==2:
            plt.savefig(folder+'/gen_0%d'%frame_counter, dpi=300)
        elif len(str(frame_counter))==3:
            plt.savefig(folder+'/gen_%d'%frame_counter, dpi=300)
        else:
            frame_counter-=1
        plt.show()
#%%
#save data to csv for later 
volume_data = dict()
volume_data['Species1']=spc_1
volume_data['Species2']=spc_2
volume_data['Species3']=spc_3
volume_data['Species1_b']=spc_1_b
volume_data['Species2_b']=spc_2_b
volume_data['Species3_b']=spc_3_b
volume_data['Species1_Density']=spc_1_d
volume_data['Species2_Density']=spc_2_d
volume_data['Species3_Density']=spc_3_d
volume_data['Time']=time_array
volume_dataframe = pd.DataFrame(data=volume_data)
volume_dataframe.to_csv(experiment_title+'.csv', index=False)
#%%
# volume_dataframe = pd.read_csv(experiment_title+'.csv')

# show snapshot of the final figure
plt.figure()
plt.axis('off')
plt.title('snap at time '+str(i))
plt.imshow(x,interpolation='nearest', origin='lower', cmap=cmap, norm=norm)
plt.savefig(graph_folder+'/ending_snap', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
    
plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
# plt.title('Stochastic Spatial Simulation')
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)
# plt.grid(alpha=0.25)
# plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig(graph_folder+'/volume.svg', dpi=300, facecolor='w', edgecolor='b',
        orientation='portrait', format='svg',
        transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()