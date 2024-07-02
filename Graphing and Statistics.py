# -*- coding: utf-8 -*-
"""
Created on Fri November 15 10:01:12 2023

@author: holtj

This script replicates the figures present in the paper using the provided datasheet. 

The figures in the paper were generated directly from .mat files or .tiff files, 
which were in turn generated using BiofilmQ. 
"""
import glob
from itertools import product
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import ternary
import seaborn as sns
import matplotlib as mpl
from ast import literal_eval
from palettable.colorbrewer.sequential import Greys_6_r
from palettable.matplotlib  import Inferno_20
#%%
"""
Figure 1A

Planktonic Sock Plot
"""

df_og = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1A')
df = df_og.groupby('time').mean().reset_index()

t = np.arange(0,len(df['time'].unique()))

width=0.9
plt.bar(t, df.pa, width=width, color='#E2E604', label=r'${\itP. aeruginosa}$', alpha=0.75)
plt.bar(t, df.ef, width=width, color='#008F8F', bottom=np.array(df.pa), label=r'${\itE. faecalis}$', alpha=0.75)
plt.bar(t, df.ec, width=width, color='#9813DE', bottom=np.array(df.pa)+np.array(df.ef), label=r'${\itE. coli}$', alpha=0.75)
plt.scatter(np.ones(4)*0,np.ones(4)*0.333, facecolors='none', edgecolors='black', alpha=1)
plt.scatter(df_og['time']/24, df_og['pa'], facecolors='none', edgecolors='black', alpha=1)
plt.ylim(0,1)
plt.xticks(ticks=t, labels=np.array(df['time'].unique()), fontsize=20)
plt.yticks(fontsize=20)
plt.grid(False)
plt.ylabel(r'R.A.')
plt.xlabel('Time (d)')
plt.title('Liquid Shaking Culture (24 hr Back Dilutions)')
# plt.savefig('Fig1A.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 1B

Biofilm Sock Plot
"""

df_og = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1B')
df = df_og.groupby('time').mean().reset_index()

t = np.arange(0,len(df['time'].unique()))

width=0.9
plt.bar(t, df.pa, width=width, color='#E2E604', label=r'${\itP. aeruginosa}$', alpha=0.75)
plt.bar(t, df.ef, width=width, color='#008F8F', bottom=np.array(df.pa), label=r'${\itE. faecalis}$', alpha=0.75)
plt.bar(t, df.ec, width=width, color='#9813DE', bottom=np.array(df.pa)+np.array(df.ef), label=r'${\itE. coli}$', alpha=0.75)
plt.scatter(np.ones(4)*0,np.ones(4)*0.333, facecolors='none', edgecolors='black', alpha=1)
plt.scatter(df_og['time'], df_og['pa'], facecolors='none', edgecolors='black', alpha=1)
plt.xticks(ticks=t, labels=np.array(df['time'].unique()), fontsize=20)
plt.ylim(0,1)
plt.yticks(fontsize=20)
plt.grid(False)
plt.ylabel(r'R.A.')
plt.xlabel('Time (d)')
plt.title('Biofilm Culture')
# plt.savefig('Fig1B.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',|
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
 #%%
"""
Figure 1C , SI Figure S2B

Tri-Culture vs Monoculture Lines
"""

ecoli_co = [0, np.mean([327.31799833195225, 74739.97154885676,  30178.86602171233, 2272,10348]), 
            np.mean([337654.3018299466, 16377.00412162256, 102800.81758497875, 65958.14699239816,  17644.44592143285, 
                     60281.3781279211, 22315.762891329276, 226175.7079400265, 87119.63391914888, 83563]),
            np.mean([174980.53583302326, 131468.29809900734, 2172.6553144305585, 291359.03185848374, 293299.52973667695, 
                     199684.7876095191, 233148.2572287676, 90213.35726961477, 34802.62569790094, 22346.36296689037, 367732.49514154455, 188472]),
            np.mean([30012.369571568208, 136101.86383347574, 12751.512534025518, 109377.53775453306, 140409.82671136397, 355895.2917991453, 387598.1942840074, 373512.4762440534, 23562])]
ecoli_co_std = [0, stats.sem([327.31799833195225, 74739.97154885676,  30178.86602171233, 2272, 10348]), 
                stats.sem([337654.3018299466, 16377.00412162256, 102800.81758497875, 65958.14699239816,  17644.44592143285, 
                           60281.3781279211, 22315.762891329276, 226175.7079400265, 87119.63391914888,83563]),
            stats.sem([174980.53583302326, 131468.29809900734, 2172.6553144305585, 291359.03185848374, 293299.52973667695, 
                       199684.7876095191, 233148.2572287676, 90213.35726961477, 34802.62569790094, 22346.36296689037, 367732.49514154455,188472]),
            stats.sem([30012.369571568208, 136101.86383347574, 12751.512534025518, 109377.53775453306, 140409.82671136397, 355895.2917991453, 387598.1942840074, 373512.4762440534, 23562])]
ecoli = [0, np.mean([429326.5281541101, 20830.831426016266, 26769.516920138325,26992.17734452068]), 
         np.mean([111347.8469921274, 183315.5309570613, 145620.88152866834, 409437.6682564925, 134370.40829945993, 277304.7236815254]),
            np.mean([387877.2050713496, 312281.15418290574, 486196.3041085362, 286454.46008238854]),
            np.mean([337383.11888024263, 348827.4483178621, 148572.00581428283, 269781.17495603,  461549.52062919])]
ecoli_std = [0, stats.sem([429326.5281541101, 20830.831426016266, 26769.516920138325,26992.17734452068]), 
             stats.sem([111347.8469921274, 183315.5309570613, 145620.88152866834, 409437.6682564925, 134370.40829945993, 277304.7236815254]),
            stats.sem([387877.2050713496, 312281.15418290574, 486196.3041085362, 286454.46008238854]),
            stats.sem([128992.59076753035, 92435.67970619304, 182715.3318089329, 89181.74673619,116083.68712655])]
ent_co = [0, np.mean([45944.71951992576, 42819.65134236022, 14419.980738746319, 63980, 132762]), 
          np.mean([109925.465309905, 79612.13354098213, 46202.46510145668, 55182.472283692354, 179308.4799678017, 72140.62085399251, 59930.49838893081,  81622.88932344121, 125376]),
          np.mean([133616.5332796591, 23329.239782307275, 38164.97239072235, 3566.42683158818, 49841.93215699958, 33499.418502635635, 32309.77994890224,  23665.44896543966, 33135.575647364, 22937.923767081345, 115016.10230890731, 187618]),
          np.mean([55616.538048072034, 10079.947927941344, 34003.70827562903, 48037.69756361926, 57930.14145661801, 26650.75047392173, 18062.719116312706, 42176.42001744546, 14012])]
ent_co_std= [0, stats.sem([45944.71951992576, 42819.65134236022, 14419.980738746319, 63980, 132762]),
             stats.sem([109925.465309905, 79612.13354098213, 46202.46510145668, 55182.472283692354, 179308.4799678017, 
                        72140.62085399251, 59930.49838893081,  81622.88932344121, 125376]),
          stats.sem([133616.5332796591, 23329.239782307275, 38164.97239072235, 35166.42683158818, 49841.93215699958, 33499.418502635635, 
                     32309.77994890224,  23665.44896543966, 33135.575647364, 22937.923767081345, 115016.10230890731, 187618]),
          stats.sem([55616.538048072034, 10079.947927941344, 34003.70827562903, 48037.69756361926, 57930.14145661801, 26650.75047392173, 18062.719116312706, 42176.42001744546, 14012])]

ent = [0, np.mean([145719.7700434657, 132853.42173086293, 49593.08968900539, 42897.058102400624, 132762]), 
       np.mean([191166.61280572205, 96122.08029061732, 133520.86963127236, 172674.3443290808, 127962.97477912421, 105634.37227330725]),
          np.mean([190393.93985167015, 82472.47297014121, 103497.2411231636, 83590.30802445539]),
          np.mean([128992.59076753035, 92435.67970619304, 182715.3318089329, 89181.74673619,116083.68712655])]

ent_std = [0, stats.sem([145719.7700434657, 132853.42173086293, 49593.08968900539, 42897.058102400624]),
           stats.sem([191166.61280572205, 96122.08029061732, 133520.86963127236, 172674.3443290808, 127962.97477912421, 105634.37227330725]),
          stats.sem([190393.93985167015, 82472.47297014121, 103497.2411231636, 83590.30802445539]),
          stats.sem([128992.59076753035, 92435.67970619304, 182715.3318089329, 89181.74673619,116083.68712655])]

pa_co = [0, np.mean([564.8839733898792, 59753.966369545655, 76084.78222377267, 1218, 1887]),
          np.mean([304214.7366569871, 314624.12837939983, 256.68741611089274, 677.4296636374746, 216.95526415455606, 495.2815851177541, 8132.448033845329, 77540.33970772212, 298198.37500867364, 1114.08660477438,1083]),
          np.mean([510447.32744726987, 232140.6416575893, 97072.44822831418, 152974.1378831626, 544.3205115672279, 359884.46961609664, 207820.33371656644, 232250.19955164785, 109811.70338191417, 49315.8829610963, 353645]),
          np.mean([282075.78929444845, 326177.65918198874, 128033.73931436759, 398119.9328066944, 103348.1002644917, 1776.8635350431903, 26057.464556832147, 365102.55372229044, 196572])]

pa_co_std = [0, stats.sem([564.8839733898792, 59753.966369545655, 76084.78222377267, 1218, 1887]), 
          stats.sem([304214.7366569871, 314624.12837939983, 256.68741611089274, 677.4296636374746, 216.95526415455606, 
                     495.2815851177541, 8132.448033845329, 77540.33970772212, 298198.37500867364, 1114.08660477438, 1083]),
          stats.sem([510447.32744726987, 232140.6416575893, 97072.44822831418, 152974.1378831626, 544.3205115672279, 
                     359884.46961609664, 207820.33371656644, 232250.19955164785, 109811.70338191417, 49315.8829610963, 353645]),
          stats.sem([282075.78929444845, 326177.65918198874, 128033.73931436759, 398119.9328066944, 103348.1002644917, 1776.8635350431903, 26057.464556832147, 365102.55372229044, 196572])]

pa = [0, np.mean([362289.6023790988, 550.316753017643, 118188.07173208568]), 
      np.mean([126345.30028869871, 80597.25082750342, 21905.3075725621, 485315.62176974496, 3045.379556830059, 2071.840862834238, 27850.31514594865, 8115.406266916172]),
          np.mean([337295.43424540985, 398535.4141324807, 242220.39670517063, 23597.645090594877, 164344.37751006486]),
          np.mean([385489.7102260408, 439545.5021568775, 227691.2968063428, 5.02094258e+05, 452658.56209157, 520242.070775, 429652.64031865])]

pa_std = [0, stats.sem([362289.6023790988, 550.316753017643, 118188.07173208568]), 
          stats.sem([126345.30028869871, 80597.25082750342, 21905.3075725621, 485315.62176974496, 3045.379556830059, 2071.840862834238, 27850.31514594865, 8115.406266916172]),
          stats.sem([337295.43424540985, 398535.4141324807, 242220.39670517063, 23597.645090594877, 164344.37751006486]),
          stats.sem([385489.7102260408, 439545.5021568775, 227691.2968063428, 5.02094258e+05, 452658.56209157, 520242.070775, 429652.64031865])]

time = [0, 24, 48, 72, 96]
plt.figure(figsize=(3.5,4))
plt.plot(time, ecoli_co, marker='o', color='#9813DE', alpha=0.75)
plt.errorbar(time, ecoli_co, yerr=ecoli_co_std, color='#9813DE', alpha=0.75)
plt.plot(time, ecoli, linestyle='--', marker='p', color='#9813DE', alpha=0.75)
plt.errorbar(time, ecoli, yerr=ecoli_std, linestyle='--', color='#9813DE', alpha=0.75)
plt.plot(time, ent_co, marker='o', color='#008F8F', alpha=0.75)
plt.errorbar(time, ent_co, yerr=ent_co_std, color='#008F8F', alpha=0.75)
plt.plot(time, ent, linestyle='--', marker='p', color='#008F8F', alpha=0.75)
plt.errorbar(time, ent, yerr=ent_std, linestyle='--', color='#008F8F', alpha=0.75)
plt.plot(time, pa_co, marker='o', color='#E2E604', alpha=0.75)
plt.errorbar(time, pa_co, yerr=pa_co_std, color='#E2E604', alpha=0.75)
plt.plot(time, pa, linestyle='--', marker='p', color='#E2E604', alpha=0.75)
plt.errorbar(time, pa, yerr=pa_std, linestyle='--', color='#E2E604', alpha=0.75)
plt.ylabel('Volume $\mu$m$^3$', fontsize=20)
plt.xlabel('Time (hr)', fontsize = 20)
plt.grid(False)
plt.xticks(fontsize=20)
plt.yticks(ticks=np.arange(0,600000,100000), labels=np.arange(0,6,1), fontsize=20)
plt.ylim(-10000,600000)
# plt.savefig('Fig1C.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 1D
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig1D')

data_list = [[df['x0'].dropna().tolist(), df['y0'].dropna().tolist()], [df['x1'].dropna().tolist(), df['y1'].dropna().tolist()], [df['x2'].dropna().tolist(), df['y2'].dropna().tolist()]]

figure, tax = ternary.figure(scale=1.0)
figure.set_size_inches(5, 5)
color_dictionary = {r"$E. coli$":'#9813DE', r'$E. faecalis$':'#008F8F',r"$P. aeruginosa$":'#E2E604'}

tax.boundary()
tax.gridlines(multiple=0.2, color="black")
tax.left_axis_label(" $E. faecalis$", fontsize=20, offset=0.2, color='#008F8F')
tax.right_axis_label("$E. coli$", fontsize=20, offset=0.2, color='#9813DE')
tax.bottom_axis_label("$P. aeruginosa$", fontsize=20, offset=0.2, color='#E2E604')
sax = tax.get_axes()

cmap = mpl.cm.inferno
norm = mpl.colors.Normalize(vmin=0, vmax=360)

for to_plot in data_list:    
    to_plot[0] = [literal_eval(x) for x in to_plot[0]]
    print(to_plot[0])
    for point in np.arange(0,len(to_plot[0])-1):
        
        pstart = ternary.project_point((to_plot[0][point]))
        print(pstart)
        pend = ternary.project_point((to_plot[0][point+1]))
        
        pstart_m = ternary.project_point((to_plot[0][point]))
        
        pend_m = ternary.project_point((to_plot[0][point+1]))    
        
        mag = 1

        dxy = (pend - pstart) * 0.99

        lw=1.75
        print(to_plot[1][point])
        color = cmap(norm(to_plot[1][point]))

        sax.arrow(*pstart, *dxy, linewidth=lw, color=color, head_width=0.02, alpha=0.75)
        if point == 0:
            sax.scatter(*pstart, marker=(5, 1), color=color)

                
tax.ticks(axis='lbr', multiple=0.2, linewidth=1, tick_formats="%.1f", offset=0.02)
tax.get_axes().axis('off')
tax.clear_matplotlib_ticks()
# tax.savefig('Fig1D.svg', dpi=300, facecolor='w', edgecolor='b',
#             orientation='portrait', format='svg',
#             transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
tax.show()

fig, ax = plt.subplots(figsize=(8,0.5))
fig.subplots_adjust(bottom=0.5)
cmap = mpl.cm.inferno
norm = mpl.colors.Normalize(vmin=0, vmax=360)

fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax, orientation='horizontal', label='Time (h)')
# plt.savefig('Fig1D_colorbar.svg', dpi=300, facecolor='w', edgecolor='b',
#             orientation='portrait', format='svg',
#             transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2A, SI Figure 2A

! -> Time is in steps  
"""
df = pd.read_csv('Stochastic_Spatial_Model/nullbiased_dispersal.csv')

time_array = df['Time']
spc_1 = df['Species1']
spc_2 = df['Species2']
spc_3 = df['Species3']
num_iter = len(df)

plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)
# plt.savefig('Fig2A.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2B, SI Figure 2B 

! -> Time is in steps  
"""
df = pd.read_csv('Stochastic_Spatial_Model/unbiased_dispersal.csv')

time_array = df['Time']
spc_1 = df['Species1']
spc_2 = df['Species2']
spc_3 = df['Species3']
num_iter = len(df)

plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)
# plt.savefig('Fig2B.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2C, SI Figure 2C

! -> Time is in steps  
"""
df = pd.read_csv('Stochastic_Spatial_Model/biased_dispersal.csv')

time_array = df['Time']
spc_1 = df['Species1']
spc_2 = df['Species2']
spc_3 = df['Species3']
num_iter = len(df)

plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)
plt.grid(False)
# plt.savefig('Fig2C.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2D

Representastive trace of three-species community dynamics in the biofilm
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2D')

plt.plot(df.Time, df.Ec_mean, color='#9813DE', label=r'$\it{E. coli}$')
plt.fill_between(df.Time, df.Ec_mean-df.Ec_std, df.Ec_mean+df.Ec_std, color='#9813DE', label=r'$\it{E. coli}$', alpha=0.25)
plt.plot(df.Time, df.Pa_mean, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.fill_between(df.Time, df.Pa_mean-df.Pa_std, df.Pa_mean+df.Pa_std, color='#E2E604', label=r'$\it{P. aeruginosa}$', alpha=0.25)
plt.plot(df.Time, df.Ef_mean, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.fill_between(df.Time, df.Ef_mean-df.Ef_std, df.Ef_mean+df.Ef_std, color='#008F8F', label=r'$\it{E. faecalis}$', alpha=0.25)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Biovolume (um^3)',fontsize=20)
plt.grid(False)
# plt.savefig('Fig2D.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2E
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2E')

x = df['v']
y = df['dv']

slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(x, y)
sns.regplot(x=x, y=y, color = 'black', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope3,intercept3, r_value3, p_value3)})
plt.ylabel(r'dv_pa14/dt', fontsize=20)
plt.xlabel(r'v_pa14', fontsize=20)
plt.xticks(fontsize=20, rotation=-90)
plt.yticks(fontsize=20)
plt.hlines(0, xmin=-10000, xmax=630000, color='red', linestyles='--')
plt.xlim(-10000,630000)
plt.ylim(-22000,22000)
plt.legend()
# plt.savefig('Fig2E.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2F
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig2F')

Pa = df['Pa']
Ec = df['Ec']
Ef = df['Ef']

slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(Pa, Pa)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(Pa, Ef)
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(Pa, Ec)

sns.regplot(x=Pa, y=Pa, color = '#E2E604', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope3,intercept3, r_value3, p_value3)})
sns.regplot(x=Pa, y=Ef, color = '#008F8F', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope2,intercept2, r_value2, p_value2)})
sns.regplot(x=Pa, y=Ec, color = '#9813DE', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope1,intercept1, r_value1, p_value1)})

plt.ylabel(r'dv_species/dt', fontsize=20)
plt.xlabel(r'dv_pa14/dt', fontsize=20)
plt.grid(False)
plt.xticks(fontsize=20, rotation=-90)
plt.yticks(fontsize=20)
plt.legend()
# plt.savefig('Fig2F.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2G
"""

plt.rcParams["figure.figsize"] = (5,6)
data=pd.read_excel('SI_Data.xlsx', sheet_name='Fig2G')
plt.plot(data.X, data.PA14/np.sum(data.PA14), color='#E2E604')
plt.plot(data.X, data.Ent/np.sum(data.Ent), color='#008F8F')
plt.plot(data.X, data.Ecoli/np.sum(data.Ecoli), color='purple')
plt.scatter(data.X, data.PA14/np.sum(data.PA14), color='#E2E604')
plt.scatter(data.X, data.Ent/np.sum(data.Ent), color='#008F8F')
plt.scatter(data.X, data.Ecoli/np.sum(data.Ecoli), color='#9813DE')

plt.ylim(0,0.06)
plt.ylabel('Frequency', fontsize=20)
plt.xlabel('Pseudomonas Cell Packing', fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid(False)
# plt.savefig('Fig2G.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 2H

Use the stochastic spatial model make_heatmap script to generate data prior to running this code. 
"""
ratio_list = list(np.flip(np.arange(0,1.1,.1)))
growth_list = list(np.arange(.03,.06003,.003))
experiment_title = 'phase_diagram'

result = []

for ratio in ratio_list:
    ratio_cake = []
    for growth in growth_list:
        to_average = []        
        for rep in np.arange(0,20):
            df = pd.read_excel('Stochastic_Spatial_Model/phase_diagram/'+experiment_title+'_p3_'+str(ratio)+'_'+str(growth)+'_'+str(rep)+'_20000.xlsx')
            if df.Species1[len(df)-1] != 0 and df.Species2[len(df)-1] != 0 and df.Species3[len(df)-1] != 0:
                to_average.append(1)
            else:
                to_average.append(0)
        ratio_cake.append(np.mean(to_average))
    result.append(ratio_cake)

ax = sns.heatmap(result, square=True, cmap=Greys_6_r.mpl_colormap)
for _, spine in ax.spines.items():
    spine.set_visible(True)
        
# plt.savefig('Fig2H.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show
#%%
"""
Figure 3A

Stochastic spatial model rep ternary plot

Use the generate_ternary_plot code in the stochastic_spatial_model directory to generate data prior to running this section
"""

figure, tax = ternary.figure(scale=1.0)
figure.set_size_inches(5, 5)
tax.boundary()
tax.gridlines(multiple=0.2, color="black")
sax = tax.get_axes()

#grab one of the runs with a 0.33,0.33,0.33 starting condition for this
volume_dataframe2 = pd.read_excel('Stochastic_Spatial_Model/biased_dispersal_vector_plot_cont_1000.xlsx')
norm = mpl.colors.Normalize(vmin=0, vmax=1000)

for i in np.arange(0,1000):    

    cmap= Inferno_20.mpl_colormap    
    color = cmap(norm(i))
    
    pstart = ternary.project_point((volume_dataframe2['Species2'][i],volume_dataframe2['Species1'][i],volume_dataframe2['Species3'][i]))
    pend = ternary.project_point((volume_dataframe2['Species2'][i+1],volume_dataframe2['Species1'][i+1],volume_dataframe2['Species3'][i+1]))
    pstart_m = ternary.project_point((volume_dataframe2['Species2'][i],volume_dataframe2['Species1'][i],volume_dataframe2['Species3'][i]))
    pend_m = ternary.project_point((volume_dataframe2['Species2'][i+1],volume_dataframe2['Species1'][i+1],volume_dataframe2['Species3'][i+1]))    
    
    mag = 1
    stepsize = np.sqrt((pend[0]-pstart[0])**2+(pend[1]-pstart[1])**2)
    dxy = (pend - pstart)*.95
    lw=1
    sax.arrow(*pstart, *dxy, linewidth=lw, color=color, head_width=0.02, alpha=.25)
    
    if i == 0:
        sax.scatter(*pstart, marker=(5, 1), color='black')

tax.get_axes().axis('off')
tax.clear_matplotlib_ticks()
# tax.savefig('single_run_trace.svg', dpi=300, facecolor='w', edgecolor='b',
#             orientation='portrait', format='svg',
#             transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
tax.show()

fig, ax = plt.subplots(figsize=(8,0.5))
fig.subplots_adjust(bottom=0.5)
cmap=Inferno_20.mpl_colormap
norm = mpl.colors.Normalize(vmin=0, vmax=1000)
# fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
#              cax=ax, orientation='horizontal', label='Time (h)')
# plt.savefig('hot_colorbar.svg', dpi=300, facecolor='w', edgecolor='b',
#             orientation='portrait', format='svg',
#             transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 3B

Stochastic spatial model ternary plot

Use the generate_ternary_plot_parallel_pool code in the stochastic_spatial_model directory to generate data
"""
file_list = glob.glob('Stochastic_Spatial_Model/ternary_plot/*.xlsx', recursive=False)

#Define the range and increments
start = 0
stop = 1
increment = 0.1

# Calculate the number of steps
num_steps = int((stop - start) / increment) + 1
# Generate all combinations
combinations = list(product(range(num_steps), repeat=3))
# Convert combinations to actual values
result = [[start + c[0] * increment, start + c[1] * increment, start + c[2] * increment] for c in combinations]

# Print the result
for combination in result:
    print(combination)

averaged = {'Species1':[], 'Species2':[], 'Species3':[]}

averaged_dispersal = {'Species1':[], 'Species2':[], 'Species3':[]}

result_final = []

for value in result:
    if round(np.sum(value), 2) == 1:
        result_final.append(value)
        
for i in np.arange(0,len(result_final)):
    
    if 1 in result_final:
        print(result_final)
    if 1 not in result_final:
            
        start_to_average_1 = []
        end_to_average_1 = []
        start_to_average_2 = []
        end_to_average_2 = []
        start_to_average_3 = []
        end_to_average_3 = []
        step = increment
        
        for file in file_list:    
            volume_dataframe= pd.read_excel(file)
            quadrant = volume_dataframe.drop(volume_dataframe[(volume_dataframe['Species1'] <= result_final[i][0]-step) |
                                               (volume_dataframe['Species1'] > result_final[i][0]) |
                                               (volume_dataframe['Species3'] <= result_final[i][1]-step) |
                                               (volume_dataframe['Species3'] > result_final[i][1]) |
                                               (volume_dataframe['Species2'] <= result_final[i][2]-step) |
                                               (volume_dataframe['Species2'] > result_final[i][2])| (volume_dataframe['Dispersal'] == 1)].index)
            
            look_length = 90
            for item in quadrant.index:
                if (item+look_length in volume_dataframe.index) and (quadrant['Iteration'][item] == volume_dataframe['Iteration'][item+look_length]) and (1 not in volume_dataframe['Dispersal'][item:item+look_length] ):
                    start_to_average_1.append(quadrant['Species1'][item])
                    end_to_average_1.append(volume_dataframe['Species1'][item+look_length])
                    start_to_average_2.append(quadrant['Species2'][item])
                    end_to_average_2.append(volume_dataframe['Species2'][item+look_length])    
                    start_to_average_3.append(quadrant['Species3'][item])
                    end_to_average_3.append(volume_dataframe['Species3'][item+look_length])    
            
        averaged['Species1'].append([np.mean(start_to_average_1), np.mean(end_to_average_1)])
        averaged['Species2'].append([np.mean(start_to_average_2), np.mean(end_to_average_2)])
        averaged['Species3'].append([np.mean(start_to_average_3), np.mean(end_to_average_3)])
        
for i in np.arange(0,len(result_final)):
    
    if 1 in result_final:
        print(result_final)
    if 1 not in result_final:
            
        start_to_average_1 = []
        end_to_average_1 = []
        start_to_average_2 = []
        end_to_average_2 = []
        start_to_average_3 = []
        end_to_average_3 = []
        step = increment
        
        for file in file_list:
            volume_dataframe= pd.read_excel(file)        
            quadrant = volume_dataframe.drop(volume_dataframe[(volume_dataframe['Species1'] <= result_final[i][0]-step) |
                                               (volume_dataframe['Species1'] > result_final[i][0]) |
                                               (volume_dataframe['Species3'] <= result_final[i][1]-step) |
                                               (volume_dataframe['Species3'] > result_final[i][1]) |
                                               (volume_dataframe['Species2'] <= result_final[i][2]-step) |
                                               (volume_dataframe['Species2'] > result_final[i][2]) | (volume_dataframe['Dispersal'] == 0)].index)
            
            look_length = 1
            for item in quadrant.index:
                if (item+look_length in volume_dataframe.index) and (quadrant['Iteration'][item] == volume_dataframe['Iteration'][item+look_length]):
                    start_to_average_1.append(quadrant['Species1'][item])
                    end_to_average_1.append(volume_dataframe['Species1'][item+look_length])
                    start_to_average_2.append(quadrant['Species2'][item])
                    end_to_average_2.append(volume_dataframe['Species2'][item+look_length])    
                    start_to_average_3.append(quadrant['Species3'][item])
                    end_to_average_3.append(volume_dataframe['Species3'][item+look_length])    
            
        averaged_dispersal['Species1'].append([np.mean(start_to_average_1), np.mean(end_to_average_1)])
        averaged_dispersal['Species2'].append([np.mean(start_to_average_2), np.mean(end_to_average_2)])
        averaged_dispersal['Species3'].append([np.mean(start_to_average_3), np.mean(end_to_average_3)])

volume_dataframe_average = pd.DataFrame(data=averaged)
volume_dataframe_average_disp = pd.DataFrame(data=averaged_dispersal)

figure, tax = ternary.figure(scale=1.0)
figure.set_size_inches(5, 5)
color_dictionary = {r"$E. coli$":'#9813DE', r'$E. faecalis$':'#008F8F',r"$P. aeruginosa$":'#E2E604'}

tax.boundary()
tax.gridlines(multiple=0.2, color="black")

for i in np.arange(0,len(volume_dataframe_average)):
    if 1 not in [volume_dataframe_average['Species1'][i][0],volume_dataframe_average['Species2'][i][0],volume_dataframe_average['Species3'][i][0]]:    
        
        pstart = ternary.project_point((volume_dataframe_average['Species2'][i][0],volume_dataframe_average['Species1'][i][0],volume_dataframe_average['Species3'][i][0]))
        pend = ternary.project_point((volume_dataframe_average['Species2'][i][1],volume_dataframe_average['Species1'][i][1],volume_dataframe_average['Species3'][i][1]))
        pstart_m = ternary.project_point((volume_dataframe_average['Species2'][i][0],volume_dataframe_average['Species1'][i][0],volume_dataframe_average['Species3'][i][0]))
        pend_m = ternary.project_point((volume_dataframe_average['Species2'][i][1],volume_dataframe_average['Species1'][i][1],volume_dataframe_average['Species3'][i][1]))    

        mag = 1
        dxy = (pend - pstart)*.25
        lw=1
        sax.arrow(*pstart, *dxy, linewidth=lw, color='black', head_width=0.025, alpha=.5)
        
for i in np.arange(0,len(volume_dataframe_average_disp)):
    if 1 not in [volume_dataframe_average_disp['Species1'][i][0],volume_dataframe_average_disp['Species2'][i][0],volume_dataframe_average_disp['Species3'][i][0]]:    

        pstart = ternary.project_point((volume_dataframe_average_disp['Species2'][i][0],volume_dataframe_average_disp['Species1'][i][0],volume_dataframe_average_disp['Species3'][i][0]))
        pend = ternary.project_point((volume_dataframe_average_disp['Species2'][i][1],volume_dataframe_average_disp['Species1'][i][1],volume_dataframe_average_disp['Species3'][i][1]))
        pstart_m = ternary.project_point((volume_dataframe_average_disp['Species2'][i][0],volume_dataframe_average_disp['Species1'][i][0],volume_dataframe_average_disp['Species3'][i][0]))
        pend_m = ternary.project_point((volume_dataframe_average_disp['Species2'][i][1],volume_dataframe_average_disp['Species1'][i][1],volume_dataframe_average_disp['Species3'][i][1]))    
        
        mag = 1       
        stepsize = np.sqrt((pend[0]-pstart[0])**2+(pend[1]-pstart[1])**2)
        dxy = (pend - pstart)
        lw=1
        sax.arrow(*pstart, *dxy, linewidth=lw, color='red', head_width=0.025, alpha=.5)

tax.get_axes().axis('off')
tax.clear_matplotlib_ticks()
# tax.savefig('Figure3B.svg', dpi=300, facecolor='w', edgecolor='b',
#             orientation='portrait', format='svg',
#             transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
tax.show()
#%%
"""
Figure 4B

Autocorrelation 
"""
plot_dict_post = {'PA14':[21.305538941558908, 5.6511008706908745, 5.0137953552405925, 4.398348976868345, 5.20582174512063, 8.403290287099342, 8.661097047023734, 4.173319267223729, 7.478522930807225, 6.939765424541318, 18.384609714330544],
                  'Ent':[1.6931792126177123, 1.8544360309360752, 1.809556537246724, 1.7076370669118321, 4.562035820904487, 5.496435664188826, 1.5519221187709673, 1.5740291047592825, 4.618614028430873, 1.7598148440561323, 1.8553295822499383],
                  'Ecoli':[1.5301727733624557, 1.8394492982430997, 1.8163981946988546, 1.6894066337125933, 1.7227547206431313, 1.7121541346500035, 1.5463180841056863, 1.5627974860670502, 1.5494989071190346, 1.6817066776283336, 0.9197616522573512]
                  }

plt.figure(figsize=(4,3))
plt.boxplot(plot_dict_post.values())
plt.scatter(np.ones(len(plot_dict_post[list(plot_dict_post.keys())[0]])), plot_dict_post[list(plot_dict_post.keys())[0]], alpha=0.75, color='#E2E604')
plt.scatter(np.ones(len(plot_dict_post[list(plot_dict_post.keys())[1]]))*2, plot_dict_post[list(plot_dict_post.keys())[1]], alpha=0.75, color='#008F8F')
plt.scatter(np.ones(len( plot_dict_post[list(plot_dict_post.keys())[2]]))*3, plot_dict_post[list(plot_dict_post.keys())[2]], alpha=0.75, color='#9813DE')
plt.xticks(ticks=[1,2, 3], labels=plot_dict_post.keys(), fontsize=18)
plt.yticks(fontsize=18)
plt.grid(False)
plt.ylabel('Autocorrelation Length ($\mu$m$^3$)', fontsize=18)
# plt.savefig('Fig4B.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(stats.mannwhitneyu(plot_dict_post['Ecoli'], plot_dict_post['Ent']))
print(stats.mannwhitneyu(plot_dict_post['PA14'], plot_dict_post['Ent']))
print(stats.mannwhitneyu(plot_dict_post['PA14'], plot_dict_post['Ecoli']))
#%%
"""
Figure 4C

Distance to Nearest Neighbor
"""

plot_dict_post = {'Pa_Ef':[16.796813975764763, 2.5752512891003794, 1.9350533052442773, 6.14428108504588, 0.8536672812931766, 6.875872974295542, 5.745158642778306, 4.37020926281347, 14.671940442099613, 3.2222448599935145, 3.8845575009545494, 0.3738369610453203],
                  'Pa_Ec':[56.27290588895764, 26.019403457337226, 21.601771031091314, 18.747747202686117, 31.374000296292976, 18.43771837207643, 18.462670814974928, 32.40908177290209, 21.604464317904053, 3.203322742605807, 2.3516076875469727, 1.1930308873578803],
                  'Ec_Ef':[10.132951476744786, 1.062532964894646, 0.7455491632806658, 4.270373677541325, 0.9280165893326213, 7.073060621291107, 3.675005094756839, 2.9271194163854073, 15.701905517496282, 2.6575859119010916, 2.846989390016894,0.4504739043519803],}

plt.figure(figsize=(4,3))
plt.boxplot(plot_dict_post.values())
plt.plot(np.ones(len(plot_dict_post[list(plot_dict_post.keys())[0]]))*1, plot_dict_post[list(plot_dict_post.keys())[0]], alpha=0.5, c='#E2E604', markerfacecoloralt='#008F8F', fillstyle='left',  marker='.', markeredgecolor='None', linestyle='', markersize=15)
plt.plot(np.ones(len( plot_dict_post[list(plot_dict_post.keys())[1]]))*2, plot_dict_post[list(plot_dict_post.keys())[1]], alpha=0.5, c='#E2E604', markerfacecoloralt='#9813DE', fillstyle='left',  marker='.', markeredgecolor='None', linestyle='', markersize=15)
plt.plot(np.ones(len(plot_dict_post[list(plot_dict_post.keys())[2]]))*3, plot_dict_post[list(plot_dict_post.keys())[2]], alpha=0.5, c='#9813DE', markerfacecoloralt='#008F8F', fillstyle='left',  marker='.', markeredgecolor='None', linestyle='', markersize=15)
plt.xticks(ticks=[1,2,3], labels=plot_dict_post.keys(), fontsize=18)
plt.yticks(fontsize=18)
plt.ylabel('Length ($\mu$m)', fontsize=18)
plt.grid(False)
# plt.savefig('Fig4C.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
print(stats.mannwhitneyu(plot_dict_post['Pa_Ec'], plot_dict_post['Pa_Ef']))
print(stats.mannwhitneyu(plot_dict_post['Ec_Ef'], plot_dict_post['Pa_Ef']))
print(stats.mannwhitneyu(plot_dict_post['Ec_Ef'], plot_dict_post['Pa_Ec']))
#%%
"""
Figure 5B 

Box and Whisker Plot, Sock Plot
"""
plt.figure(figsize=(2,3))
plot_dict = {r'WT':[0.732701392,0.855511269,0.605011441,0.8603337031017434,0.674319833], 
             '$\Delta bqsS$':[0.8604058699519027, 0.9975449725269114, 0.9744276316075082, 0.975276188580307, 0.9935968885138542]}

plt.boxplot(plot_dict.values(), widths=0.5)
plt.scatter(np.ones(len(plot_dict['WT'])), plot_dict['WT'], alpha=0.75, color='black', s=100)
plt.scatter(np.ones(len(plot_dict['$\Delta bqsS$']))*2, plot_dict['$\Delta bqsS$'], marker='D', alpha=0.75, color='green', s=100)
plt.ylim(0.5,1)

plt.grid(False)
plt.xticks(ticks=[1,2,], labels=plot_dict.keys(), fontsize=18)
plt.yticks(ticks = np.arange(0.5, 1.01, 0.1), fontsize=18)
plt.ylabel(r'Relative Abundance', fontsize=18)
# plt.savefig('Fig3B_1.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

print(stats.mannwhitneyu(plot_dict['WT'], plot_dict['$\Delta bqsS$']))

plt.figure(figsize=(3,4))

plot_dict = {r'WT':[[0.732701392, 0.000642044, 0.266656565], [0.855511269,0.064441976,0.080046755], [0.605011441, 0.393689205, 0.001299354], [0.674319833, 0.000281737, 0.32539843]], 
             '$\Delta bqsS$':[[0.8604058699519027, 0.09367418792941473, 0.04591994211868235], [0.9975449725269114, 0.0004435961750302079, 0.0020114312980582445], [0.9744276316075082, 0.02450629927006042, 0.0010660691224313199], [0.975276188580307, 0.016741727878660777, 0.007982083541032315], [0.9935968885138542, 0.00354877813774713, 0.002854333348398737]]}

data = [[],[],[]]
for key in list(plot_dict.keys()):
    pa = []
    ef = []
    ec = []
    for i in plot_dict[key]:
        pa.append(i[0])
        ef.append(i[1])
        ec.append(i[2])

    pa = np.mean(pa)
    ef = np.mean(ef)
    ec = np.mean(ec)
    total = pa+ef+ec
    pa = pa/total
    ef = ef/total
    ec = ec/total
    if key == 'WT':
        data[0].append(pa)
        data[1].append(ef)
        data[2].append(ec)
    elif key == '$\Delta bqsS$':
        data[0].append(pa)
        data[1].append(ef)
        data[2].append(ec)
        
        
plt.bar([0.25,0.75], data[0][:], color='#E2E604', label=r'${\itP. aeruginosa}$', alpha=0.75, width=0.2)
plt.bar([0.25,0.75], data[1][:], color='#008F8F', bottom=np.array(data[0][:]), label=r'${\itE. faecalis}$', alpha=0.75,  width=0.2)
plt.bar([0.25,0.75], data[2][:], color='#9813DE', bottom=np.array(data[0][:])+np.array(data[1][:]), label=r'${\itE. coli}$', alpha=0.75,  width=0.2)
plt.ylim(0,1)
plt.xlim(0,1)
plt.grid(False)
plt.ylabel(r'Relative Abundance', fontsize=18)
# plt.savefig('Fig3B_2.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 5C

Time average biovolume
"""
plt.figure(figsize=(4,3))
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig5C_1')
df_wt =  pd.read_excel('SI_Data.xlsx', sheet_name='Fig5C_0') 

df = df.sort_values(by='Time', ascending=True)
df_wt = df_wt.sort_values(by='Time', ascending=True)
            
times = df['Time']
times= np.arange(0,372,24)

def time_average(df, subject, ranger = 72):
    time_averaged_m = []
    time_averaged_std = []
    for t2 in times:
        df2 = df.drop(df[(df.Time > t2+ranger)|(df.Time < t2-ranger)].index, inplace = False)
        time_averaged_m.append(np.mean(df2[subject]))
        time_averaged_std.append(np.std(df2[subject], ddof=1))
    return((time_averaged_m, time_averaged_std))

pa14_mean, pa14_std = time_average(df_wt, 'PA_WT')
pa14_mean_del, pa14_std_del = time_average(df, 'PA_bqsS')

c_pa='#E2E604'
c_ec = '#9813DE'
c_ent = '#008F8F'

plt.scatter(df_wt['Time'],df_wt['PA_WT'], color='black', alpha=0.25)
plt.plot(times, pa14_mean, 'black', linestyle='--')
plt.fill_between(times, np.asarray(pa14_mean) - np.asarray(pa14_std), np.asarray(pa14_mean) + np.asarray(pa14_std), color = 'black', alpha=0.15) 

plt.scatter(df['Time'],df['PA_bqsS'], color='green', alpha=0.25, marker='D')
plt.plot(times, pa14_mean_del, 'green', linestyle='--')
plt.fill_between(times, np.asarray(pa14_mean_del) - np.asarray(pa14_std_del), np.asarray(pa14_mean_del) + np.asarray(pa14_std_del), color = 'green', alpha=0.15) 
 
plt.xlim(0,360)
plt.xticks(fontsize=15)
plt.ylim(-50000,700000)
plt.yticks(fontsize=15)
plt.grid(False)
# plt.savefig('Fig3C.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 5D
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig5D').T.reset_index()
new_header = df.iloc[0] 
df = df[1:] 
df.columns = new_header 

df.drop(columns=['Index'])
my_pal = {'Pa_mutant':'#E2E604', 'Pa_wt':'#E2E604', 
          'Ec_mutant':'#9813DE', 'Ec_wt':'#9813DE', 
          'Ef_mutant':'#008F8F', 'Ef_wt':'#008F8F'}

#loop for running one sided test
for inp in list(df.columns)[1:]:
    arr = df[inp].dropna()
    print(inp)
    print(stats.wilcoxon(arr, alternative='less', zero_method = 'zsplit'))

sns.boxplot(df, palette =my_pal)
sns.stripplot(df, color='black')
plt.xticks(rotation=45)
plt.hlines(y=0, xmin=-0.5, xmax=5.5, linewidth=2, color='r', linestyle='--')
plt.ylabel('Growth of Invaders', fontsize=16)
plt.grid(False)
# plt.savefig('Fig3D.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 6A
"""
to_plot = []
to_plot2=[]

plot_dict_pre = {r'WT':[216041.71467526804, 376725.7601454683, 285949.8442536336,309289.94785010454, 367276.81454676663, 415122.3440446011], 
                 '$\Delta bqsS$':[397859.0728222483, 341144.0556416029, 521094.96936566086, 613590.5673621687, 668634.2710488527, 501324.46009506844, 473167.6972423181]}


plt.figure(figsize=(3,4))
plt.boxplot(plot_dict_pre.values(), widths=0.3)
plt.xticks(ticks=[1,2], labels=plot_dict_pre.keys(), fontsize=20)

plt.scatter(np.ones(len(plot_dict_pre['WT'])), plot_dict_pre['WT'], alpha=0.75, color='black', s=60)
plt.scatter(np.ones(len(plot_dict_pre['$\Delta bqsS$']))*2, plot_dict_pre['$\Delta bqsS$'], marker='D', alpha=0.75, color='green', s=60)
plt.yticks(ticks=np.arange(200000, 800000, 100000), labels=np.arange(2,8,1), fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel(r'Volume ($um^3$)', fontsize=20)
plt.xlabel(r'$P. aeruginosa$ Strain', fontsize = 20)
plt.title('Parent Chamber', fontsize = 20)
plt.grid(False)
# plt.savefig('Fig4A.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
print(stats.mannwhitneyu(plot_dict_pre['WT'], plot_dict_pre['$\Delta bqsS$']))
#%%
"""
Figure 6B
"""
plot_dict_post = {r'WT':[42167.54187893376, 138331.075653775, 197788.37685636582, 55864.185694406304, 129865.69168123294,8784.902380816098], 
                  '$\Delta bqsS$':[29323.739031003217, 97789.70671581058, 18374.87530246589, 10775.614615360977, 1439.9204144279065, 3611.4793872326204, 4521.823017211099]}


plt.figure(figsize=(3,4))
plt.boxplot(plot_dict_post.values(), widths=0.3)
plt.xticks(ticks=[1,2], labels=plot_dict_pre.keys(), fontsize=20)

plt.scatter(np.ones(len(plot_dict_post['WT'])), plot_dict_post['WT'], alpha=0.75, color='black', s=60)
plt.scatter(np.ones(len(plot_dict_post['$\Delta bqsS$']))*2, plot_dict_post['$\Delta bqsS$'], marker='D', alpha=0.75, color='green', s=60)

plt.yticks(ticks=np.arange(0, 250000, 50000), labels=np.arange(0,2.5,0.5), fontsize=20)
plt.ylabel(r'Volume ($um^3$)', fontsize=20)
plt.xlabel(r'$P. aeruginosa$ Strain', fontsize = 20)
plt.title('Progeny Chamber', fontsize = 20)
plt.grid(False)
# plt.savefig('Fig4B.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
print(stats.mannwhitneyu(plot_dict_post['WT'], plot_dict_post['$\Delta bqsS$']))
#%%
"""
Figure 6C
"""
plot_dict_post = {'WT':np.log10([1e6,1.2e6,1.6e5,0.2e6,1.2e6,0.4e5,0.6e6,0.6e5]),
                  '$\Delta bqsS$':np.log10([3.2e5,0.4e5,0.4e5,0.4e5,0.4e5,0.2e5,0.8e5,0.4e5])
                  }

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict_post.values(), widths=0.3)
plt.scatter(np.ones(len(plot_dict_post[list(plot_dict_post.keys())[0]])), plot_dict_post[list(plot_dict_post.keys())[0]], alpha=0.75, color='black', s=60)
plt.scatter(np.ones(len(plot_dict_post[list(plot_dict_post.keys())[1]]))*2, plot_dict_post[list(plot_dict_post.keys())[1]], alpha=0.75, marker='D', color='green', s=60)
plt.xticks(ticks=[1,2], labels=plot_dict_post.keys(), fontsize=20)
plt.yticks(np.arange(4,7.5,0.75), fontsize=20)
plt.ylabel('log$_{10}$ CFUs / ($\mu$l$\cdot$hr)', fontsize=20)
plt.xlabel('$P. aeruginosa$ strain', fontsize=20)
plt.title('Dispersal Rate', fontsize = 20)
plt.grid(False)
# plt.savefig('Fig4C.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
print(stats.mannwhitneyu(plot_dict_post['WT'], plot_dict_post['$\Delta bqsS$']))
#%%
"""
Figure 6E
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig6E')
my_pal = {'WT_Pre':'black', 'WT_Post':'black', 'bqsS_Pre':'green', 'bqsS_Post':'green'}

#this one is reported in paper
#colonized chamber
print(stats.mannwhitneyu(df['WT_Post'].dropna(), df['bqsS_Post']))


plot_dict_post = {'WT':df['WT_Post'].dropna(),
                  '$\Delta bqsS$':df['bqsS_Post'].dropna()
                  }

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict_post.values(), widths=0.3)
plt.scatter(np.ones(len(plot_dict_post[list(plot_dict_post.keys())[0]])), plot_dict_post[list(plot_dict_post.keys())[0]], alpha=0.75, color='black', s=60)
plt.scatter(np.ones(len(plot_dict_post[list(plot_dict_post.keys())[1]]))*2, plot_dict_post[list(plot_dict_post.keys())[1]], alpha=0.75, marker='D', color='green', s=60)
plt.xticks(ticks=[1,2], labels=plot_dict_post.keys(), fontsize=20)
plt.yticks(np.arange(0,.9,.2), fontsize=20)
plt.xlabel('$P. aeruginosa$ strain', fontsize=20)
plt.ylabel('Biofilm Joint Local Density', fontsize=20)
plt.grid(False)
# plt.savefig('Fig4E.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
Figure 6F
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='Fig6F')
my_pal = {'WT_Pre':'black', 'WT_Post':'black', 'bqsS_Pre':'green', 'bqsS_Post':'green'}

#this one is reported in paper
#colonized chamber
print(stats.mannwhitneyu(df['WT_Post'].dropna(), df['bqsS_Post']))

plot_dict_post = {'WT':df['WT_Post'].dropna(),
                  '$\Delta bqsS$':df['bqsS_Post'].dropna()
                  }

plt.figure(figsize=(3,4))
plt.boxplot(plot_dict_post.values(), widths=0.3)
plt.scatter(np.ones(len(plot_dict_post[list(plot_dict_post.keys())[0]])), plot_dict_post[list(plot_dict_post.keys())[0]], alpha=0.75, color='black', s=60)
plt.scatter(np.ones(len(plot_dict_post[list(plot_dict_post.keys())[1]]))*2, plot_dict_post[list(plot_dict_post.keys())[1]], alpha=0.75, marker='D', color='green', s=60)
plt.xticks(ticks=[1,2], labels=plot_dict_post.keys(), fontsize=20)
plt.yticks(np.arange(3,19,3), fontsize=20)
plt.ylabel('Biofilm Height', fontsize=20)
plt.xlabel('$P. aeruginosa$ strain', fontsize=20)
plt.grid(False)
# plt.savefig('Fig4D.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S2A

Biofilm Sock Plots form co culture
"""

df_og = pd.read_excel('SI_Data.xlsx', sheet_name='FigS2A_1')
df = df_og.groupby('time').mean().reset_index()

t = np.arange(0,len(df['time'].unique()))

width=0.9
plt.bar(t, df.pa, width=width, color='#E2E604', label=r'${\itP. aeruginosa}$', alpha=0.75)
plt.bar(t, df.ec, width=width, color='#9813DE', bottom=np.array(df.pa), label=r'${\itE. coli}$', alpha=0.75)
plt.scatter(np.ones(4)*0,np.ones(4)*0.5, facecolors='none', edgecolors='black', alpha=1)
plt.scatter(df_og['time'], df_og['pa'], facecolors='none', edgecolors='black', alpha=1)
plt.xticks(ticks=t, labels=np.array(df['time'].unique()), fontsize=20)
plt.ylim(0,1)
plt.yticks(fontsize=20)
plt.grid(False)
plt.ylabel(r'R.A.')
plt.xlabel('Time (d)')
plt.title('Biofilm Culture')
# plt.savefig('FigS2A_1.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

df_og = pd.read_excel('SI_Data.xlsx', sheet_name='FigS2A_2')
df = df_og.groupby('time').mean().reset_index()

t = np.arange(0,len(df['time'].unique()))

width=0.9
plt.bar(t, df.pa, width=width, color='#E2E604', label=r'${\itP. aeruginosa}$', alpha=0.75)
plt.bar(t, df.ef, width=width, color='#008F8F', bottom=np.array(df.pa), label=r'${\itE. faecalis}$', alpha=0.75)
plt.scatter(np.ones(4)*0,np.ones(4)*0.5, facecolors='none', edgecolors='black', alpha=1)
plt.scatter(df_og['time'], df_og['pa'], facecolors='none', edgecolors='black', alpha=1)
plt.xticks(ticks=t, labels=np.array(df['time'].unique()), fontsize=20)
plt.ylim(0,1)
plt.yticks(fontsize=20)
plt.grid(False)
plt.ylabel(r'R.A.')
plt.xlabel('Time (d)')
plt.title('Biofilm Culture')
# plt.savefig('FigS2A_2.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

df_og = pd.read_excel('SI_Data.xlsx', sheet_name='FigS2A_3')
df = df_og.groupby('time').mean().reset_index()

t = np.arange(0,len(df['time'].unique()))

width=0.9
plt.bar(t, df.ef, width=width, color='#008F8F', label=r'${\itE. faecalis}$', alpha=0.75)
plt.bar(t, df.ec, width=width, color='#9813DE', bottom=np.array(df.ef), label=r'${\itE. coli}$', alpha=0.75)
plt.scatter(np.ones(4)*0,np.ones(4)*0.5, facecolors='none', edgecolors='black', alpha=1)
plt.scatter(df_og['time'], df_og['ef'], facecolors='none', edgecolors='black', alpha=1)
plt.xticks(ticks=t, labels=np.array(df['time'].unique()), fontsize=20)
plt.ylim(0,1)
plt.yticks(fontsize=20)
plt.grid(False)
plt.ylabel(r'R.A.')
plt.xlabel('Time (d)')
plt.title('Biofilm Culture')
# plt.savefig('FigS2A_3.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S2B

Planktonic Sock Plots from co culture 
"""

df_og = pd.read_excel('SI_Data.xlsx', sheet_name='FigS2B_1')
df = df_og.groupby('time').mean().reset_index()

t = np.arange(0,len(df['time'].unique()))

width=0.9
plt.bar(t, df.pa, width=width, color='#E2E604', label=r'${\itP. aeruginosa}$', alpha=0.75)
plt.bar(t, df.ec, width=width, color='#9813DE', bottom=np.array(df.pa), label=r'${\itE. coli}$', alpha=0.75)
plt.scatter(np.ones(4)*0,np.ones(4)*0.5, facecolors='none', edgecolors='black', alpha=1)
plt.scatter(df_og['time']/24, df_og['pa'], facecolors='none', edgecolors='black', alpha=1)
plt.xticks(ticks=t, labels=np.array(df['time'].unique()), fontsize=20)
plt.ylim(0,1)
plt.yticks(fontsize=20)
plt.grid(False)
plt.ylabel(r'R.A.')
plt.xlabel('Time (d)')
plt.title('Biofilm Culture')
# plt.savefig('FigS2B_1.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

df_og = pd.read_excel('SI_Data.xlsx', sheet_name='FigS2B_2')
df = df_og.groupby('time').mean().reset_index()

t = np.arange(0,len(df['time'].unique()))

width=0.9
plt.bar(t, df.pa, width=width, color='#E2E604', label=r'${\itP. aeruginosa}$', alpha=0.75)
plt.bar(t, df.ef, width=width, color='#008F8F', bottom=np.array(df.pa), label=r'${\itE. faecalis}$', alpha=0.75)
plt.scatter(np.ones(4)*0,np.ones(4)*0.5, facecolors='none', edgecolors='black', alpha=1)
plt.scatter(df_og['time']/24, df_og['pa'], facecolors='none', edgecolors='black', alpha=1)
plt.xticks(ticks=t, labels=np.array(df['time'].unique()), fontsize=20)
plt.ylim(0,1)
plt.yticks(fontsize=20)
plt.grid(False)
plt.ylabel(r'R.A.')
plt.xlabel('Time (d)')
plt.title('Biofilm Culture')
# plt.savefig('FigS2B_1.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

df_og = pd.read_excel('SI_Data.xlsx', sheet_name='FigS2B_3')
df = df_og.groupby('time').mean().reset_index()

t = np.arange(0,len(df['time'].unique()))

width=0.9
plt.bar(t, df.ef, width=width, color='#008F8F', label=r'${\itE. faecalis}$', alpha=0.75)
plt.bar(t, df.ec, width=width, color='#9813DE', bottom=np.array(df.ef), label=r'${\itE. coli}$', alpha=0.75)
plt.scatter(np.ones(4)*0,np.ones(4)*0.5, facecolors='none', edgecolors='black', alpha=1)
plt.scatter(df_og['time']/24, df_og['ef'], facecolors='none', edgecolors='black', alpha=1)
plt.xticks(ticks=t, labels=np.array(df['time'].unique()), fontsize=20)
plt.ylim(0,1)
plt.yticks(fontsize=20)
plt.grid(False)
plt.ylabel(r'R.A.')
plt.xlabel('Time (d)')
plt.title('Biofilm Culture')
# plt.savefig('FigS2B_1.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S2C

SI Figure S2D is taken from Fig1C
"""

ecoli_co_pa = [0, np.mean([521, 1162, 1506, 6241, 14282, 3599]), 
            np.mean([7696, 58535, 180358, 238798,32897]),
            np.mean([238670, 37391, 227485, 179853,47073]),
            np.mean([292677, 131934,133782, 46141])]
ecoli_co_pa_std = [0, stats.sem([521, 1162, 1506, 6241, 14282, 3599]), 
                stats.sem([7696, 58535, 180358, 238798, 32897]),
            stats.sem([238670, 37391, 227485, 179853, 47073]),
            stats.sem([292677, 131934,133782, 46141])]
ecoli = [0, np.mean([429326.5281541101, 20830.831426016266, 26769.516920138325,26992.17734452068]), 
         np.mean([111347.8469921274, 183315.5309570613, 145620.88152866834, 409437.6682564925, 134370.40829945993, 277304.7236815254]),
            np.mean([387877.2050713496, 312281.15418290574, 486196.3041085362, 286454.46008238854]),
            np.mean([337383.11888024263, 348827.4483178621, 148572.00581428283, 269781.17495603,  461549.52062919])]
ecoli_std = [0, stats.sem([429326.5281541101, 20830.831426016266, 26769.516920138325,26992.17734452068]), 
             stats.sem([111347.8469921274, 183315.5309570613, 145620.88152866834, 409437.6682564925, 134370.40829945993, 277304.7236815254]),
            stats.sem([387877.2050713496, 312281.15418290574, 486196.3041085362, 286454.46008238854]),
            stats.sem([128992.59076753035, 92435.67970619304, 182715.3318089329, 89181.74673619,116083.68712655])]

ec_co_ent = [0, np.mean([4055,20499, 5187, 6332, 1957, 8876]), 
          np.mean([29622, 84081, 6241, 104382, 148197, 147312, 123539]),
          np.mean([229821, 8775, 340910, 144991, 233076, 194977]),
          np.mean([91793, 325309,158899])]
ec_co_ent_std= [0, stats.sem([4055, 20499, 5187, 6332, 1957,8876]),
             stats.sem([29622, 84081, 6241, 104382, 148197, 147312,123539]),
          stats.sem([229821, 8775, 340910, 144991, 233076, 194977]),
          stats.sem([91793, 325309, 158899])]

ent_co_ec = [0, np.mean([95967, 64789, 139213, 103088, 57349, 33617]), 
          np.mean([65875, 97997, 21398, 14864, 50231, 95568, 97224]),
          np.mean([57668, 26571, 40493, 34867, 46649, 66529]),
          np.mean([5863, 25565, 83554])]
ent_co_ec_std= [0, stats.sem([95967, 64789, 139213, 103088, 57349, 33617]),
             stats.sem([65875, 97997, 21398, 14864, 50231, 95568,97224]),
          stats.sem([57668, 26571, 40493, 34867, 46649,66529]),
          stats.sem([5863, 25565, 83554])]

ent_co_pa = [0, np.mean([9527, 16130, 40180, 19456, 8663, 29653]), 
          np.mean([51496, 109658, 29879, 14297, 26301, 13789]),
          np.mean([52988, 60321, 24687,34963, 115294]),
          np.mean([92694, 34963,116788])]
ent_co_pa_std= [0, stats.sem([9527,16130, 40180,38727,19456, 8663,29653]),
             stats.sem([51496, 109658, 29879, 14297, 26301, 13789]),
          stats.sem([52988, 60321, 24687, 34963,115294]),
          stats.sem([92694, 34963,116788])]

pa_co_ent = [0, np.mean([4553,110590, 4249, 342978,2067, 210, 1051, 2606]), 
          np.mean([806,87617, 340569,389, 1011,7313, 13789]),
          np.mean([6793, 249928, 437070,11604, 78296, 115294]),
          np.mean([166180, 78296, 116788])]
pa_co_ent_std= [0, stats.sem([4553, 110590, 4249, 342978, 2067, 210,1051, 2606]),
             stats.sem([806, 87617, 340569, 389, 1011, 7313, 13789]),
          stats.sem([6793, 249928, 437070, 11604, 78296, 115294]),
          stats.sem([166180,78296, 116788])]

ent = [0, np.mean([145719.7700434657, 132853.42173086293, 49593.08968900539, 42897.058102400624]), 
       np.mean([191166.61280572205, 96122.08029061732, 133520.86963127236, 172674.3443290808, 127962.97477912421, 105634.37227330725]),
          np.mean([190393.93985167015, 82472.47297014121, 103497.2411231636, 83590.30802445539]),
          np.mean([128992.59076753035, 92435.67970619304, 182715.3318089329, 89181.74673619,116083.68712655])]

ent_std = [0, stats.sem([145719.7700434657, 132853.42173086293, 49593.08968900539, 42897.058102400624]),
           stats.sem([191166.61280572205, 96122.08029061732, 133520.86963127236, 172674.3443290808, 127962.97477912421, 105634.37227330725]),
          stats.sem([190393.93985167015, 82472.47297014121, 103497.2411231636, 83590.30802445539]),
          stats.sem([128992.59076753035, 92435.67970619304, 182715.3318089329, 89181.74673619,116083.68712655])]

pa_co_ec = [0, np.mean([24205, 60648,2129, 571,1363, 14791, 1735]),
          np.mean([15080, 133186, 432725, 41, 13178, 8369]),
          np.mean([12266,302415,304050, 13450, 419061, 98712]),
          np.mean([181554, 79787, 255382,368915])]

pa_co_ec_std = [0, stats.sem([24205, 60648, 2129, 571, 1363,14791, 1735]), 
          stats.sem([15080, 133186, 432725, 41,13178,8369]),
          stats.sem([12266,302415, 304050, 13450, 419061,98712]),
          stats.sem([181554, 79787, 255382, 368915])]

pa = [0, np.mean([362289.6023790988, 550.316753017643, 118188.07173208568]), 
      np.mean([126345.30028869871, 80597.25082750342, 21905.3075725621, 485315.62176974496, 3045.379556830059, 2071.840862834238, 27850.31514594865, 8115.406266916172]),
          np.mean([337295.43424540985, 398535.4141324807, 242220.39670517063, 23597.645090594877, 164344.37751006486]),
          np.mean([385489.7102260408, 439545.5021568775, 227691.2968063428, 5.02094258e+05, 452658.56209157, 520242.070775, 429652.64031865])]

pa_std = [0, stats.sem([362289.6023790988, 550.316753017643, 118188.07173208568]), 
          stats.sem([126345.30028869871, 80597.25082750342, 21905.3075725621, 485315.62176974496, 3045.379556830059, 2071.840862834238, 27850.31514594865, 8115.406266916172]),
          stats.sem([337295.43424540985, 398535.4141324807, 242220.39670517063, 23597.645090594877, 164344.37751006486]),
          stats.sem([385489.7102260408, 439545.5021568775, 227691.2968063428, 5.02094258e+05, 452658.56209157, 520242.070775, 429652.64031865])]



time = [0, 24, 48, 72, 96]
plt.figure(figsize=(3.5,4))
plt.plot(time, ecoli_co_pa, alpha=0.75, c='#9813DE', markerfacecoloralt='#E2E604', fillstyle='left',  marker='.', markeredgecolor='None', linestyle='', markersize=15)
plt.errorbar(time, ecoli_co_pa, yerr=ecoli_co_pa_std, color='#9813DE', alpha=0.75)
plt.plot(time, ec_co_ent, alpha=0.75, c='#9813DE', markerfacecoloralt='#008F8F', fillstyle='left',  marker='.', markeredgecolor='None', linestyle='', markersize=15)
plt.errorbar(time, ec_co_ent, yerr=ec_co_ent_std, color='#9813DE', alpha=0.75)
plt.plot(time, ecoli, linestyle='--', marker='p', color='#9813DE', alpha=0.75)
plt.errorbar(time, ecoli, yerr=ecoli_std, linestyle='--', color='#9813DE', alpha=0.75)
plt.plot(time, ent_co_pa, alpha=0.75, c='#008F8F', markerfacecoloralt='#E2E604', fillstyle='left',  marker='.', markeredgecolor='None', linestyle='', markersize=15)
plt.errorbar(time, ent_co_pa, yerr=ent_co_pa_std, color='#008F8F', alpha=0.75)
plt.plot(time, ent_co_ec, alpha=0.75, c='#008F8F', markerfacecoloralt='#9813DE', fillstyle='left',  marker='.', markeredgecolor='None', linestyle='', markersize=15)
plt.errorbar(time, ent_co_ec, yerr=ent_co_ec_std, color='#008F8F', alpha=0.75)
plt.plot(time, ent, linestyle='--', marker='p', color='#008F8F', alpha=0.75)
plt.errorbar(time, ent, yerr=ent_std, linestyle='--', color='#008F8F', alpha=0.75)
plt.plot(time, pa_co_ec, alpha=0.75, c='#E2E604', markerfacecoloralt='#9813DE', fillstyle='left',  marker='.', markeredgecolor='None', linestyle='', markersize=15)
plt.errorbar(time, pa_co_ec, yerr=pa_co_ec_std, color='#E2E604', alpha=0.75)
plt.plot(time, pa_co_ent, alpha=0.75, c='#E2E604', markerfacecoloralt='#008F8F', fillstyle='left',  marker='.', markeredgecolor='None', linestyle='', markersize=15)
plt.errorbar(time, pa_co_ent, yerr=pa_co_ent_std, color='#E2E604', alpha=0.75)
plt.plot(time, pa, linestyle='--', marker='p', color='#E2E604', alpha=0.75)
plt.errorbar(time, pa, yerr=pa_std, linestyle='--', color='#E2E604', alpha=0.75)
plt.ylabel('Volume $\mu$m$^3$', fontsize=20)
plt.xlabel('Time (hr)', fontsize = 20)
plt.xticks(fontsize=20)
plt.grid(False)
plt.yticks(ticks=np.arange(0,600000,100000), labels=np.arange(0,6,1), fontsize=20)
plt.ylim(-10000,500000)
# plt.savefig('FigS2C.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figures S3 and S4 are composed of .tiff images exported from ZEN Blue  
"""
#%%
"""
SI Figure S5D

Made paper plot in BiofilmQ
local Density Histograms from Monoculture 
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS5A')

plt.hist(df['Pa2'].dropna(), color='#E2E604', bins=50)
plt.hist(df['Pa'].dropna(), color='#E2E604', bins=50)
plt.hist(df['Ec'].dropna(), color='#9813DE', bins=50)
plt.hist(df['Ef'].dropna(), color='#008F8F', bins=50)
plt.grid(False)
plt.ylabel('Local Density (r=10um)')
plt.xlabel('Species')
# plt.savefig('FigS5D.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S5E
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS5B')

my_pal = {'Pa_8d':'#E2E604', 'Pa_5d':'#E2E604', 'Ec_8d':'#9813DE', 'Ef_8d':'#008F8F'}
sns.boxplot(df, palette =my_pal)
sns.stripplot(df, color='black')
plt.grid(False)
plt.ylim(0,1)
plt.ylabel('Local Density (r=10um)')
plt.xlabel('Species')
# plt.savefig('FigS5E.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

#print mean to get carrying capacity
print(df['Pa_8d'].mean())
print(df['Pa_5d'].mean())
print(df['Ec_8d'].mean())
print(df['Ef_8d'].mean())
#%%
"""
SI Figure S5F
"""
Pa = pd.read_excel('SI_Data.xlsx', sheet_name='FigS5C_1')
Ec = pd.read_excel('SI_Data.xlsx', sheet_name='FigS5C_2')
Ef = pd.read_excel('SI_Data.xlsx', sheet_name='FigS5C_3')

#fit curve to the natural log to find intrinsic rate of increase
Pa['Pa'] = np.log(Pa['Pa'])
Ec['Ec'] = np.log(Ec['Ec'])
Ef['Ef'] = np.log(Ef['Ef'])

slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(Pa['Time'], Pa['Pa'])
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(Ec['Time'], Ec['Ec'])
slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(Ef['Time'], Ef['Ef'])

plt.figure(figsize=(6,3))
sns.regplot(x=Pa['Time'], y=Pa['Pa'], color = '#E2E604', marker='p', line_kws={'label':"y={0:.2f}x+{1:0.1f}".format(np.round(slope1, decimals=2), np.round(intercept1, decimals=0))},ci=None)
sns.regplot(x=Ec['Time'], y=Ec['Ec'], color = '#9813DE', marker='p', line_kws={'label':"y={0:.2f}x+{1:0.1f}".format(np.round(slope2, decimals=2), np.round(intercept2, decimals=0))}, ci=None)
sns.regplot(x=Ef['Time'], y=Ef['Ef'], color = '#008F8F', marker='p', line_kws={'label':"y={0:.2f}x+{1:0.1f}".format(np.round(slope3, decimals=2), np.round(intercept3, decimals=0))}, ci=None)
plt.xlabel(r'Time', fontsize=12)
plt.ylabel(r'ln(biovolume (um3))', fontsize=12)
plt.xticks(fontsize=12, rotation=-90)
plt.yticks(fontsize=12)
plt.legend()
plt.grid(False)
# plt.savefig('FigS5F.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S6B

PA14 Monoculture Growth Traces
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS6B')

plt.plot(df.x1, df.y1, color='#E2E604')
plt.fill_between(df.x1, df.y1-df.y_err1, df.y1+df.y_err1, color='#E2E604', alpha=0.25)
plt.plot(df.x2, df.y2, color='#E2E604')
plt.fill_between(df.x2, df.y2-df.y_err2, df.y2+df.y_err2, color='#E2E604', alpha=0.25)
plt.plot(df.x3, df.y3, color='#E2E604')
plt.fill_between(df.x3, df.y3-df.y_err3, df.y3+df.y_err3, color='#E2E604', alpha=0.25)

plt.title('PA14 Monoculture')
# plt.savefig('FigS6B.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S7A

Stochastic-Spatial Model
"""

df = pd.read_excel('Stochastic_Spatial_Model/no_dispersal_cont_20000.xlsx')

time_array = list(df.index)
spc_1 = df['Species1']
spc_2 = df['Species2']
spc_3 = df['Species3']
num_iter = len(df)

plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)
# plt.savefig('FigS7A.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S7B
"""

df = pd.read_excel('Stochastic_Spatial_Model/unbiased_dispersal_cont_20000.xlsx')

time_array = list(df.index)
spc_1 = df['Species1']
spc_2 = df['Species2']
spc_3 = df['Species3']
num_iter = len(df)

plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)
# plt.savefig('FigS7B.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S7C
"""

df = pd.read_excel('Stochastic_Spatial_Model/biased_dispersal_final_cont_20000.xlsx')

time_array = list(df.index)
spc_1 = df['Species1']
spc_2 = df['Species2']
spc_3 = df['Species3']
num_iter = len(df)

plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)
# plt.savefig('FigS7C.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S7D

Rapid Mixing Stochastic-Spatial Model
"""

df = pd.read_excel('Stochastic_Spatial_Model/biased_dispersal_rapid_mixing_cont_20000.xlsx')

time_array = list(df.index)
spc_1 = df['Species1']
spc_2 = df['Species2']
spc_3 = df['Species3']
num_iter = len(df)

plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)

# plt.savefig('FigS7D.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S7E

Pa invasion into the stochastic spatial model

Run the invasion script in stochastic_spatial_models to generate data
"""
df = pd.read_csv('Stochastic_Spatial_Model/invasion_data/pa_invasion_biased_dispersal.csv')

time_array = list(df.index)
spc_1 = df['Species1']
spc_2 = df['Species2']
spc_3 = df['Species3']
num_iter = len(df)

plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)

# plt.savefig('FigS7E.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S7F

Ec invasion into the stochastic spatial model
"""
df = pd.read_csv('Stochastic_Spatial_Model/invasion_data/ec_invasion_biased_dispersal.csv')

time_array = list(df.index)
spc_1 = df['Species1']
spc_2 = df['Species2']
spc_3 = df['Species3']
num_iter = len(df)

plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)

# plt.savefig('FigS7F.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S7G

Ef invasion into the stochastic spatial model
"""
df = pd.read_csv('Stochastic_Spatial_Model/invasion_data/ef_invasion_biased_dispersal.csv')

time_array = list(df.index)
spc_1 = df['Species1']
spc_2 = df['Species2']
spc_3 = df['Species3']
num_iter = len(df)

plt.plot(time_array, spc_1, color='#9813DE', label=r'$\it{E. coli}$')
plt.plot(time_array, spc_2, color='#E2E604', label=r'$\it{P. aeruginosa}$')
plt.plot(time_array, spc_3, color='#008F8F', label=r'$\it{E. faecalis}$')
plt.ylim(-5,2100)
plt.xlim(-5,num_iter)
plt.xticks(fontsize=20, rotation=-45)
plt.yticks(fontsize=20)
plt.xlabel('Time',fontsize=20)
plt.ylabel('Individuals',fontsize=20)

# plt.savefig('FigS7G.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S8 was made in MATLAB using the scripts found in the Mean_Field_Model directory
"""
#%%
"""
SI Figure S9A

PA14 change predicted by abundance
Monoculture PA14 change predicted by abundance
Triculture data taken from Fig2E
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS9A_1')
x = df['v']
y = df['dv']

df2 = pd.read_excel('SI_Data.xlsx', sheet_name='Fig9A_2')
x2 = df2['v']
y2 = df2['dv']

slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(x, y)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x2, y2)

sns.regplot(x=x2, y=y2, color = 'black', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope2,intercept2, r_value2, p_value2)})
sns.regplot(x=x, y=y, color = '#E2E604', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope3,intercept3, r_value3, p_value3)})
plt.ylabel(r'dv_pa14/dt', fontsize=20)
plt.xlabel(r'v_pa14', fontsize=20)
plt.xticks(fontsize=20, rotation=-90)
plt.yticks(fontsize=20)
plt.legend()
plt.hlines(0, 0, 600000, color='red', linestyles='--', alpha=0.5)
# plt.savefig('FigS9A.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()

#%%
"""
SI Figure S9B

Species change as a fraction. V1-V0 over V0
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS9B')

Pa = df['Pa']
Ec = df['Ec']
Ef = df['Ef']


slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(Pa, Pa)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(Pa, Ef)
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(Pa, Ec)

sns.regplot(x=Pa, y=Pa, color = '#E2E604', marker='D', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope3,intercept3, r_value3, p_value3)})
sns.regplot(x=Pa, y=Ef, color = '#008F8F', marker='D', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope2,intercept2, r_value2, p_value2)})
sns.regplot(x=Pa, y=Ec, color = '#9813DE', marker='D', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope1,intercept1, r_value1, p_value1)})

plt.ylabel(r'dv_species', fontsize=20)
plt.xlabel(r'dv_pa14', fontsize=20)
plt.xticks(fontsize=20, rotation=-90)
plt.yticks(fontsize=20)
plt.legend()
# plt.savefig('FigS9B.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S9C

Species change predicts density change
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS9C')

x = df['x']
y = df['y']

slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(x, y)
sns.regplot(x=x, y=y, color = 'black', marker='D', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope3,intercept3, r_value3, p_value3)})
plt.ylabel(r'dv_pa14/dt', fontsize=20)
plt.xlabel(r'v_pa14', fontsize=20)
plt.xticks(fontsize=20, rotation=-90)
plt.yticks(fontsize=20)
plt.legend()
# plt.savefig('FigS9C.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S10 is composed of .tiff files exported form ZEN Blue
"""
#%%
"""
SI Figure S11A

Histogram of Local Density. Paper figure generated in BiofilmQ.
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS11A')

plt.hist(df['bqsS'].dropna(), bins=50, color='green', alpha=0.5)
plt.hist(df['WT'].dropna(), bins=50, color='black', alpha=0.5)
# plt.savefig('FigS11A.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S11B

Histogram of Joint Local Density

Paper has the tops of the bars summing to one. This plot has the area under the curve equaling 1. 
"""

df1 = pd.read_excel('SI_Data.xlsx', sheet_name='FigS19B_2')
df2 = pd.read_excel('SI_Data.xlsx', sheet_name='FigS19B_1')

bqsS = df1['Pa_bqsS_Density10']
wt = df2['Pa_WT_Density10um']

plt.hist(bqsS, bins=20, density=True, color='green', alpha=0.5)
plt.hist(wt, bins=20, density = True, color='black', alpha=0.5)
plt.xlabel('Joint Local Density')
plt.ylabel('Normalized Frequency')
plt.title('Tri-Culture')
# plt.savefig('FigS11B.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S11C
"""

df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS11C')
x = df['pa']
y = df['dv_pa']

slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(x, y)
sns.regplot(x=x, y=y, color = 'green', marker='D', line_kws={'label':"y={0:.1f}x+{1:.1f}, r={2:.1f},  p={3:.4f}".format(slope3,intercept3, r_value3, p_value3)})
plt.ylabel(r'dv_pa14/dt', fontsize=20)
plt.xlabel(r'v_pa14', fontsize=20)
plt.xticks(fontsize=20, rotation=-90)
plt.yticks(fontsize=20)
plt.hlines(0, xmin=-10000, xmax=630000, color='red', linestyles='--')
plt.xlim(-10000,630000)
plt.ylim(-22000,22000)
plt.legend()
# plt.savefig('FigS11C.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S12A

Invasion First Time Point
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS12A').T.reset_index()

new_header = df.iloc[0] 
df = df[1:] 
df.columns = new_header 

df = df.drop(columns=['Index'])

plt.figure(figsize=(4,3))
sns.boxplot(df)
plt.grid(False)
plt.xticks(rotation=90)
plt.ylabel(r'Biovolume', fontsize=18)
# plt.savefig('FigS12A.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S12B

Invasion Second Time Point
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS12B').T.reset_index()

new_header = df.iloc[0] 
df = df[1:] 
df.columns = new_header 
df = df.drop(columns=['Index'])

plt.figure(figsize=(4,3))
sns.boxplot(df)
plt.grid(False)
plt.xticks(rotation=90)
plt.ylabel(r'Biovolume', fontsize=18)
# plt.savefig('FigS12B.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S12C_1
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS12C_1').T.reset_index()

new_header = df.iloc[0] 
df = df[1:] 
df.columns = new_header 
df = df.drop(columns=['Index'])


for i in np.arange(1,len(df)+1):
    if df['Ec_0_WT'][i] != 'NaN':
        plt.plot([90,96,112], [0,df['Ec_0_WT'][i],df['Ec_16_WT'][i]], color='black')
    if df['Ec_0_bqsS'][i] != 'NaN':
        plt.plot([90,96,112], [0,df['Ec_0_bqsS'][i],df['Ec_16_bqsS'][i]], color='green', linestyle='dashed')
plt.xlim(85, 117)
plt.hlines(y=40000, xmin=80, xmax=170, color='#9813DE', linestyle='--', alpha=0.5)
# plt.savefig('FigS12C_1.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure 12C_2
"""
df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS12C_2').T.reset_index()

new_header = df.iloc[0] 
df = df[1:] 
df.columns = new_header
df = df.drop(columns=['Index'])
 

for i in np.arange(1,len(df)+1):
    if df['Ef_0_WT'][i] != 'NaN':
        plt.plot([90,96,112], [0,df['Ef_0_WT'][i],df['Ef_16_WT'][i]], color='black')
    if df['Ef_0_bqsS'][i] != 'NaN':
        plt.plot([90,96,112], [0,df['Ef_0_bqsS'][i],df['Ef_16_bqsS'][i]], color='green', linestyle='dashed')
plt.hlines(y=80000, xmin=80, xmax=170, color='#008F8F', linestyle='--', alpha=0.5)
plt.xlim(85, 117)
# plt.savefig('FigS12C_2.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()
#%%
"""
SI Figure S12C_3
""" 
df = pd.read_excel('SI_Data.xlsx', sheet_name='FigS12C_3').T.reset_index()

new_header = df.iloc[0] 
df = df[1:] 
df.columns = new_header 
df = df.drop(columns=['Index'])


for i in np.arange(1,len(df)+1):
    if df['Pa_0_WT'][i] != 'NaN':
        plt.plot([90,96,112], [0,df['Pa_0_WT'][i],df['Pa_16_WT'][i]], color='black')
    if df['Pa_0_bqsS'][i] != 'NaN':
        plt.plot([90,96,112], [0,df['Pa_0_bqsS'][i],df['Pa_16_bqsS'][i]], color='green', linestyle='dashed')
plt.hlines(y=300000, xmin=80, xmax=170, color='#E2E604', linestyle='--', alpha=0.5)
plt.xlim(85, 117)
# plt.savefig('FigS12C_2.svg', dpi=300, facecolor='w', edgecolor='b',
#         orientation='portrait', format='svg',
#         transparent=False, bbox_inches='tight', pad_inches=.05, metadata=None)
plt.show()