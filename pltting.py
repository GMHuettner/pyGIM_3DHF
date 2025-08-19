"""
Georg HÃ¼ttner - October 2024

Load the results of both actual scripts to do some nice plotting
also the  2D/3D comparison is in here - dont ask why, it doesnt like it and will start to scream
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import pathlib
import cartopy.crs as ccrs
import pickle
import glob
import os.path
from scipy.interpolate import RegularGridInterpolator
import matplotlib.gridspec as gridspec
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter)
import shapely
from shapely.geometry import LineString, Point

linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),
     ('long dash with offset', (5, (10, 3))),
     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]

#%% load data

# input data for nice overview plots

# 2d data if used

name_2D = '2024-11-15_18-22-15_2d_ASB_1'
with open('output/'+name_2D+'.pkl','rb') as file:
    data_plot_2D = pickle.load(file)

# 3d data if used

name_3D = '2024-11-08_01-04-22_3d_DC_full'
with open('output/'+name_3D+'.pkl','rb') as file:
    data_plot_3D = pickle.load(file)
    
name_3D_2 = '2024-11-07_16-38-11_3d_LV_full_2'
with open('output/'+name_3D_2+'.pkl','rb') as file:
    data_plot_3D_2 = pickle.load(file)

#%% 4 with feature

profile_coord_LV_3=([1325000,1365000],[-275000,-375000])

# formatter = ticker.ScalarFormatter(useOffset=False)
# # formatter.set_scientific(True)
# formatter.set_powerlimits((0, 0))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(8,8),sharex=True,sharey=True,dpi=500,layout='compressed')

a1 = ax1.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['shf_3D'][1:-1,1:-1],vmin=45,vmax=60,cmap="afmhot")
pr1 = ax1.plot([1325,1365],[-275,-375],'black',linestyle=(0, (5, 5)),linewidth=1)

ax1.set_xlabel('Easting in km')
ax1.set_ylabel('Northing in km')
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top')
# ax1.yaxis.set_major_formatter(formatter)
ax1.set_aspect('equal')
ax1.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
# ax1.set_xlim(data_plot_3D['x_int_t'][1]/1000000,data_plot_3D['x_int_t'][-1]/1000000)
ax1.set_ylim(data_plot_3D['y_int_t'][1]/1000,data_plot_3D['y_int_t'][-1]/1000)
# ax1.colorbar(label='mW/m$^2$')

a2 = ax2.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['shf_1D'][1:-1,1:-1],vmin=45,vmax=60,cmap="afmhot")
ax2.set_xlabel('Easting in km')
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')
# ax2.set_ylabel('y in m')
fig.colorbar(a2,ax=[ax1,ax2],label='GHF in mW/m$^2$',location='bottom',shrink=0.7)
ax2.set_aspect('equal')
ax2.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

a3 = ax3.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['shf_3D'][1:-1,1:-1]-data_plot_3D['shf_1D'][1:-1,1:-1],vmin=-10,vmax=10,cmap="bwr")
# ax3.set_xlabel('x in m')
ax3.set_ylabel('Northing in km')
fig.colorbar(a3,ax=ax3,label='difference in mW/m$^2$',location='bottom')
ax3.set_aspect('equal')
ax3.tick_params(top=True, labeltop=False, bottom=False, labelbottom=False)

# a4 = ax4.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['grid_topo'][1:-1,1:-1]/1000,cmap="gist_earth")
# a4 = ax4.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['grid_sed'][1:-1,1:-1]/1000,cmap="YlGn_r")
# a4 = ax4.pcolormesh(data_plot_3D_2['x_int_t'][1:-1]/1000,data_plot_3D_2['y_int_t'][1:-1]/1000,data_plot_3D_2['grid_moho_for_topo'][1:-1,1:-1]/1000,vmin=40,vmax=42,cmap="gist_earth_r")
# a4 = ax4.pcolormesh(data_plot_3D_2['x_int_t'][1:-1]/1000,data_plot_3D_2['y_int_t'][1:-1]/1000,data_plot_3D_2['grid_lab_for_topo'][1:-1,1:-1]/1000,vmin=178,vmax=188,cmap="gist_earth_r")
a4 = ax4.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['grid_A_for_topo'][1:-1,1:-1]/1e-6,cmap="magma")
# ax4.set_xlabel('x in m')
# ax4.set_ylabel('y in m')
# fig.colorbar(a4,ax=ax4,label='Height in km',location='bottom')
# fig.colorbar(a4,ax=ax4,label='Sediment thickness in km',location='bottom')
# fig.colorbar(a4,ax=ax4,label='Depth in km',location='bottom')
fig.colorbar(a4,ax=ax4,label='HP in µW/m$^3$',location='bottom')
ax4.set_aspect('equal')
ax4.tick_params(top=True, labeltop=False, bottom=False, labelbottom=False)

fig.text(0.165, 0.985, "A", ha='left', va='top', fontsize=14)
fig.text(0.9, 0.985, "B", ha='left', va='top', fontsize=14)
fig.text(0.165, 0.5, "C", ha='left', va='top', fontsize=14)
fig.text(0.9, 0.5, "D", ha='left', va='top', fontsize=14)


# fig.savefig('abb/3d_LV_HP.pdf',bbox_inches='tight')
# fig.savefig('abb/3d_LV_HP.png',bbox_inches='tight')


#%% 4 features with profiles test 2 DC

area_coord_LV = [1300000,1400000,-275000,-375000] # Lake Vostok "south-western" edge
area_coord_DC = [1355000,1435000,-855000,-985000] # Dome C
area_coord_ASB = [1800000,1900000,-650000,-750000] # ASB

profile_coord_LV_1=([1300000,1400000],[-320000,-300000]) # LV prof 1
profile_coord_LV_2=([1300000,1400000],[-340000,-320000]) # LV prof 2
profile_coord_LV_3=([1325000,1365000],[-275000,-375000]) # LV prof 3

profile_coord_DC_1=([1405000,1435000],[-985000,-855000]) # DC prof 1
profile_coord_DC_2=([1360000,1415000],[-985000,-855000]) # DC prof 2
profile_coord_DC_3=([1355000,1420000],[-920000,-985000]) # DC prof 3
profile_coord_DC_4=([1355000,1435000],[-905000,-970000]) # DC prof 3

profile_coord_ASB_1=([1900000,1800000],[-650000,-750000]) # ASB prof 1
profile_coord_ASB_2=([1800000,1900000],[-700000,-700000]) # ASB prof 2

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(8,10),dpi=500,layout='compressed', subplot_kw={'projection': None}) #LV (8,8), DC ()
#,sharex=True,sharey=True
ax1.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax1 = plt.subplot(2,2,1,projection=ccrs.AzimuthalEquidistant(central_latitude=-90))
pr1 = ax1.plot(profile_coord_DC_1[0],profile_coord_DC_1[1],'gold',transform=ccrs.AzimuthalEquidistant(central_latitude=-90))
pr2 = ax1.plot(profile_coord_DC_2[0],profile_coord_DC_2[1],'firebrick',transform=ccrs.AzimuthalEquidistant(central_latitude=-90))
pr3 = ax1.plot(profile_coord_DC_3[0],profile_coord_DC_3[1],'darkorange')
pr4 = ax1.plot(profile_coord_DC_4[0],profile_coord_DC_4[1],'darkviolet')
# ax1.set_extent([area_coord_LV[0],area_coord_LV[1],area_coord_LV[2],area_coord_LV[3]],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
ax1.set_extent([area_coord_DC[0],area_coord_DC[1],area_coord_DC[2],area_coord_DC[3]],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
# ax1.set_extent([area_coord_ASB[0],area_coord_ASB[1],area_coord_ASB[2],area_coord_ASB[3]],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, rotate_labels = False)
gl.top_labels = False
gl.right_labels = False
gl.ylabel_style = {'rotation': 45, 'rotation_mode': 'anchor'}
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()
ax1.set_title('Profile locations')

a2 = ax2.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['grid_topo'][1:-1,1:-1]/1000,cmap="gist_earth")
ax2.set_xlabel('Easting in km')
ax2.set_ylabel('Northing in km')
# ax2.xaxis.tick_top()
# ax2.xaxis.set_label_position('top')
ax2.yaxis.set_label_position('right')
fig.colorbar(a2,ax=ax2,label='Height in km',location='bottom')#,pad=-0.02
ax2.set_aspect('equal')
ax2.tick_params(left=False, labelleft=False, right=True, labelright=True)
ax2.set_title('Topography')

a3 = ax3.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['grid_sed'][1:-1,1:-1]/1000,cmap="YlGn_r")
ax3.set_xlabel('Easting in km')
ax3.set_ylabel('Northing in km')
fig.colorbar(a3,ax=ax3,label='km',location='bottom')
ax3.set_aspect('equal')
# ax3.tick_params(top=True, labeltop=False, bottom=False, labelbottom=False)
ax3.set_title('Sediment thickness')

a4 = ax4.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['grid_A_for_topo'][1:-1,1:-1]/1e-6,cmap="magma")
ax4.set_xlabel('Easting in km')
ax4.set_ylabel('Northing in km')
ax4.yaxis.set_label_position('right')
fig.colorbar(a4,ax=ax4,label='µW/m$^3$',location='bottom')
ax4.set_aspect('equal')
ax4.tick_params(left=False, labelleft=False, right=True, labelright=True)
ax4.set_title('Heat production')

fig.text(0.165, 0.995, "A", ha='left', va='top', fontsize=14)
fig.text(0.815, 0.995, "B", ha='left', va='top', fontsize=14)
fig.text(0.165, 0.5, "C", ha='left', va='top', fontsize=14)
fig.text(0.815, 0.5, "D", ha='left', va='top', fontsize=14)


#%% 3 without feature DC

fig = plt.figure(figsize=(6, 10), dpi=500,layout='constrained')
# ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(8,8),sharex=True,sharey=True,dpi=300,layout='compressed')
gs = fig.add_gridspec(nrows=2, ncols=2)

ax1 = fig.add_subplot(gs[0, 0])
a1 = ax1.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['shf_3D'][1:-1,1:-1],vmin=40,vmax=70,cmap="afmhot")
ax1.set_xlabel('Easting in km')
ax1.set_ylabel('Northing in km')
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top')
# ax1.yaxis.set_major_formatter(formatter)
ax1.set_aspect('equal')
ax1.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
# ax1.set_xlim(data_plot_3D['x_int_t'][1]/1000000,data_plot_3D['x_int_t'][-1]/1000000)
ax1.set_ylim(data_plot_3D['y_int_t'][1]/1000,data_plot_3D['y_int_t'][-1]/1000)
# ax1.colorbar(label='mW/m$^2$')

ax2 = fig.add_subplot(gs[0, 1])
a2 = ax2.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['shf_1D'][1:-1,1:-1],vmin=40,vmax=70,cmap="afmhot")
ax2.set_xlabel('Easting in km')
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')
# ax2.set_ylabel('y in m')
fig.colorbar(a2,ax=[ax1,ax2],label='GHF in mW/m$^2$',location='bottom',shrink=0.65)#0.65
ax2.set_aspect('equal')
ax2.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, labelleft=False)

ax3 = fig.add_subplot(gs[1, :])
a3 = ax3.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['shf_3D'][1:-1,1:-1]-data_plot_3D['shf_1D'][1:-1,1:-1],vmin=-10,vmax=10,cmap="bwr")
# ax3.set_xlabel('x in m')
ax3.set_ylabel('Northing in km')
fig.colorbar(a3,ax=ax3,label='difference in mW/m$^2$',location='bottom',shrink=0.49)
ax3.set_aspect('equal')
ax3.tick_params(top=True, labeltop=False, bottom=False, labelbottom=False)

fig.text(0.075, 0.985, "A", ha='left', va='top', fontsize=14)
fig.text(0.965, 0.985, "B", ha='left', va='top', fontsize=14)
fig.text(0.3, 0.495, "C", ha='left', va='top', fontsize=14)

#%% 4 features with profiles test 2

area_coord_LV = [1300000,1400000,-275000,-375000] # Lake Vostok "south-western" edge
area_coord_DC = [1355000,1435000,-855000,-985000] # Dome C
area_coord_ASB = [1800000,1900000,-650000,-750000] # ASB

profile_coord_LV_1=([1300000,1400000],[-320000,-300000]) # LV prof 1
profile_coord_LV_2=([1300000,1400000],[-340000,-320000]) # LV prof 2
profile_coord_LV_3=([1325000,1365000],[-275000,-375000]) # LV prof 3

profile_coord_DC_1=([1405000,1435000],[-985000,-855000]) # DC prof 1
profile_coord_DC_2=([1360000,1415000],[-985000,-855000]) # DC prof 2
profile_coord_DC_3=([1355000,1420000],[-920000,-985000]) # DC prof 3
profile_coord_DC_4=([1355000,1435000],[-905000,-970000]) # DC prof 3

profile_coord_ASB_1=([1900000,1800000],[-650000,-750000]) # ASB prof 1
profile_coord_ASB_2=([1800000,1900000],[-700000,-700000]) # ASB prof 2

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(8,8),dpi=500,layout='compressed', subplot_kw={'projection': None}) #LV (8,8), DC ()
#,sharex=True,sharey=True
ax1.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
ax1 = plt.subplot(2,2,1,projection=ccrs.AzimuthalEquidistant(central_latitude=-90))
pr1 = ax1.plot(profile_coord_ASB_1[0],profile_coord_ASB_1[1],'gold',transform=ccrs.AzimuthalEquidistant(central_latitude=-90))
pr2 = ax1.plot(profile_coord_ASB_2[0],profile_coord_ASB_2[1],'firebrick',transform=ccrs.AzimuthalEquidistant(central_latitude=-90))
# pr3 = ax1.plot(profile_coord_LV_3[0],profile_coord_LV_3[1],'darkorange')
# pr4 = ax1.plot(profile_coord_DC_4[0],profile_coord_DC_4[1],'darkviolet')
# ax1.set_extent([area_coord_LV[0],area_coord_LV[1],area_coord_LV[2],area_coord_LV[3]],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
# ax1.set_extent([area_coord_DC[0],area_coord_DC[1],area_coord_DC[2],area_coord_DC[3]],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
ax1.set_extent([area_coord_ASB[0],area_coord_ASB[1],area_coord_ASB[2],area_coord_ASB[3]],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
gl = ax1.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, rotate_labels = False)
gl.top_labels = False
gl.right_labels = False
gl.ylabel_style = {'rotation': 45, 'rotation_mode': 'anchor'}
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()
ax1.set_title('Profile locations')

a2 = ax2.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['grid_topo'][1:-1,1:-1]/1000,cmap="gist_earth")
ax2.set_xlabel('Easting in km')
ax2.set_ylabel('Northing in km')
# ax2.xaxis.tick_top()
# ax2.xaxis.set_label_position('top')
ax2.yaxis.set_label_position('right')
fig.colorbar(a2,ax=ax2,label='Height in km',location='bottom')#,pad=-0.02
ax2.set_aspect('equal')
ax2.tick_params(left=False, labelleft=False, right=True, labelright=True)
ax2.set_title('Topography')

a3 = ax3.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['grid_sed'][1:-1,1:-1]/1000,cmap="YlGn_r")
ax3.set_xlabel('Easting in km')
ax3.set_ylabel('Northing in km')
fig.colorbar(a3,ax=ax3,label='km',location='bottom')
ax3.set_aspect('equal')
# ax3.tick_params(top=True, labeltop=False, bottom=False, labelbottom=False)
ax3.set_title('Sediment thickness')

a4 = ax4.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['grid_A_for_topo'][1:-1,1:-1]/1e-6,cmap="magma")
ax4.set_xlabel('Easting in km')
ax4.set_ylabel('Northing in km')
ax4.yaxis.set_label_position('right')
fig.colorbar(a4,ax=ax4,label='µW/m$^3$',location='bottom')
ax4.set_aspect('equal')
ax4.tick_params(left=False, labelleft=False, right=True, labelright=True)
ax4.set_title('Heat production')

# fig.text(0.165, 0.995, "A", ha='left', va='top', fontsize=14)
# fig.text(0.815, 0.995, "B", ha='left', va='top', fontsize=14)
# fig.text(0.165, 0.5, "C", ha='left', va='top', fontsize=14)
# fig.text(0.815, 0.5, "D", ha='left', va='top', fontsize=14)

fig.text(0.175, 0.995, "A", ha='left', va='top', fontsize=14)
fig.text(0.805, 0.995, "B", ha='left', va='top', fontsize=14)
fig.text(0.175, 0.5, "C", ha='left', va='top', fontsize=14)
fig.text(0.805, 0.5, "D", ha='left', va='top', fontsize=14)


#%% 3 without feature

fig = plt.figure(figsize=(7, 8), dpi=500,layout='constrained')
# ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(8,8),sharex=True,sharey=True,dpi=300,layout='compressed')
gs = fig.add_gridspec(nrows=2, ncols=2)

ax1 = fig.add_subplot(gs[0, 0])
a1 = ax1.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['shf_3D'][1:-1,1:-1],vmin=40,vmax=70,cmap="afmhot")
ax1.set_xlabel('Easting in km')
ax1.set_ylabel('Northing in km')
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top')
# ax1.yaxis.set_major_formatter(formatter)
ax1.set_aspect('equal')
ax1.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
# ax1.set_xlim(data_plot_3D['x_int_t'][1]/1000000,data_plot_3D['x_int_t'][-1]/1000000)
ax1.set_ylim(data_plot_3D['y_int_t'][1]/1000,data_plot_3D['y_int_t'][-1]/1000)
# ax1.colorbar(label='mW/m$^2$')

ax2 = fig.add_subplot(gs[0, 1])
a2 = ax2.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['shf_1D'][1:-1,1:-1],vmin=40,vmax=70,cmap="afmhot")
ax2.set_xlabel('Easting in km')
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position('top')
# ax2.set_ylabel('y in m')
fig.colorbar(a2,ax=[ax1,ax2],label='GHF in mW/m$^2$',location='bottom',shrink=0.65)#0.65
ax2.set_aspect('equal')
ax2.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, labelleft=False)

ax3 = fig.add_subplot(gs[1, :])
a3 = ax3.pcolormesh(data_plot_3D['x_int_t'][1:-1]/1000,data_plot_3D['y_int_t'][1:-1]/1000,data_plot_3D['shf_3D'][1:-1,1:-1]-data_plot_3D['shf_1D'][1:-1,1:-1],vmin=-10,vmax=10,cmap="bwr")
# ax3.set_xlabel('x in m')
ax3.set_ylabel('Northing in km')
fig.colorbar(a3,ax=ax3,label='difference in mW/m$^2$',location='bottom',shrink=0.49)
ax3.set_aspect('equal')
ax3.tick_params(top=True, labeltop=False, bottom=False, labelbottom=False)

fig.text(0.115, 0.98, "A", ha='left', va='top', fontsize=14)
fig.text(0.965, 0.98, "B", ha='left', va='top', fontsize=14)
fig.text(0.35, 0.495, "C", ha='left', va='top', fontsize=14)

# fig.subplots_adjust(wspace=0.05, hspace=0.1)
# 
# fig.savefig('abb/3d_ASB.pdf',bbox_inches='tight')
# fig.savefig('abb/3d_ASB.png',bbox_inches='tight')


#%% intersection points of profiles

# area_coord_LV = [1300000,1400000,-275000,-375000] # Lake Vostok "south-western" edge
# area_coord_DC = [1355000,1435000,-855000,-985000] # Dome C
area_coord_ASB = [1800000,1900000,-650000,-750000] # ASB

# profile_coord_LV_1=([1300000,1400000],[-320000,-300000]) # LV prof 1
# profile_coord_LV_2=([1300000,1400000],[-340000,-320000]) # LV prof 2
# profile_coord_LV_3=([1325000,1365000],[-275000,-375000]) # LV prof 3

# profile_coord_DC_1=([1405000,1435000],[-985000,-855000]) # DC prof 1
# profile_coord_DC_2=([1360000,1415000],[-985000,-855000]) # DC prof 2
# profile_coord_DC_3=([1355000,1435000],[-905000,-970000]) # DC prof 3
# profile_coord_DC_4=([1355000,1420000],[-920000,-985000]) # DC prof 4

profile_coord_ASB_1=([1900000,1800000],[-650000,-750000]) # ASB prof 1
profile_coord_ASB_2=([1800000,1900000],[-700000,-700000]) # ASB prof 2

# A = (profile_coord_LV_1[0][0],profile_coord_LV_1[1][0])
# B = (profile_coord_LV_1[0][1],profile_coord_LV_1[1][1])
# C = (profile_coord_LV_2[0][0],profile_coord_LV_2[1][0])
# D = (profile_coord_LV_2[0][1],profile_coord_LV_2[1][1])
# A = (profile_coord_DC_1[0][0],profile_coord_DC_1[1][0])
# B = (profile_coord_DC_1[0][1],profile_coord_DC_1[1][1])
# C = (profile_coord_DC_2[0][0],profile_coord_DC_1[1][0])
# D = (profile_coord_DC_2[0][1],profile_coord_DC_1[1][1])
A = (profile_coord_ASB_1[0][0],profile_coord_ASB_1[1][0])
B = (profile_coord_ASB_1[0][1],profile_coord_ASB_1[1][1])
C = (profile_coord_ASB_2[0][0],profile_coord_ASB_1[1][0])
D = (profile_coord_ASB_2[0][1],profile_coord_ASB_1[1][1])

# E = (profile_coord_LV_3[0][0],profile_coord_LV_3[1][0])
# F = (profile_coord_LV_3[0][1],profile_coord_LV_3[1][1])
# E = (profile_coord_DC_3[0][0],profile_coord_DC_3[1][0])
# F = (profile_coord_DC_3[0][1],profile_coord_DC_3[1][1])

# G = (profile_coord_DC_4[0][0],profile_coord_DC_4[1][0])
# H = (profile_coord_DC_4[0][1],profile_coord_DC_4[1][1])

line1 = LineString([A, B])
line2 = LineString([C, D])

# line3 = LineString([A, B])
# line4 = LineString([G, H])

int_pt = line1.intersection(line2)
POI = int_pt.x, int_pt.y

# int_pt_2 = line3.intersection(line4)
# POI_2 = int_pt_2.x, int_pt_2.y

print(POI)
# print(POI_2)

#%% 2d/3d with funky pcolormesh in 2nd plot for A

x_shf_3D = data_plot_3D['x_int_t']
y_shf_3D = data_plot_3D['y_int_t']
x_shf_2D = data_plot_2D['x_int_t']
y_shf_2D = data_plot_2D['y_int_t']

shf_3D = data_plot_3D['shf_3D']
profile_2D = data_plot_2D['shf_2D']
profile_1D = data_plot_2D['shf_1D']
dist_prof_t = data_plot_2D['dist_prof_t']

interp_func = RegularGridInterpolator((y_shf_3D,x_shf_3D),shf_3D)
profile_3D = interp_func((y_shf_2D,x_shf_2D),method="linear")

idx,idy = np.argmin(np.abs(x_shf_2D-POI[0])),np.argmin(np.abs(y_shf_2D-POI[1]))
if idx == idy:
    print(profile_2D[idx])
    print(dist_prof_t[idx])
    poi_p = dist_prof_t[idx]
else:
    print('idk you fucked somethign up')
    print(profile_2D[idx])

# idx_2,idy_2 = np.argmin(np.abs(x_shf_2D-POI_2[0])),np.argmin(np.abs(y_shf_2D-POI_2[1]))
# if idx_2 == idy_2:
#     print(profile_2D[idx_2])
#     print(dist_prof_t[idx_2])
#     poi_p_2 = dist_prof_t[idx_2]
# else:
#     print('idk you fucked somethign up')
#     print(profile_2D[idx_2])

fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(8,4),sharex=True,dpi=300,layout='constrained')#,layout='constrained'
# fig.set_constrained_layout_pads(h_pad=0.1, hspace=0.3)

ax1.plot(dist_prof_t[1:-1]/1000,profile_1D[1:-1][::-1],color='cornflowerblue',label="1D") 
ax1.plot(dist_prof_t[1:-1]/1000,profile_2D[1:-1][::-1],linestyle=(0, (1, 1)),color='crimson',label="2D")
ax1.plot(dist_prof_t[2:-1]/1000,profile_3D[2:-1][::-1],linestyle=(0, (5, 1)),color='seagreen',label="3D")
ax1.set_ylabel('GHF in mW/m$^2$')
ax1.set_xlim(0,np.max(data_plot_2D['dist_prof_t']/1000))
ax1.grid(alpha=0.5)
ax1.axvline(dist_prof_t[idx]/1000,color='k',linestyle='--',linewidth=.5)
ax1.annotate('ASB_C1', xy=(dist_prof_t[idx]/1000,ax1.get_ylim()[1]+0.3), xycoords='data', ha='center',annotation_clip=False)
# ax1.axvline(dist_prof_t[idx_2]/1000,color='k',linestyle='--',linewidth=.5)
# ax1.annotate('DC_C2', xy=(dist_prof_t[idx_2]/1000,ax1.get_ylim()[1]+0.3), xycoords='data', ha='center',annotation_clip=False)
# ax1.legend(loc='center right', bbox_to_anchor=(0.93, 0.72),bbox_transform=fig.transFigure)
legend = fig.legend(ax1.get_lines(),  # Get lines from ax1 for legend
                    [line.get_label() for line in ax1.get_lines()],  # Labels
                    loc='center right', bbox_to_anchor=(1, 0.74),  # Adjust position to the right of the plot
                    fontsize=10)

lns1 = ax2.plot(dist_prof_t/1000, data_plot_2D['profile_topo'][::-1]-data_plot_2D['profile_sed'][::-1], color='gold', label='Sediment depth')  # Second plot on left y-axis
lns2 = ax2.plot(dist_prof_t/1000, data_plot_2D['profile_topo'][::-1], color='slategrey', label='Topography')
ax2.set_ylabel("Height in m")
ax2.set_xlim(0,np.max(data_plot_2D['dist_prof_t']/1000))
ax2.grid(alpha=0.5)
ax3 = ax2.twinx()


profile_A_for_topo = np.interp(data_plot_2D['dist_prof_t'],data_plot_2D['dist_prof_A'],[data_plot_2D['profile_A'][i]/1e-6 for i in range(len(data_plot_2D['profile_A']))])
profile_A_for_topo = profile_A_for_topo[::-1]

ax2up = ax2.get_ylim()[1]
ax2low = ax2.get_ylim()[0]-1000

ax2.set_ylim(ax2.get_ylim()[0]-1000,ax2.get_ylim()[1])

ax3y = np.arange(np.round(ax2low,2),np.round(ax2up,3),1)
ax3meshx,ax3meshy = np.meshgrid(data_plot_2D['dist_prof_t']/1000,ax3y)
ax3data = np.zeros((ax3meshx.shape[0],ax3meshx.shape[1]))
for i in range(len(data_plot_2D['dist_prof_t'])):
    for j in range(len(ax3y)):
        if ax3y[j] <= (data_plot_2D['profile_topo'][::-1][i]-data_plot_2D['profile_sed'][::-1][i]):
            ax3data[j,i] = profile_A_for_topo[i]
        else:
            ax3data[j,i] = np.nan


# lns3 = ax3.plot(dist_prof_t/1000, profile_A_for_topo[::-1], color='darkorchid', label='Heat production')
lns3 = ax3.pcolormesh(ax3meshx,ax3meshy,ax3data,cmap='magma',alpha=1)
ax3.tick_params(right=False, labelright=False)
fig.colorbar(lns3,ax=ax3,label='µW/m$^3$',location='right',pad=0.01,panchor=False)
# ax3.set_ylabel("µW/m$^3$")
# ax2.set_zorder(ax3.get_zorder()+1)
lns = lns1+lns2#+lns3
labs = [l.get_label() for l in lns]
ax2.legend(lns,labs,loc="lower left",bbox_to_anchor=(0,0))#,bbox_to_anchor=(0.84,0)
ax2.set_xlabel("Profile length in km")

ax2.set_zorder(1)  # default zorder is 0 for ax1 and ax2
ax2.patch.set_visible(False)

fig.text(0.045, 0.98, "E", ha='left', va='top', fontsize=12)
fig.text(0.95, 0.98, "W", ha='right', va='top', fontsize=12)



#%% 2d/3d comparison with single layer

x_shf_3D = data_plot_3D['x_int_t']
y_shf_3D = data_plot_3D['y_int_t']
x_shf_2D = data_plot_2D['x_int_t']
y_shf_2D = data_plot_2D['y_int_t']

shf_3D = data_plot_3D['shf_3D']
profile_2D = data_plot_2D['shf_2D']
profile_1D = data_plot_2D['shf_1D']
dist_prof_t = data_plot_2D['dist_prof_t']

interp_func = RegularGridInterpolator((y_shf_3D,x_shf_3D),shf_3D)
profile_3D = interp_func((y_shf_2D,x_shf_2D),method="linear")

# idx_2,idy_2 = np.argmin(np.abs(x_shf_2D-POI_2[0])),np.argmin(np.abs(y_shf_2D-POI_2[1]))
# if idx_2 == idy_2:
#     print(profile_2D[idx_2])
#     print(dist_prof_t[idx_2])
#     poi_p_2 = dist_prof_t[idx_2]
# else:
#     print('idk you fucked somethign up')

fig, (ax1,ax2) = plt.subplots(2, 1,figsize=(8,4),dpi=300,sharex=True)#,layout='constrained'
# fig.set_constrained_layout_pads(h_pad=0.1, hspace=0.3)

ax1.plot(dist_prof_t[1:-1]/1000,profile_1D[1:-1],color='cornflowerblue',label="1D") 
ax1.plot(dist_prof_t[1:-1]/1000,profile_2D[1:-1],linestyle=(0, (1, 1)),color='crimson',label="2D")
ax1.plot(dist_prof_t[1:-1]/1000,profile_3D[1:-1],linestyle=(0, (5, 1)),color='seagreen',label="3D")
# ax1.set_xlabel("Profile length in km")
ax1.set_ylabel('GHF in mW/m$^2$')
ax1.set_xlim(0,np.max(data_plot_2D['dist_prof_t']/1000))
ax1.grid(alpha=0.5)
# ax1.axvline(dist_prof_t[idx]/1000,color='k',linestyle='--',linewidth=.5)
# ax1.annotate('ASB_C1', xy=(dist_prof_t[idx]/1000,ax1.get_ylim()[1]+0.4), xycoords='data', ha='center',annotation_clip=False)
# ax1.axvline(dist_prof_t[idx_2]/1000,color='k',linestyle='--',linewidth=.5)
# ax1.annotate('DC_C3', xy=(dist_prof_t[idx_2]/1000,ax1.get_ylim()[1]+0.4), xycoords='data', ha='center',annotation_clip=False)
# ax1.legend(loc='lower right')
legend = fig.legend(ax1.get_lines(),  # Get lines from ax1 for legend
                    [line.get_label() for line in ax1.get_lines()],  # Labels
                    loc='center right', bbox_to_anchor=(1, 0.74),  # Adjust position to the right of the plot
                    fontsize=10)


# lns1 = ax2.plot(dist_prof_t/1000, data_plot_2D['profile_topo']/1000-data_plot_2D['profile_sed']/1000, color='gold', label='Sediment depth')  # Second plot on left y-axis
# lns1 = ax2.plot(data_plot_2D['dist_prod_l']/1000, -data_plot_2D['profile_lab']/1000, color='orangered', label='Sediment depth')  # Second plot on left y-axis
# lns1 = ax2.plot(data_plot_2D['dist_prof_m']/1000, -data_plot_2D['profile_moho']/1000, color='darkorange', label='Sediment depth')  # Second plot on left y-axis
profile_A_for_topo = np.interp(data_plot_2D['dist_prof_t'],data_plot_2D['dist_prof_A'],[data_plot_2D['profile_A'][i]/1e-6 for i in range(len(data_plot_2D['profile_A']))])
lns1 = ax2.plot(dist_prof_t/1000, profile_A_for_topo, color='darkorchid', label='Sediment depth')  # Second plot on left y-axis
# lns2 = ax2.plot(dist_prof_t/1000, data_plot_2D['profile_topo']/1000, color='slategrey', label='Topography')
# ax2.set_ylabel("Height in km")
# ax2.set_ylabel("Depth in km")
ax2.set_ylabel("µM/m$^3$")
ax2.set_xlim(0,np.max(data_plot_2D['dist_prof_t']/1000))
ax2.grid(alpha=0.5)
# ax3 = ax2.twinx()
# profile_A_for_topo = np.interp(data_plot_2D['dist_prof_t'],data_plot_2D['dist_prof_A'],[data_plot_2D['profile_A'][i]/1e-6 for i in range(len(data_plot_2D['profile_A']))])
# lns3 = ax3.plot(dist_prof_t/1000, profile_A_for_topo[::-1], color='darkorchid', label='Heat production')
# ax3.set_ylabel("µW/m$^3$")
# lns = lns1+lns2+lns3
# labs = [l.get_label() for l in lns]
# ax2.legend(lns,labs,loc="center right")#,bbox_to_anchor=(0.84,0)
ax2.set_xlabel("Profile length in km")

fig.text(0.08, 0.95, "SW", ha='left', va='top', fontsize=12)
fig.text(0.94, 0.95, "NE", ha='right', va='top', fontsize=12)

# fig.savefig('abb/ASB_P1.pdf',bbox_inches='tight')
# fig.savefig('abb/ASB_P1.png',bbox_inches='tight')


#%% 2d/3d comparison

x_shf_3D = data_plot_3D['x_int_t']
y_shf_3D = data_plot_3D['y_int_t']
x_shf_2D = data_plot_2D['x_int_t']
y_shf_2D = data_plot_2D['y_int_t']

shf_3D = data_plot_3D['shf_3D']
profile_2D = data_plot_2D['shf_2D']
profile_1D = data_plot_2D['shf_1D']
dist_prof_t = data_plot_2D['dist_prof_t']

interp_func = RegularGridInterpolator((y_shf_3D,x_shf_3D),shf_3D)
profile_3D = interp_func((y_shf_2D,x_shf_2D),method="linear")

idx,idy = np.argmin(np.abs(x_shf_2D-POI[0])),np.argmin(np.abs(y_shf_2D-POI[1]))
if idx == idy:
    print(profile_2D[idx])
    print(dist_prof_t[idx])
    poi_p = dist_prof_t[idx]
else:
    print('idk you fucked somethign up')
    print(profile_2D[idx])

# idx_2,idy_2 = np.argmin(np.abs(x_shf_2D-POI_2[0])),np.argmin(np.abs(y_shf_2D-POI_2[1]))
# if idx_2 == idy_2:
#     print(profile_2D[idx_2])
#     print(dist_prof_t[idx_2])
#     poi_p_2 = dist_prof_t[idx_2]
# else:
#     print('idk you fucked somethign up')

fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(8,4),sharex=True,dpi=300)#,layout='constrained'
# fig.set_constrained_layout_pads(h_pad=0.1, hspace=0.3)

ax1.plot(dist_prof_t[1:-1]/1000,profile_1D[1:-1][::-1],label="1D") 
ax1.plot(dist_prof_t[1:-1]/1000,profile_2D[1:-1][::-1],label="2D")
ax1.plot(dist_prof_t[3:-1]/1000,profile_3D[3:-1][::-1],label="3D")
ax1.set_ylabel('GHF in mW/m$^2$')
ax1.set_xlim(0,np.max(data_plot_2D['dist_prof_t']/1000))
ax1.grid(alpha=0.5)
ax1.axvline(dist_prof_t[idx]/1000,color='k',linestyle='--',linewidth=.5)
ax1.annotate('LV_C1', xy=(dist_prof_t[idx]/1000,ax1.get_ylim()[1]+0.4), xycoords='data', ha='center',annotation_clip=False)
# ax1.axvline(dist_prof_t[idx_2]/1000,color='k',linestyle='--',linewidth=.5)
# ax1.annotate('DC_C3', xy=(dist_prof_t[idx_2]/1000,ax1.get_ylim()[1]+0.4), xycoords='data', ha='center',annotation_clip=False)
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

lns1 = ax2.plot(dist_prof_t/1000, data_plot_2D['profile_topo'][::-1]-data_plot_2D['profile_sed'][::-1], color='gold', label='Sediment depth')  # Second plot on left y-axis
lns2 = ax2.plot(dist_prof_t/1000, data_plot_2D['profile_topo'][::-1], color='slategrey', label='Topography')
ax2.set_ylabel("Height in m")
ax2.set_xlim(0,np.max(data_plot_2D['dist_prof_t']/1000))
ax2.grid(alpha=0.5)
ax3 = ax2.twinx()
profile_A_for_topo = np.interp(data_plot_2D['dist_prof_t'],data_plot_2D['dist_prof_A'],[data_plot_2D['profile_A'][i]/1e-6 for i in range(len(data_plot_2D['profile_A']))])
lns3 = ax3.plot(dist_prof_t/1000, profile_A_for_topo[::-1], color='darkorchid', label='Heat production')
ax3.set_ylabel("µW/m$^3$")
lns = lns1+lns2+lns3
labs = [l.get_label() for l in lns]
ax2.legend(lns,labs,loc="center right")#,bbox_to_anchor=(0.84,0)
ax2.set_xlabel("Profile length in km")

fig.text(0.05, 0.95, "E", ha='left', va='top', fontsize=12)
fig.text(0.95, 0.95, "W", ha='right', va='top', fontsize=12)

# fig.savefig('abb/ASB_P1.pdf',bbox_inches='tight')
# fig.savefig('abb/ASB_P1.png',bbox_inches='tight')

#%% 2d/3d with funky pcolormesh in 2nd plot for A

x_shf_3D = data_plot_3D['x_int_t']
y_shf_3D = data_plot_3D['y_int_t']
x_shf_2D = data_plot_2D['x_int_t']
y_shf_2D = data_plot_2D['y_int_t']

shf_3D = data_plot_3D['shf_3D']
profile_2D = data_plot_2D['shf_2D']
profile_1D = data_plot_2D['shf_1D']
dist_prof_t = data_plot_2D['dist_prof_t']

interp_func = RegularGridInterpolator((y_shf_3D,x_shf_3D),shf_3D)
profile_3D = interp_func((y_shf_2D,x_shf_2D),method="linear")

idx,idy = np.argmin(np.abs(x_shf_2D-POI[0])),np.argmin(np.abs(y_shf_2D-POI[1]))
if idx == idy:
    print(profile_2D[idx])
    print(dist_prof_t[idx])
    poi_p = dist_prof_t[idx]
else:
    print('idk you fucked somethign up')
    print(profile_2D[idx])

# idx_2,idy_2 = np.argmin(np.abs(x_shf_2D-POI_2[0])),np.argmin(np.abs(y_shf_2D-POI_2[1]))
# if idx_2 == idy_2:
#     print(profile_2D[idx_2])
#     print(dist_prof_t[idx_2])
#     poi_p_2 = dist_prof_t[idx_2]
# else:
#     print('idk you fucked somethign up')

fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(8,4),sharex=True,dpi=300,layout='constrained')#,layout='constrained'
# fig.set_constrained_layout_pads(h_pad=0.1, hspace=0.3)

ax1.plot(dist_prof_t[1:-1]/1000,profile_1D[1:-1][::-1],label="1D") 
ax1.plot(dist_prof_t[1:-1]/1000,profile_2D[1:-1][::-1],label="2D")
ax1.plot(dist_prof_t[3:-1]/1000,profile_3D[3:-1][::-1],label="3D")
ax1.set_ylabel('GHF in mW/m$^2$')
ax1.set_xlim(0,np.max(data_plot_2D['dist_prof_t']/1000))
ax1.grid(alpha=0.5)
ax1.axvline(dist_prof_t[idx]/1000,color='k',linestyle='--',linewidth=.5)
ax1.annotate('ASB_C1', xy=(dist_prof_t[idx]/1000,ax1.get_ylim()[1]+0.3), xycoords='data', ha='center',annotation_clip=False)
# ax1.axvline(dist_prof_t[idx_2]/1000,color='k',linestyle='--',linewidth=.5)
# ax1.annotate('DC_C3', xy=(dist_prof_t[idx_2]/1000,ax1.get_ylim()[1]+0.4), xycoords='data', ha='center',annotation_clip=False)
# ax1.legend(loc='center right', bbox_to_anchor=(0.93, 0.72),bbox_transform=fig.transFigure)
legend = fig.legend(ax1.get_lines(),  # Get lines from ax1 for legend
                    [line.get_label() for line in ax1.get_lines()],  # Labels
                    loc='center right', bbox_to_anchor=(0.99, 0.74),  # Adjust position to the right of the plot
                    fontsize=10)

lns1 = ax2.plot(dist_prof_t/1000, data_plot_2D['profile_topo'][::-1]-data_plot_2D['profile_sed'][::-1], color='gold', label='Sediment depth')  # Second plot on left y-axis
lns2 = ax2.plot(dist_prof_t/1000, data_plot_2D['profile_topo'][::-1], color='slategrey', label='Topography')
ax2.set_ylabel("Height in m")
ax2.set_xlim(0,np.max(data_plot_2D['dist_prof_t']/1000))
ax2.grid(alpha=0.5)
ax3 = ax2.twinx()


profile_A_for_topo = np.interp(data_plot_2D['dist_prof_t'],data_plot_2D['dist_prof_A'],[data_plot_2D['profile_A'][i]/1e-6 for i in range(len(data_plot_2D['profile_A']))])

ax2up = ax2.get_ylim()[1]
ax2low = ax2.get_ylim()[0]-1000

ax2.set_ylim(ax2.get_ylim()[0]-1000,ax2.get_ylim()[1])

ax3y = np.arange(np.round(ax2low,2),np.round(ax2up,3),1)
ax3meshx,ax3meshy = np.meshgrid(data_plot_2D['dist_prof_t']/1000,ax3y)
ax3data = np.zeros((ax3meshx.shape[0],ax3meshx.shape[1]))
for i in range(len(data_plot_2D['dist_prof_t'])):
    for j in range(len(ax3y)):
        if ax3y[j] <= (data_plot_2D['profile_topo'][::-1][i]-data_plot_2D['profile_sed'][::-1][i]):
            ax3data[j,i] = profile_A_for_topo[::-1][i]
        else:
            ax3data[j,i] = np.nan


# lns3 = ax3.plot(dist_prof_t/1000, profile_A_for_topo[::-1], color='darkorchid', label='Heat production')
lns3 = ax3.pcolormesh(ax3meshx,ax3meshy,ax3data,cmap='magma',alpha=1)
ax3.tick_params(right=False, labelright=False)
fig.colorbar(lns3,ax=ax3,label='µW/m$^3$',location='right',pad=0.01,panchor=False)
# ax3.set_ylabel("µW/m$^3$")
# ax2.set_zorder(ax3.get_zorder()+1)
lns = lns1+lns2#+lns3
labs = [l.get_label() for l in lns]
ax2.legend(lns,labs,loc="center left",bbox_to_anchor=(0,0.6))#,bbox_to_anchor=(0.84,0)
ax2.set_xlabel("Profile length in km")

ax2.set_zorder(1)  # default zorder is 0 for ax1 and ax2
ax2.patch.set_visible(False)

fig.text(0.045, 0.97, "E", ha='left', va='top', fontsize=12)
fig.text(0.93, 0.97, "W", ha='right', va='top', fontsize=12)

#%% 2d/3d comparison without anything else (only for across LV)

x_shf_3D = data_plot_3D['x_int_t']
y_shf_3D = data_plot_3D['y_int_t']
x_shf_2D = data_plot_2D['x_int_t']
y_shf_2D = data_plot_2D['y_int_t']

shf_3D = data_plot_3D['shf_3D']
profile_2D = data_plot_2D['shf_2D']
profile_1D = data_plot_2D['shf_1D']
dist_prof_t = data_plot_2D['dist_prof_t']

interp_func = RegularGridInterpolator((y_shf_3D,x_shf_3D),shf_3D)
profile_3D = interp_func((y_shf_2D,x_shf_2D),method="linear")

# idx_2,idy_2 = np.argmin(np.abs(x_shf_2D-POI_2[0])),np.argmin(np.abs(y_shf_2D-POI_2[1]))
# if idx_2 == idy_2:
#     print(profile_2D[idx_2])
#     print(dist_prof_t[idx_2])
#     poi_p_2 = dist_prof_t[idx_2]
# else:
#     print('idk you fucked somethign up')

fig, ax1 = plt.subplots(1, 1,figsize=(8,3),dpi=300)#,layout='constrained'
# fig.set_constrained_layout_pads(h_pad=0.1, hspace=0.3)

ax1.plot(dist_prof_t[1:-1]/1000,profile_1D[1:-1],label="1D") 
ax1.plot(dist_prof_t[1:-1]/1000,profile_2D[1:-1],label="2D")
ax1.plot(dist_prof_t[3:-1]/1000,profile_3D[3:-1],label="3D")
ax1.set_xlabel("Profile length in km")
ax1.set_ylabel('GHF in mW/m$^2$')
ax1.set_xlim(0,np.max(data_plot_2D['dist_prof_t']/1000))
ax1.grid(alpha=0.5)
# ax1.axvline(dist_prof_t[idx]/1000,color='k',linestyle='--',linewidth=.5)
# ax1.annotate('ASB_C1', xy=(dist_prof_t[idx]/1000,ax1.get_ylim()[1]+0.4), xycoords='data', ha='center',annotation_clip=False)
# ax1.axvline(dist_prof_t[idx_2]/1000,color='k',linestyle='--',linewidth=.5)
# ax1.annotate('DC_C3', xy=(dist_prof_t[idx_2]/1000,ax1.get_ylim()[1]+0.4), xycoords='data', ha='center',annotation_clip=False)
ax1.legend(loc='upper right')

# lns1 = ax2.plot(dist_prof_t/1000, data_plot_2D['profile_topo'][::-1]-data_plot_2D['profile_sed'][::-1], color='gold', label='Sediment depth')  # Second plot on left y-axis
# lns2 = ax2.plot(dist_prof_t/1000, data_plot_2D['profile_topo'][::-1], color='slategrey', label='Topography')
# ax2.set_ylabel("Height in m")
# ax2.set_xlim(0,np.max(data_plot_2D['dist_prof_t']/1000))
# ax2.grid(alpha=0.5)
# ax3 = ax2.twinx()
# profile_A_for_topo = np.interp(data_plot_2D['dist_prof_t'],data_plot_2D['dist_prof_A'],[data_plot_2D['profile_A'][i]/1e-6 for i in range(len(data_plot_2D['profile_A']))])
# lns3 = ax3.plot(dist_prof_t/1000, profile_A_for_topo[::-1], color='darkorchid', label='Heat production')
# ax3.set_ylabel("µW/m$^3$")
# lns = lns1+lns2+lns3
# labs = [l.get_label() for l in lns]
# ax2.legend(lns,labs,loc="center right")#,bbox_to_anchor=(0.84,0)
# ax2.set_xlabel("Profile length in km")

fig.text(0.05, 0.95, "S", ha='left', va='top', fontsize=12)
fig.text(0.95, 0.95, "N", ha='right', va='top', fontsize=12)

# fig.savefig('abb/ASB_P1.pdf',bbox_inches='tight')
# fig.savefig('abb/ASB_P1.png',bbox_inches='tight')

#%% east antarctica

area_coord_LV = [1300000,1400000,-275000,-375000] # Lake Vostok "south-western" edge
area_coord_DC = [1355000,1435000,-855000,-985000] # Dome C
area_coord_ASB = [1800000,1900000,-650000,-750000] # ASB

profile_coord=([1300000,1400000],[-275000,-375000])

plt.figure(dpi=300)
ax = plt.axes(projection=ccrs.AzimuthalEquidistant(central_latitude=-90))
# sc = ax.pcolormesh(data_calc[1][0],data_calc[1][1],data_calc[1][2]/1000,cmap='YlGn_r',vmin=0,vmax=6)
# rect = patches.Rectangle((area_coord_LV[0], area_coord_LV[3]), area_coord_LV[1]-area_coord_LV[0], area_coord_LV[2]-area_coord_LV[3], linewidth=1, edgecolor='firebrick', facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((area_coord_DC[0], area_coord_DC[3]), area_coord_DC[1]-area_coord_DC[0], area_coord_DC[2]-area_coord_DC[3], linewidth=1, edgecolor='gold', facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((area_coord_ASB[0], area_coord_ASB[3]), area_coord_ASB[1]-area_coord_ASB[0], area_coord_ASB[2]-area_coord_ASB[3], linewidth=1, edgecolor='darkorange', facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((area_coord_LV[0], area_coord_LV[3]), area_coord_LV[1]-area_coord_LV[0], area_coord_LV[2]-area_coord_LV[3], linewidth=1, edgecolor='firebrick', facecolor='none')
# rect = patches.Rectangle((area_coord_ASB[0], area_coord_ASB[3]), area_coord_ASB[1]-area_coord_ASB[0], area_coord_ASB[2]-area_coord_ASB[3], linewidth=1, edgecolor='darkorange', facecolor='none')
rect = patches.Rectangle((area_coord_DC[0], area_coord_DC[3]), area_coord_DC[1]-area_coord_DC[0], area_coord_DC[2]-area_coord_DC[3], linewidth=1, edgecolor='gold', facecolor='none')
ax.add_patch(rect)
# pr = ax.plot(profile_coord[0],profile_coord[1],'m')
ax.set_extent([250000,2750000,-2200000,-38000],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, rotate_labels = False)
gl.top_labels = False
gl.right_labels = False
gl.bottom_labels = False
gl.left_labels = False
# gl.xlabel_style = {'size': 8}
# gl.ylabel_style = {'rotation': 45, 'rotation_mode': 'anchor','size': 5}
# gl.xformatter = LongitudeFormatter()
# gl.yformatter = LatitudeFormatter()
# ax.contour(data_calc[1][0],data_calc[1][1],data_calc[1][2],levels=[0,1500,3000,4500],colors='k',linewidths=0.5) # sed
# ax.contour(data_calc[0][0],data_calc[0][1],data_calc[0][2],levels=[-1500,-1000,-500,0,500,1000,1500],colors='k',linewidths=0.5) # topo
# cbar = plt.colorbar(cm.ScalarMappable(norm=norm,cmap='YlGn_r'), ax=ax,fraction=0.038, pad=0.04)
# cbar = plt.colorbar(sc, ax=ax,fraction=0.038, pad=0.04)
# cbar.set_label("Sediment depth in km")
# cbar.formatter.set_powerlimits((-2, 2))
# cbar.formatter.set_useMathText(True)
# cbar.locator = ticker.MaxNLocator(nbins=7)
# cbar.update_ticks()
ax.coastlines(resolution='10m',linewidth=0.5)
plt.show()

#%% small areas

# area_coord_LV = [1300000,1400000,-275000,-375000] # Lake Vostok "south-western" edge
# area_coord_DC = [1355000,1435000,-855000,-985000] # Dome C
area_coord_ASB = [1800000,1900000,-650000,-750000] # ASB

# profile_coord_LV_1=([1300000,1400000],[-320000,-300000]) # LV prof 1
# profile_coord_LV_2=([1300000,1400000],[-340000,-320000]) # LV prof 2
# profile_coord_LV_3=([1325000,1365000],[-275000,-375000]) # LV prof 3

# profile_coord_DC_1=([1405000,1435000],[-985000,-855000]) # DC prof 1
# profile_coord_DC_2=([1360000,1415000],[-985000,-855000]) # DC prof 2
# profile_coord_DC_3=([1355000,1420000],[-920000,-985000]) # DC prof 3
# profile_coord_DC_4=([1355000,1435000],[-905000,-970000]) # DC prof 3

profile_coord_ASB_1=([1900000,1800000],[-650000,-750000]) # ASB prof 1
profile_coord_ASB_2=([1800000,1900000],[-700000,-700000]) # ASB prof 2

plt.figure(dpi=300)
ax = plt.axes(projection=ccrs.AzimuthalEquidistant(central_latitude=-90))
sc = ax.pcolormesh(data_plot_3D['x_int_t'],data_plot_3D['y_int_t'],data_plot_3D['grid_sed']/1000,cmap='YlGn_r')
pr = ax.plot(profile_coord_ASB_1[0],profile_coord_ASB_1[1],'gold')
pr2 = ax.plot(profile_coord_ASB_2[0],profile_coord_ASB_2[1],'firebrick')
# pr3 = ax.plot(profile_coord_DC_3[0],profile_coord_DC_3[1],'darkorange')
# pr4 = ax.plot(profile_coord_DC_4[0],profile_coord_DC_4[1],'darkviolet')
# ax.set_extent([area_coord_LV[0],area_coord_LV[1],area_coord_LV[2],area_coord_LV[3]],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
# ax.set_extent([area_coord_DC[0],area_coord_DC[1],area_coord_DC[2],area_coord_DC[3]],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
ax.set_extent([area_coord_ASB[0],area_coord_ASB[1],area_coord_ASB[2],area_coord_ASB[3]],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, rotate_labels = False)
gl.top_labels = False
gl.right_labels = False
gl.ylabel_style = {'rotation': 45, 'rotation_mode': 'anchor'}
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()
cbar = plt.colorbar(sc,ax=ax,fraction=0.038, pad=0.04)#cm.ScalarMappable(norm=norm,cmap='YlGn_r'),
cbar.set_label("Sediment thickness in km")
cbar.formatter.set_powerlimits((-2, 2))
cbar.formatter.set_useMathText(True)
cbar.locator = ticker.MaxNLocator(nbins=7)
cbar.update_ticks()
ax.coastlines(resolution='10m',linewidth=0.5)
plt.show()
#%% 2d results

def shf_2D_with_feature(data,a,b,c,x2):
    x1 = 'dist_{}'.format(c[0])
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1,sharex=True,dpi=300)
    ax1.plot(data[x1], [a/1000 for a in data[c]])
    if c == 'A':
        ax1.set_ylabel(c+' in µW/m^3')
    else:
        ax1.set_ylabel(c+' in km')
    ax1.set_title('Influence of '+c)
    ax1.grid(True)

    ax2.plot(data[x2], data[a],label='2d')
    ax2.plot(data['dist_1d'], data[b],label='1d')
    ax2.set_ylabel('ghf in mW/m^2')
    ax2.legend()
    ax2.grid(True)
    
    diff = [data[a][i]-data[b][i] for i in range(len(data[b]))]
    ax3.plot(data[x1],diff,label='diff')
    ax3.set_xlabel('x in m')
    ax3.set_ylabel('2D-1D')
    ax3.grid(True)
    return data[a],diff

shf_2d, diff = shf_2D_with_feature(data_plot_2D,'shf','f_hf','topo','mspkt')

plt.figure(figsize=(8,3),dpi=300)
plt.plot(data_plot_2D['dist_prof_t'],data_plot_2D['shf_2D'],color='orangered',label="2D")
plt.plot(data_plot_2D['dist_prof_t'],data_plot_2D['shf_1D'],color='forestgreen',label="1D")  
plt.xlabel('x in m')
plt.xlim(0,141000)
plt.ylabel('GHF in mW/m$^2$')
plt.grid(alpha=0.7)
plt.legend()
plt.savefig('abb/full2d.pdf',bbox_inches='tight')


fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(8,4),sharex=True,dpi=300,layout='constrained')

ax1.plot(data_plot_2D['dist_prof_t'],data_plot_2D['shf_2D'],color='orangered',label="2D")
ax1.plot(data_plot_2D['dist_prof_t'],data_plot_2D['shf_1D'],color='forestgreen',label="1D")
# ax1.set_xlim(0,141000)
# ax1.set_ylim(50,65)
ax1.set_ylabel('GHF in mW/m$^2$')
ax1.grid(alpha=0.7)
ax1.legend(loc='upper right')

ax2.plot(data_plot_2D['dist_prof_t'],data_plot_2D['profile_topo']/1000,color='slategrey')
ax2.plot(data_plot_2D['dist_prof_t'],(data_plot_2D['profile_topo']-data_plot_2D['profile_sed'])/1000,color='gold')
# ax2.plot(data_plot_2D['dist_prof_t'],data_plot_2D['profile_sed']/1000,color='gold')
# ax2.plot(data_plot_2D['dist_prof_m'],-data_plot_2D['profile_moho']/1000,color='darkorange')
# ax2.plot(data_plot_2D['dist_prod_l'],-data_plot_2D['profile_lab']/1000,color='orangered')
# profile_A_for_topo = np.interp(data_plot_2D['dist_prof_t'],data_plot_2D['dist_prof_A'],data_plot_2D['profile_A'])
# ax2.plot(data_plot_2D['dist_prof_t'],profile_A_for_topo,color='darkorchid')
ax2.set_xlabel('x in m')
# ax2.set_xlim(0,141000)
ax2.set_ylabel('height in km')
ax2.grid(alpha=0.7)
# ax2.set_ylabel('thickness in km')
# ax2.set_ylabel('depth in km')
# ax2.set_ylabel('HP in W/m$^3$')
fig.savefig('abb/hp2d.pdf',bbox_inches='tight')
