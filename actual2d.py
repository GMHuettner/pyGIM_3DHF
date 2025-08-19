"""
Georg HÃ¼ttner - September 2024

This script performs two dimensional thermal modelling from provided lithospheric
interfaces. Based on this a surface geothermal heat flow can be calculated.
"""

import pygimli as pg
import pygimli.meshtools as mt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib import ticker
import cartopy.crs as ccrs
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter)
from scipy.interpolate import griddata,RegularGridInterpolator
import datetime
import pickle
from scipy.stats import norm
import gepy

#%% research area and params
area_coord_s = [1300000,1400000,-275000,-375000]
# profile_coord=([1300000,1400000],[-275000,-375000])

# profile_coord=([1300000,1400000],[-320000,-300000]) # LV prof 1
# profile_coord=([1300000,1400000],[-340000,-320000]) # LV prof 2
profile_coord=([1325000,1365000],[-275000,-375000]) # LV prof 3

# area_coord_s = [1355000,1435000,-855000,-985000] # Dome C

# profile_coord=([1405000,1435000],[-985000,-855000]) # DC prof 1
# profile_coord=([1360000,1415000],[-985000,-855000]) # DC prof 2
# profile_coord=([1355000,1420000],[-920000,-985000]) # DC prof 3
# profile_coord=([1355000,1435000],[-905000,-970000]) # DC prof 3
# 
# area_coord_s = [1800000,1900000,-650000,-750000] # ASB

# profile_coord=([1900000,1800000],[-650000,-750000]) # ASB prof 1
# profile_coord=([1800000,1900000],[-750000,-650000]) # ASB prof 1
# profile_coord=([1800000,1900000],[-700000,-700000]) # ASB prof 2

k0 = 1.5
k1 = 2.7
k2 = 3.5

km = 1000

notes = '2d_LV_3_Tprof'
print(notes)

#%% load data and interpolate on profiles

data_calc,data_resolution = gepy.load_data()

x_profile,y_profile = (profile_coord[0],profile_coord[1])

profile_length = np.round(np.sqrt((x_profile[1]-x_profile[0])**2 + (y_profile[1]-y_profile[0])**2))

x_int_t, y_int_t, dist_prof_t = gepy.define_profile(x_profile, y_profile, data_resolution[0])
x_int_s, y_int_s, dist_prof_s = gepy.define_profile(x_profile, y_profile, data_resolution[1])
x_int_m, y_int_m, dist_prof_m = gepy.define_profile(x_profile, y_profile, data_resolution[2])
x_int_l, y_int_l, dist_prof_l = gepy.define_profile(x_profile, y_profile, data_resolution[3])

dist_prof_t = np.round(dist_prof_t)
dist_prof_s = np.round(dist_prof_s)
dist_prof_m = np.round(dist_prof_m)
dist_prof_l = np.round(dist_prof_l)

xvals_t = np.unique(np.asarray(data_calc[0][0]))
yvals_t = np.unique(np.asarray(data_calc[0][1]))

xvals_s = np.unique(np.asarray(data_calc[1][0]))
yvals_s = np.unique(np.asarray(data_calc[1][1]))

xvals_m = np.unique(np.asarray(data_calc[2][0]))
yvals_m = np.unique(np.asarray(data_calc[2][1]))

xvals_l = np.unique(np.asarray(data_calc[3][0]))
yvals_l = np.unique(np.asarray(data_calc[3][1]))

interp_func = RegularGridInterpolator((yvals_t,xvals_t),np.asarray(data_calc[0][2]))
profile_topo = interp_func((y_int_t,x_int_t),method="linear").tolist()
profile_topo = np.round(profile_topo,decimals=3)

interp_func = RegularGridInterpolator((yvals_s,xvals_s),np.asarray(data_calc[1][2]))
profile_sed = interp_func((y_int_s,x_int_s),method="linear").tolist()
profile_sed = np.round(profile_sed,decimals=3)

interp_func = RegularGridInterpolator((yvals_m,xvals_m),np.asarray(data_calc[2][2]))
profile_moho = interp_func((y_int_m,x_int_m),method="linear").tolist()
profile_moho = np.round(profile_moho,decimals=3)

interp_func = RegularGridInterpolator((yvals_l,xvals_l),np.asarray(data_calc[3][2]))
profile_lab = interp_func((y_int_l,x_int_l),method="linear").tolist()
profile_lab = np.round(profile_lab,decimals=3)

profile_topo_for_sed = np.interp(dist_prof_s,dist_prof_t,profile_topo)
profile_sed_for_topo = np.interp(dist_prof_t,dist_prof_s,profile_sed)
profile_sed_for_topo = np.round(profile_sed_for_topo,decimals=-1)

profile_moho_for_topo = np.interp(dist_prof_t,dist_prof_m,profile_moho)
profile_lab_for_topo = np.interp(dist_prof_t,dist_prof_l,profile_lab)

# heat production:
data_HP = np.load("data/SMOS_hp.npz")
grid_A_ws = data_HP['A_ws']
grid_A_wo = data_HP['A_wo']

xi_smos,yi_smos = data_HP['xi_smos'],data_HP['yi_smos']
x_int_A, y_int_A, dist_prof_A = gepy.define_profile(x_profile, y_profile, 12500)
dist_prof_A = np.round(dist_prof_A)
xvals_A = np.unique(xi_smos)
yvals_A = np.unique(yi_smos)
interp_func = RegularGridInterpolator((yvals_A,xvals_A),np.asarray(grid_A_ws))
profile_A = interp_func((y_int_A,x_int_A),method="linear").tolist()
# profile_A = np.round(profile_A,decimals=3)
profile_A_for_topo = np.interp(dist_prof_t,dist_prof_A,profile_A)

# profile_A_norm = profile_A_for_topo/np.max(profile_A_for_topo)
# profile_A_for_topo = profile_A_for_topo*profile_A_norm*1.5


# single layer scheiße

xi_topo_s_a2, yi_topo_s_a2, grid_topo_s_a2 = gepy.in_area_s(area_coord_s,data_calc[0][0],data_calc[0][1],data_calc[0][2])
xi_sed_s_a2, yi_sed_s_a2, grid_sed_s_a2    = gepy.in_area_s(area_coord_s,data_calc[1][0],data_calc[1][1],data_calc[1][2])
xi_moho_s_a2, yi_moho_s_a2, grid_moho_s_a2 = gepy.in_area_s(area_coord_s,data_calc[2][0],data_calc[2][1],data_calc[2][2])
xi_lab_s_a2, yi_lab_s_a2, grid_lab_s_a2    = gepy.in_area_s(area_coord_s,data_calc[3][0],data_calc[3][1],data_calc[3][2])
xii_topo_s_a2,yii_topo_s_a2 = np.meshgrid(xi_topo_s_a2,yi_topo_s_a2)
xii_sed_s_a2,yii_sed_s_a2   = np.meshgrid(xi_sed_s_a2,yi_sed_s_a2)
xii_moho_s_a2,yii_moho_s_a2 = np.meshgrid(xi_moho_s_a2,yi_moho_s_a2)
xii_lab_s_a2,yii_lab_s_a2   = np.meshgrid(xi_lab_s_a2,yi_lab_s_a2)
res_topo,res_moho,res_lab = (1000,10000,10000)
x_int_t_a2,y_int_t_a2 = np.arange(area_coord_s[0],area_coord_s[1]+res_topo,res_topo),np.arange(area_coord_s[3],area_coord_s[2]+res_topo,res_topo)
x_int_m_a2,y_int_m_a2 = np.arange(area_coord_s[0],area_coord_s[1]+res_moho,res_moho),np.arange(area_coord_s[3],area_coord_s[2]+res_moho,res_moho)
x_int_l_a2,y_int_l_a2 = np.arange(area_coord_s[0],area_coord_s[1]+res_lab,res_lab),np.arange(area_coord_s[3],area_coord_s[2]+res_lab,res_lab)
xi_int_t_a2,yi_int_t_a2 = np.meshgrid(x_int_t_a2,y_int_t_a2)
xi_int_m_a2,yi_int_m_a2 = np.meshgrid(x_int_m_a2,y_int_m_a2)
xi_int_l_a2,yi_int_l_a2 = np.meshgrid(x_int_l_a2,y_int_l_a2)
interp_topo_a2 = RegularGridInterpolator((np.unique(xi_topo_s_a2),np.unique(yi_topo_s_a2)),grid_topo_s_a2)
grid_topo_a2 = interp_topo_a2((xi_int_t_a2,yi_int_t_a2))
interp_sed_a2 = RegularGridInterpolator((np.unique(xi_sed_s_a2),np.unique(yi_sed_s_a2)),grid_sed_s_a2)
grid_sed_a2 = interp_sed_a2((xi_int_t_a2,yi_int_t_a2),method='linear')
interp_moho_a2 = RegularGridInterpolator((np.unique(xi_moho_s_a2),np.unique(yi_moho_s_a2)),grid_moho_s_a2)
grid_moho_a2 = interp_moho_a2((xi_int_m_a2,yi_int_m_a2))
grid_moho_for_topo_a2 = interp_moho_a2((xi_int_t_a2,yi_int_t_a2))
interp_lab_a2 = RegularGridInterpolator((np.unique(xi_lab_s_a2),np.unique(yi_lab_s_a2)),grid_lab_s_a2)
grid_lab_a2 = interp_lab_a2((xi_int_l_a2,yi_int_l_a2))
grid_lab_for_topo_a2 = interp_lab_a2((xi_int_t_a2,yi_int_t_a2))
xi_smos_s_a2, yi_smos_s_a2, A_s_a2 = gepy.in_area_s(area_coord_s,xi_smos,yi_smos,grid_A_ws)
xii_smos_s_a2,yii_smos_s_a2   = np.meshgrid(xi_smos_s_a2,yi_smos_s_a2)    
res_HP = 12500
x_int_HP_a2,y_int_HP_a2 = np.arange(area_coord_s[0],area_coord_s[1]+res_HP,res_HP),  np.arange(area_coord_s[3],area_coord_s[2]+res_HP,res_HP)
xi_int_HP_a2,yi_int_HP_a2 = np.meshgrid(x_int_HP_a2,y_int_HP_a2)
grid_A_a2  = griddata((xii_smos_s_a2.flatten(), yii_smos_s_a2.flatten()), np.transpose(A_s_a2).flatten(), (xi_int_HP_a2,yi_int_HP_a2),fill_value=0.000001)
# interp_A = RegularGridInterpolator((x_int_HP,y_int_HP),np.transpose(grid_A))
interp_A_a2 = RegularGridInterpolator((x_int_HP_a2,y_int_HP_a2),np.transpose(grid_A_a2))
grid_A_for_topo_a2 = interp_A_a2((xi_int_t_a2,yi_int_t_a2))

# select which interfaces are part of the world

# profile_topo = np.full(len(profile_topo),np.mean(grid_topo_a2))
# profile_sed = np.full(len(profile_sed),np.mean(grid_sed_a2))
# profile_sed_for_topo = np.full(len(profile_sed_for_topo),np.mean(grid_sed_a2))
# profile_moho = np.full(len(profile_moho),np.mean(grid_moho_a2))
# profile_moho_for_topo = np.full(len(profile_topo),np.mean(grid_moho_for_topo_a2))
# profile_lab = np.full(len(profile_lab),np.mean(grid_lab_a2))
# profile_lab_for_topo = np.full(len(profile_topo),np.mean(grid_lab_for_topo_a2))
# profile_A = np.full(len(profile_A),np.mean(grid_A_a2))
# profile_A_for_topo = np.full(len(profile_A_for_topo),np.mean(grid_A_for_topo_a2))

profile_sed_depth = profile_topo-profile_sed_for_topo

#%% plot that shit

vmin = np.nanmin(data_calc[0][2])
vmax = np.nanmax(data_calc[0][2])

plt.figure(dpi=300)
ax = plt.axes(projection=ccrs.AzimuthalEquidistant(central_latitude=-90))
norm = Normalize(vmin=vmin, vmax=vmax)
# sc = ax.pcolormesh(data_calc[0][0],data_calc[0][1],data_calc[0][2],cmap='jet',vmin=vmin,vmax=vmax)
ax.set_extent([250000,2750000,-2200000,-38000],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, rotate_labels = False)
gl.top_labels = False
gl.right_labels = False
gl.ylabel_style = {'rotation': 45, 'rotation_mode': 'anchor'}
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()
pr = ax.plot(x_int_s,y_int_s,'m')
ax.contour(data_calc[1][0],data_calc[1][1],data_calc[1][2],levels=[0,1500,3000,4500],colors='k',linewidths=0.5)
# cbar = plt.colorbar(cm.ScalarMappable(norm=norm,cmap='jet'), ax=ax,fraction=0.038, pad=0.04)
# cbar.set_label("Sediment depth in m")
# cbar.formatter.set_powerlimits((-2, 2))
# cbar.formatter.set_useMathText(True)
# cbar.locator = ticker.MaxNLocator(nbins=7)
# cbar.update_ticks()
ax.coastlines(resolution='10m',linewidth=0.5)
plt.show()

#%% build the world

# create 2d mesh from the input data

world_start = [0,4*km] 
world_end = [profile_length,-220*km]

world = mt.createWorld(start=world_start,end=world_end, worldMarker=False)

topo_line = mt.createPolygon([[a,b] for a,b in zip(dist_prof_t,profile_topo)],marker = 5,boundaryMarker = 5,isClosed=False)
sed_line = mt.createPolygon([[a,b] for a,b in zip(dist_prof_t,profile_topo-profile_sed_for_topo)],marker = 6,boundaryMarker = 6,isClosed=False)
moho_line = mt.createPolygon([[a,b] for a,b in zip(dist_prof_m,-profile_moho)],marker = 7,boundaryMarker = 7,isClosed=False)
lab_line = mt.createPolygon([[a,b] for a,b in zip(dist_prof_l,-profile_lab)],marker = 8,boundaryMarker = 8,isClosed=False)

world = world + sed_line + moho_line + lab_line + topo_line

mesh = mt.createMesh(world,quality=34,area=profile_length*4)#,area=500000

# create lists for temperature and heat production that are to be applied to the nodes
force = np.zeros(mesh.nodeCount())
force_marker = np.zeros(mesh.nodeCount())
Tnode = []
for i, node in enumerate(mesh.nodes()):
    x,y,z = node.pos()
    
    idx_t_0,idx_t_1 = np.argsort(np.abs(dist_prof_t-x))[0],np.argsort(np.abs(dist_prof_t-x))[1]
    
    m_s = (profile_sed_depth[idx_t_0]-profile_sed_depth[idx_t_1])/(dist_prof_t[idx_t_0]-dist_prof_t[idx_t_1])
    y_s = m_s * (x-dist_prof_t[idx_t_1]) + profile_sed_depth[idx_t_1] 
    
    m_t = (profile_topo[idx_t_0]-profile_topo[idx_t_1])/(dist_prof_t[idx_t_0]-dist_prof_t[idx_t_1])
    y_t = m_t * (x-dist_prof_t[idx_t_1]) + profile_topo[idx_t_1] 
    
    idx_tA = np.argmin(np.abs(dist_prof_t-x))
    idx_bA = np.argmin(np.abs(dist_prof_m-x))
    
    idx_t = np.argmin(np.abs(dist_prof_t-x))
    idx_l = np.argmin(np.abs(dist_prof_l-x))
    
    idx_A = np.argmin(np.abs(dist_prof_t-x))
    
    if y >= -profile_moho[idx_bA] and y <= y_s:
        force[i] = profile_A_for_topo[idx_A]
        # force[i] = 0.000001
        force_marker[i] = node.marker()
    else:
        force[i] = 0
    
    if y >= y_t:
        T0fornode = [node.id(),0]
        Tnode.append(T0fornode)
    if y <= -profile_lab[idx_l]:
        TLABfornode = [node.id(),1315]
        Tnode.append(TLABfornode)

# apply correct marker to the cells
crust_cells = []
crust_nodes = []
for cell in mesh.cells():
    center = cell.center()
    x, y, z = center.x(), center.y(), center.z()
    
    idx_t_0,idx_t_1 = np.argsort(np.abs(dist_prof_t-x))[0],np.argsort(np.abs(dist_prof_t-x))[1]
    
    m_t = (profile_topo[idx_t_0]-profile_topo[idx_t_1])/(dist_prof_t[idx_t_0]-dist_prof_t[idx_t_1])
    m_s = (profile_sed_depth[idx_t_0]-profile_sed_depth[idx_t_1])/(dist_prof_t[idx_t_0]-dist_prof_t[idx_t_1])
    y_t = m_t * (x-dist_prof_t[idx_t_1]) + profile_topo[idx_t_1] 
    y_s = m_s * (x-dist_prof_t[idx_t_1]) + profile_sed_depth[idx_t_1] 
    
    idx_m = np.argmin(np.abs(dist_prof_m-x))
    idx_l = np.argmin(np.abs(dist_prof_l-x))
    
    if y < y_t and y > y_s:
        cell.setMarker(5)
    if y < y_s:
        cell.setMarker(6)
    if y < -profile_moho[idx_m]:
        cell.setMarker(7)
    if y < -profile_lab[idx_l]:
        cell.setMarker(8)
    
# #%% quick world plot to check

# fig, ax = plt.subplots(figsize=(9,9),dpi=300)
# pg.show(mesh,showMesh=True,markers=True,ax=ax)
# # pg.show(mesh,data=T,showMesh=True, label='Temperature in °C', cMap="inferno",nCols=265,ax=ax)
# # # pg.show(mesh,data=gradient[:,1],showMesh=True, ax=ax,fitView=False,vmin=-0.04,vmax=-0.03)
# # # ax.set_aspect('equal','box')
# # ax.set_xlim([0,20000])
# # ax.set_ylim([-1000,-250])
# ax.set_xlabel('x in m')
# ax.set_ylabel('depth in m')

#%% run calc and get shf

T = pg.solver.solveFiniteElements(mesh,
                                a={0: 1.5, 5: 1.5, 6: 2.7, 7: 3.5, 8: 4.0},#0: 1.5, 1: 1.5
                                f=force,
                                bc={'Dirichlet': {8: 1315, 5: 0}},verbose=True)# {'Node': Tnode} }, 
# pg.show(mesh,data=T,showMesh=True, label='Temperature in °C', cMap="inferno",nCols=256)

T_list = [T[i] for i in range(len(T))]

# similarly to the 3D case, I'll probably want to select other locations for the measurement points
topcoord = []
for i in range(len(dist_prof_t)):
    topcoord.append([dist_prof_t[i],profile_topo[i]-1,0])

gradientTop = pg.solver.grad(mesh, T, topcoord)

geology_list = []
shf = np.zeros(len(gradientTop))
for i in range(len(gradientTop)):
    point = pg.RVector3(topcoord[i][0],topcoord[i][1],topcoord[i][2])
    cell = mesh.findCell(point)
    geology = cell.marker()
    geology_list.append(geology)
    if geology == 6:   
        shf[i] = np.linalg.norm(gradientTop[i,:])*2.7
    elif geology == 5:
        shf[i] = np.linalg.norm(gradientTop[i,:])*1.5

shf_2D = shf*1000        

#%%

# fig,ax = plt.subplots(figsize=(9,9),dpi=300)
# pg.show(mesh,showMesh=True,markers=True,ax=ax,clipBoundaryMarkers=True)
# ax.set_xlim([40000,60000])
# ax.set_ylim([-2000,2000])
# ax.set_xlabel('x in m')
# ax.set_ylabel('depth in m')

nodes_pos = [mesh.nodes()[i].pos() for i in range(len(mesh.nodes()))]

depth_bins = np.arange(np.min(-profile_lab),np.max(profile_topo)+100,100)
depth_mean = np.zeros(len(depth_bins))
T_arr = np.array(T)
n_pos = np.array(nodes_pos)
for i,d in enumerate(depth_bins[:-1]):
    sel = (n_pos[:,1] >= depth_bins[i]) & (n_pos[:,1] < depth_bins[i+1])
    depth_mean[i] = T_arr[sel].mean()

pred_T = np.interp(n_pos[:,1],depth_bins,depth_mean)
T_anom = T-pred_T

for i, node in enumerate(mesh.nodes()):
    x,y,z = node.pos()
    
    idx_t_0,idx_t_1 = np.argsort(np.abs(dist_prof_t-x))[0],np.argsort(np.abs(dist_prof_t-x))[1]
    
    m_t = (profile_topo[idx_t_0]-profile_topo[idx_t_1])/(dist_prof_t[idx_t_0]-dist_prof_t[idx_t_1])
    y_t = m_t * (x-dist_prof_t[idx_t_1]) + profile_topo[idx_t_1] 
    
    if y >= y_t:
        T_anom[i] = 0


# gradient_full = pg.solver.grad(mesh,T)
# gradient_nodes = pg.solver.grad(mesh, T, nodes_pos)
# gradient_nodes_norm = np.linalg.norm(gradient_nodes,axis=1)

fig,ax = plt.subplots(figsize=(9,4),dpi=300)
pg.show(mesh,data=T_anom,cMin=-25,cMax=25,showMesh=False,showBoundary=True,label='$\Delta$T in °C',orientation='vertical', cMap="RdYlBu_r",nCols=50,ax=ax,shading='flat')#,cMin=0,cMax=10
# pg.show(mesh,data=gradient_nodes_norm,showMesh=True,showBoundary=True,label='idk', cMap="plasma",nCols=50,ax=ax,shading='gouraud')#,cMin=0,cMax=10
# pg.show(mesh,linewidth=.5,data=np.linalg.norm(gradient_full,axis=1),cMin=0.01,cMax=0.04,showMesh=True,showBoundary=True,label='Tempterature gradient in K/m', cMap="plasma",nCols=50,ax=ax)#,cMin=0,cMax=10
# pg.show(mesh,data=gradient_full[:,1],showMesh=True,showBoundary=True,label='idk', cMap="plasma",nCols=50,ax=ax,shading='gouraud')#,cMin=0,cMax=10
# pg.show(mesh,data=T,showMesh=True,showBoundary=True,label='Temperature in °C',cMin=0,cMax=10, cMap="plasma",nCols=50,ax=ax)#,cMin=0,cMax=10
# ax.set_xlim([40000,65000])
ax.set_ylim([-12500,1000])
ax.set_xlabel('Profile length in m')
ax.set_ylabel('Depth in m')
ax.set_aspect('auto')

# from pygimli.viewer.mpl import drawStreams
# fig,ax = plt.subplots(figsize=(9,9),dpi=300)
# pg.show(mesh,data=T,showMesh=True,showBoundary=True,label='Temperature in °C',cMin=0,cMax=50, cMap="plasma",nCols=50,ax=ax,shading='gouraud')#,cMin=0,cMax=10
# drawStreams(ax,mesh,gradient_full*100000000,color='green',quiver=True)
# ax.set_xlim([40000,50000])
# ax.set_ylim([-500,500])
# ax.set_xlabel('x in m')
# ax.set_ylabel('depth in m')
# ax.set_aspect('auto')

#%% Plot data and compare to 1D result

# plt.figure(dpi=300)
# plt.plot(dist_prof_t,gradientTop[:,1])   
# plt.plot(dist_prof_t,shf_2D)      
# plt.plot(dist_prof_A,profile_A)   

profile_moho_t = profile_moho_for_topo - profile_sed_for_topo + profile_topo
profile_lab_t = profile_lab_for_topo - profile_moho_for_topo

A = profile_A_for_topo
# A = 0.000001

shf_1D = ( 1315 - 0 + (A*profile_moho_t*profile_lab_t)/k2 + (1/2*A*profile_moho_t**2)/k1 )/( profile_lab_t/k2 + profile_moho_t/k1 + profile_sed_for_topo/k0 )*1000

plt.figure(figsize=(8,3),dpi=300)
plt.plot(dist_prof_t,shf_2D,label="2D")
# plt.plot(dist_prof_t,shf_2D_interp,label="2D interp")
plt.plot(dist_prof_t,shf_1D,label="1D")  
plt.xlabel('x in m')
plt.ylabel('GHF in mW/m$^2$')
plt.legend(loc='upper right')

# plt.figure(figsize=(8,3),dpi=300)
# plt.plot(dist_prof_t,shf_2D,label="2D")
# plt.plot(collect_cells_unique[:,1],collect_cells_unique[:,3]*-1.5*1000,label="2D interp")
# plt.plot(dist_prof_t,shf_1D,label="1D")  
# plt.xlabel('x in m')
# plt.ylabel('GHF in mW/m$^2$')
# plt.legend()

# plt.figure(dpi=300) 
# plt.plot(dist_prof_t,shf_2D-shf_1D)  

# plt.figure(dpi=300) 
# plt.plot(dist_prof_t,profile_topo,label='topo')  
# plt.plot(dist_prof_t,profile_topo-profile_sed_for_topo,label='sed depth')
# plt.xlabel('x in m')
# plt.ylabel('height in m')
# plt.legend() 


# interface (change acordingly):

# fig, (ax1, ax2) = plt.subplots(2, 1,sharex=True,dpi=300)

# ax1.plot(dist_prof_t,shf_2D,label="2D")
# ax1.plot(dist_prof_t,shf_1D,label="1D")  
# ax1.set_ylabel('GHF in mW/m$^2$')
# ax1.legend(loc='lower right')

# # ax2.plot(dist_prof_t,profile_topo/1000,color='slategrey')
# # ax2.plot(dist_prof_t,(profile_topo-profile_sed_for_topo)/1000,color='gold')
# # ax2.plot(dist_prof_m,-profile_moho/1000,color='darkorange')
# # ax2.plot(dist_prof_l,-profile_lab/1000,color='orangered')
# ax2.plot(dist_prof_t,profile_A_for_topo,color='darkorchid')
# ax2.set_xlabel('x in m')
# # ax2.set_ylabel('height in km')
# # ax2.set_ylabel('thickness in km')
# # ax2.set_ylabel('depth in km')
# ax2.set_ylabel('HP in W/m$^3$')


# all interfaces:

# fig, (ax1, ax2, ax3) = plt.subplots(3, 1,sharex=True,dpi=300)

# ax1.plot(dist_prof_t,(profile_topo-profile_sed_for_topo)/1000,color='yellow',label='Sediment depth')
# ax1.plot(dist_prof_t,profile_topo/1000,label='Topography')
# # ax1.set_xlabel('x in m')
# ax1.set_ylabel('height in km')
# # ax1.colorbar(label='mW/m$^2$')

# ax2.plot(dist_prof_m,-profile_moho/1000,color='orange',label='Moho')
# # ax2.set_xlabel('x in m')
# ax2.set_ylabel('depth in km')

# ax3.plot(dist_prof_l,-profile_lab/1000,color='red',label='LAB')
# ax3.set_xlabel('x in m')
# ax3.set_ylabel('depth in km')

# fig.legend()

#%% save
print(notes)

# still need to add seperate msp
save_dict = {'shf_2D': shf_2D,
             'topcoord': topcoord,
             'T': T_list,
             'shf_1D': shf_1D,
             'profile_topo': profile_topo,
             'dist_prof_t': dist_prof_t,
             'x_int_t': x_int_t,
             'y_int_t': y_int_t,
             'profile_sed': profile_sed_for_topo,
             'profile_A': profile_A,
             'dist_prof_A': dist_prof_A,
             'x_int_A': x_int_A,
             'y_int_A': y_int_A,
             'profile_moho': profile_moho,
             'dist_prof_m': dist_prof_m,
             'x_int_m': x_int_m,
             'y_int_m': y_int_m,
             'profile_lab': profile_lab,
             'dist_prod_l': dist_prof_l,
             'x_int_l': x_int_l,
             'y_int_l': y_int_l,}

now = datetime.datetime.now()
mesh.exportVTK('outputVTK/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes))
world.exportPLC('outputVTK/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes))
with open('output/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes)+'.pkl','wb+') as file:
    pickle.dump(save_dict,file)

