"""
Georg HÃ¼ttner - September 2024

This script performs three dimensional thermal modelling from provided lithospheric
interfaces. Based on this a surface geothermal heat flow can be calculated.

Look on my works, ye Mighty, and despair!

There are a bunch of really stupid things in this script and the entire method needs improvement.
Some points to keep in mind: 
 - world size and extent of data should align with data resolution, as
   otherwise the gaps at the edges could lead "misapplied" markers
 - there are a lot of transposed 2d arrays in here, really make sure that your data is oriented correctly
 - mesh creation should probably be done with an external program, as the pygimli internal stuff
   takes unbearably long for large regions
 - marker application could probably be overhauled to allow for faster result (who would have thunk that
   a for loop over a couple of million cells takes long)
 - the solver itself is incredibly ram hungry and for larger models (cell count > 3 mil) a normal pc is
   gonna die
"""

import pygimli as pg
import pygimli.meshtools as mt
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib import ticker
import cartopy.crs as ccrs
from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter)
#import polygon_setup as ps
from scipy.interpolate import griddata,LinearNDInterpolator, RegularGridInterpolator
from scipy.spatial import KDTree
import datetime
import pickle
from scipy.stats import norm
import gepy
from pyproj import CRS, Transformer
# from skimage.util.shape import view_as_windows

data_calc,data_resolution = gepy.load_data()

#%% research area and params

# like profile coord, but take x0,x1 and y0,y1 of area, area width/length should be divisible by 10000
area_coord_s = [1300000,1400000,-275000,-375000] # Lake Vostok "south-western" edge
# area_coord_s = [1355000,1435000,-855000,-985000] # Dome C
# area_coord_s = [1800000,1900000,-650000,-750000] # ASB

k0 = 1.5
k1 = 2.7
k2 = 3.5

km = 1000

#%% quick plot for yall

vmin = np.nanmin(data_calc[1][2])
vmax = np.nanmax(6000)

plt.figure(dpi=300)
ax = plt.axes(projection=ccrs.AzimuthalEquidistant(central_latitude=-90))
norm = Normalize(vmin=vmin, vmax=vmax)
sc = ax.pcolormesh(data_calc[1][0],data_calc[1][1],data_calc[1][2],cmap='YlGnBu_r',vmin=vmin,vmax=vmax)
# sc = ax.pcolormesh(xi_smos,yi_smos,grid_A_ws,cmap='afmhot')
rect = patches.Rectangle((area_coord_s[0], area_coord_s[3]), area_coord_s[1]-area_coord_s[0], area_coord_s[2]-area_coord_s[3], linewidth=1, edgecolor='r', facecolor='none')
ax.add_patch(rect)
ax.set_extent([250000,2750000,-2200000,-38000],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, rotate_labels = False)
gl.top_labels = False
gl.right_labels = False
gl.ylabel_style = {'rotation': 45, 'rotation_mode': 'anchor'}
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()
ax.contour(data_calc[1][0],data_calc[1][1],data_calc[1][2],levels=[0,1500,3000,4500],colors='k',linewidths=0.5) # sed
# ax.contour(data_calc[0][0],data_calc[0][1],data_calc[0][2],levels=[-1500,-1000,-500,0,500,1000,1500],colors='k',linewidths=0.5) # topo
# cbar = plt.colorbar(cm.ScalarMappable(norm=norm,cmap='jet'), ax=ax,fraction=0.038, pad=0.04)
# cbar.set_label("Sediment depth in m")
# cbar.formatter.set_powerlimits((-2, 2))
# cbar.formatter.set_useMathText(True)
# cbar.locator = ticker.MaxNLocator(nbins=7)
# cbar.update_ticks()
ax.coastlines(resolution='10m',linewidth=0.5)
plt.show()

#%% load data and cut it to size

xi_topo_s, yi_topo_s, grid_topo_s = gepy.in_area_s(area_coord_s,data_calc[0][0],data_calc[0][1],data_calc[0][2])
xi_sed_s, yi_sed_s, grid_sed_s    = gepy.in_area_s(area_coord_s,data_calc[1][0],data_calc[1][1],data_calc[1][2])
xi_moho_s, yi_moho_s, grid_moho_s = gepy.in_area_s(area_coord_s,data_calc[2][0],data_calc[2][1],data_calc[2][2])
xi_lab_s, yi_lab_s, grid_lab_s    = gepy.in_area_s(area_coord_s,data_calc[3][0],data_calc[3][1],data_calc[3][2])

xii_topo_s,yii_topo_s = np.meshgrid(xi_topo_s,yi_topo_s)
xii_sed_s,yii_sed_s   = np.meshgrid(xi_sed_s,yi_sed_s)
xii_moho_s,yii_moho_s = np.meshgrid(xi_moho_s,yi_moho_s)
xii_lab_s,yii_lab_s   = np.meshgrid(xi_lab_s,yi_lab_s)

# 

res_topo,res_moho,res_lab = (1000,10000,10000)

x_int_t,y_int_t = np.arange(area_coord_s[0],area_coord_s[1]+res_topo,res_topo),np.arange(area_coord_s[3],area_coord_s[2]+res_topo,res_topo)
# x_int_m,y_int_m = np.arange(area_coord_s[0],area_coord_s[1]+res_moho,res_moho),np.arange(area_coord_s[3],area_coord_s[2]+res_moho,res_moho)
# x_int_l,y_int_l = np.arange(area_coord_s[0],area_coord_s[1]+res_lab,res_lab),  np.arange(area_coord_s[3],area_coord_s[2]+res_lab,res_lab)
x_int_m,y_int_m = np.arange(area_coord_s[0],area_coord_s[1]+res_moho,res_moho),np.arange(area_coord_s[3],area_coord_s[2]+res_moho,res_moho)
x_int_l,y_int_l = np.arange(area_coord_s[0],area_coord_s[1]+res_lab,res_lab),np.arange(area_coord_s[3],area_coord_s[2]+res_lab,res_lab)

xi_int_t,yi_int_t = np.meshgrid(x_int_t,y_int_t)
xi_int_m,yi_int_m = np.meshgrid(x_int_m,y_int_m)
xi_int_l,yi_int_l = np.meshgrid(x_int_l,y_int_l)

interp_topo = RegularGridInterpolator((np.unique(xi_topo_s),np.unique(yi_topo_s)),grid_topo_s)
grid_topo = interp_topo((xi_int_t,yi_int_t))

interp_sed = RegularGridInterpolator((np.unique(xi_sed_s),np.unique(yi_sed_s)),grid_sed_s)
grid_sed = interp_sed((xi_int_t,yi_int_t),method='linear')

for i in range(np.shape(grid_sed)[0]):
    for j in range(np.shape(grid_sed)[1]):
        grid_sed[i,j] = round(grid_sed[i,j]/35)*35

interp_moho = RegularGridInterpolator((np.unique(xi_moho_s),np.unique(yi_moho_s)),grid_moho_s)
grid_moho = interp_moho((xi_int_m,yi_int_m))
grid_moho_for_topo = interp_moho((xi_int_t,yi_int_t))

interp_lab = RegularGridInterpolator((np.unique(xi_lab_s),np.unique(yi_lab_s)),grid_lab_s)
grid_lab = interp_lab((xi_int_l,yi_int_l))
grid_lab_for_topo = interp_lab((xi_int_t,yi_int_t))

# heat production

data_HP = np.load("data/SMOS_hp.npz")
grid_A_ws = data_HP['A_ws']
grid_A_wo = data_HP['A_wo']

xi_smos,yi_smos = data_HP['xi_smos'],data_HP['yi_smos']

A = grid_A_ws # heat production with (ws) or without (wo) sediments? 

xi_smos_s, yi_smos_s, A_s = gepy.in_area_s(area_coord_s,xi_smos,yi_smos,A)
xii_smos_s,yii_smos_s   = np.meshgrid(xi_smos_s,yi_smos_s)
    
res_HP = 12500
    
x_int_HP,y_int_HP = np.arange(area_coord_s[0],area_coord_s[1]+res_HP,res_HP),  np.arange(area_coord_s[3],area_coord_s[2]+res_HP,res_HP)
xi_int_HP,yi_int_HP = np.meshgrid(x_int_HP,y_int_HP)

grid_A  = griddata((xii_smos_s.flatten(), yii_smos_s.flatten()), np.transpose(A_s).flatten(), (xi_int_HP,yi_int_HP),fill_value=0.000001)
# interp_A = RegularGridInterpolator((x_int_HP,y_int_HP),np.transpose(grid_A))
interp_A = RegularGridInterpolator((x_int_HP,y_int_HP),np.transpose(grid_A))
grid_A_for_topo = interp_A((xi_int_t,yi_int_t))
# grid_A_for_topo = griddata((xii_smos_s.flatten(), yii_smos_s.flatten()), A_s.flatten(), (xi_int_t,yi_int_t),method='linear')
    

# # heat production part 2, electric bogaloo

# data_HP_Lös = np.loadtxt("data/Predicted_HP_Ant_Lös.txt",skiprows=1)

# spae_crs = CRS(proj="aeqd", lat_0=-90, lon_0=0, datum="WGS84")

# ll_to_spae = Transformer.from_crs(CRS("EPSG:4326"),spae_crs, always_xy=True)

# x_HP_Lös, y_HP_Lös = ll_to_spae.transform(data_HP_Lös[:,0],data_HP_Lös[:,1])

# plt.figure(dpi=300)
# plt.scatter(x_HP_Lös,y_HP_Lös,c=data_HP_Lös[:,2]/np.e,s=5)
# plt.colorbar()
# plt.show()


del A_s,xi_smos,yi_smos,xi_smos_s,yi_smos_s,xii_smos_s,yii_smos_s    
del grid_topo_s,xi_topo_s,yi_topo_s,xii_topo_s,yii_topo_s
del grid_sed_s,xi_sed_s,yi_sed_s,xii_sed_s,yii_sed_s
del grid_moho_s,xi_moho_s,yi_moho_s,xii_moho_s,yii_moho_s
del grid_lab_s,xi_lab_s,yi_lab_s,xii_lab_s,yii_lab_s

#%% ja/nein/vielleicht?

notes = '3d_LV_A_1.5_norm'
print(notes)

grid_topo = np.full(np.shape(grid_topo),np.mean(grid_topo))
# grid_topo = grid_topo/1000
grid_sed = np.full(np.shape(grid_sed),np.mean(grid_sed))
grid_moho = np.full(np.shape(grid_moho),np.mean(grid_moho))
grid_lab = np.full(np.shape(grid_lab),np.mean(grid_lab))
# grid_A_for_topo = np.full(np.shape(grid_A_for_topo),np.mean(grid_A_for_topo))

grid_A_norm = grid_A_for_topo/np.max(grid_A_for_topo)
grid_A_for_topo = grid_A_for_topo*grid_A_norm*1.5


#%% build the world

# create 3d mesh from the input data

# the +/- 10 are to create a slight buffer around the input data and world edge, there are issues
# in the mesh creation when edges of boundaries touch other boundaries
border = 10
start = [(area_coord_s[0]-border)/km,(area_coord_s[2]+border)/km,4] 
end = [(area_coord_s[1]+border)/km,(area_coord_s[3]-border)/km,-220]

world = mt.createWorld(start=start,end=end, worldMarker=False)

# create layer polygons and shift them to the correct height
mesh_moho = mt.createMesh2D(x_int_m/km,y_int_m/km)
mesh_lab  = mt.createMesh2D(x_int_l/km,y_int_l/km)

surface_moho = mt.createSurface(mesh_moho)
surface_lab  = mt.createSurface(mesh_lab)

surface_moho = gepy.fix_surface_height(surface_moho, x_int_m/km, y_int_m/km, -grid_moho/km)
surface_lab  = gepy.fix_surface_height(surface_lab, x_int_l/km, y_int_l/km, -grid_lab/km)

for boundary in surface_moho.boundaries():
    boundary.setMarker(6)
for boundary in surface_lab.boundaries():
    boundary.setMarker(7)

surface_topo,surface_sed = gepy.create_sed_interface(x_int_t, y_int_t, grid_topo, grid_sed)

geometry = world + surface_topo + surface_sed + surface_moho + surface_lab

mesh = mt.createMesh(geometry,quality=34,area=2.5)#,area=2.5

# print(geometry)
# pg.show(geometry,showMesh=True,alpha=0.7)

# print(mesh)
# pg.show(mesh,showMesh=True,alpha=0.7)

#%% apply correct markers to cells and nodes, and create temp and force vectors

grid_xy_ts = np.array(list(zip(xi_int_t.flatten()/km,yi_int_t.flatten()/km)))
tree_ts = KDTree(grid_xy_ts)
grid_xy_m = np.array(list(zip(xi_int_m.flatten()/km,yi_int_m.flatten()/km)))
tree_m = KDTree(grid_xy_m)
grid_xy_l = np.array(list(zip(xi_int_l.flatten()/km,yi_int_l.flatten()/km)))
tree_l = KDTree(grid_xy_l)

# create lists for temperature and heat production that are to be applied to the nodes

force = np.zeros(mesh.nodeCount())
# Tnode = []
for i, node in enumerate(mesh.nodes()):
    x,y,z = node.pos()
    z_topo = gepy.compare_to_plane(x,y,grid_topo,grid_xy_ts,tree_ts)
    z_sed = gepy.compare_to_plane(x,y,(grid_topo-grid_sed),grid_xy_ts,tree_ts)
    z_moho = gepy.compare_to_plane(x,y,-grid_moho,grid_xy_m,tree_m)
    z_lab = gepy.compare_to_plane(x,y,-grid_lab,grid_xy_l,tree_l)
    
    idx_A,idy_A = np.argmin(np.abs(x_int_t/km-x)),np.argmin(np.abs(y_int_t/km-y))
    
    if z >= z_moho and z <= z_sed:
        force[i] = grid_A_for_topo[idy_A,idx_A]
        # force[i] = 0.000001
    else:
        force[i] = 0
    
# apply correct marker to the cells
for i,cell in enumerate(mesh.cells()):
    center = cell.center()
    x, y, z = center.x(), center.y(), center.z()
    
    z_topo = gepy.compare_to_plane(x,y,grid_topo,grid_xy_ts,tree_ts)
    z_sed = gepy.compare_to_plane(x,y,(grid_topo-grid_sed),grid_xy_ts,tree_ts)
    z_moho = gepy.compare_to_plane(x,y,-grid_moho,grid_xy_m,tree_m)
    z_lab = gepy.compare_to_plane(x,y,-grid_lab,grid_xy_l,tree_l)
    
    if z < z_topo:
        cell.setMarker(4)
    if z < z_sed:
        cell.setMarker(5)
    if z < z_moho:
        cell.setMarker(6)
    if z < z_lab:
        cell.setMarker(7)

# pg.show(mesh,showMesh=True)

#%% run calc and get shf

T = pg.solver.solveFiniteElements(mesh,
                                a={1: 1.0*km, 4: 1.5*km, 5: 2.7*km, 6: 3.5*km, 7: 4.0*km},
                                f=force*km*km*km,
                                bc={'Dirichlet': {7: 1315, 4: 0}},verbose=True)#{'Node': Tnode}

T_list = [T[i] for i in range(len(T))]

# pg.show(mesh,data=T,showMesh=True, label='Temperature in °C', cMap="inferno")

# gradient = pg.solver.grad(mesh, T)

# pg.show(mesh,data=gradient[:,2],filter={'clip':{'origin':(1050, 0, 0)},})

topcoord = []
for i in range(len(x_int_t)):
    for j in range(len(y_int_t)):
        topcoord.append([x_int_t[i]/km,y_int_t[j]/km,(grid_topo[j,i]-5)/km])

gradientTop = pg.solver.grad(mesh, T, topcoord)

geology_list = []

shf = np.zeros(len(gradientTop))
for i in range(len(gradientTop)):
    point = pg.RVector3(topcoord[i][0],topcoord[i][1],topcoord[i][2])
    cell = mesh.findCell(point)
    geology = cell.marker()
    geology_list.append(geology)
    if geology == 5:   
        shf[i] = np.linalg.norm(gradientTop[i,:])*2.7
    elif geology == 4:
        shf[i] = np.linalg.norm(gradientTop[i,:])*1.5

shf_3D = np.transpose(np.reshape(shf,(len(x_int_t),len(y_int_t))))
gradientTop_rs = np.reshape(gradientTop[:,2],(len(x_int_t),len(y_int_t)))

#%% Plot data and compare to 1D result

plt.figure(dpi=300)
plt.pcolormesh(x_int_t,y_int_t,shf_3D,vmin=40,vmax=80)
# plt.pcolormesh(x_int_t,y_int_t,shf_3D_norm,vmin=30,vmax=100)
# plt.pcolormesh(gradientTop_rs)
# plt.pcolormesh(x_int_t,y_int_t,geology_rs)
# plt.pcolormesh(x_int_t,y_int_t,grid_topo)
# plt.pcolormesh(x_int_t,y_int_t,grid_sed,vmin=0,vmax=70)
# plt.pcolormesh(x_int_HP,y_int_HP,np.transpose(grid_A))
# plt.pcolormesh(grid_A_for_topo)
plt.xlabel('x in m')
plt.ylabel('y in m')
plt.colorbar(label='mW/m$^2$')

grid_moho_t = grid_moho_for_topo - grid_sed + grid_topo
grid_lab_t = grid_lab_for_topo - grid_moho_for_topo

A_1D = grid_A_for_topo
# A_1D = 0.000001

shf_1D = (1315-0+(A_1D*grid_moho_t*grid_lab_t)/k2+(1/2*A_1D*grid_moho_t**2)/k1)/(grid_lab_t/k2+grid_moho_t/k1+grid_sed/k0)*1000

plt.figure(dpi=300)
plt.pcolormesh(x_int_t,y_int_t,shf_1D,vmin=40,vmax=80)#
plt.xlabel('x in m')
plt.ylabel('y in m')
plt.colorbar(label='mW/m$^2$')

plt.figure(dpi=300)
plt.pcolormesh(x_int_t,y_int_t,shf_3D-shf_1D,cmap='bwr',vmin=-15,vmax=15)
plt.xlabel('x in m')
plt.ylabel('y in m')
plt.colorbar(label='mW/m$^2$')


#%% save 

print(notes)
save_dict = {'shf_3D': shf_3D,
             'topcoord': topcoord,
             'T': T_list,
             'shf_1D': shf_1D,
             'grid_topo': grid_topo,
             'x_int_t': x_int_t,
             'y_int_t': y_int_t,
             'grid_sed': grid_sed,
             'grid_A': grid_A,
             'grid_A_for_topo': grid_A_for_topo,
             'x_int_HP': x_int_HP,
             'y_int_HP': y_int_HP,
             'grid_moho': grid_moho,
             'grid_moho_for_topo': grid_moho_for_topo,
             'x_int_m': x_int_m,
             'y_int_m': y_int_m,
             'grid_lab': grid_lab,
             'grid_lab_for_topo': grid_lab_for_topo,
             'x_int_l': x_int_l,
             'y_int_l': y_int_l
             }

now = datetime.datetime.now()
mesh.exportVTK('outputVTK/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes))
geometry.exportPLC('outputVTK/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes))
# geometry.exportSTL('outputVTK/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes))
with open('output/'+str(now.strftime('%Y-%m-%d_%H-%M-%S'))+'_'+str(notes)+'.pkl','wb+') as file:
    pickle.dump(save_dict,file)

#%% code dump

# area_coord_LV = [1300000,1400000,-275000,-375000] # Lake Vostok "south-western" edge
# area_coord_DC = [1355000,1435000,-855000,-985000] # Dome C
# area_coord_ASB = [1800000,1900000,-650000,-750000] # ASB

# profile_coord=([1300000,1400000],[-275000,-375000])

# plt.figure(dpi=300)
# ax = plt.axes(projection=ccrs.AzimuthalEquidistant(central_latitude=-90))
# sc = ax.pcolormesh(data_calc[1][0],data_calc[1][1],data_calc[1][2]/1000,cmap='YlGn_r',vmin=0,vmax=6)
# rect = patches.Rectangle((area_coord_LV[0], area_coord_LV[3]), area_coord_LV[1]-area_coord_LV[0], area_coord_LV[2]-area_coord_LV[3], linewidth=1, edgecolor='firebrick', facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((area_coord_DC[0], area_coord_DC[3]), area_coord_DC[1]-area_coord_DC[0], area_coord_DC[2]-area_coord_DC[3], linewidth=1, edgecolor='gold', facecolor='none')
# ax.add_patch(rect)
# rect = patches.Rectangle((area_coord_ASB[0], area_coord_ASB[3]), area_coord_ASB[1]-area_coord_ASB[0], area_coord_ASB[2]-area_coord_ASB[3], linewidth=1, edgecolor='darkorange', facecolor='none')
# ax.add_patch(rect)
# pr = ax.plot(profile_coord[0],profile_coord[1],'m',linewidth=1)
# ax.set_extent([250000,2750000,-2200000,-38000],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
# gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, rotate_labels = False)
# gl.top_labels = False
# gl.right_labels = False
# gl.ylabel_style = {'rotation': 45, 'rotation_mode': 'anchor'}
# gl.xformatter = LongitudeFormatter()
# gl.yformatter = LatitudeFormatter()
# # ax.contour(data_calc[1][0],data_calc[1][1],data_calc[1][2],levels=[0,1500,3000,4500],colors='k',linewidths=0.5) # sed
# # ax.contour(data_calc[0][0],data_calc[0][1],data_calc[0][2],levels=[-1500,-1000,-500,0,500,1000,1500],colors='k',linewidths=0.5) # topo
# # cbar = plt.colorbar(cm.ScalarMappable(norm=norm,cmap='YlGn_r'), ax=ax,fraction=0.038, pad=0.04)
# cbar = plt.colorbar(sc,ax=ax,fraction=0.038, pad=0.04)
# cbar.set_label("Sediment thickness in km")
# cbar.formatter.set_powerlimits((-2, 2))
# cbar.formatter.set_useMathText(True)
# cbar.locator = ticker.MaxNLocator(nbins=7)
# cbar.update_ticks()
# ax.coastlines(resolution='10m',linewidth=0.5)
# plt.show()

# topoplot

from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import cmcrameri.cm as cmc

vmin = np.nanmin(-3000)
vmax = np.nanmax(3000)

cmap_terrain = plt.get_cmap('terrain')
cmap_below_sea = plt.get_cmap('Blues_r')  # Reverse 'Blues' so that lighter blues are closer to zero

# Create a combined colormap
colors_below = cmap_below_sea(np.linspace(0, 1, 128))   # First half for below zero
colors_above = cmap_terrain(np.linspace(0, 1, 128))     # Second half for above zero
colors = np.vstack((colors_below, colors_above))
custom_cmap = LinearSegmentedColormap.from_list('custom_terrain', colors)

# Normalize the data with 0 as the midpoint
norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

plt.figure(dpi=300)
ax = plt.axes(projection=ccrs.AzimuthalEquidistant(central_latitude=-90))
norm = Normalize(vmin=vmin, vmax=vmax)
sc = ax.pcolormesh(data_calc[0][0],data_calc[0][1],data_calc[0][2],cmap=cmc.bukavu,vmin=vmin,vmax=vmax)
# sc = ax.pcolormesh(xi_smos,yi_smos,grid_A_ws,cmap='afmhot')
# rect = patches.Rectangle((area_coord_s[0], area_coord_s[3]), area_coord_s[1]-area_coord_s[0], area_coord_s[2]-area_coord_s[3], linewidth=1, edgecolor='r', facecolor='none')
[250000,2750000,-2200000,-38000]
rect = patches.Rectangle((250000, -2200000), 2750000-250000, -38000--2200000, linewidth=1, edgecolor='r', facecolor='none')
ax.add_patch(rect)
ax.set_extent([-3000000,3000000,-2600000,2600000],crs=ccrs.AzimuthalEquidistant(central_latitude=-90))
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, alpha=0.5, rotate_labels = False)
gl.top_labels = False
gl.right_labels = False
gl.ylabel_style = {'rotation': 45, 'rotation_mode': 'anchor'}
gl.xformatter = LongitudeFormatter()
gl.yformatter = LatitudeFormatter()
# ax.contour(data_calc[1][0],data_calc[1][1],data_calc[1][2],levels=[0,1500,3000,4500],colors='k',linewidths=0.5) # sed
# ax.contour(data_calc[0][0],data_calc[0][1],data_calc[0][2],levels=[-1500,-1000,-500,0,500,1000,1500],colors='k',linewidths=0.5) # topo
cbar = plt.colorbar(cm.ScalarMappable(norm=norm,cmap=cmc.bukavu), ax=ax,fraction=0.038, pad=0.04)
cbar.set_label("Topography in m")
# cbar.formatter.set_powerlimits((-2, 2))
cbar.formatter.set_useMathText(True)
cbar.locator = ticker.MaxNLocator(nbins=7)
cbar.update_ticks()
ax.coastlines(resolution='10m',linewidth=0.5)
plt.show()


# I'll let this coat stand for now, as I'll probably choose other measurement point resolutions for the
# surface GHF than the topo locations

# msp_abst = 1000
# msp_x = np.arange(area_coord_s[0]+msp_abst,area_coord_s[1],msp_abst)
# msp_y = np.arange(area_coord_s[3]+msp_abst,area_coord_s[2],msp_abst)
# mspi_x,mspi_y = np.meshgrid(msp_x,msp_y)

# msp_grid = griddata((xii_topo_s.flatten(),yii_topo_s.flatten()),grid_topo_s.flatten(),(mspi_x,mspi_y),method='linear')
# msp_grid_fill = griddata((xii_topo_s.flatten(),yii_topo_s.flatten()),grid_topo_s.flatten(),(mspi_x,mspi_y),method='nearest')
# msp_grid[np.isnan(msp_grid)] = msp_grid_fill[np.isnan(msp_grid)]
# msp_grid = (msp_grid)/1000


# # interpolation:

# x_int_10m = np.arange(x_int_t[0],x_int_t[-1]+100,100)
# y_int_10m = np.arange(y_int_t[0],y_int_t[-1]+100,100)

# xi_int_10m,yi_int_10m = np.meshgrid(x_int_10m,y_int_10m)

# grid_topo_10m = griddata((xi_int_t.flatten(),yi_int_t.flatten()),grid_topo.flatten(),(xi_int_10m,yi_int_10m))

# gradient = pg.solver.grad(mesh, T)

# # collect_cells = np.arange(0,profile_length+5,5)
# # collect_cells_d = np.interp(collect_cells,dist_prof_t,profile_topo)

# collect_cells_all = []
# for i in range(len(x_int_10m)):
#     for j in range(len(y_int_10m)):
#         cell = mesh.findCell(pg.RVector3(x_int_10m[i]/km,y_int_10m[j]/km,(grid_topo_10m[i,j]-5)/km))
#         x = cell.center().x()
#         y = cell.center().y()
#         ID = cell.id()
#         marker = cell.marker()
#         vertgrad = np.linalg.norm(gradient[ID,:])
    
#         celldata = [ID,x,y,marker,vertgrad]
#         collect_cells_all.append(celldata)
    
# collect_cells_all = np.asarray(collect_cells_all)    
# unique_keys, indices = np.unique(collect_cells_all[:,0], return_index=True)
# collect_cells_unique = collect_cells_all[np.sort(indices)]


# # auf 100 m

# shf_interp = griddata((collect_cells_unique[:,1],collect_cells_unique[:,2]),collect_cells_unique[:,4],(xi_int_10m/km,yi_int_10m/km))
# shf_interp = np.transpose(shf_interp)

# geology_list = []
# shf_3D_interp = np.zeros(np.shape(shf_interp))
# for i in range(len(x_int_10m)):
#     for j in range(len(y_int_10m)):
#         point = pg.RVector3(x_int_10m[i]/km,y_int_10m[j]/km,(grid_topo_10m[i,j]-5)/km)
#         cell = mesh.findCell(point)
#         geology = cell.marker()
#         geology_list.append(geology)
#         if geology == 5:
#             shf_3D_interp[i,j] = shf_interp[i,j]*2.7
#         elif geology == 4:
#             shf_3D_interp[i,j] = shf_interp[i,j]*1.5

# plt.figure(dpi=300)
# plt.pcolormesh(x_int_10m,y_int_10m,shf_3D_interp,vmin=30,vmax=100)#,vmin=30,vmax=100
# # plt.pcolormesh(x_int_10m,y_int_10m,shf_interp)#,vmin=30,vmax=100
# # plt.pcolormesh(x_int_t,y_int_t,shf_3D_interp_t,vmin=30,vmax=100)#,vmin=30,vmax=100
# # plt.pcolormesh(x_int_10m,y_int_10m,shf_interp,vmin=-150,vmax=0)
# # plt.pcolormesh(x_int_10m,y_int_10m,grid_topo_10m)
# plt.xlabel('x in m')
# plt.ylabel('y in m')
# plt.colorbar(label='mW/m$^2$')


###

# x_int_t_tc = x_int_t[:-1]+500
# y_int_t_tc = y_int_t[:-1]+500
# xi_int_t_tc, yi_int_t_tc = np.meshgrid(x_int_t_tc,y_int_t_tc)
# grid_topo_tc = interp_topo((xi_int_t_tc,yi_int_t_tc))


# topcoord_tc = list(zip(xi_int_t_tc.flatten()/km,yi_int_t_tc.flatten()/km,(grid_topo_tc.flatten()-50)/km))
# gradientTop_tc = pg.solver.grad(mesh, T, topcoord_tc)
# geology_list_tc = []
# shf_tc = np.zeros(len(gradientTop_tc))
# for i in range(len(gradientTop_tc)):
#     point = pg.RVector3(topcoord_tc[i][0],topcoord_tc[i][1],topcoord_tc[i][2])
#     cell = mesh.findCell(point)
#     geology = cell.marker()
#     geology_list_tc.append(geology)
#     if geology == 5:   
#         shf_tc[i] = np.linalg.norm(gradientTop_tc[i,:])*2.7
#         # shf[i] = -gradientTop[i,2]*2.7
#     elif geology == 4:
#         shf_tc[i] = np.linalg.norm(gradientTop_tc[i,:])*1.5
#         # shf[i] = -gradientTop[i,2]*1.5
# shf_3D_tc = np.reshape(shf_tc,(len(x_int_t_tc),len(y_int_t_tc)))
# gradientTop_tc_rs = np.reshape(gradientTop_tc[:,2],(len(x_int_t_tc),len(y_int_t_tc)))

# plt.figure(dpi=300)
# # plt.pcolormesh(x_int_t_tc,y_int_t_tc,shf_3D_tc,vmin=40,vmax=80)
# plt.pcolormesh(x_int_t,y_int_t,shf_3D,vmin=40,vmax=80)
# # plt.pcolormesh(x_int_t,y_int_t,shf_3D_norm,vmin=30,vmax=100)
# # plt.pcolormesh(gradientTop_rs)
# # plt.pcolormesh(x_int_t,y_int_t,geology_rs)
# # plt.pcolormesh(x_int_t,y_int_t,grid_topo)
# # plt.pcolormesh(x_int_t,y_int_t,grid_sed,vmin=0,vmax=70)
# # plt.pcolormesh(grid_A_for_topo)
# plt.xlabel('x in m')
# plt.ylabel('y in m')
# plt.colorbar(label='mW/m$^2$')

# grid_moho_t = grid_moho_for_topo - grid_sed + grid_topo
# grid_lab_t = grid_lab_for_topo - grid_moho_for_topo

# # A_1D = grid_A_for_topo
# A_1D = 0.000001

# shf_1D = (1315-0+(A_1D*grid_moho_t*grid_lab_t)/k2+(1/2*A_1D*grid_moho_t**2)/k1)/(grid_lab_t/k2+grid_moho_t/k1+grid_sed/k0)*1000

# plt.figure(dpi=300)
# plt.pcolormesh(x_int_t_tc,y_int_t_tc,shf_1D,vmin=30,vmax=100)
# plt.xlabel('x in m')
# plt.ylabel('y in m')
# plt.colorbar(label='mW/m$^2$')

# plt.figure(dpi=300)
# plt.pcolormesh(x_int_t_tc,y_int_t_tc,shf_3D-shf_1D,cmap='bwr',vmin=-15,vmax=15)
# plt.xlabel('x in m')
# plt.ylabel('y in m')
# plt.colorbar(label='mW/m$^2$')

# check in between topo points

# xi_int_t_mid,yi_int_t_mid = np.meshgrid(x_int_t[:-1]+500,y_int_t[:-1]+500)

# grid_topo_mid = griddata((xi_int_t.flatten(),yi_int_t.flatten()),grid_topo.flatten(),(xi_int_t_mid,yi_int_t_mid))

# topcoord_mid = []
# for i in range(len(x_int_t)-1):
#     for j in range(len(y_int_t)-1):
#         topcoord_mid.append([x_int_t[i]/km+0.5,y_int_t[j]/km+0.5,(grid_topo_mid[i,j]-5)/km])

# gradientTop_mid = pg.solver.grad(mesh, T, topcoord_mid)

# geology_list_mid = []

# shf_mid = np.zeros(len(gradientTop_mid))
# for i in range(len(gradientTop_mid)):
#     point = pg.RVector3(topcoord_mid[i][0],topcoord_mid[i][1],topcoord_mid[i][2])
#     cell = mesh.findCell(point)
#     geology = cell.marker()
#     geology_list_mid.append(geology)
#     if geology == 5:   
#         shf_mid[i] = np.linalg.norm(gradientTop_mid[i,:])*2.7
#         # shf[i] = -gradientTop[i,2]*2.7
#     elif geology == 4:
#         shf_mid[i] = np.linalg.norm(gradientTop_mid[i,:])*1.5
#         # shf[i] = -gradientTop[i,2]*1.5

# shf_3D_mid = np.reshape(shf_mid,(len(x_int_t)-1,len(y_int_t)-1))
# gradientTop_rs_mid = np.reshape(gradientTop_mid[:,2],(len(x_int_t)-1,len(y_int_t)-1))
# plt.pcolormesh(shf_3D_mid,vmin=30,vmax=100)

