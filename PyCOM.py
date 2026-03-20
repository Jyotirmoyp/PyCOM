#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 17:19:05 2026

@author: jyotirmp
"""


import os
import gzip
import numpy as np


model = 'SMEAN2_high_vel_removed_0_1_100_129_200_km_craton'
nproc = 96

#define empty list
coord = []
t= []
e = []
v = []

# extract data from each node
for i in range(0,nproc): 
    #print(i)
    t_path = os.path.join(model, model, '0', f't.{i}.0.gz')
    e_path = os.path.join(model, model, '0', f'visc.{i}.0.gz')
    v_path = os.path.join(model, model, '0', f'vtk_v.{i}.0.gz')
    coord_path = os.path.join(model, model, f'coord.{i}.gz')
    print(t_path)

    with gzip.open(coord_path, 'rt') as input:
        next(input)     #skip first row
        coord.append(np.loadtxt(input))
        
    with gzip.open(t_path, 'rt') as input: 
        next(input)
        next(input)
        t.append(np.loadtxt(input).reshape(-1,1))
        
    with gzip.open(e_path, 'rt') as input: 
        next(input)
        e.append(np.loadtxt(input).reshape(-1,1))
        
    with gzip.open(v_path, 'rt') as input: 
        v.append(np.loadtxt(input))
        
#### ---- vertical stack coord and find depth

coord_full = np.vstack(coord)

theta = coord_full[:, 0]   # colatitude (rad)
phi   = coord_full[:, 1]   # longitude (rad)
r     = coord_full[:, 2]   # radius (normalized)

lon = np.degrees(phi)
lat = 90 - np.degrees(theta)

# Convert normalized radius to depth (km)
depth_km = (1 - r) * 6371

# Reassemble
coord_full = np.column_stack((lon, lat, depth_km))

depth = np.unique(coord_full[:,2])

# stack temperature, viscosity, velocity
t_full = np.vstack(t)
e_full = np.vstack(e)
v_full = np.vstack(v)

#### --- horizontal stack and extract from each depth

depth_layers = []

model_output = np.hstack((coord_full, t_full, e_full, v_full))
for d in depth: 
    print('Working on depth', d, 'km')
    mask = np.isclose(model_output[:,2], d)
    depth_layer = model_output[mask]
    depth_layers.append(depth_layer)
    
    
#%% plot scalar
import pygmt
import numpy as np

layer_index = 16 
layer = depth_layers[layer_index]
lon = layer[:, 0]
lat = layer[:, 1]
temperature = layer[:,3]
viscosity = np.log10(1e21*layer[:,4])
fig = pygmt.Figure()

#region = '0/359/-89.5/89.5'
proj = 'H12c'

temperature_grid = pygmt.surface(
        x = lon, 
        y = lat, 
        z=temperature, 
        region = 'd', 
        spacing = 1/1,
        lower="d",
        upper="d",
        )

viscosity_grid = pygmt.surface(
        x = lon, 
        y = lat, 
        z=viscosity, 
        region = 'd', 
        spacing = 1/1,
        lower="d",
        upper="d",
        )


fig.grdimage(grid = viscosity_grid, 
             region = 'd', 
             projection = proj, 
             cmap = True)

fig.coast(shorelines = True, 
          region = 'd', 
          projection = proj, 
          #frame = True
          )


fig.colorbar(
    frame='af+lTemperature'
)
fig.show()



#%% plot velocity
import pygmt
import numpy as np

layer_index = 1     #check which layer depth to plot
layer = depth_layers[layer_index]
lon = layer[:, 0]
lat = layer[:, 1]
temperature = layer[:,3]
viscosity = np.log10(1e21*layer[:,4])
vel = layer[:,[5,6,7]]

kappa = 1e-6
R = 6371e3
vscale = kappa / R
#vscale=1
mps2cmpy = 3.15576e9

# velocity components
vx = layer[:,5]
vy = layer[:,6]
vz = layer[:,7]

# convert lon/lat to radians
phi = np.radians(lon)
theta = np.radians(90 - lat)   # same as awk

# trig
ct = np.cos(theta)
st = np.sin(theta)
cp = np.cos(phi)
sp = np.sin(phi)

# spherical basis vectors
er_x = st * cp
er_y = st * sp
er_z = ct

etheta_x = ct * cp
etheta_y = ct * sp
etheta_z = -st

ephi_x = -sp
ephi_y = cp
ephi_z = 0.0

# projection (dot products)
v_r = er_x*vx + er_y*vy + er_z*vz
v_theta = etheta_x*vx + etheta_y*vy + etheta_z*vz
v_phi = ephi_x*vx + ephi_y*vy + ephi_z*vz

# combine results
vel_spherical = np.column_stack((v_r, v_theta, v_phi))


# dimentional velocity
vp = v_phi * vscale*mps2cmpy
vt = v_theta * vscale*mps2cmpy
vr = v_r * vscale*mps2cmpy

vp_grid = pygmt.surface(
        x=lon, 
        y = lat, 
        z = vp, 
        region = 'd', 
        spacing = 1/1,
        lower="d",
        upper="d",
        )

vt_grid = pygmt.surface(
        x=lon, 
        y = lat, 
        z = vt, 
        region = 'd', 
        spacing = 1/1,
        lower="d",
        upper="d",
        )

vr_grid = pygmt.surface(
        x=lon, 
        y = lat, 
        z = vr, 
        region = 'd', 
        spacing = 1/1,
        lower="d",
        upper="d",
        )


vel_mag = np.sqrt(vt**2 + vp**2)
vel_dir = np.degrees(np.arctan2(vp, -vt))

fig = pygmt.Figure()

fig.grdimage(
    grid=vt_grid,
    region="d",
    projection="W15c",
    cmap=True
)

fig.coast(
    shorelines=True,
    region="d",
    projection="W15c", 
    frame = True
)

scale = 0.1
step = 43

fig.plot(
    x=lon[::step],
    y=lat[::step],
    direction=[vel_dir[::step], vel_mag[::step]*scale],
    style="V0.1c+e",
    fill="black",
    pen="0.5p,black"
)
fig.colorbar(
    frame='af+lvr'
)

fig.show()
#%% verifying from pmodel_vtk

import pygmt
import numpy as np

#model = 'OSL_0_1_100_nocrat_LLNL_129'
model = 'SMEAN2_high_vel_removed_0_1_100_129'

fig = pygmt.Figure()
#region = '0/359/-89.5/89.5'
region='d'
proj = 'H10c'


vt = pygmt.grd2xyz(f"{model}/{model}/0/vt.0.128.grd")
vp = pygmt.grd2xyz(f"{model}/{model}/0/vp.0.128.grd")






velocity_ref = np.loadtxt(f'{model}/{model}/0/vel.0.128.gbl')
subsampled_vel_ref = velocity_ref[::10]
subsampled_vel_ref[:, 3] *= 0.1


fig.grdimage(
    grid = f'{model}/{model}/0/vt.0.128.grd',
    region = region,
    projection = proj,
    cmap = True)
  
fig.coast(
    shorelines = True, 
    region = region, 
    projection = proj, 
    frame = True)
fig.colorbar(
    frame='af+lvr'
)


fig.plot(
    data=subsampled_vel_ref,
    pen = '0.3p, black',
    fill = 'black',
    style = 'V0.15c+e')

fig.show()