# -*- coding: utf-8 -*-
"""
Created on Tue May 13 14:55:54 2025

@author: middl
"""
## --- Libraries ---
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import configparser
from scipy.interpolate import RegularGridInterpolator
import os
from scipy.interpolate import interp1d
import subprocess

## --- Config Setup ---
config = configparser.ConfigParser()
config.read('C:/Users/middl/OneDrive/Documents/config.txt')

# --- Defining Paths ---
bathy_path = config['paths']['bathy']
output_folder = config['paths']['output']

os.makedirs(output_folder, exist_ok=True)


### ----- Make sure this is that this reads the right columns!!! -----

# --- Reading the Bathymetry ---
x = pd.read_csv(bathy_path, usecols=[2], dtype=np.float32).squeeze()
y = pd.read_csv(bathy_path, usecols=[3], dtype=np.float32).squeeze()
z = pd.read_csv(bathy_path, usecols=[1], dtype=np.float32).squeeze()

# --- Change to arrays ---
x = x.to_numpy()
y = y.to_numpy()
z = z.to_numpy()

x_unique = np.unique(x)
y_unique = np.unique(y)

# --- Make into a grid ---
X, Y = np.meshgrid(x_unique, y_unique)

# --- Interpolate z-values onto grid with defined configuration ---
method = config['settings'].get('interpolation_method', 'cubic') # 'cubic' is the default and fallback
Z = griddata((x, y), z, (X, Y), method='cubic')  # 'linear' and 'nearest' also available

# --- Sanity check ---
plt.figure(figsize=(10, 6))
plt.pcolormesh(X, Y, Z, shading='auto', cmap='viridis')
plt.colorbar(label='Elevation')
plt.title('Interpolated Bathymetry')
plt.xlabel('X')
plt.ylabel('Y')
plt.axis('equal')
plt.show()

# --- Define endpoints of the transect  ---
x1, y1 = int(config['transect']['xstart']), int(config['transect']['ystart'])
x2, y2 = int(config['transect']['xend']), int(config['transect']['yend'])
number_of_points = int(config['transect']['number_of_points'])

x_points = np.linspace(x1, x2, number_of_points)
y_points = np.linspace(y1, y2, number_of_points)

# --- Plot original bathymetry and transect line ---
plt.figure(figsize=(10, 6))
plt.pcolormesh(X, Y, Z, shading='auto', cmap='viridis')
plt.plot([x1, x2], [y1, y2], 'r')
plt.colorbar(label='Elevation [m]')
plt.title('Measured Bathymetry [m]')
plt.xlabel('X')
plt.ylabel('Y')
plt.axis('equal')
plt.show()

# --- Interpolate bathymetry along the transect ---
interp_func = RegularGridInterpolator((y_unique, x_unique), Z, bounds_error=False, fill_value=np.nan)

# --- Stack x/y into (N, 2) array of points for interpolation ---
points = np.column_stack((y_points, x_points))
bathymetry_values = interp_func(points)

# --- Plot bathymetry profile along transect ---
plt.figure(figsize=(10, 4))
plt.plot(np.linspace(0, 1, number_of_points), bathymetry_values)
plt.xlabel('Relative Distance Along Transect')
plt.ylabel('Bathymetry (m)')
plt.title('Transect of Bathymetry')
plt.grid(True)
plt.show()

# --- Wave Parameters ---
Hm0 = float(config['waves']['Hm0'])   # Significant wave height
Tp = float(config['waves']['Tp'])    # Wave period
mainang = float(config['waves']['mainang'])  # Wave angle

# --- Storm Surge ---
zs0 = float(config['waves']['zs0'])  # Surge level in meters

# --- Grid Size Parameters ---
dxmin = float(config['grid']['dxmin'])
nx = int(config['grid']['nx'])
ny = int(config['grid']['ny'])
dx = float(config['grid']['dx'])
dy = float(config['grid']['dy'])

x_grid = np.arange(0, nx) * dx
y_grid = np.arange(0, ny) * dy
zgr = bathymetry_values  # Use for 1D transect simulation




# --- Compute transect length ---
L = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

# --- Create x.grd using dxmin spacing ---
dxmin = float(config['grid'].get('dxmin', 1.0))
num_points = int(np.round(L / dxmin)) + 1
xgr = np.linspace(0, L, num_points)

# --- Interpolate z values (bathymetry) onto xgr ---
# Assume zgr is the original elevation profile along the transect
x_original = np.linspace(0, L, len(zgr))
z_interp = interp1d(x_original, zgr, kind='linear', fill_value='extrapolate')
zgr_resampled = z_interp(xgr)

# --- Flatten both arrays to 1D ---
xgr = np.asarray(xgr).flatten()
zgr_resampled = np.asarray(zgr_resampled).flatten()

# --- Validate ---
if len(xgr) != len(zgr_resampled):
    raise ValueError(f"Mismatch: x.grd has {len(xgr)} points, bed.dep has {len(zgr_resampled)}.")

if np.isnan(zgr_resampled).any() or np.isinf(zgr_resampled).any():
    raise ValueError("bed.dep contains NaN or Inf values.")

# --- Save both files ---
output_folder = config['paths']['output']
os.makedirs(output_folder, exist_ok=True)

np.savetxt(os.path.join(output_folder, 'x.grd'), xgr, fmt='%.6f')
np.savetxt(os.path.join(output_folder, 'bed.dep'), zgr_resampled, fmt='%.6f')

print(f"Saved x.grd ({len(xgr)} points) and bed.dep ({len(zgr_resampled)} points) to {output_folder}")





# --- Load from config ---
tintm = int(config['model']['tintm'])
tintp = float(config['model']['tintp'])
tintg = int(config['model']['tintg'])
tstop = int(config['model'].get('tstop', 3600))
rt = int(config['waves'].get('rt', tstop))  # fallback to tstop

def get_bool(config, section, key, default=0):
    try:
        value = config.get(section, key, fallback=str(default)).strip()
        return value.lower() in ['1', 'true', 'yes']
    except Exception:
        return bool(default)

vegetation = get_bool(config, 'options', 'vegetation')
sedtrans = get_bool(config, 'options', 'sedtrans')
morphology = get_bool(config, 'options', 'morphology')

wavemodel = config['model']['wavemodel']




def parse_output_vars(config, section, default):
    raw = config.get(section, 'vars', fallback=default)
    vars_clean = [v.strip() for v in raw.split(',') if v.strip()]
    return vars_clean

# --- Parse each output section ---
global_vars = parse_output_vars(config, 'global_output', 'zs, u, hh, H, zb')
mean_vars = parse_output_vars(config, 'mean_output', 'zs, u')
point_vars = parse_output_vars(config, 'point_output', 'zs, u')

# Output path
params_path = os.path.join(output_folder, 'params.txt')

with open(params_path, 'w') as f:
    f.write(f"order         = 2\n")
    f.write(f"wavemodel     = {wavemodel}\n")
    f.write(f"wbctype       = parametric\n")
    f.write(f"depfile       = bed.dep\n")
    f.write(f"posdwn        = 1\n")
    f.write(f"nx            = {nx}\n")
    f.write(f"ny            = {ny}\n")
    f.write(f"vardx         = 1\n")
    f.write(f"xfile         = x.grd\n")
    f.write(f"xori          = 0\n")
    f.write(f"yori          = 0\n")
    f.write(f"thetamin      = -90\n")
    f.write(f"thetamax      = 90\n")
    f.write(f"dtheta        = 10\n")
    f.write(f"zs0           = {zs0:.6f}\n")
    f.write(f"tstop         = {tstop}\n")
    f.write(f"vegetation    = {int(vegetation)}\n")
    if vegetation:
        f.write(f"veggiefile    = veggiefile.txt\n")
        f.write(f"veggiemapfile = vegmap.txt\n")

    f.write(f"sedtrans      = {int(sedtrans)}\n")
    f.write(f"morphology    = {int(morphology)}\n")

    f.write(f"bcfile        = jonswap.txt\n")
    f.write(f"rt            = {rt}\n")
    f.write(f"outputformat  = netcdf\n")
    f.write(f"tintm         = {tintm}\n")
    f.write(f"tintp         = {tintp:.6f}\n")
    f.write(f"tintg         = {tintg}\n")
    f.write(f"tstart        = 0\n")
    f.write(f"nglobal       = {len(global_vars)}\n")
    for var in global_vars:
        f.write(f"{var}\n")
    f.write("\n")

    f.write(f"nmeanvar      = {len(mean_vars)}\n")
    for var in mean_vars:
        f.write(f"{var}\n")
    f.write("\n")

    f.write(f"npointvar     = {len(point_vars)}\n")
    for var in point_vars:
        f.write(f"{var}\n")
    f.write("\n")
    
    f.write(f"nmeanvar      = 2\nzs\nu\n")
   
    f.write(f"npointvar     = {point_vars}\n")
    for var in point_vars:
        f.write(f"{var}\n")

        f.write(f"npoints       = 2\n100 1\n550 1\n")

probe_raw = config.get('probes', 'locations', fallback='').strip()

if probe_raw:
    # Split and clean the list of locations
    probe_lines = [line.strip() for line in probe_raw.split(',') if line.strip()]
    npoints = len(probe_lines)

    f.write(f"npoints       = {npoints}\n")
    for line in probe_lines:
        f.write(f"{line}\n")
else:
    print("No probe locations defined — skipping npoints section.")


print(f"params.txt saved to {params_path}")

# Run XBeach
xbeach_executable = config.get('paths', 'output', fallback='xbeach.exe')  # or full path if needed
try:
    subprocess.run([xbeach_executable], cwd=output_folder, check=True)
    print("XBeach simulation started successfully.")
except Exception as e:
    print(f"Error running XBeach: {e}")

input("Press Enter to exit...")


