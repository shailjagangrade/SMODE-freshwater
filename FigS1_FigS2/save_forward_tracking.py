#!/python/bash

## simple forward tracking advection model
import numpy as np
import xarray as xr

def forward_tracking(xds,lat_seed,lon_seed,tstart,tend):
    # ------ INPUT ------
    # xds:       xarray containing lat,lon,u,v
    # lat_seed: latitude starting point, ex: 36.8, can also be multiple
    # lon_seed: longitude starting point, ex:-122
    # tstart:   date start ex:'2023-04-15'
    # tend:     date end ex:'2023-05-01'

    # ------ OUTPUT ------
    # lat_frwd: trajectory latitude array
    # lon_frwd: trajectory longitude array

    xds= xds.sel(time=slice(tstart,tend))
    tt=12; print(f'dt divided by {tt}')
    u = xds.u[::tt,...]
    v = xds.v[::tt,...]
    time = xds.time[::tt,...]
    T = len(time)
    dt = float(np.diff(time)[0]*1e-9) # time resolution in [seconds]

    N = np.size(lat_seed) # number of particles
    lat_fwrd = np.empty([T,N])*np.nan
    lon_fwrd = np.empty([T,N])*np.nan

    # initialize particle positions
    lat_fwrd[0,:] = lat_seed
    lon_fwrd[0,:] = lon_seed
    
    # forward velocities: u(t) = dx/dt
    for t in range(T-1):
        # velocity u(t) at particle location
        u_fwrd = u[t,...].sel(lat=lat_fwrd[t,:], lon=lon_fwrd[t,:], method='nearest')
        v_fwrd = v[t,...].sel(lat=lat_fwrd[t,:], lon=lon_fwrd[t,:], method='nearest')

        u_fwrd = np.diag(u_fwrd)
        v_fwrd = np.diag(v_fwrd)
        
        # compute displacement [meters]
        dx = u_fwrd*dt
        dy = v_fwrd*dt
    
        # convert into degrees
        dlat = dy / 111320.
        dlon = dx / (np.cos(np.deg2rad(lat_fwrd[t,:])) * 111320.)

        
        # new position
        lat_fwrd[t+1,:] = lat_fwrd[t,:]+dlat
        lon_fwrd[t+1,:] = lon_fwrd[t,:]+dlon
            
    return lat_fwrd, lon_fwrd

# Find lat and lon pairs for seeding
def find_lat_lon(u,lat,lon):
    u_mean = u[0:5,...].mean(dim='time') # mean u field
    
    ilats = np.where((lat>35) & (lat<40))[0]
    ilons = np.empty(len(ilats),dtype=np.int64); n=0
    for ilat in ilats:
        idx = int(np.argmax(~np.isnan(u_mean[ilat,::-1]).values))
        ilons[n]=len(lon)-idx-4  ; n+=1
    lat_seed = lat[ilats]
    lon_seed = lon[ilons]
    return lat_seed, lon_seed

# Load data
url='http://hfrnet-tds.ucsd.edu/thredds/dodsC/HFR/USWC/6km/hourly/RTV/HFRADAR_US_West_Coast_6km_Resolution_Hourly_RTV_best.ncd'
ds= xr.open_dataset(url).sel(time=slice('2023-02-01','2023-06-01'))
u = ds.u
v = ds.v

# regular spacing
lat=np.linspace(ds.lat[0].values,ds.lat[-1].values,len(ds.lat))
lon=np.linspace(ds.lon[0].values,ds.lon[-1].values,len(ds.lon))

# Forward tracking multiple pairs
lat_seed, lon_seed = find_lat_lon(u,lat,lon)
tstart   = '2023-04-01T00'
tend     = '2023-04-15T00'
lat_fwrd, lon_fwrd = forward_tracking(ds,lat_seed[0:2],lon_seed[0:2],tstart,tend)

# save as netcdf
forward_tracking = xr.Dataset({
        'lat_seed': (('seed'), lat_seed[0:2]),
        'lon_seed': (('seed'), lon_seed[0:2]),
        'lat_fwrd': (('seed','time'), lat_fwrd.T),
        'lon_fwrd': (('seed','time'), lon_fwrd.T),
        'u':        (('time','lat','lon'), u.sel(time=slice(tstart,tend)).values[::12,...]),
        'v':        (('time','lat','lon'), v.sel(time=slice(tstart,tend)).values[::12,...])
    }, coords={
        'seed': np.arange(len(lat_seed[0:2])),
        'time': u.time.sel(time=slice(tstart,tend))[::12],
        'lat': lat,
        'lon': lon
    })
ncname=f'forward_tracking_{tstart}_{tend}.nc'
forward_tracking.to_netcdf(ncname)

# remove tt=12 in function
# remove [0:2] from lat_seed[0:2],lon_seed & 'seed' coords
# remove [::12] from 'time' coords, u & v variable
