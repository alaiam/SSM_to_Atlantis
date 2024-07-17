import numpy as np
import xarray as xr
import dask.array as da
import netCDF4
import datetime as dt
import os
import dask
import pyinterp
from dask.distributed import Client
from multiprocessing import freeze_support #might need to make this windows-compatible


if __name__ == '__main__':

    # Create a Dask client for parallel computation
    client = Client()

    # Create an RTree instance for spatial indexing using pyinterp
    mesh = pyinterp.RTree()

    def kriging_universal(values, lons, lats):
        # Pack the original data into the RTree for spatial indexing
        mesh.packing(np.vstack((lons, lats)).T, values)
        # Perform universal kriging to interpolate new_lon, new_lat points
        kriging, neighbors = mesh.universal_kriging(np.vstack((lats.ravel(), lons.ravel())).T, within=True, k=11,
                                                    covariance='matern_12', alpha=1_000_000, num_threads=0)
        return kriging.reshape(lons.shape)

    # Define the file name of the original Salish Sea Model NetCDF data
    filename = "//nfsdata/SSM_example/HYD/ssm_00115.nc"

    # Open the NetCDF file and read the original Salish Sea Model data with Dask
    ssm_solution = xr.open_dataset(
        filename, 
        decode_cf=True, decode_times=False, 
        chunks={'time': 100, 'lat': 100, 'lon': 100}
        )

    print('NetCDF file read!')


    lats = ssm_solution.lat.values
    lons = ssm_solution.lon.values

    # mx, my = np.meshgrid(lats, lons, indexing='ij', copy=False)
    # print(f'mx shape: {mx.shape} my shape: {my.shape}')


    original_siglay = ssm_solution.siglay.values
    original_siglev = ssm_solution.siglev.values
    original_time = ssm_solution.time.values

    siglay_size = len(original_siglay)
    time_size = len(original_time)

    # Create Dask arrays for interpolated data
    new_regular_temp = da.empty((time_size, siglay_size, len(lats), len(lons)), chunks=(1, siglay_size, 100, 100))
    new_regular_salt = da.empty((time_size, siglay_size, len(lats), len(lons)), chunks=(1, siglay_size, 100, 100))
    new_regular_sigmalay = da.empty((siglay_size, len(lats), len(lons)), chunks=(siglay_size, 100, 100))
    new_regular_u = da.empty((time_size, siglay_size, len(lats), len(lons)), chunks=(1, siglay_size, 100, 100))
    new_regular_v = da.empty((time_size, siglay_size, len(lats), len(lons)), chunks=(1, siglay_size, 100, 100))

    # Interpolate sigma layer values in parallel
    def interpolate_sigma(d):
        org_sigma = ssm_solution.siglay[d].values
        return kriging_universal(org_sigma, lons, lats)

    new_regular_sigmalay = da.map_blocks(interpolate_sigma, da.arange(siglay_size, chunks=1), dtype=np.float64)

    print('Interpolation sigma layer done!')
    # Interpolate temperature and salinity values in parallel
    def interpolate_temp(t, d):
        new_temp = kriging_universal(ssm_solution.temp[t, d].values, lons, lats)
        return new_temp

    def interpolate_salt(t, d):
        new_salt = kriging_universal(ssm_solution.salinity[t, d].values, lons, lats)
        return new_salt

    print('Interpolation variables done!')

    new_regular_temp  = da.map_blocks(interpolate_temp, da.arange(time_size, chunks=1), da.arange(siglay_size, chunks=1), dtype=np.float64)

    new_regular_salt = da.map_blocks(interpolate_salt, da.arange(time_size, chunks=1), da.arange(siglay_size, chunks=1), dtype=np.float64)

    # Interpolate velocity fields in parallel
    def interpolate_u(t, d):
        new_u = kriging_universal(ssm_solution.u[t, d].values, lons, lats)
        return new_u
    def interpolate_v(t, d):
        new_v = kriging_universal(ssm_solution.v[t, d].values, lons, lats)
        return new_v

    new_regular_u = da.map_blocks(
        interpolate_u, da.arange(time_size, chunks=1), da.arange(siglay_size, chunks=1), dtype=np.float64)
    new_regular_v = da.map_blocks(
        interpolate_v, da.arange(time_size, chunks=1), da.arange(siglay_size, chunks=1), dtype=np.float64)

    # Compute the results in parallel
    new_regular_temp = dask.compute(new_regular_temp)
    print("temp done")
    new_regular_salt = dask.compute(new_regular_salt)
    print("salt done")
    new_regular_sigmalay = dask.compute(new_regular_sigmalay)
    print("sigmalay done")
    new_regular_u = dask.compute(new_regular_u)
    print("u done")
    new_regular_v = dask.compute(new_regular_v)



    # new_regular_temp, new_regular_salt, new_regular_sigmalay, new_regular_u, new_regular_v = dask.compute(
    #     new_regular_temp, 
    #     new_regular_salt, 
    #     new_regular_sigmalay, 
    #     new_regular_u, 
    #     new_regular_v
    #     )

    print("all done!")

    # Create a new NetCDF file with the interpolated data
    data_path = r'.'
    file_name_output = 'regular_grid_' + filename[-12:]
    nc = netCDF4.Dataset(os.path.join(data_path, file_name_output), 'w')

    # Define NetCDF global attributes
    nc.title = 'Regular grid'
    nc.Conventions = 'CF-1.6'
    nc.history = '{0} creation of regular grid NetCDF file'.format(
        dt.datetime.now().strftime("%Y-%m-%d"))

    # Create NetCDF variables and set attributes
    lat_dim = nc.createDimension('latitude', len(lats))
    lon_dim = nc.createDimension('longitude', len(lons))
    siglev_dim = nc.createDimension('sigma_layer', siglay_size)
    time_dim = nc.createDimension('time', time_size)

    lat_var = nc.createVariable('latitude', np.single, ('latitude'))
    lat_var.units = 'degrees_north'
    lat_var.standard_name = 'latitude'
    lat_var.axis = 'Y'
    lat_var[:] = lats.astype('float')

    lon_var = nc.createVariable('longitude', np.single, ('longitude'))
    lon_var.units = 'degrees_east'
    lon_var.standard_name = 'longitude'
    lon_var.axis = 'X'
    lon_var[:] = lons.astype('float')

    time_var = nc.createVariable('time_vector', np.intc, ('time'))
    time_var.units = 'days since 1858-11-17 00:00:00'
    time_var.standard_name = 'time'
    time_var.format = 'modified julian day (MJD)'
    time_var.time_zone = 'UTC'
    time_var[:] = original_time.astype('int')

    siglay_var = nc.createVariable('siglay', np.single, ('sigma_layer'))
    siglay_var.units = 'sigma_layers'
    siglay_var.standard_name = 'ocean_sigma/general_coordinate'
    siglay_var[:] = original_siglay.astype('float')

    siglev_var = nc.createVariable(
        'siglev', np.single, ('sigma_layer', 'latitude', 'longitude'))
    siglev_var.units = 'sigma_level'
    siglev_var.standard_name = 'ocean_sigma/general_coordinate'
    siglev_var[:] = new_regular_sigmalay.astype('float')

    temp_var = nc.createVariable(
        'temp', np.single, ('time', 'sigma_layer', 'latitude', 'longitude'))
    temp_var.units = 'degrees_C'
    temp_var.standard_name = 'sea_water_temperature'
    temp_var[:] = new_regular_temp.astype('float')

    salt_var = nc.createVariable(
        'salinity', np.single, ('time', 'sigma_layer', 'latitude', 'longitude'))
    salt_var.units = '1e-3'
    salt_var.standard_name = 'sea_water_salinity'
    salt_var[:] = new_regular_salt.astype('float')

    u_var = nc.createVariable(
        'u', np.single, ('time', 'sigma_layer', 'latitude', 'longitude'))
    u_var.units = 'm s-1'
    u_var.standard_name = 'eastward_sea_water_velocity'
    u_var[:] = new_regular_u.astype('float')

    v_var = nc.createVariable(
        'v', np.single, ('time', 'sigma_layer', 'latitude', 'longitude'))
    v_var.units = 'm s-1'
    v_var.standard_name = 'northward_sea_water_velocity'
    v_var[:] = new_regular_v.astype('float')

    nc.close()
    print('NetCDF file created!')


    # Close the Dask client
    client.close()



