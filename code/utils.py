import numpy as np
import xarray as xr
import rioxarray as rio
from rasterio.features import shapes
from rasterio.enums import Resampling
import geopandas as gpd
import pandas as pd
from geocube.api.core import make_geocube
import dask

# define function to calculate surface area of each pixel
def calc_pixel_area(raster:xr.DataArray) -> xr.DataArray:
    '''
    Calculate the area of each pixel in a raster

    Parameters:
    raster (xarray.DataArray): raster to calculate pixel area for

    Returns:
    xarray.DataArray: raster with pixel area as values
    '''

    # get the resolution of the raster
    res = raster.rio.resolution()

    l1 = np.radians(raster['y']- np.abs(res[1])/2)
    l2 = np.radians(raster['y']+ np.abs(res[1])/2)
    dx = np.radians(np.abs(res[0]))    
    _R = 6371e3  # Radius of earth in m. Use 3956e3 for miles

    # calculate the area of each pixel
    area = _R**2 * dx * (np.sin(l2) - np.sin(l1))

    # create a new xarray with the pixel area as values
    result = ((raster-raster+1)*area)

    # set the nodata value    
    if raster.rio.nodata is None:
        result.rio.set_nodata(np.nan,inplace=True)
    else:
        result.rio.set_nodata(raster.rio.nodata,inplace=True)
    
    return result

def calc_area(raster: xr.DataArray, surf_area_file = '../results/00_preprocessing/land_surface_area.nc') -> xr.DataArray:
    '''
    Calculate the area of each raster cell

    Parameters:
    raster: xarray.DataArray
        The raster to calculate the area for
    surf_area_file: str
        The path to the terrestrial area raster file

    Returns:
    xarray.DataArray
        The area of each raster cell
    '''
    # load the terrestrial area raster file
    surf_area = rio.open_rasterio(surf_area_file, masked=True).drop_vars('band').squeeze().rio.set_nodata(0)
    
    # set the nodata value to np.nan if it is None
    _nodata = np.nan if raster.rio.nodata is None else surf_area.rio.nodata
    _nodata_surf = np.nan if surf_area.rio.nodata is None else surf_area.rio.nodata
    
    # calculate the area of each raster cell
    result = surf_area.rio.set_nodata(_nodata_surf).rio.reproject_match(raster,resampling=Resampling.sum).rio.set_nodata(_nodata)

    result = result.where(result>0) #TODO - test the effect
    return result

def raster_vector_zonal_stats(vector: gpd.GeoDataFrame,raster: xr.DataArray,stat:str,interp=False) -> pd.DataFrame:
    '''
    Calculate zonal statistics for a raster and a vector

    Parameters:
    vector: geopandas.GeoDataFrame
        The vector to calculate the zonal statistics for
    raster: xarray.DataArray
        The raster to calculate the zonal statistics for
    stat: str    
        The statistic to calculate. Either 'sum' or 'mean'
    interp: bool
        Whether to interpolate the index variable before calculating the statistics
    
    Returns:
    pandas.DataFrame
        The zonal statistics for the vector

    '''

    # create a xarray.Dataset with two variables: 'index' and 'raster' - the index is the data from the vector and the raster is the data from the raster file
    vector_ras= make_geocube(
        vector_data = vector.reset_index(),
        measurements = ['index'],
        like = raster
    )
    vector_ras['raster'] = (raster.dims, raster.values, raster.attrs, raster.encoding)

    # interpolate NA values in index variable if interp is True
    if interp:
        vector_ras['index'] = vector_ras['index'].rio.interpolate_na()
    
    # calculate the zonal statistics - group by according to index and calculate the sum or mean of the raster values
    if stat == 'sum':
        res = vector_ras.drop_vars('spatial_ref').set_coords('index').groupby('index').sum()['raster']
    elif stat == 'mean':
        res = vector_ras.drop_vars('spatial_ref').set_coords('index').groupby('index').mean()['raster']
    
    # copy the dimensions of the raster to the result
    for d in res.dims[1:]:
        res[d] = raster[d]

    # return the result as a pandas DataFrame
    return res.to_dataframe()['raster']

def down_sample(src: xr.DataArray,x_factor:int,y_factor:int,stat:str) -> xr.DataArray:
    '''
    Down sample a raster by a given factor

    Parameters:
    src: xarray.DataArray
        The raster to down sample
    x_factor: int
        The factor to down sample in x-direction
    y_factor: int
        The factor to down sample in y-direction
    stat: str
        The statistic to calculate. Either 'sum' or 'mean'

    Returns:
    xarray.DataArray
        The down sampled raster
    '''
    with dask.config.set(**{'array.slicing.split_large_chunks': False}): # set dask config to avoid large chunk splitting
        # get source resolution
        src_res = np.array(src.rio.resolution())

        # set the nodata value to np.nan if it is None
        src.rio.write_nodata(np.nan if src.rio.nodata is None else src.rio.nodata,inplace=True)

        # down sample the raster by the given factors
        res = src.coarsen(x=x_factor,y=y_factor,boundary='pad')

        # calculate the sum or mean of the raster values
        if stat == 'sum':
            result = res.sum().rio.write_nodata(src.rio.nodata)
        elif stat == 'mean':
            result = res.mean().rio.write_nodata(src.rio.nodata)

        # calculate the new transform and set it to the result
        new_transform,new_width,new_height = rio.rioxarray.raster_array._make_dst_affine(src, src.rio.crs,src.rio.crs,src_res*x_factor)
        result.rio.write_transform(new_transform,inplace=True)

        return result

def resample_match(src: xr.DataArray, dst: xr.DataArray, weighted=True) -> xr.DataArray:
    '''
    Resample a raster to match the resolution of another raster

    Parameters:
    src: xarray.DataArray
        The raster to resample
    dst: xarray.DataArray
        The raster to match the resolution to
    weighted: bool
        Whether to use weighted resampling

    Returns:
    xarray.DataArray
        The resampled raster
    '''

    with dask.config.set(**{'array.slicing.split_large_chunks': False}):
            src.rio.write_nodata(np.nan if src.rio.nodata is None else src.rio.nodata,inplace=True)
            dst.rio.write_nodata(np.nan if dst.rio.nodata is None else src.rio.nodata,inplace=True)
            src_res = np.array(src.rio.resolution())
            dst_res =  np.array(dst.rio.resolution())
            factor = np.round(dst_res/src_res,0).astype(int)
            area = calc_area(src).rio.write_nodata(np.nan if src.rio.nodata is None else src.rio.nodata)
            if weighted:
                area_down_sample = area.coarsen(x=factor[0],y=factor[1],boundary='pad').sum().rio.write_nodata(src.rio.nodata)
                down_sample = ((src*area).coarsen(x=factor[0], y=factor[1],boundary='pad').sum()/area_down_sample).rio.write_nodata(src.rio.nodata)
            else:
                down_sample = src.coarsen(x=factor[0], y=factor[1]).mean().rio.write_nodata(src.rio.nodata)
            
            new_transform,new_width,new_height = rio.rioxarray.raster_array._make_dst_affine(down_sample, down_sample.rio.crs,down_sample.rio.crs,src_res*factor)
            down_sample.rio.write_transform(new_transform,inplace=True)

            if weighted:
                area_reproj = area_down_sample.rio.reproject_match(dst,resampling=Resampling.sum)
                result = ((down_sample*area_down_sample).transpose(*(src.dims)).rio.reproject_match(dst,resampling=Resampling.sum)/area_reproj).rio.write_nodata(src.rio.nodata)
                        
            else:
                result = down_sample.rio.reproject_match(dst,resampling=Resampling.average).rio.write_nodata(src.rio.nodata)
            
            return result

def polygonize(raster:xr.DataArray) -> gpd.GeoDataFrame:
    '''
    Convert a raster to a polygon

    Parameters:
    raster: xarray.DataArray
        The raster to convert to a polygon

    Returns:
    geopandas.GeoDataFrame
        The polygonized raster
    '''

    # get the shapes of the raster and store them in a list of dictionaries with keys id and geometry
    results = ({'properties': {'id': v}, 'geometry': s} for i, (s, v) in enumerate(shapes(raster.astype(np.int16), transform=raster.rio.transform())))

    # convert the list of dictionaries to a GeoDataFrame and dissolve by id
    out = gpd.GeoDataFrame.from_features(list(results)).dissolve(by='id')

    # set the crs of the GeoDataFrame to the crs of the raster
    out.crs = raster.rio.crs
    
    return out
