import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import geopandas as gpd
import mplleaflet
import pandas as pd
import subprocess
from pysheds.grid import Grid
import os
import shutil
import pyproj
from pyproj import CRS
from IPython.display import IFrame
from polycircles import polycircles
from functools import partial
from shapely.ops import transform
from shapely.geometry import Point
import argparse
# %matplotlib inline


ap = argparse.ArgumentParser()
ap.add_argument("-d", "--dem", required=True,
	help="path to input dem")
ap.add_argument("-f", "--farmer", required=True,
	help="path to farmer data")
ap.add_argument("-lat", "--latitude", required=True, type=float,
	help="Enter the latitude")
ap.add_argument("-lon", "--longitude", required=True, type=float,
	help="Enter the longitude")
ap.add_argument("-b", "--buffer", required=True, type=float,
	help="Enter the buffer distance")
ap.add_argument("-t","--threshold", type=int, default=1000,
	help="give threshold value for extracting the river network")
args = vars(ap.parse_args())


def create_folder(foldername):
    pathaddress = os.getcwd()
    newfolder = pathaddress +"/"+foldername
    foldercheck = os.path.exists(newfolder)
    if foldercheck==False:
        os.mkdir(foldername)
        return "new folder, {} created.".format(foldername)
    else:
        shutil.rmtree(newfolder)
        os.mkdir(foldername)
        return "folder already exist"

def flow_acc(path_dem):
  grid = Grid.from_raster(path_dem, data_name='dem')
  print("The raster has been read")

  # plt.figure(figsize=(12,10))
  # plt.imshow(grid.dem, extent=grid.extent, cmap='Blues')
  # plt.imshow(data, cmap=cmap)

  # plt.colorbar(label='Elevation (m)')
  # plt.grid()

  # plotFigure(grid.dem, 'Elevation (m)')
  
  depressions = grid.detect_depressions('dem')
   
  # Plot pits
  # fig, ax = plt.subplots(figsize=(8,8))
  # ax.imshow(depressions, cmap='cubehelix', zorder=1)
  # ax.set_yticklabels([])
  # ax.set_xticklabels([])
  # plt.title('Depressions', size=14)

  grid.fill_depressions(data='dem', out_name='flooded_dem')
  flats = grid.detect_flats('flooded_dem')

  # Plot flats
  # fig, ax = plt.subplots(figsize=(8,8))
  # plt.imshow(flats, cmap='cubehelix', zorder=1)
  # ax.set_yticklabels([])
  # ax.set_xticklabels([])
  # plt.title('Flats', size=14)

  grid.resolve_flats(data='flooded_dem', out_name='inflated_dem')

  # plt.figure(figsize=(12,10))
  # plt.imshow(grid.inflated_dem, extent=grid.extent, cmap='cubehelix')
  # plt.imshow(data, cmap=cmap)

  # plt.colorbar(label='Flats')
  # plt.grid()


  # plotFigure(grid.inflated_dem, 'flats', 'cubehelix')

  flats_test = grid.detect_flats('inflated_dem')
  print("The Depressions have been filled and the Flats have been resolved")
  # Plot flats
  # fig, ax = plt.subplots(figsize=(8,8))
  # plt.imshow(flats_test, cmap='cubehelix', zorder=1)
  # ax.set_yticklabels([])
  # ax.set_xticklabels([])
  # plt.title('Flats', size=14)

  #N    NE    E    SE    S    SW    W    NW
  dirmap = (64,  128,  1,   2,    4,   8,    16,  32)

# Compute flow direction based on corrected DEM
  grid.flowdir(data='inflated_dem', out_name='dir', dirmap=dirmap)
  print("Flow direction is computed")
  # plt.figure(figsize=(12,10))
  # plt.imshow(grid.dir, extent=grid.extent, cmap='viridis')
  # plt.imshow(data, cmap=cmap)

  # plt.colorbar(label='flow_dir')
  # plt.grid()

  # plotFigure(grid.dir , 'flow dir', 'viridis')

  # Compute flow accumulation based on computed flow direction
  grid.accumulation(data='dir', out_name='acc', dirmap=dirmap)
  accView = grid.view('acc', nodata=np.nan)
  print("Flow accumulation is computed")
  # plt.figure(figsize=(12,10))
  # plt.imshow(accView, extent=grid.extent, cmap='PuRd')
  # plt.imshow(data, cmap=cmap)

  # plt.colorbar(label='Cell_number')
  # plt.grid()

  # plotFigure(accView,"Cell Number",'PuRd')

  # Delineate catchment at point of high accumulation
  y, x = np.unravel_index(np.argsort(grid.acc.ravel())[-2], grid.acc.shape)
  grid.catchment(x, y, data='dir', out_name='catch', dirmap=dirmap, xytype='index')
  print("Catchment Delineated")
  
  streams = grid.extract_river_network('catch', 'acc', args["threshold"], dirmap=dirmap)
  print(streams["features"][:2])
  print("River network has been extracted")
  f = open('OUTPUT/streams.geojson','w')
  f.write(str(streams))
  f.close()
    
  # streamNet = gpd.read_file('streams.geojson')

  # The polygonize argument defaults to the grid mask when no arguments are supplied
  # shapes = grid.polygonize()

  # Plot catchment boundaries
  # fig, ax = plt.subplots(figsize=(6.5, 6.5))

  # for shape in shapes:
  #     coords = np.asarray(shape[0]['coordinates'][0])
  #     ax.plot(coords[:,0], coords[:,1], color='cyan')
      
  # ax.set_xlim(grid.bbox[0], grid.bbox[2])
  # ax.set_ylim(grid.bbox[1], grid.bbox[3])
  # ax.set_title('Catchment boundary (vector)')
  # gpd.plotting.plot_dataframe(streamNet, None, cmap='Blues', ax=ax)

#   df = pd.read_csv(path_loc)
#   df1 = df[['lattitude','longitude']]
#   geometry = [Point(xy) for xy in zip(df1['longitude'] , df1['lattitude'])]

#   # crs = {'init': 'epsg:4326'}
#   crs = CRS('epsg:4326')

#   gdf = gpd.GeoDataFrame(df1 , crs=crs , geometry=geometry)
#   gdf.to_file(filename = 'output/pointfile.geojson' , driver = 'GeoJSON')
#   df_new = gpd.read_file('output/pointfile.geojson')
# # ax = df_new.plot(figsize =(10,10))
  # mplleaflet.show(fig=ax.figure, epsg=4326)



def farmer_loc_buffer_zone(path, lat, lon, km):
    df = pd.read_csv(path)
    df1 = df[['lattitude','longitude']]
    geometry = [Point(xy) for xy in zip(df1['longitude'] , df1['lattitude'])]

    # crs = {'init': 'epsg:4326'}
    crs = CRS('epsg:4326')

    gdf = gpd.GeoDataFrame(df1 , crs=crs , geometry=geometry)
    gdf.to_file(filename = 'OUTPUT/loc.geojson', driver='GeoJSON')
#     df_new = gpd.read_file('.geojson')
    print("The farmer locations are plotted")
    proj_wgs84 = pyproj.Proj('+proj=longlat +datum=WGS84')


    def geodesic_point_buffer():
        # Azimuthal equidistant projection
        aeqd_proj = '+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0'
        project = partial(
            pyproj.transform,
            pyproj.Proj(aeqd_proj.format(lat=lat, lon=lon)),
            proj_wgs84)
        buf = Point(0, 0).buffer(km * 1000)  # distance in metres
        return transform(project, buf).exterior.coords[:]
    
    b = geodesic_point_buffer()
    dff = pd.DataFrame(b,columns=['LAT','LONG'])
#     print (dff)
    geometry1 = [Point(xy) for xy in zip(dff['LAT'] , dff['LONG'])]

    crs = CRS('epsg:4326')
    gdf_buffer = gpd.GeoDataFrame(dff , crs=crs , geometry=geometry1)
    gdf_buffer.to_file(filename = 'OUTPUT/buffer.geojson', driver='GeoJSON')
    print("Buffer created")
    print("\n\nThe results have been added onto the OUTPUT folder")

#     df_buffer = gpd.read_file('buffer.geojson')


create_folder('OUTPUT')
flow_acc(args["dem"])
farmer_loc_buffer_zone(args["farmer"],args["latitude"],args["longitude"],args["buffer"])