import geopandas as gpd
import pandas as pd
import fiona
import rasterio
from rasterio.mask import mask
from rasterio.merge import merge
from rasterio import shutil
import os
import shutil
import glob
from osgeo import ogr
from pathlib import Path
from json import dump
import shapefile as shp
import xml.etree.ElementTree as ET
import warnings


def XMLparser(deliveryMetadata):
    tree = ET.parse(deliveryMetadata)
    root = tree.getroot()
    project = {}
    metadata = {}
    for product in root.findall('{http://xsd.digitalglobe.com/xsd/dm}product'):
        bands = product.findall('{http://xsd.digitalglobe.com/xsd/dm}band')
        ph = product.find('{http://xsd.digitalglobe.com/xsd/dm}pixelHeight').text
        pw = product.find('{http://xsd.digitalglobe.com/xsd/dm}pixelWidth').text
        bpp = product.find('{http://xsd.digitalglobe.com/xsd/dm}bitsPerPixel').text
        offNadir = product.find('{http://xsd.digitalglobe.com/xsd/dm}offNadirAngle').text
        for band in bands:
            metadata[band.text] = {'PixelHeight' : ph, 'PixelWidth': pw, 'OffNadirAngle' : offNadir, 'BitsPerPixel' : bpp}
    return metadata

def multi2single(gpdf):
    gpdf_singlepoly = gpdf[gpdf.geometry.type == 'Polygon']
    gpdf_multipoly = gpdf[gpdf.geometry.type == 'MultiPolygon']
    for i, row in gpdf_multipoly.iterrows():
        Series_geometries = pd.Series(row.geometry)
        df = pd.concat([gpd.GeoDataFrame(row, crs = gpdf_multipoly.crs).T] * len(Series_geometries), ignore_index = True)
        df['geometry']  = Series_geometries
        gpdf_singlepoly = pd.concat([gpdf_singlepoly, df])
    gpdf_singlepoly.reset_index(inplace = True, drop = True)
    return gpdf_singlepoly

def getFeatures(gds):
    """Function to parse features from GeoDataFrame in such a manner that rasterio wants them"""
    import json
    return [json.loads(gds.to_json())['features'][0]['geometry']]

def getTile(tileName, directory):
    tiles = Path(directory).iterdir()
    for tile in tiles:
        if os.path.isfile(tile) and tileName in tile.__str__() :
            return tile.__str__()

def createChip(data_row, band, folder, path):
    fileName = ''
    tileName = str(data_row['tileName'])
    fileName = getTile(tileName, folder)
    inputRaster = fileName
    try:
        with rasterio.open(inputRaster, 'r') as src:
            out_image, out_transform = mask(src, getFeatures(gpd.GeoSeries(data_row['geometry'])), crop = True)
            out_meta = src.meta
            out_meta.update({"driver": "GTiff",
                                 "height": out_image.shape[1],
                                 "width": out_image.shape[2],
                                 "transform": out_transform})
        src.close()
        out = os.path.join(path, str(data_row['UID']), str(data_row['index_right']) + band + ".tif")
        with rasterio.open(out, "w", **out_meta) as dest:
            dest.write(out_image)
        dest.close()
    except TypeError:
        print (inputRaster)
    return out
def mergeChips(image, chips, band):
    for uid in chips:
        rasters = []
        for chip in chips[uid]:
            src = rasterio.open(chip)
            rasters.append(src)

        mosaic, out_trans = merge(rasters)
        out_meta = src.meta.copy()
        out_meta.update({"driver": "GTiff",
                         "height": mosaic.shape[1],
                         "width": mosaic.shape[2],
                         "transform": out_trans})
        out = os.path.join("CHIPS", str(uid), image + '_' + str(uid) + '_' + band + ".tif")
        with rasterio.open(out, "w", **out_meta) as dest:
            dest.write(mosaic)
        dest.close()
        for raster in rasters:
            raster.close()
        for chip in chips[uid]:
            rasterio.shutil.delete(chip)
    
    
def createTrainingData(image, features, tileShape, panFolder, mulFolder, deliveryXML):
    g1 = gpd.GeoDataFrame.from_file(features)
    tile_shp = gpd.read_file(tileShape)    
    buffer = g1.copy()
    buffer.geometry = buffer['geometry'].buffer(50)
    u = buffer.unary_union
    un = gpd.GeoDataFrame(crs = g1.crs, geometry = [u])
    un_single = multi2single(un)
    mbr = un_single.envelope
    mbrf = gpd.GeoDataFrame(geometry = gpd.GeoSeries(mbr), columns=['UID'])
    mbrf.crs = mbr.crs
    

    #probably better defined elsewhere
    path = "CHIPS"
    try:
        os.makedirs(path)
    except OSError:
        shutil.rmtree(path)
        os.makedirs(path)

    for i in mbrf.index:
        mbrf.at[i, 'UID'] = i
        #s = gpd.GeoSeries(mbrf.at[i, 'geometry'])
        #s.crs = mbr.crs
        a = g1.within(mbrf.at[i, 'geometry'])
        b = g1[a]
        just_geom = b[["geometry"]]
        c = just_geom.copy()
        c.crs = mbr.crs
        try:
            chip_path = os.path.join(path, str(i))
            os.makedirs(chip_path)
            #c.to_file(chip_path + '\shapes.shp')
            c.to_file(chip_path + '\\' + image + '_' + str(i) + '_anno.geojson', driver='GeoJSON')
            with open(chip_path + '\\' + image + '_' + str(i) + '_metadata.json', 'w') as outfile:
                dump(deliveryXML, outfile)
            
        except OSError:
            shutil.rmtree(path)
            os.makedirs(path)   
    mbrf.to_file('mbr.shp')
    tile = tile_shp[['geometry', 'tileName']]
    tile = tile.to_crs(mbr.crs)
    join = gpd.sjoin(mbrf, tile, how = 'left', op = 'intersects')

    panChips = {}
    mulChips = {}

    for index, row in join.iterrows():
        if row['UID'] in panChips:
            panChips[row['UID']].append(createChip(row, 'PAN', panFolder, path))
        else:
            panChips[row['UID']] = [createChip(row, 'PAN', panFolder, path)]
        if row['UID'] in mulChips:
            mulChips[row['UID']].append(createChip(row, 'MUL', mulFolder, path))
        else:
            mulChips[row['UID']] = [createChip(row, 'MUL', mulFolder, path)]

    mergeChips(image, panChips, 'PAN')
    mergeChips(image, mulChips, 'MUL')

def treeDir(image):
    features =''
    tileShape = ''
    panFolder = ''
    mulFolder = ''
    os.chdir(image) #we switch to the main image folder to change our context
    print(os.getcwd())
    if os.path.isfile(image.__str__() + ".shp"):
        features = image.__str__() + ".shp" #see if our features exists       
        subs = Path('.').iterdir() #need to find the data folder which matches this pattern 012232499010_01_003\012232499010_01
        for sub in subs:
            if os.path.isdir(sub):
                x = sub.__str__().split('_')
                if len(x) > 2: 
                    rootFolder = sub.__str__() #this is our root data folder, which contains the GIS_FILEs directory and MUL/PAN folders
                    folders = Path(sub).iterdir()
                    for folder in folders:
                        if os.path.isfile(folder) and "DeliveryMetadata.xml" in folder.__str__():
                            deliveryMetadata = XMLparser(folder.__str__())
                        if os.path.isdir(folder):                                   
                            tifFolders = Path(folder).iterdir()
                            for tifFolder in tifFolders:
                                if "PAN" in tifFolder.__str__():
                                    panFolder = tifFolder
                                if "MUL" in tifFolder.__str__():
                                    mulFolder = tifFolder
                                if "GIS_FILES" in tifFolder.__str__():
                                    shps = os.listdir(tifFolder)
                                    for shp in shps:
                                        if "TILE_SHAPE.shp" in shp.__str__():
                                            t = os.path.join(tifFolder, shp)
                                            if os.path.isfile(t):
                                                tileShape = t.__str__()
        createTrainingData(image, features, tileShape, panFolder, mulFolder, deliveryMetadata)                         
                                            
    else:
        print("Couldn't get the features file for " + image.__str__())
warnings.filterwarnings("error")
DIRNAMES=1
images = next(os.walk('.'))[DIRNAMES]
parent = os.getcwd()
for image in images:
    treeDir(image)
    os.chdir(parent)