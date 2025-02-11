import numpy as np
import glob
import wrf
from wrf import (to_np, getvar, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)
from datetime import datetime, timedelta
from netCDF4 import Dataset,num2date,date2num
import netCDF4
from scipy.ndimage import gaussian_filter
import argparse
import os 
import cartopy.io.shapereader as shp
from cartopy.feature import ShapelyFeature
import cartopy.crs as crs
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib

from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Polygon
import pandas as pd

plt.rcParams['pcolor.shading'] ='nearest'
#############################################################################################################
#############################################################################################################

# Ejemplo de uso:
# python regridding.py './TEST/*.nc4' '/home/msuarez/Documents/Regridding/geo_em.d03.nc' './TEST/output'

# Modificar las opciones de ejecucion y el directorio de los shapefiles

#############################################################################################################
#############################################################################################################

# Opciones de Ejecucion

# True para generar una imagen del regrillado
plot_regridded = True
# True para generar una imagen del producto imerg (sirve solo para centro de Argentina!)
plot_imerg = True      
# True para mostrar detalle de la ejecucion
debug = False 
# True para calcular la precipitación media en el dominio de regrillado y exportar en un archivo csv
mean_precipitation = True          

#############################################################################################################
#############################################################################################################


def get_plot_element(infile):
    import netCDF4 as nc
    rootgroup = nc.Dataset(infile, 'r')
    p = wrf.getvar(rootgroup, 'HGT_M')
    #lats, lons = wrf.latlon_coords(p)
    cart_proj = wrf.get_cartopy(p)
    xlim = wrf.cartopy_xlim(p)
    ylim = wrf.cartopy_ylim(p)
    rootgroup.close()
    return cart_proj, xlim, ylim

def draw_screen_poly( lats, lons, m, edgecolor):
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( list(xy), facecolor='none', edgecolor=edgecolor, lw=3 )
    plt.gca().add_patch(poly)


parser = argparse.ArgumentParser()

parser.add_argument("imerg_files", type=str, help='Path a los archivos IMERG. Ej: /home/user/IMERG-Files/*.nc4')
parser.add_argument("geo_file"   , type=str, help='Path al archivo geogrid. Ej: /home/user/geo_em.d01.nc')
parser.add_argument("output_path", type=str, help='Path al directorio donde se guardaran los archivos. Ej: /home/user/output_dir')

args = parser.parse_args()

imerg_path = args.imerg_files
geo_path = args.geo_file
output_dir = args.output_path

# If folder doesn't exist, then create it.
check_output_dir = os.path.isdir(output_dir)
if not check_output_dir:
    os.makedirs(output_dir)
    print("created folder : ", output_dir)
else:
    print(output_dir, "folder already exists.")


imergFilelist = glob.glob(imerg_path)
geofile = geo_path

imergFilelist.sort()

geof = Dataset(geofile)
var = wrf.getvar(geof, 'HGT_M')

nws_precip_colors = [
    "#fdfdfd",
    "#04e9e7",  # 0.01 - 0.10 inches
    "#019ff4",  # 0.10 - 0.25 inches
    "#0300f4",  # 0.25 - 0.50 inches
    "#02fd02",  # 0.50 - 0.75 inches
    "#01c501",  # 0.75 - 1.00 inches
    "#008e00",  # 1.00 - 1.50 inches
    "#fdf802",  # 1.50 - 2.00 inches
    "#e5bc00",  # 2.00 - 2.50 inches
    "#fd9500",  # 2.50 - 3.00 inches
    "#fd0000",  # 3.00 - 4.00 inches
    "#d40000",  # 4.00 - 5.00 inches
    "#bc0000",  # 5.00 - 6.00 inches
    "#f800fd",  # 6.00 - 8.00 inches
    "#9854c6"  # 8.00 - 10.00 inches  # 10.00+
]
precip_colormap = matplotlib.colors.ListedColormap(nws_precip_colors)

if plot_regridded or plot_imerg: 
    cart_proj, xlim_d01, ylim_d01 = get_plot_element(geofile)

    shapes_dir = '/home/msuarez/Documents/ohmc_workspace/Shapes'
    dptos = shp.Reader(shapes_dir+'/departamentos.shp' )
    pcias = shp.Reader(shapes_dir+'/provincias.shp' )

    shape_dp = ShapelyFeature(dptos.geometries(),
                              crs.PlateCarree(), edgecolor='black')
    shape_pcias = ShapelyFeature(pcias.geometries(),
                              crs.PlateCarree(), edgecolor='black')

lats_geofile, lons_geofile = latlon_coords(var)
geof.close()
#############################################################################################################
#############################################################################################################

def find_nearest(array, value):
    '''
    Busca el valor mas cercano
    array: vector 1D
    value: valor
    Retorna el valor (array[idx]) y el indice (idx) (ubicacion en el array del valor)
    '''

    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx], idx 

def imerg_bounds(lat_imerg,lon_imerg,lat_geofile,lon_geofile):
    '''
    Funcion para encontrar los limites del dominio de imerg
    en base a los limites del dominio del geofile
    Retorna xlim e ylim
    '''

    xlim = np.empty(2,dtype=int)
    ylim = np.empty(2,dtype=int)
    
    xlim[0] = find_nearest(lon_imerg, float(np.min(lon_geofile)))[1]-4
    xlim[1] = find_nearest(lon_imerg, float(np.max(lon_geofile)))[1]+4
    ylim[0] = find_nearest(lat_imerg, float(np.min(lat_geofile)))[1]-4
    ylim[1] = find_nearest(lat_imerg, float(np.max(lat_geofile)))[1]+4
    
    return xlim,ylim

def crop_imerg(xlim,ylim,lat_imerg,lon_imerg,imerg_precip):
    '''
    Funcion para recortar el dominio de IMERG acorde al tamaño del dominio del geofile
    Regresa un array recortado (lat,lon,precip)
    '''

    xdim = xlim[1]-xlim[0]
    ydim = ylim[1]-ylim[0]
    
    # Latitud
    cropped_lat_imerg = lat_imerg[ylim[0]:ylim[1]]
    # Longitud
    cropped_lon_imerg = lon_imerg[xlim[0]:xlim[1]]
    # Precipitacion
    cropped_precip_imerg = imerg_precip[0,xlim[0]:xlim[1],ylim[0]:ylim[1]]
    
    return cropped_lat_imerg,cropped_lon_imerg,cropped_precip_imerg


def find_nearest_corners(cropped_lat_imerg,cropped_lon_imerg,
                         lat_target:float,lon_target:float):
    '''
    Funcion que busca las coordenadas de los 4 puntos de imerg-f
    mas cercanos al punto target (objetivo)
    
    Entrada: 
        cropped_lat_imerg: latitud imerg-f
        cropped_lon_imerg: longitud imerg-f
        lat_target: latitud objetivo
        lon_target: longitud objetivo
        
    Salida:
        min_lats: indices de las dos latitudes mas cercanas
        min_lons: indices de las dos longitudes mas cercanas
        lat_array[min_lats]: coordenadas de las dos latitudes mas cercanas
        lon_array[min_lons]: coordenadas de las dos longitudes mas cercanas
    '''

    # Paso las latitudes y longitudes de imerg a un formato para trabajar
    lat_array = np.asarray(cropped_lat_imerg)
    lon_array = np.asarray(cropped_lon_imerg)

    # Creo dos arrays vacios de dim 2
    min_lats = np.empty(2,dtype=int)
    min_lons = np.empty(2,dtype=int)

    # Busco y guardo los indices de las latitudes mas cercanas
    min_lats[0] = list(np.abs(lat_array - lat_target)).index(sorted(np.abs(lat_array - lat_target))[:2][0])
    min_lats[1] = list(np.abs(lat_array - lat_target)).index(sorted(np.abs(lat_array - lat_target))[:2][1])

    # Busco y guardo los indices de las longitudes mas cercanas
    min_lons[0] = list(np.abs(lon_array - lon_target)).index(sorted(np.abs(lon_array - lon_target))[:2][0])
    min_lons[1] = list(np.abs(lon_array - lon_target)).index(sorted(np.abs(lon_array - lon_target))[:2][1])
    
    return min_lats, min_lons, lat_array[min_lats], lon_array[min_lons]

def distance(lat0,lon0,lat1,lon1):
    '''
    Funcion que calcula la distancia entre dos puntos con coordenadas
    lat0,lon0 y lat1,lon1
    '''

    dist = np.sqrt( (lat0-lat1)**2 + (lon0-lon1)**2 )
    return dist

#############################################################################################################
#############################################################################################################

df_mean_precip = pd.DataFrame()

for i in range(len(imergFilelist)):
# for i in range(13,27):

    imerg_file = Dataset(imergFilelist[i], mode='r')
    print('Regrillando archivo: ',imergFilelist[i].split('/')[-1])

    try:
        lons_imerg = imerg_file.variables['lon'][:]
        lats_imerg = imerg_file.variables['lat'][:]
        imerg_precip = imerg_file.variables['precipitationCal'][:]
    except:
        raise KeyError('No se encontró alguna de las variables (lat,lon,precip) en el archivo IMERG')


    if np.min(lats_geofile) < np.min(lats_imerg) or \
       np.max(lats_geofile) > np.max(lats_imerg) or \
       np.min(lons_geofile) < np.min(lons_imerg) or \
       np.max(lons_geofile) > np.max(lons_imerg):
        raise RuntimeError('El dominio del archivo IMERG no abarca completamente al dominio del geofile')


    xlim,ylim = imerg_bounds(lats_imerg,lons_imerg,lats_geofile,lons_geofile)
    cropped_lat_imerg,cropped_lon_imerg,cropped_precip_imerg = crop_imerg(xlim,ylim,lats_imerg,
                                                                          lons_imerg,imerg_precip)


    temp_array = np.empty(var.shape)


    for row in range(len(np.asarray(lats_geofile[:,0]))):
        for col in range(len(np.asarray(lons_geofile[0,:]))):

            lat_punto_objetivo = float(lats_geofile[row,col])
            lon_punto_objetivo = float(lons_geofile[row,col])

            if debug:
                print('lat_punto_objetivo: ',lat_punto_objetivo,'lon_punto_objetivo: ',lon_punto_objetivo)

            min_lat_i, min_lon_i, lat_array, lon_array = find_nearest_corners(cropped_lat_imerg,
                                                                              cropped_lon_imerg,
                                                                             lat_punto_objetivo,
                                                                             lon_punto_objetivo)
            if debug:
                print('Lats mas cercanas:',lat_array)
                print('Lons mas cercanas:',lon_array)

            peso1  = 1/distance(lat_punto_objetivo,lon_punto_objetivo,
                                lat_array[0],lon_array[0])
            peso2  = 1/distance(lat_punto_objetivo,lon_punto_objetivo,
                                lat_array[0],lon_array[1])
            peso3  = 1/distance(lat_punto_objetivo,lon_punto_objetivo,
                                lat_array[1],lon_array[0])
            peso4  = 1/distance(lat_punto_objetivo,lon_punto_objetivo,
                                lat_array[1],lon_array[1])
            suma_pesos = peso1 + peso2 + peso3 + peso4


            precip_pto1 = cropped_precip_imerg[min_lon_i[0],min_lat_i[0]]
            precip_pto2 = cropped_precip_imerg[min_lon_i[1],min_lat_i[0]]
            precip_pto3 = cropped_precip_imerg[min_lon_i[0],min_lat_i[1]]
            precip_pto4 = cropped_precip_imerg[min_lon_i[1],min_lat_i[1]]

            if debug:
                print('Precipitacion de los puntos cercanos')
                print('1=',precip_pto1,'2=',precip_pto2,'3=',precip_pto3,'4=',precip_pto4)

            temp_array[row,col] = (precip_pto1*peso1+precip_pto2*peso2+precip_pto3*peso3+precip_pto4*peso4)/suma_pesos

            if debug:
                print('Precipitacion calculada:',temp_array[row,col])
                print('\n')
            
            # Aca termina el regrillado (temp_array).

    #############################################################################################################
    #############################################################################################################
            
    # Suavizo el array con un filtro gausiano
    temp_array_filtered = gaussian_filter(temp_array, sigma=1)
    # Selecciono el string datetime
    datetime_str = imergFilelist[i].split('3IMERG')[-1][1:25]
    # Lo convierto a datetime object
    datetime_object = datetime.strptime(datetime_str[:8]+' '+datetime_str[10:14], '%Y%m%d %H%M')

    unout = 'days since 2000-01-01 00:00:00'

    # using netCDF3 for output format 
    ncout = Dataset(output_dir+'/'+imergFilelist[i].split('3IMERG')[-1][1:9]+datetime_object.strftime("%H%M")+'.PRECIP_FORCING.nc','w','NETCDF4_CLASSIC'); 
    
    ncout.createDimension('Time', 1);
    ncout.createDimension('DateStrLen', 20);
    ncout.createDimension('south_north', temp_array_filtered.shape[0]);
    ncout.createDimension('west_east', temp_array_filtered.shape[1]);

    var_time = ncout.createVariable('Times', 'S1', ('DateStrLen'))
    str_out = netCDF4.stringtochar(np.array([datetime_object.strftime("%Y-%m-%d_%H:%M:%S")], 'S20'))
    var_time[:] = str_out

    rainrate = ncout.createVariable('precip_rate','float32',('Time','south_north','west_east'))
    rainrate.setncattr('units','mm/s'); rainrate[:] = temp_array_filtered/86400;

    #Add global attributes
    rainrate.remap = "regridded via inverse-distance-weighting method"
    today = datetime.today()
    ncout.history = "Created " + today.strftime("%d/%m/%y")
    #Add local attributes to variable instances
    rainrate.units = 'mm s^-1'
    rainrate.description = 'RAINRATE'
    rainrate.long_name = 'RAINRATE'

    ncout.close();
            
    #############################################################################################################
    #############################################################################################################
    if plot_regridded:  

        # Se inicia el plot
        fig_regridded = plt.figure(figsize=[8,8], tight_layout=True)
        # Utilizar la proyección de WRF
        ax = plt.axes(projection=cart_proj)  

        ax.coastlines('10m', linewidth=1)
        ax.add_feature(shape_dp, linewidth=.5, facecolor='none')
        ax.add_feature(shape_pcias, linewidth=1., facecolor='none')

        # d01
        ax.set_xlim([xlim_d01[0]-(xlim_d01[1]-xlim_d01[0])/2, xlim_d01[1]+(xlim_d01[1]-xlim_d01[0])/2])
        ax.set_ylim([ylim_d01[0]-(ylim_d01[1]-ylim_d01[0])/2, ylim_d01[1]+(ylim_d01[1]-ylim_d01[0])/2])

        vmin = 0
        vmax = imerg_precip.max()
        clevs = np.linspace(vmin,vmax,10) #levels / mm
        norm = matplotlib.colors.BoundaryNorm(clevs, len(clevs)) # set boundary of data by normalizing (0,1)

        plt.pcolor(to_np(lons_geofile), to_np(lats_geofile), temp_array_filtered, vmin=vmin, vmax=vmax,
                   transform=crs.PlateCarree(),
                   cmap=precip_colormap)

        ax.add_patch(mpl.patches.Rectangle((xlim_d01[0], ylim_d01[0]), xlim_d01[1]-xlim_d01[0],
                                            ylim_d01[1]-ylim_d01[0],
                                            fill=None, lw=1.5, edgecolor='black', zorder=10))

        cbar = plt.colorbar(ax=ax, shrink=.5)

        plt.title('Regridded IMERG-F '+imergFilelist[i].split('3IMERG')[-1][1:25]+' Z', fontsize=15)
        plt.tight_layout()
        plt.savefig(output_dir+'/Regridded_Plot_'+imergFilelist[i].split('3IMERG')[-1][1:9]+datetime_object.strftime("%H%M")+'.png')

        fig_regridded.clear()
        plt.close()
    #############################################################################################################
    #############################################################################################################
    if mean_precipitation:
        mean = temp_array_filtered.mean()
        df_temp = pd.DataFrame({'mean_rainrate':[mean]}, index=[datetime_object])
        df_mean_precip = pd.concat([df_mean_precip,df_temp])
    #############################################################################################################
    #############################################################################################################

    if plot_imerg:

        fig_imerg = plt.figure(figsize=(10, 10))

        units = imerg_file.variables['precipitationCal'].units

        lats_red = [ np.min(lats_geofile).values, np.max(lats_geofile).values, np.max(lats_geofile).values, np.min(lats_geofile).values ]
        lons_red = [ np.min(lons_geofile).values, np.min(lons_geofile).values, np.max(lons_geofile).values, np.max(lons_geofile).values ]

        # Get some parameters for the Stereographic Projection
        lon_0 = lons_imerg.mean()
        lat_0 = lats_imerg.mean()

        m = Basemap(llcrnrlon=np.min(to_np(lons_geofile))-4.5,
                    llcrnrlat=np.min(to_np(lats_geofile))-5.5,
                    urcrnrlon=np.max(to_np(lons_geofile))+4.5,
                    urcrnrlat=np.max(to_np(lats_geofile))+3,
                    projection='lcc', 
                    lat_0=np.mean(to_np(lats_geofile)),
                    lon_0=np.mean(to_np(lons_geofile)),
                    resolution ='l',area_thresh=1000.)

        # Because our lon and lat variables are 1D,
        # use meshgrid to create 2D arrays
        # Not necessary if coordinates are already in 2D arrays.
        lon, lat = np.meshgrid(lons_imerg, lats_imerg)
        xi, yi = m(lon, lat)

        # Plot Data
        vmin = 0
        vmax = imerg_precip.max()
        clevs = np.linspace(vmin,vmax,10) #levels / mm

        cs = m.pcolor(xi, yi, np.transpose(np.squeeze(imerg_precip)),
                      cmap=precip_colormap, vmin=vmin, vmax=vmax)

        # Add Grid Lines
        m.drawstates()
        m.drawcountries()

        draw_screen_poly( lats_red, lons_red, m, 'black' )

        # Add Colorbar
        cbar = m.colorbar(cs, location='bottom')
        cbar.set_label(units)

        # Add Title
        plt.title('IMERG-F '+imergFilelist[i].split('3IMERG')[-1][1:25]+' Z', fontsize=15)

        plt.tight_layout()
        plt.savefig(output_dir+'/IMERG_Plot_'+imergFilelist[i].split('3IMERG')[-1][1:9]+datetime_object.strftime("%H%M")+'.png')

        fig_imerg.clear()
        plt.close()

    print('Regrillado del archivo exitoso')
    imerg_file.close()

if mean_precipitation:
    df_mean_precip.to_csv(output_dir+'/mean_rainrate.csv')
print('Proceso de regrillado finalizado')
