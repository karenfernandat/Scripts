# Meteorologista Karen Teixeira
# karenbrazteixeira@gmail.com
# Script para calcular a velocidade do vento mensal a partir das componentes u e v
# ERA5 2003 a 2020
# Grade do BRAZIL  ax.set_extent([-74, -34.7, 5.30, -32.3])
# Grade do SUDESTE -56.0, -37, -13, -25.5

# Packages
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Read data 
u10m = xr.open_dataset('~/ERA5_u10m_mean_mon_2003_2020.nc')
v10m = xr.open_dataset('~/ERA5_v10m_mean_mon_2003_2020.nc')

# Monthly mean
med_u10m = u10m.groupby('time.month').mean()
med_v10m = v10m.groupby('time.month').mean()

# Loop para calcular a média mensal da velocidade do vento a partir das componentes U e V
variables = {'janeiro': (med_u10m.u10[0]**2 + med_v10m.v10[0]**2)**(1/2),
             'fevereiro': (med_u10m.u10[1]**2 + med_v10m.v10[1]**2)**(1/2),
             'março': (med_u10m.u10[2]**2 + med_v10m.v10[2]**2)**(1/2),
             'abril': (med_u10m.u10[3]**2 + med_v10m.v10[3]**2)**(1/2),
             'maio': (med_u10m.u10[4]**2 + med_v10m.v10[4]**2)**(1/2),
             'junho': (med_u10m.u10[5]**2 + med_v10m.v10[5]**2)**(1/2),
             'julho': (med_u10m.u10[6]**2 + med_v10m.v10[6]**2)**(1/2),
             'agosto': (med_u10m.u10[7]**2 + med_v10m.v10[7]**2)**(1/2),
             'setembro': (med_u10m.u10[8]**2 + med_v10m.v10[8]**2)**(1/2),
             'outubro': (med_u10m.u10[9]**2 + med_v10m.v10[9]**2)**(1/2),
             'novembro': (med_u10m.u10[10]**2 + med_v10m.v10[10]**2)**(1/2),
             'dezembro': (med_u10m.u10[11]**2 + med_v10m.v10[11]**2)**(1/2)}

for season, variable in variables.items():

    # Plot velocidade do média do vento
    fig = plt.figure(dpi=180)
    ax = fig.add_subplot(projection=ccrs.PlateCarree())
    ax.set_extent([-56.0, -37, -13, -25.5])
    ax.add_feature(cfeature.LAND, color='white', facecolor='0.9')
    ax.add_feature(cfeature.LAKES, alpha=0.2)
    ax.add_feature(cfeature.BORDERS)
    ax.coastlines()
    states_provinces = cfeature.NaturalEarthFeature(category='cultural',
                                                    name='admin_1_states_provinces_lines',
                                                    scale='50m',
                                                    facecolor='none')
    ax.add_feature(states_provinces, edgecolor='k', linewidth=0.5)
    grid = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5,
                        color='gray', alpha=0.5, linestyle='--')
    grid.xlabels_top = False
    grid.ylabels_right = False
    grid.xlines = False
    grid.ylines = False
    grid.xformatter = LONGITUDE_FORMATTER
    grid.yformatter = LATITUDE_FORMATTER
    grid.xlabel_style = {'size': 14, 'color': 'black'}
    grid.ylabel_style = {'size': 14, 'color': 'black'}
    var = ax.contourf(variable.longitude[:], variable.latitude[:], variable[:], cmap='turbo',
    levels=np.arange(1.0, 10.0, 1.0), extend='both', transform=ccrs.PlateCarree())
    cb_var = fig.colorbar(var, ax=ax, shrink=0.9, aspect=12)
    cb_var.ax.tick_params(labelsize=14)
    ax.set_title(f'{season.upper()} - MENSAL', fontdict={'family': 'serif', 'weight': 'normal',
                                    'color': 'k', 'size': 14})
    plt.savefig(f'~/Imagens/vento_sudeste_{season}.png')
    plt.close()
