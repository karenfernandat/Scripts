# Meteorologista Karen Teixeira
# karenbrazteixeira@gmail.com
# Script para calcular a climatologia da temperatura mensal para região SUDESTE
# Dados do CAMS e ERA5 de 2003 a 2020.
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
ds = xr.open_dataset('~/dados/CAMS_temp_mean_mon_2003_2020.nc')

# Agrupar mensalmente os dados
med = ds.groupby('time.month').mean()

# Loop para calcular a média mensal da temperatura na região Sudeste
variables = {'janeiro': med.t[0] - 273.15, 'fevereiro': med.t[1] - 273.15, 'março': med.t[2] - 273.15,
              'abril': med.t[3] - 273.15, 'maio': med.t[4] - 273.15, 'junho': med.t[5] - 273.15,
              'julho': med.t[6] - 273.15, 'agosto': med.t[7] - 273.15, 'setembro': med.t[8] - 273.15,
              'outubro': med.t[9] - 273.15, 'novembro': med.t[10] - 273.15, 'dezembro': med.t[11] - 273.15}
for season, variable in variables.items():
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
    grid.top_labels = False
    grid.right_labels = False
    grid.xlines = False
    grid.ylines = False
    grid.xformatter = LONGITUDE_FORMATTER
    grid.yformatter = LATITUDE_FORMATTER
    grid.xlabel_style = {'size': 14, 'color': 'black'}
    grid.ylabel_style = {'size': 14, 'color': 'black'}
    var = ax.contourf(variable.longitude[:], variable.latitude[:], variable[:], cmap='RdBu_r',
    levels=np.arange(15, 30, 1.0), extend='both', transform=ccrs.PlateCarree())
    cb_var = fig.colorbar(var, ax=ax, shrink=0.9, aspect=12)
    cb_var.ax.tick_params(labelsize=14)
    ax.set_title(f'TEMPERATURA - {season.upper()}', fontdict={'family': 'serif', 'weight': 'normal',
                                    'color': 'k', 'size': 14})
    plt.savefig(f'~/Imagens/temp_mensal_sudeste_{season}.png')
    plt.close()