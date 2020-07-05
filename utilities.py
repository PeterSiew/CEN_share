from scipy.stats.stats import pearsonr
import matplotlib.path as mpath
import numpy as np
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import pickle as pickle
from mpl_toolkits.basemap import Basemap
import datetime as dt
from datetime import datetime
from scipy import stats
import scipy
import xarray as xr
import ipdb
import timeseries_treatment as tt
import pandas as pd
import cartopy.crs as ccrs

def read_timeseries(timescales):

	var_f = {}
	if timescales in ['monthly', 'submonthly']: # monthly_timescale
		ts_folder = '/home/ysi082/codes/cen_1/ts_anomaly/'
		var_f['U10'] = ts_folder + 'U10_ts.p'
		var_f['SPV'] = ts_folder + 'SPV_ts.p'
		var_f['ICE'] = ts_folder + 'ICE_ts.p'
		var_f['THF'] = ts_folder + 'THF_ts.p' # Original THF\
		var_f['THFnew'] = ts_folder + 'THFnew_ts.p' # Original THF\
		var_f['URALS'] = ts_folder + 'URALS_ts.p'
		var_f['V*T*'] = ts_folder + 'VTstar_ts.p'
		var_f['VTstar'] = ts_folder + 'VTstar_ts.p'
		var_f['VTstardaily'] = ts_folder + 'VTstar_ts.p'
		var_f['VTstar_Urals'] = ts_folder + 'VTstar_Urals_ts.p'
		var_f['VTstar_0to180'] = ts_folder + 'VTstar_0to180_ts.p'
		var_f['VTstar_60to190'] = ts_folder + 'VTstar_60to190_ts.p'
		var_f['V*T*_half'] = ts_folder + 'VTstar0to180_ts.p'
		var_f['V*T*1'] = ts_folder + 'VTstar_wave1_monthly_ts.p'
		var_f['V*T*2'] = ts_folder + 'VTstar_wave2_monthly_ts.p'
		var_f['V*T*123'] = ts_folder + 'VTstar_wave123_monthly_ts.p'
		var_f['Z10'] = ts_folder + 'Z10_ts.p'
		var_f['Z50'] = ts_folder + 'Z50_ts.p'
		var_f['Z500'] = ts_folder + 'Z500_ts.p'
		var_f['IR'] =  ts_folder + 'IR_ts.p'
		var_f['QBO_U10_raw'] = ts_folder + 'QBO_U10_raw_ts.p'
		var_f['Nino34'] = ts_folder + 'nino34_ts.p'
		#var_f['NAO'] = ts_folder + 'NAO_ts.p'
		var_f['NAO_NOAA'] = ts_folder + 'NAO_ts.p'
		var_f['NAO_EOFbased'] = ts_folder + 'NAO_EOFbased_ts.p' # Only monthly timeseries
		var_f['NAO_SLPbased'] = ts_folder + 'NAO_SLPbased_ts.p' # Only monthly timeseries

		# Read to delete
		var_f['Z10e1'] = './ts_anomaly/Z10eof1_ts.p'
		var_f['Z10e2'] = './ts_anomaly/Z10eof2_ts.p'
		var_f['Z10e3'] = './ts_anomaly/Z10eof3_ts.p'
		var_f['SLPA'] = './ts_anomaly/SLPA_ts.p'
		var_f['SLPL'] = './ts_anomaly/SLPL_ts.p'
		var_f['VT0to60'] = './ts_anomaly/VTstar0to60_ts.p'
		var_f['ICE_new'] = './ts_anomaly/ICE_new_ts.p'
		var_f['BKS_ice_anom'] = '/home/ysi082/codes/cen_1/ts_anomaly_new/lingling/BKS_ice_ts_anom.p'
		var_f['BKS_ice_raw'] = '/home/ysi082/codes/cen_1/ts_anomaly_new/lingling/BKS_ice_ts_raw.p'
		var_f['Arctic_ice_raw'] = '/home/ysi082/codes/cen_1/ts_anomaly_new/lingling/panArctic_ice_ts_raw.p'
		var_f['Arctic_ice_anom'] = '/home/ysi082/codes/cen_1/ts_anomaly_new/lingling/panArctic_ice_ts_anom.p'

	return var_f

	
######### Read files

def read_xarray_cache_new(vars, vartype='anom', slicing=False,
						resolution='monthly', cache_folder = '/Data/skd/scratch/ysi082/EraIn_daily/', unit=True, mon=False, year_control=True, dt_start=None, dt_end=None):
	# Note that cache might have different length in time. e.g., CI is just up to 31 Dec 2017

	vards = xr.Dataset()
	for var in vars:
		if vartype == 'raw':
			var3d = xr.open_dataset(cache_folder+'%s_raw.nc'%var, chunks={'time':500})
		elif vartype == 'anom':
			var3d = xr.open_dataset(cache_folder+'%s_anom.nc'%var, chunks={'time':500})
		var_key = list(var3d.keys())[0]
		var3d = var3d[var_key]
		if var == 'ICE':
			var3d = var3d.where(var3d != 0)
		if var  == 'THF':
			var3d = var3d * -1.0
		# By this way, it will use its longest time length if they have different time
		vards[var] = var3d

	if resolution == 'daily':
		pass
	elif resolution == 'monthly':
		# To monthly_data
		vards = vards.resample(time = '1MS').mean()
	mon_abbre_dict = mon_num_abbre()[0]
	if (resolution == 'monthly') & (mon in mon_abbre_dict):
		mask = (vards.time.dt.month == mon_abbre_dict[mon])
		vards = vards.where(mask, drop=True)
	if slicing: # It seems casuing some unknow error
		lons = vards.longitude
		lats = vards.latitude
		vards = vards.sel(longitude=vards.longitude[::3], latitude=vards.latitude[::3])
	if unit:
		vards = standard_unit(vards)
	if year_control: #Control the date. The ful llength is from 1979 to 2018
		vards = vards.sel(time=slice('1979-01-01', '2016-12-31'))
	else: 
		vards = vards.sel(time=slice(dt_start, dt_end))

	return vards

def read_ts_cache(vars, resolution = 'monthly', year_control=True ,dt_start=None, dt_end=None, mon=False, unit=True):
	# Report
	# in vars only contain NAO, it will return the 29 Feb since NAO file has the 29Feb
	# if vars ['ICE', 'NAO'], it follow the shape of ICE. If vars is ['NAO', 'SLP'], it follows the shape of NAO(29 Feb)

	# Read the timeseries 
	ts_paths = read_timeseries('monthly')
	ts_daily, dt_daily = {}, {}
	for var in vars:
		ts_daily[var], dt_daily[var] = ts_extraction(ts_paths[var])
	
	# Xarray will automatically fit the data
	# Put the timeseries into Xarray Dataset
	ts_ds = xr.Dataset()
	for var in vars:
		dt64 = pd.to_datetime(dt_daily[var])
		ts_ds[var] = xr.DataArray(ts_daily[var], dims='time', coords = {'time': dt64})

	if resolution == 'monthly':
		# Read the daily into monthly data
		monthly_ts_ds = ts_ds.resample(time = '1MS').mean()
		print(('You have read the %s monthly anomaly timeseries' %vars))
		ts_ds = monthly_ts_ds
	
	elif resolution == 'halfmonthly':
		#halfmonthly_ts_ds = ts_ds.resample(time = '1SMS').mean()
		# Change it to pandas
		ts_ds_pd = ts_ds.to_dataframe()
		# Resample, but starting on 16
		pd_offset = pd.tseries.offsets.SemiMonthBegin(n=1, normalize=False, day_of_month=16)
		halfmonthly_ts_ds = ts_ds_pd.resample(pd_offset).mean().to_xarray()
		ts_ds = halfmonthly_ts_ds
	elif resolution=='pentad':
		# 1: 1st-5st, 2: 6th-10th, 3:11th-15th, 4:16th-20th, 5:21th-25th, 6:26th-28/30/31th
		# Loop the datetime_index. Extract 1-5, 6-10, 11-15...
		# Group them and then mean
		pentad_time = []
		for i in ts_ds.time.to_index():
			if i.day in [1,2,3,4,5]:
				day=1
			elif i.day in [6,7,8,9,10]:
				day=6
			elif i.day in [11,12,13,14,15]:
				day=11
			elif i.day in [16,17,18,19,20]:
				day=16
			elif i.day in [21,22,23,24,25]:
				day=21
			elif i.day in [26,27,28,29,30,31]:
				day=26
			dt_day = dt.date(i.year, i.month, day)
			pentad_time.append(dt_day)
		pentad_time = pd.to_datetime(pentad_time)
		# Assign the pentad_group to be the time corrindate
		ts_ds = ts_ds.assign_coords(time=pentad_time)
		# Group by the pentad group. And then then mean
		ts_ds = ts_ds.groupby('time').mean(dim='time')

	elif resolution == 'daily':
		print(('You have read the %s daily anomaly timeseries' %vars))
		pass # Don't need to do any change

	else:
		print('Wrong resolution')

	mon_abbre_dict = mon_num_abbre()[0]
	if (resolution == 'monthly') & (mon in mon_abbre_dict):
		mask = (monthly_ts_ds.time.dt.month == mon_abbre_dict[mon])
		ts_ds = ts_ds.where(mask, drop=True)
	if unit:
		ts_ds = standard_unit(ts_ds)
	if 'SPV' in ts_ds:
		ts_ds['SPV'] = ts_ds['SPV'] * -1.0
	# Should do the selction after aggrgation
	if year_control==True:
		ts_ds = ts_ds.sel(time=slice('1979-01-01', '2016-12-31'))
	else: 
		# To select all
		ts_ds = ts_ds.sel(time=slice(dt_start, dt_end))

	return ts_ds

def standard_unit(vards):
	# for only ECMWF Era-In data

	vars = list(vards.keys())
	for var in vars:
		if var in ['ICE']:
			vards[var] = vards[var] * 100
		elif var in ['Z10', 'Z50', 'Z100', 'Z200', 'Z300', 'Z300star', 'Z300starwave1', 'Z300starwave2', 'Z300starwave3', 'Z500', 'SPV', 'All_Z']:
			vards[var] = vards[var] / 9.80665
		#elif var in ['SLP', 'URALS']:
		#	vards[var] = vards[var] / 100.0

	return vards


def ts_extraction(f_path, numpy_file=False):
	# f_path is the file path of the numpy file
	if numpy_file:
		uts, anom_ts = np.loadtxt(f_path)
		daily_dt = [datetime.fromtimestamp(i) for i in uts]

	else: # new pickle way than numpy. Numpy way is above 
		print('Loading pickle timeseries')
		anom_ts, daily_dt = pickle.load(open(f_path, 'rb'), encoding='bytes')

	return anom_ts, daily_dt

def mon_num_abbre():

	mon_abbre_dict = {'May':5, 'Jun':6, 'Jul':7, 'Aug':8, 'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12, 'Jan': 1, 'Feb':2, 'Mar':3}	
	mon_num_dict= {5:'May', 6:'Jun', 7:'Jul', 8:'Aug', 9:'Sep', 10:'Oct', 11:'Nov', 12:'Dec', 1: 'Jan',  2:'Feb', 3:'Mar'}

	return mon_abbre_dict, mon_num_dict


##################### Plotting

def grid_plotting_height_lat(map_grids, rows, cols, mapcolors, clevels_rows, clabels_rows, top_titles=[], left_titles=[], pval_map=None, sig=0.05, contour_w_shading=False, ind_titles=None):

	plt.close()
	fig = plt.figure(figsize=(cols*2.5, rows*2.5))
	
	for i in range(rows):
		for j in range(cols):
			# If the numbers inside the map_grids is less than rows*cols
			if len(map_grids) <= i*cols+j:
				continue
			ax = fig.add_subplot(rows, cols, i*cols+j+1)
			lats = map_grids[i*cols+j].latitude.values
			heights = map_grids[i*cols+j].lev_hPa.values
			log_heights = np.log(heights)
			cs = ax.contourf(lats, log_heights, map_grids[i*cols+j], clevels_rows[i], cmap=mapcolors, extend='both')
			if True:
				idx_0 = (len(clevels_rows[i]) -1) / 2.0
				csf = ax.contour(lats, log_heights, map_grids[i*cols+j], np.delete(clevels_rows[i], [idx_0]), colors='k', linewidths=0.5)
			if pval_map!=None:
				ax.contourf(lats, log_heights, pval_map[i*cols+j], [0, sig, 1000], hatches=['..', None], colors='none', extend='neither')
			if i==0:
				ax.set_title(top_titles[j], size=30)
			if j==0:
				ax.text(-0.65, 0.35, '%s'%left_titles[i], rotation='vertical', transform=ax.transAxes, size=30)

			if j==(cols-1): # the last column
				ax_pos = ax.get_position()
				cax = fig.add_axes([ax_pos.x1 + 0.1*(ax_pos.x1-ax_pos.x0), ax_pos.y0, 0.015, ax_pos.y1-ax_pos.y0]) 
				cbar = fig.colorbar(cs, cax=cax, orientation='vertical')
				cbar.set_label(clabels_rows[i])
			height_yticks = [10,50,200,500,1000]
			ax.set_yticks(np.log(height_yticks))
			ax.set_yticklabels(height_yticks)
			ax.set_ylim(np.log(height_yticks[0]), np.log(height_yticks[-1]))
			ax.invert_yaxis()
			ax.set_xticks([10,45,80])
			if ind_titles is not None:
				ax.set_title(ind_titles[i*cols+j])

	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.3)

def map_grid_plotting(map_grids, rows, cols, mapcolors, clevels_rows, clabels_rows, top_titles=[], left_titles=[], region_boxes=None, pval_map=None, sig=0.05, low_lat=45, contour_w_shading=False, ind_titles=None, projection=ccrs.NorthPolarStereo(), pval_hatch_small=True, colorbar=True, pltf=None, ax_all=None, map_extents=None, xsize=2.5, ysize=2.5, transpose=False, subplot_adjust=False, contour_map_grids=None, contour_clevels=None, gridline=True, leftcorner_text=None, circle=False):

	gidx= np.empty((rows,cols)).astype('int')
	for i in range(rows):
		for j in range(cols):
			gidx[i,j] = i*cols+j
	if transpose:
		rows, cols = cols, rows
		gidx = gidx.T
		left_titles, top_titles = top_titles, left_titles

	# use the Catropy
	plt.close() if pltf is None else ''
	fig = plt.figure(figsize=(cols*xsize, rows*ysize)) if pltf is None else pltf
	map_extents = {'left':-180, 'right':180, 'bottom':low_lat, 'top':90} if map_extents is None else map_extents
	k=0
	for i in range(rows):
		for j in range(cols):
			if len(map_grids) <= gidx[i,j]: # If the numbers inside the map_grids is less than rows*cols
				continue
			ax = fig.add_subplot(rows, cols, k+1, projection=projection) if ax_all is None else ax_all[gidx[i,j]]
			lons = map_grids[gidx[i,j]].longitude.values
			lats = map_grids[gidx[i,j]].latitude.values
			cs = ax.contourf(lons, lats, map_grids[gidx[i,j]], clevels_rows[gidx[i,j]], cmap=mapcolors, transform=ccrs.PlateCarree(), extend='both')
			ax.coastlines(color='darkgray')
			ax.set_extent([map_extents['left'], map_extents['right'], map_extents['bottom'], map_extents['top']], ccrs.PlateCarree())
			if circle:
				theta = np.linspace(0, 2*np.pi, 100); center, radius = [0.5, 0.5], 0.5; verts = np.vstack([np.sin(theta), np.cos(theta)]).T
				circle = mpath.Path(verts * radius + center)
				ax.set_boundary(circle, transform=ax.transAxes)
			if contour_w_shading:
				if clevels_rows[gidx[i,j]][0]*-1 == clevels_rows[gidx[i,j]][-1]:
					idx_0 = (len(clevels_rows[gidx[i,j]]) -1) / 2.0
					clevels_new = np.delete(clevels_rows[gidx[i,j]], idx_0)
					csf = ax.contour(lons, lats, map_grids[gidx[i,j]], clevels_new, colors='k', linewidths=0.5, transform=ccrs.PlateCarree())
			if not contour_map_grids is None:
				if not contour_map_grids[gidx[i,j]] is None:
					idx_0 = (len(contour_clevels[gidx[i,j]]) -1) / 2.0
					clevels_new = np.delete(contour_clevels[gidx[i,j]], idx_0)
					csf = ax.contour(lons, lats, contour_map_grids[gidx[i,j]], clevels_new, colors='k', linewidths=1, transform=ccrs.PlateCarree())
			if not region_boxes is None:
				lons_l, lats_l = region_boxes[gidx[i,j]][0], region_boxes[gidx[i,j]][1]
				ax.plot(lons_l, lats_l, transform=ccrs.PlateCarree(), color='k', linewidth=2)
			if not pval_map is None:
				if pval_hatch_small[gidx[i,j]]: # Especially for the model plotting
					ax.contourf(lons, lats, pval_map[gidx[i,j]], [0, sig, 100], hatches=['.', None], colors='none', extend='neither', transform=ccrs.PlateCarree())
				else:
					#ax.contourf(lons, lats, pval_map[gidx[i,j]], [0, 0.80, 0.90, 1], hatches=[None, 'x', '.'], colors='none', extend='neither', transform=ccrs.PlateCarree())
					ax.contourf(lons, lats, pval_map[gidx[i,j]], [0, 0.30, 1], hatches=[None, 'x'], colors='none', extend='neither', transform=ccrs.PlateCarree())
			if not leftcorner_text is None:
				if not leftcorner_text[gidx[i,j]] is None:
					ax.annotate(leftcorner_text[gidx[i,j]], xy=(0.01, 0.99), xycoords='axes fraction', fontsize=16, verticalalignment='top')
			if gridline:
				g1 = ax.gridlines(draw_labels=False, linewidth=1, color='lightgray', linestyle='-', alpha=0.5)
				g1.xlocator = matplotlib.ticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
				g1.ylocator = matplotlib.ticker.FixedLocator([0, 45, 90])
			if i==0:
				ax.set_title(top_titles[j], size=20)
			if j==0:
				ax.text(-0.30, 0.35, '%s'%left_titles[i], rotation='vertical', transform=ax.transAxes, size=20)
			if colorbar:
				ax_pos = ax.get_position()
				if not transpose:
					if j==(cols-1): # the last column
						cax = fig.add_axes([ax_pos.x1 + 0.1*(ax_pos.x1-ax_pos.x0), ax_pos.y0+(ax_pos.y1-ax_pos.y0)*0.1, 0.015, (ax_pos.y1-ax_pos.y0)-(ax_pos.y1-ax_pos.y0)*0.1]) 
						#cax = fig.add_axes([ax_pos.x1 + 0.1*(ax_pos.x1-ax_pos.x0), ax_pos.y0 , 0.015, ax_pos.y1-ax_pos.y0]) 
						cbar = fig.colorbar(cs, cax=cax, orientation='vertical')
						cbar.set_label(clabels_rows[i])
				else:
					if i==(rows-1): # the last column
						cax = fig.add_axes([ax_pos.x0, ax_pos.y0-0.2*(ax_pos.y1-ax_pos.y0), ax_pos.x1-ax_pos.x0, 0.015])
						cbar = fig.colorbar(cs, cax=cax, orientation='horizontal')
						cbar.set_label(clabels_rows[j])
						cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
			if ind_titles is not None:
				ax.set_title(ind_titles[gidx[i,j]])
			k=k+1
	if subplot_adjust:
		plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.02, hspace=0.02)

def quick_map_plotting(map_2d, lats, lons, fn_name, clevels=10, **aargus):

	plt.close()
	fig, ax1 = plt.subplots(1, 1, figsize=(10,10))
	#clevels = 10
	cmap = 'bwr'
	m = Basemap(projection='npaeqd', boundinglat=25, lon_0=0, resolution='l', ax=ax1)
	m.drawcoastlines(color='lightgray')
	x, y = m(*np.meshgrid(lons,lats))
	cs = m.contourf(x, y, map_2d , clevels, cmap=cmap)
	
	from mpl_toolkits.axes_grid1 import make_axes_locatable
	divider = make_axes_locatable(ax1)
	cax = divider.append_axes("bottom", size="3%", pad=0.1)
	cbar = fig.colorbar(cs, cax=cax, orientation='horizontal')

	if 'box' in aargus: # plotting a box
		lat1, lat2, lon1, lon2 = aargus['box']
		box_lon, box_lat = [lon1, lon2, lon2, lon1, lon1], [lat2, lat2, lat1, lat1, lat2]
		box_lon=np.array([np.linspace(j, box_lon[i+1], 30) for i, j  in enumerate(box_lon[0:-1])]).reshape(-1)
		box_lat=np.array([np.linspace(j, box_lat[i+1], 30) for i, j  in enumerate(box_lat[0:-1])]).reshape(-1)
		m.plot(box_lon, box_lat, latlon=True, color='brown', linewidth=1)

	plt.savefig('../graphs/%s'%fn_name, bbox_inches='tight')

def map_plotting(map_2d, titles, lats, lons, clevels, mapcolors, clabel, ax, sig_lev=0.05, cb=True, low_lat_bound=-10, **aargus):
	
	# a more detailed version of quick_map_plotting
	import matplotlib

	cmap = matplotlib.colors.ListedColormap(mapcolors)
	if 'proj' in aargus:
		m = Basemap(projection=aargus['proj'],llcrnrlat=-15,urcrnrlat=90,\
	            llcrnrlon=-180,urcrnrlon=180,resolution='c', ax=ax)
	else:
		m = Basemap(projection='npaeqd', boundinglat=low_lat_bound, lon_0=0, resolution='l', round = True, ax=ax)
		#m = Basemap(projection='npstere', boundinglat=0, lon_0=0, resolution='l', round =False, ax=ax)
		#m.fillcontinents(color='gainsboro', lake_color='gainsboro', ax=None, zorder=None, alpha=None)
	x, y = m(*np.meshgrid(lons,lats))
	cs = m.contourf(x, y, map_2d , clevels, cmap=cmap, extend = 'both')
	if 'pval_map' in aargus:
		m.contourf(x, y, aargus['pval_map'], [0, sig_lev], hatches=['..', None], colors='none', extend='neither')
		#m.contourf(x, y, aargus['pval_map'], [0, sig_lev], hatches=[None, '..'], colors='none', extend='neither')
		
	if 'contour_map' in aargus:
		 css = m.contour(x, y, aargus['contour_map'], 6, colors='k') 
		 plt.clabel(css, fontsize=9, inline=1)
	m.drawcoastlines(color='gray')
	m.drawparallels(np.arange(-90,91,45), color='gray')
	m.drawmeridians(np.arange(0,360,60), color='gray')
	ax.set_title(titles)

	if 'box' in aargus: # plotting a box
		lat1, lat2, lon1, lon2 = aargus['box']
		box_lon, box_lat = [lon1, lon2, lon2, lon1, lon1], [lat2, lat2, lat1, lat1, lat2]
		box_lon=np.array([np.linspace(j, box_lon[i+1], 30) for i, j  in enumerate(box_lon[0:-1])]).reshape(-1)
		box_lat=np.array([np.linspace(j, box_lat[i+1], 30) for i, j  in enumerate(box_lat[0:-1])]).reshape(-1)
		m.plot(box_lon, box_lat, latlon=True, color='brown', linewidth=1)

	if cb:
		from mpl_toolkits.axes_grid1 import make_axes_locatable
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("bottom", size="3%", pad=0.1)
		cbar = plt.colorbar(cs, cax=cax, orientation='horizontal')
		cbar.set_label(clabel)
		cbar.ax.set_xticklabels(cbar.ax.get_xticklabels(), rotation='vertical')
		cbar.ax.tick_params(labelsize=6)
	
	return cs

def map_fill_white_gap(map2d, lons):

	map2d_filled = np.column_stack((map2d, map2d[: , -1]))
	lon_diff = lons[1] - lons[0]
	lons_filled = np.insert(lons, lons.size, lons[-1] + lon_diff)
	
	return map2d_filled, lons_filled
	
def scatter_plot(x, y, xlabel, ylabel, ax, clabel='', mcolor=False, color_normalize=True, regress_line=False, anom=False, texts=None, legend_no = '', tsize=12, single_color='k', markersize=10, s_alpha=1, xlim=None, ylim=None, text_xadjust=0):

	if isinstance(mcolor, np.ndarray): # Colors are some values. Plot the color bar
		mapcolors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
		cmap= matplotlib.colors.ListedColormap(mapcolors)
		if color_normalize:
			mcolor = (mcolor - mcolor.mean()) / mcolor.std()
		cs = ax.scatter(x, y, c=mcolor, s=10, vmin=-2, vmax=2, cmap=cmap)
		from mpl_toolkits.axes_grid1 import make_axes_locatable
		divider = make_axes_locatable(ax)
		cax = divider.append_axes("right", size="3%", pad=0.1)
		cbar = plt.colorbar(cs, cax=cax, orientation='vertical', extend='both')
		cbar.set_label(clabel, size=tsize)

	elif not mcolor: 
		ax.scatter(x, y, s=markersize, c=single_color, alpha=s_alpha) 

	if not texts is None:
		for i in range(len(x)):
			ax.annotate(texts[i], (x[i]+text_xadjust, y[i]), fontsize=8)

	if regress_line:
		m, __, __, __, __ = stats.linregress(x,y)
		slope = round(m, 3)
		fit = np.polyfit(x, y ,1)
		fit_fn = np.poly1d(fit) 
		corr = round(scipy.stats.pearsonr(x, y)[0], 3)
		xx = [np.min(x), np.max(x)]
		ax.plot(xx, fit_fn(xx), '-.', color=single_color, label=legend_no + 'corr: ' + str(corr) + '/ slope: '+str(slope))
		ax.legend()
	
	if anom:
		ax.axvline(x=0, color='lightgray', linestyle='--', linewidth=1)
		ax.axhline(y=0, color='lightgray', linestyle='--', linewidth=1)
	
	if xlim != None:
		ax.set_xlim(xlim[0],xlim[1])
	if ylim != None:
		ax.set_ylim(ylim[0],ylim[1])
	ax.set_xlabel(xlabel, size=tsize)
	ax.set_ylabel(ylabel, size=tsize)

def timeseries_plotting(times, timeseries, legend_label, ylabel, twin_axis=False, vertical_lines=[None], ylims=None, fname='', pltf=None, ax=None, grid=False, xsize=8, ysize=3, xticklabels=None):

	fig = plt.figure(figsize=(xsize, ysize)) if pltf is None else pltf
	ax = fig.add_subplot(111) if ax is None else ax
	# Plotting the timeseries
	for i, ts in enumerate(timeseries):
		ax.plot(times[i], ts, marker='o', markersize=5, lw=1, label=legend_label[i])
	if not ylims is None:
		ax.set_ylim(ylims[0], ylims[1])
	if twin_axis==True:
		axn = ax[1].twinx()
		axn.plot(time, ts_twin, color='red', lw=0.3)
		axn.set_ylabel(ylabel_twin, color='red')
	if not vertical_lines[0] is None:
		for vl in vertical_lines:
			ax.axvline(x=vl, color='k', ls='--', lw=0.2)
	if not legend_label[0] is None:
		ax.legend(loc='upper right')
	ax.grid() if grid else ''
	ax.axhline(y=0, color='gray', ls='--', lw=1)
	ax.set_ylabel(ylabel)
	ax.set_xticklabels(xticklabels, rotation='vertical') if not xticklabels is None else ""

	plt.savefig('../graphs/%s_%s_ts.png' %(dt.date.today(), fname), bbox_inches='tight') if pltf is None else ""

def regression_two_fields_not_used(var3d, ts, corr=True):
	# Regress between a 3d-field and one timeseries
	# please detrend the data grid by grid before doing the regression
	# var3d and ts should have the same dimensions
	
	# Standardize the ts
	ts_std = (ts - ts.mean()) / ts.std()
	
	# Doing the regession starting from here
	var3d_timedim = var3d.shape[0]
	var3d_2d = var3d.reshape(var3d_timedim, -1) # have a dimension of time x grids
	grid_dim = var3d_2d.shape[1]
	
	# create two empty arrays to store regression values
	regress = np.ones(grid_dim); p_vals = np.ones(grid_dim)
	for grid in range(0, grid_dim, 1):
		# standerdize the data
		grid_ts = var3d_2d[:, grid]
		if corr:
			grid_ts_std = (grid_ts - grid_ts.mean()) / grid_ts.std()
			regress[grid], _, _, p_vals[grid], _ = stats.linregress(ts_std, grid_ts_std)
		else:
			regress[grid], _, _, p_vals[grid], _ = stats.linregress(ts, grid_ts)
			#m, c = np.linalg.lstsq(grid_ts_std, ts_std, rcond=None)[0]
		
	regress = regress.reshape(var3d.shape[1], var3d.shape[2])
	p_vals = p_vals.reshape(var3d.shape[1], var3d.shape[2])

	return regress, p_vals		


def defined_regions():
	# All data in ECMWF go from -180 to +180 longitude right now
	regions = {}
	regions['ICE'] = (80, 70, 30, 105) # lat1, lat2, lon1, lon2
	regions['THF'] = (80, 70, 30, 105) # lat1, lat2, lon1, lon2
	regions['URALS'] = (70, 45, 40, 85) # lat1, lat2, lon1, lon2
	regions['SLP'] = (70, 45, 40, 85) # lat1, lat2, lon1, lon2
	regions['VTstar'] = (75, 45, -180, 180) # Fit the era-in data in skdcyclone 
	regions['V*T*'] = (75, 45, -180, 180) # Fit the era-in data in skdcyclone 
	regions['Z500'] = (70, 45, 40, 85)
	regions['Z10'] = (90, 65, -180, 180) # fit the data directly from ECMWF
	regions['U10'] = (90, 65, -180, 180) # fit the data directly from ECMWF
	regions['SPV'] = (90, 65, -180, 180) # fit the data directly from ECMWF
	regions['DLR'] = (90, 70, -180, 180) # fit the data directly from ECMWF
	regions['Sibcoast'] = (80, 65, 40, 130) # For Inoue et al. 2009. The positive high pressure anomaly
	return regions

	
def clevels_var(ratio=1):

	# Raw values
	clevels_raw = {}
	clevels_raw['ICE'] = np.linspace(0,100,6)
	clevels_raw['THF'] = np.linspace(-50,200,6)
	clevels_raw['SPV'] = np.linspace(21000,22500,6)
	clevels_raw['URALS'] = np.linspace(100000,104000,6)
	clevels_raw['SLP'] = np.linspace(100000,104000,6)
	clevels_raw['DLR'] = np.linspace(100,300,6)
	clevels_raw['VTstar'] = np.linspace(-20,80,6)
	clevels_raw['VTstardaily'] = np.linspace(-20,80,6)
	clevels_raw['Z500'] = np.linspace(5100,5600, 6)
	clevels_raw['U10'] = np.linspace(-120,120,13)
	clevels_raw['Z50'] = np.linspace(19000,21000,6)
	clevels_raw['U850'] = np.linspace(-5,15,6)
	clevels_raw['SST'] = np.linspace(-8,32,6)
	clevels_raw['Z300'] = np.linspace(8300,10000, 6)
	clevels_raw['Z300star'] = np.linspace(-200,300,6)

	# Anomaly composites
	clevels_anom = {}
	clevels_anom['ICE'] = np.linspace(-24, 24, 13)
	clevels_anom['THF'] = np.linspace(-72,72,13)
	clevels_anom['DLR'] = np.linspace(-12,12,13)
	clevels_anom['SPV'] = np.linspace(-630,630,13)
	clevels_anom['Z50'] = np.linspace(-720,720,13)
	clevels_anom['URALS'] = np.linspace(-2400, 2400, 13)
	clevels_anom['SLP'] = np.linspace(-2400, 2400, 13)
	clevels_anom['VTstar'] = np.linspace(-72,72,13)
	clevels_anom['VTstardaily'] = np.linspace(-72,72,13)
	clevels_anom['Z500'] = np.linspace(-200,200,13)
	clevels_anom['Z300'] = np.linspace(-200,200,13)
	clevels_anom['Z300star'] = np.linspace(-200,200,13)
	clevels_anom['U10'] = np.linspace(-12,12,13)
	clevels_anom['Z10'] = np.linspace(-200,200,13)
	clevels_anom['U850'] = np.linspace(-3.6, 3.6, 13)
	clevels_anom['SST'] = np.linspace(-12,12,13)

	# Regress with sea ice
	clevels_regress= {}
	clevels_regress['ICE'] = np.linspace(-1.5, 1.5,13)
	clevels_regress['THF'] = np.linspace(-3,3,13)
	clevels_regress['SPV'] = np.linspace(-9,9,13)
	clevels_regress['URALS'] = np.linspace(-48,48,13)
	clevels_regress['SLP'] = np.linspace(-54,54,13)
	clevels_regress['Z500'] = np.linspace(-5.4, 5.4, 13)
	clevels_regress['VTstar'] = np.linspace(-1.8,1.8,13)
	clevels_regress['VTstardaily'] = np.linspace(-1.8,1.8,13)
	clevels_regress['V*T*'] = np.linspace(-1.8,1.8,13)
	clevels_regress['U10'] = np.linspace(-0.60,0.60,13)
	clevels_regress['Z50'] = np.linspace(-18,18,13)
	clevels_regress['U850'] = np.linspace(-0.30,0.30,13)
	clevels_regress['T850'] = np.linspace(-0.24,0.24,13)
	clevels_regress['DLR'] = np.linspace(-1.2,1.2,13)
	clevels_regress['snowd'] = np.linspace(-1.2,1.2,13)
	clevels_regress['SST'] = np.linspace(-0.12,0.12,13)
	clevels_regress['NAO'] = np.linspace(-0.05,0.05,13)
	clevels_regress['Z300'] = np.linspace(-0.05,0.05,13)
	clevels_regress['Z300star'] = np.linspace(-0.05,0.05,13)
	
	for var in clevels_anom:
		clevels_anom[var] = clevels_anom[var] * ratio
	for var in clevels_raw:
		clevels_raw[var] = clevels_raw[var] * ratio
	for var in clevels_raw:
		clevels_regress[var] = clevels_regress[var] * ratio

	return clevels_raw, clevels_anom, clevels_regress

def map_raw_anom(raw_colors=6):
	if raw_colors==6:
		mapcolor_raw = ['#f0f9e8', '#bae4bc', '#7bccc4', '#43a2ca', '#0868ac']
	elif raw_colors==8:
		mapcolor_raw = ['#f4a582', '#fddbc7', '#fff6e5', '#a6bddb', '#74a9cf', '#3690c0', '#0570b0', '#034e7b'] # the fist two has different color
	elif raw_colors=='8_allblue':
		mapcolor_raw = ['#ece7f2', '#d0d1e6', '#a6bddb', '#74a9cf', '#3690c0', '#0570b0', '#045a8d', '#023858']
	mapcolor_anom  = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#eef7fa', '#74c476', '#88419d', '#fff6e5', '#fddbc7',  '#f4a582', '#d6604d', '#b2182b']
	mapcolor_anom  = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#eef7fa', '#eef7fa', '#fff6e5', '#fff6e5', '#fddbc7',  '#f4a582', '#d6604d', '#b2182b']
	#mapcolor_anom  = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#eef7fa', '#ffffff', '#ffffff', '#fff6e5', '#fddbc7',  '#f4a582', '#d6604d', '#b2182b']

	return mapcolor_raw, mapcolor_anom

def units_var():

	cblab = {}
	cblab['ICE'] = '%'
	cblab['siconc'] = '%'
	cblab['THF'] = r'Wm$^{-2}$'
	cblab['SPV'] = 'm'
	cblab['URALS'] = 'Pa'
	cblab['SLP'] = 'Pa'
	cblab['msl'] = 'Pa'
	cblab['IR'] = r'Wm$^{-2}$'
	cblab['VTstar'] = r'Kms$^{-1}$'
	cblab['VTstardaily'] = r'Kms$^{-1}$'
	cblab['V*T*'] = r'Kms$^{-1}$'
	cblab['V100'] = r'ms$^{-1}$'
	cblab['T100'] = r'K'
	cblab['TS'] = r'K'
	cblab['Z500'] = r'm'
	cblab['Z300'] = r'm'
	cblab['Z300star'] = r'm'
	cblab['Z300starraw'] = r'm'
	cblab['Z300starwave1'] = r'm'
	cblab['Z300starwave2'] = r'm'
	cblab['Z300starwave3'] = r'm'
	cblab['Z10'] = r'm'
	cblab['Z50'] = r'm'
	cblab['U850'] = r'm/s'
	cblab['U500'] = r'm/s'
	cblab['U300'] = r'm/s'
	cblab['U'] = r'm/s'
	cblab['U30'] = r'm/s'
	cblab['U10'] = r'm/s'
	cblab['U10'] = r'm/s'
	cblab['T850'] = r'K'
	cblab['DLR'] = r'W/m2'
	cblab['SST'] = r'K'
	cblab['snowd'] = r'kg$^{3}$'
	cblab['SAT'] = r'K'
	return cblab

def linregress_xarray(y, x, null_hypo=0):
	# y is usually the 3darray
	# x is the timeseries
	# x and y should have the same time dimension
	# null_hypo is the null hypothsis of the slope 
	"""
	Input: Two xr.Datarrays of any dimensions with the first dim being time. 
	Thus the input data could be a 1D time series, or for example, have three dimensions (time,lat,lon). 
	Datasets can be provied in any order, but note that the regression slope and intercept will be calculated
	for y with respect to x.
	Output: Covariance, correlation, regression slope and intercept, p-value, and standard error on regression
	between the two datasets along their aligned time dimension.  
	Lag values can be assigned to either of the data, with lagx shifting x, and lagy shifting y, with the specified lag amount. 
	""" 
	#3. Compute data length, mean and standard deviation along time axis for further use: 
	size  = x.time.shape[0]
	xmean = x.mean(dim='time')
	ymean = y.mean(dim='time')
	xstd  = x.std(dim='time')
	ystd  = y.std(dim='time')

	#4. Compute covariance along time axis
	cov = ((x-xmean)*(y-ymean)).sum(dim='time', skipna=True)/size

	#5. Compute correlation along time axis
	cor = cov/(xstd*ystd)

	#6. Compute regression slope and intercept:
	slope     = cov/(xstd**2)
	intercept = ymean - xmean*slope  

	#7. Compute P-value and standard error
	#Compute t-statistics
	tstats = cor*np.sqrt(size-2)/np.sqrt(1-cor**2)
	stderr = slope/tstats

	if null_hypo!=0:
		# Calculate the standard error manually
		predicted_y = x*slope+intercept
		stderr_new = np.sqrt(((((predicted_y-y)**2).sum(dim='time'))/(size-2)) / (((x-xmean)**2).sum(dim='time')))
		tstats = (slope - null_hypo) / stderr_new

	pval   = scipy.stats.t.sf(np.abs(tstats), size-2)*2
	pval   = xr.DataArray(pval, dims=cor.dims, coords=cor.coords)

	results_ds = xr.Dataset()
	results_ds['covariance'] = cov
	results_ds['correlation'] = cor
	results_ds['slope'] = slope
	results_ds['intercept'] = intercept
	results_ds['pvalues'] = pval
	results_ds['standard_error'] = stderr

	#return cov,cor,slope,intercept,pval,stderr
	return results_ds

def xarray_expand_time_var_dims(ds, expand_time_key, expand_time_sel, ts, compact_var=False):

	# In order to carry out regression via several dimension in 1-step without looping
	# ds should be the raw ds containing all the time
	# ts is the xarray time for regression
	# expand_time_sel is the 

	ts = ts.assign_coords(time=ts.time.dt.year)

	### the 3d-array
	ds_append = []
	for k in expand_time_key:
		ds_temp=ds.sel(time=expand_time_sel[k]) # It also captures the year at this point
		# Always assign the time axis of the time series to it.
		ds_temp = ds_temp.assign_coords(time=ts.time)
		ds_append.append(ds_temp)

	ds_time_expand = xr.concat(ds_append, dim='expand_time_sel').assign_coords(expand_time_sel=expand_time_key)

	if compact_var: # Expand the variable
		ds_time_expand = ds_time_expand.to_array().rename({'variable':'vars'})

	return ds_time_expand, ts
	# Load the results
	#results_ds = linregress_xarray(ds_final, ts)
	#correlation = results_ds['correlation'] 
	#regression = results_ds['slope'] 
	#pvals = results_ds['pvalues'] 
	#return correlation, regression, pvals
