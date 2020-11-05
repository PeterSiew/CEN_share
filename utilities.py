import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from datetime import datetime
import xarray as xr
import ipdb
import pandas as pd
import pickle

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
		var_f['NAO'] = ts_folder + 'NAO_ts.p'
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


