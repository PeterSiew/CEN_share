import numpy as np
import datetime as dt
import copy
from scipy import stats
import scipy
import itertools
import ipdb
import statsmodels.api as sm

import arrows_plotting_temp as apt
import utilities as ut
import CEN as CEN



def CEN_mon_var_lag_initiate(timescales, vars, lags, step3_px=None):

	if timescales == 'monthly':
		months = ['Nov', 'Dec', 'Jan', 'Feb', 'Mar']
		monthly_ts_ds = ut.read_ts_cache(vars, year_control=False, dt_start = '1979-09-01', dt_end = '2018-03-01', resolution='monthly')

	elif timescales == 'submonthly':
		months = ['Nov-1', 'Nov-2', 'Dec-1', 'Dec-2', 'Jan-1', 'Jan-2', 'Feb-1', 'Feb-2', 'Mar-1', 'Mar-2']
		# It has to be 16 Mar in order to include the 2nd part
		monthly_ts_ds = ut.read_ts_cache(vars, year_control=False, dt_start = '1979-09-01', dt_end = '2018-03-16', resolution = 'halfmonthly')

	elif timescales == 'pentad':
		# For the IR, URALS and ICE
		# For the ICE, THF, URALS, V*T*, SPV, NAO
		months = ['Oct-1', 'Oct-2', 'Oct-3', 'Oct-4', 'Oct-5', 'Oct-6', 'Nov-1', 'Nov-2', 'Nov-3', 'Nov-4', 'Nov-5', 'Nov-6', 'Dec-1', 'Dec-2', 'Dec-3', 'Dec-4', 'Dec-5', 'Dec-6',
		          'Jan-1', 'Jan-2', 'Jan-3', 'Jan-4', 'Jan-5', 'Jan-6', 'Feb-1', 'Feb-2', 'Feb-3', 'Feb-4', 'Feb-5', 'Feb-6',
				  'Mar-1', 'Mar-2', 'Mar-3', 'Mar-4', 'Mar-5', 'Mar-6']
		months = ['Nov-1', 'Nov-2', 'Nov-3', 'Nov-4', 'Nov-5', 'Nov-6', 'Dec-1', 'Dec-2', 'Dec-3', 'Dec-4', 'Dec-5', 'Dec-6',
		          'Jan-1', 'Jan-2', 'Jan-3', 'Jan-4', 'Jan-5', 'Jan-6', 'Feb-1', 'Feb-2', 'Feb-3', 'Feb-4', 'Feb-5', 'Feb-6',
				  'Mar-1', 'Mar-2', 'Mar-3', 'Mar-4', 'Mar-5', 'Mar-6']
		monthly_ts_ds = ut.read_ts_cache(vars, year_control=False, dt_start = '1979-09-01', dt_end = '2018-03-26', resolution='pentad')

	if (step3_px!=None) & (timescales=='monthly'): # Just for temporary. Please delete in the future
		months = ['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar']
		monthly_ts_ds = ut.read_ts_cache(vars, year_control=False, dt_start = '1979-07-01', dt_end = '2018-03-01', resolution='monthly')

	mon_var_lag_ts, mon_var_lag_dt, var_lag_ds, var_lag_dt = CEN.read_CEN_format_ts(vars, months, monthly_ts_ds, lags, timescales)

	return mon_var_lag_ts

def CEN_step1to4_exe(mon_var_lag_ts, step12_sigs, step34_sig, no_autoc=False, skip_step3=False, step3_fdrc=False, step3_px=None):

	# Do the AIC test for each significant level and for each variable
	dri_dict_condition = {}
	AIC = {}
	step1_log, step2_log = {}, {}
	for sig in step12_sigs:
		poten_dri_dict, step1_log[sig] = CEN.CEN_step1(mon_var_lag_ts, sig_lev=sig, scheume_autoc=no_autoc)
		dri_dict_condition[sig], step2_log[sig] = CEN.CEN_step2(mon_var_lag_ts, poten_dri_dict, sig_lev=sig, parents_combination=False)
		AIC[sig] = AIC_calculation(mon_var_lag_ts, dri_dict_condition[sig])
	dri_dict_condition_sel, AIC_log = AIC_selection(mon_var_lag_ts, AIC, dri_dict_condition, step12_sigs)

	if skip_step3:
		dri_dict = dri_dict_condition_sel
		step3_log = []
	else:
		dri_dict, step3_log = CEN_step3(mon_var_lag_ts, dri_dict_condition_sel, sig_lev=step34_sig, px=step3_px, scheume_autoc=no_autoc, fdrc=step3_fdrc)

	org_dri_beta_dict, add_dri_beta_dict, step4_log = CEN.CEN_step4(mon_var_lag_ts, dri_dict, no_add_dri=False, sig_lev=step34_sig)

	return org_dri_beta_dict, add_dri_beta_dict, (step1_log, step2_log, AIC_log, step3_log, step4_log)

def AIC_calculation(mvl_ts, dri_dict):
	mons = [*mvl_ts]
	vars = [*mvl_ts[mons[0]]]

	AIC = {mon:{var:None for var in vars} for mon in mons}
	for mon in mons:
		for var in vars:
			var_ts = mvl_ts[mon][var][0]
			var_ts_s = (var_ts - var_ts.mean()) / var_ts.std()
			dris = [i[0] for i in dri_dict[mon][var]]
			dri_lags = [i[1] for i in dri_dict[mon][var]]
			
			if len(dris)==0: # It has no driver.
				# Create an virtual driver with all zero values. It won't affect the results
				#dris_ts_s = np.zeros(len(var_ts_s)).T
				AIC[mon][var] = np.nan

			else: # It has other drivers
				dris_ts = [mvl_ts[mon][dri][dri_lag] for dri, dri_lag in zip(dris, dri_lags)]
				dris_ts = np.array(dris_ts).T
				dris_ts_s = (dris_ts - dris_ts.mean(axis=0)) / dris_ts.std(axis=0)
				dris_ts_matrix = sm.add_constant(dris_ts_s)
				regress_results = sm.OLS(endog = var_ts_s, exog = dris_ts_matrix).fit()
				RSS = ((regress_results.fittedvalues  - var_ts_s)**2).sum()
				AIC[mon][var] = 2*len(dris) + len(var_ts) * np.log(RSS)

	return AIC

def AIC_selection(mvl_ts, AIC, dri_dict_condition_all, sig_levs):

	mons = [*mvl_ts]
	vars = [*mvl_ts[mons[0]]]

	# Choose the best significant level for each driver
	dri_dict_condition_sel = {mon:{var:{} for var in vars} for mon in mons}
	AIC_log = {mon:{var:{} for var in vars} for mon in mons}
	for mon in mons:
		for var in vars:
			AIC_sig = [AIC[sig][mon][var] for sig in sig_levs]
			# Select the minimum AIC and the significant level
			if not np.isnan(AIC_sig).all(): 
				min_sig_idx = np.nanargmin(AIC_sig)
			else: # All are nan, np.isnan(AIC_sig).all() return True
				min_sig_idx = 0
			# In the case of all np.nan (no driver)
			chosen_sig = sig_levs[min_sig_idx]
			dri_dict_condition_sel[mon][var] = dri_dict_condition_all[chosen_sig][mon][var]
			AIC_log[mon][var] = chosen_sig if not np.isnan(AIC_sig).all() else None 
	
	return dri_dict_condition_sel, AIC_log


def read_CEN_format_ts(vars, months, ts_ds, lags, timescales):
	
	# 1: 1st-5st, 2: 6th-10th, 3:11th-15th, 4:16th-20th, 5:21th-25th, 6:26th-28/30/31th
	months_no = {'Sep':9, 'Oct':10, 'Nov':11, 'Dec':12, 'Jan':1, 'Feb':2, 'Mar':3}
	halfmonthly_day = {'1':1, '2':16}
	pentad_day = {'1':1, '2':6, '3':11, '4':16, '5':21, '6':26}
	
	mon_var_lag = {}
	for mon in months:
		mon_var_lag[mon] = {}
		if timescales == 'monthly':
			mask = ts_ds.time.dt.month == months_no[mon]
		elif timescales == 'submonthly':
			mask = (ts_ds.time.dt.month == months_no[mon[0:3]]) & (ts_ds.time.dt.day == halfmonthly_day[mon[-1]])
		elif timescales == 'pentad':
			mask = (ts_ds.time.dt.month == months_no[mon[0:3]]) & (ts_ds.time.dt.day == pentad_day[mon[-1]])
		elif timescales == 'DJF_pentad':
			mask = ts_ds.time.dt.month.isin([12,1,2])
		else:
			raise ValueError('Please check')
	
		for var in vars:
			mon_var_lag[mon][var] = {}
			mon_var_lag[mon][var][0] = ts_ds[var].where(mask, drop=True).values
			
			for lag in lags:
				mask_lag = mask.values.nonzero()[0] - lag # mask lag indeed is an array of index
				if (mask_lag < 0).any():
					raise ValueError('the mask has negative value which is invalid')
				# the mask is based on monthly_ts_ds. Monthly is already selected via indexing
				mon_var_lag[mon][var][lag] =  ts_ds[var][mask_lag]
	
	# Extract data out of xarray format
	all_lag_m_ts = copy.deepcopy(mon_var_lag)
	all_lag_m_dt = copy.deepcopy(mon_var_lag)
	for mon in months:
		for var in vars:
			for lag in lags:
				 all_lag_m_ts[mon][var][lag] = mon_var_lag[mon][var][lag].values
				 all_lag_m_dt[mon][var][lag] = mon_var_lag[mon][var][lag].time.values

	# Create an all timeseries(old format)
	all_m_ts, all_m_dt = {}, {}
	for var in vars:
		all_m_ts[var] = ts_ds[var].values
		all_m_dt[var] = ts_ds[var].time
		
	return all_lag_m_ts, all_lag_m_dt, all_m_ts, all_m_dt

def CEN_step1(mvl_ts, scheume_autoc=False, sig_lev=0.05):
	# mvl_ts refers to mon_var_lag ts
	mons = list(mvl_ts.keys())
	vars = list(mvl_ts[mons[0]].keys())
	lags = np.sort(list(mvl_ts[mons[0]][vars[0]].keys()))
	
	corrs_log = []
	
	pdris_dict = {}
	for mon in mons:
		pdris_dict[mon] = {}
		for var in vars:
			pdris_dict[mon][var] = {} # It will be empty like [] if there is no potential driver

	for mon in mons:
		for var in vars:
			# 1st step: Look for potentail drivers by normal correlation
			var_ts = mvl_ts[mon][var][0]
			for pdri in vars:
				if scheume_autoc: # Not doing the autocoorelation
					if pdri == var:
						continue # Skip this loop and look for other potential driver
				for lag in lags: # Not doing lag 0 in step 1
					if (pdri == var) & ~(lag == 1): # Autocorrelation just do for one month lag. Not doing 0 lag (which gives 1 correlation)
						continue
	
					pdri_ts = mvl_ts[mon][pdri][lag]
					corr, pval = stats.pearsonr(var_ts, pdri_ts)
					corrs_log.append([mon, var, pdri, lag, corr, pval])
					
					if (pval <= sig_lev) & (lag > 0):
						# Save the potential driver for each variable
						pdris_dict[mon][var][(pdri, lag)] = corr
						
	return pdris_dict, corrs_log

def CEN_step2(mvl_ts, pdri_dict, sig_lev=0.05, parents_combination=False):

	mons = list(mvl_ts.keys())
	vars = list(mvl_ts[mons[0]].keys())
	lags = np.sort(list(mvl_ts[mons[0]][vars[0]].keys()))
	
	step2_log = []
	dri_dict = {}
	for mon in mons:
		dri_dict[mon] = {}
		for var in vars:
			dri_dict[mon][var] = {}
	
	for mon in mons:
		for var in vars:
			# Read all the drivers and do the sorted. To loop for the conditions
			
			pdri = pdri_dict[mon][var]
			
			if len(list(pdri.keys())) == 0: # this potential driver might be empty
				continue
			
			pdri_keys = [i[0] for i in list(pdri.items())]
			pdri_abs_coeff = [abs(i[1]) for i in list(pdri.items())]
			pdri_sort = [j for i, j  in sorted(zip(pdri_abs_coeff, pdri_keys), reverse = True)]
			pdri_valid = np.ones((len(pdri_sort)), dtype=bool) # this determines if this driver is still valid. Initially all are true
	
			# Loop for the pdri and find the condition first
			for i, dri_lag in enumerate(pdri_sort): # dri_list is sorted
				
				if pdri_valid[i] == False: # if the dri is stil a ponentrial dirver
					continue
	
				cvar = dri_lag[0]
				clag = dri_lag[1]
				cts = mvl_ts[mon][cvar][clag]

				# Loop for other drivers
				for i, dri_lag1 in enumerate(pdri_sort):

					if pdri_valid[i] == False: # this this is invalid, skip this loop
						continue
						
					pdri = dri_lag1[0]
					pdri_lag = dri_lag1[1]  
					pdri_ts = mvl_ts[mon][pdri][pdri_lag]
		
					if (pdri == cvar) & (pdri_lag == clag): # if it is equal to the conditions, skip this loop
						continue
				
					partial_corr_data = np.column_stack((mvl_ts[mon][var][0], pdri_ts, cts))
					p_corr =  partial_corr(partial_corr_data)[0][1] # correlation of the first column and 2nd column given the 3rd column
					K = partial_corr_data.shape[1]
					N = len(mvl_ts[mon][var][0])
					t_stat = p_corr * ((N-2-K) / (1-p_corr**2))**0.5
					pval = stats.t.sf(np.abs(t_stat), N-2-K)*2
					step2_log.append([mon, var, pdri, pdri_lag, cvar, clag, p_corr, pval])
					if pval <= sig_lev:
						pass
						
					else:
						pdri_valid[i] = False

			pdri_sort = [pdri_sort[i] for i in pdri_valid.nonzero()[0]] # To remove the False value item
			pdri_length = len(pdri_sort)
									
			####### Combination. We find out two conditions first. And then do partial correlation test based on those two condtions
			if (pdri_length >= 3) & parents_combination: #there has to be at least 3 potential drivers
				
				pdri_valid = np.ones((len(pdri_sort)), dtype=bool)
				z_ncr = [i for i in  itertools.combinations(list(range(pdri_length)), 2)] # Conditions
				for i, j in z_ncr: #z_ncr size can be different from pdri_sort. E.g., pdri_length = 5, zcr can be 10 (5C2 =10)
					if (pdri_valid[i] == False) | (pdri_valid[j] == False):
						continue

					cvar1 = pdri_sort[i][0]
					clag1 = pdri_sort[i][1]
					cts1 = mvl_ts[mon][cvar1][clag1]
					cvar2 = pdri_sort[j][0]
					clag2 = pdri_sort[j][1]
					cts2 = mvl_ts[mon][cvar2][clag2]
			
					for k, dri_lag in enumerate(pdri_sort):
						
						if pdri_valid[k] == False:
							continue			
						pdri = dri_lag[0]
						pdri_lag = dri_lag[1]
						pdri_ts = mvl_ts[mon][pdri][pdri_lag]
					
						if ((pdri == cvar1) & (pdri_lag == clag1)) | ((pdri == cvar2) & (pdri_lag == clag2)): # either in one case Finish the loop
							continue

						partial_corr_data = np.column_stack((mvl_ts[mon][var][0], pdri_ts, cts1, cts2))
						p_corr =  partial_corr(partial_corr_data)[0][1]
						K = partial_corr_data.shape[1]
						N = len(mvl_ts[mon][var][0])
						t_stat = p_corr * ((N-2-K) / (1-p_corr**2))**0.5
						pval = stats.t.sf(np.abs(t_stat), N-2-K)*2
						step2_log.append([mon, var, pdri, pdri_lag, cvar1, clag1, cvar2, clag2, p_corr, pval])

						if pval <= sig_lev:
							pass
						else:
							pdri_valid[k] = False

				# The new pdri_sort after combination test
				pdri_sort = [pdri_sort[i] for i in pdri_valid.nonzero()[0]] # To remove the False value item
			
			dri_dict[mon][var] = pdri_sort # just save the driver name and lag
			
	return dri_dict, step2_log

def CEN_step3(mvl_ts, dri_dict, sig_lev=0.05, px=None, scheume_autoc=False, fdrc=False):

	# dri_dict is used as conditions
	conditions_dict = dri_dict
	mons = list(mvl_ts.keys())
	vars = list(mvl_ts[mons[0]].keys())
	lags = np.sort(list(mvl_ts[mons[0]][vars[0]].keys()))

	if not px is None: # Just for temporary. Please delete in the future
		mons = ['Nov', 'Dec', 'Jan', 'Feb', 'Mar']

	dri_dict = {mon:{var:[] for var in vars} for mon in mons}
	step3_log = []
	mvpl_pvals = []
	# Test all links again with 
	for mon in mons:
		for var in vars:
			var_ts = mvl_ts[mon][var][0]
			for pdri in vars:

				if (pdri==var) & scheume_autoc:
					continue

				for lag in np.delete(lags,0): # Delete the lag 0 from lags
					if (pdri == var) & ~(lag == 1): # Autocorrelation just do for one month lag. Not doing 0 lag (which gives 1 correlation)
						continue

					condition_vars_ts = [mvl_ts[mon][v][l] for v,l in conditions_dict[mon][var] if (v,l)!=(pdri,lag)] # Eliminate the X from the codition_ts if X is within the condition
					pdri_ts = mvl_ts[mon][pdri][lag] # sane as mvl_ts[mon_lag][pdri][0]

					if not px is None:
						months_all = list(mvl_ts.keys())
						if (months_all.index(mon)-lag) < 0:
							ipdb.set_trace(); raise ValueError('check this line')
						mon_lag = months_all[months_all.index(mon) - lag]
						condition_pdri_ts = [mvl_ts[mon_lag][v][l] for v,l in conditions_dict[mon_lag][pdri]]
					else:
						condition_pdri_ts = []
					conditions_combine = condition_vars_ts + condition_pdri_ts
					conditions_combine = np.zeros(var_ts.shape) if len(conditions_combine)==0 else np.column_stack(conditions_combine)

					if False:
						print('Mon, Var: ',mon,var)
						print('Lag, Pdri: ',lag,pdri)
						print('Vars condition: ', conditions_dict[mon][var])
						#print('Pdri condition: ', conditions_dict[mon_lag][pdri])

					# Carry out the partial correlation test
					partial_corr_data = np.column_stack((var_ts, pdri_ts, conditions_combine))
					partial_corr_data = partial_corr_data[:, ~np.all(partial_corr_data==0, axis=0)] # Delete the column with 0
					p_corr =  CEN.partial_corr(partial_corr_data)[0][1]
					K = partial_corr_data.shape[1]
					N = len(var_ts)
					t_stat = p_corr * ((N-2-K) / (1-p_corr**2))**0.5
					pval = stats.t.sf(np.abs(t_stat), N-2-K)*2
					step3_log.append([mon, var, pdri, lag, conditions_dict[mon][var], p_corr, pval])
					mvpl_pvals.append([mon, var, pdri, lag, pval])
					if pval <= sig_lev:
						dri_dict[mon][var].append((pdri, lag))

	if fdrc: # calculate the adjusted q or p value
		mvpl_pvals = np.array(mvpl_pvals)
		sort_idx = np.array(mvpl_pvals[:,4], dtype='float').argsort()
		mvpl_pvals_sort = mvpl_pvals[sort_idx]
		pvals_sort = np.array(mvpl_pvals_sort[:,4], dtype='float')
		# Adjust
		ranks = np.arange(1,len(mvpl_pvals)+1)
		pvals_adjust = (pvals_sort * len(mvpl_pvals)) / ranks
		#sig_adjust = 0.05*ranks/len(mvpl_pvals)
		ipdb.set_trace()
		
		mask_idx = (pvals_adjust < sig_lev).nonzero()[0]
		ipdb.set_trace()
		dri_dict = {mon:{var:[] for var in vars} for mon in mons}
		for m in mask_idx:
			mon = mvpl_pvals_sort[m][0]
			var = mvpl_pvals_sort[m][1]
			pdri = mvpl_pvals_sort[m][2]
			lag = int(mvpl_pvals_sort[m][3])
			dri_dict[mon][var].append((pdri, lag))
		ipdb.set_trace()


	return dri_dict, step3_log
			
def CEN_step4(mvl_ts, dri_dict, no_add_dri=False, sig_lev=0.05):
	
	mons = list(dri_dict.keys())
	vars = list(mvl_ts[mons[0]].keys())
	if no_add_dri:
		lags = [0]
	else:
		lags = np.sort(list(mvl_ts[mons[0]][vars[0]].keys()))

	dris_beta_dict, add_dris_beta_dict = {}, {}
	for mon in mons:
		dris_beta_dict[mon], add_dris_beta_dict[mon] = {}, {}
		for var in vars:
			dris_beta_dict[mon][var], add_dris_beta_dict[mon][var] = {}, {}

	step4_log = []
	for mon in mons:
		for var in vars:
			# The ts of the variable
			
			var_ts = mvl_ts[mon][var][0]
			var_ts_s = (var_ts - var_ts.mean()) / var_ts.std()
			dris = [i[0] for i in dri_dict[mon][var]]
			dri_lags = [i[1] for i in dri_dict[mon][var]]
			
			if len(dris) == 0: # It has no driver.
				# Create an virtual driver with all zero values. It won't affect the results
				dris_ts_s = np.zeros(len(var_ts_s)).T

			else:  # It has other drivers
				
				dris_ts = [mvl_ts[mon][dri][dri_lag] for dri, dri_lag in zip(dris, dri_lags)]
				dris_ts = np.array(dris_ts).T
				dris_ts_s = (dris_ts - dris_ts.mean(axis=0)) / dris_ts.std(axis=0)
				dris_ts_matrix = sm.add_constant(dris_ts_s)
				regress_results = sm.OLS(endog = var_ts_s, exog = dris_ts_matrix).fit()
				b_coeff = regress_results.params
				p_vals = regress_results.pvalues
				
				# The first one is beta-coeff of the intercept
				for i, dri in enumerate(dris):
					# the 0 index of b_coeff and p_coeff is the Beta value for the constant
					dris_beta_dict[mon][var][(dri, dri_lags[i])] = b_coeff[i+1]
					step4_log.append([mon, var, dri, dri_lags[i], b_coeff[i+1], p_vals[i+1]])
				
			# For additional variable. Including simultaneous relationship
			for add_dri in vars:
				for lag in lags: # including all the lags
					if (add_dri == var) & (lag == 0): #itself. Then skip this
						continue
					if (add_dri, lag) in list(zip(dris, dri_lags)): #if it is already in the driver list. 
						continue
					add_ts = mvl_ts[mon][add_dri][lag]
					add_ts_s = (add_ts - add_ts.mean()) / add_ts.std()
					dris_add_ts_s = np.column_stack((dris_ts_s, add_ts_s))
					#dris_ts_matrix = np.column_stack((np.ones(len(var_ts_s)), dris_add_ts_s))
					dris_ts_matrix = sm.add_constant(dris_add_ts_s)
					regress_results = sm.OLS(endog = var_ts_s, exog = dris_ts_matrix).fit()
					b_coeff = regress_results.params
					p_vals = regress_results.pvalues
					
					if p_vals[-1] <= sig_lev:
						# the additional driver is the last value
						#add_dris_beta_dict[mon][var].append([add_dri, lag, b_coeff[-1]])
						add_dris_beta_dict[mon][var][(add_dri, lag)] = b_coeff[-1]
						step4_log.append([mon, var, add_dri, lag, b_coeff[-1], p_vals[-1]])
	
	return dris_beta_dict, add_dris_beta_dict, step4_log


def partial_corr(C):

	if False:
		# The direct equation. Answer is the same as this function
		p_corr1 = (rxy-(rxz*ryz)) / (((1-ryz**2)**0.5)*((1-rxz**2)**0.5))

	C = np.asarray(C)
	p = C.shape[1]
	P_corr = np.zeros((p, p), dtype=np.float)
	for i in range(p):
		P_corr[i, i] = 1
		for j in range(i+1, p):
			idx = np.ones(p, dtype=np.bool)
			idx[i] = False
			idx[j] = False
			beta_i = scipy.linalg.lstsq(C[:, idx], C[:, j])[0]
			beta_j = scipy.linalg.lstsq(C[:, idx], C[:, i])[0]
			res_j = C[:, j] - C[:, idx].dot( beta_i)
			res_i = C[:, i] - C[:, idx].dot(beta_j)
			corr = stats.pearsonr(res_i, res_j)[0]
			P_corr[i, j] = corr
			P_corr[j, i] = corr
	return P_corr

def aggregated_submonthly(org_dri_dict, add_dri_dict):

	# This is a function to change the submonthly result into monthly format

	def submonthly_to_monthly_step4_dict(dri_dict_submonthly):

		# These two full mons list are only for calculating the lag. Not for looping
		full_monthly_mons = ['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar']
		full_halfmonthly_mons = ['Sep-1', 'Sep-2', 'Oct-1', 'Oct-2', 'Nov-1', 'Nov-2', 'Dec-1', 'Dec-2', 'Jan-1', 'Jan-2', 'Feb-1', 'Feb-2', 'Mar-1', 'Mar-2']

		halfmonthly_mons = list(dri_dict_submonthly.keys())
		vars = list(dri_dict_submonthly[halfmonthly_mons[-1]].keys())
		monthly_mons = ['Nov', 'Dec', 'Jan', 'Feb', 'Mar']

		# Initiate the dict
		dri_dict_monthly = {mon:{var:{} for var in vars} for mon in monthly_mons}

		for mon in halfmonthly_mons:
			for var in vars:

				dris_lags = list(dri_dict_submonthly[mon][var].keys())
				
				if len(dris_lags) == 0: # It has not driver
					continue

				for dri, lag in dris_lags:
					beta = dri_dict_submonthly[mon][var][(dri, lag)]
					# Put this beta into a new dict that resemeble the monthly resolution

					new_mon = mon[0:3]
					dri_new_mon = full_halfmonthly_mons[full_halfmonthly_mons.index(mon) - lag][0:3]
					lag_new = full_monthly_mons.index(new_mon) - full_monthly_mons.index(dri_new_mon)

					if (dri, lag_new) not in dri_dict_monthly[new_mon][var]:
						dri_dict_monthly[new_mon][var][(dri, lag_new)] = []
					dri_dict_monthly[new_mon][var][(dri, lag_new)].append(beta)

					
		# Do the average if it has two
		# This situation is really rare
		for mon in monthly_mons:
			for var in vars:
				dris_lags = list(dri_dict_monthly[mon][var].keys())
				
				if len(dris_lags) == 0: # It has not driver
					continue
				for dri, lag in dris_lags:
					if ((mon, var, dri, lag) == ('Dec', 'V*T*', 'THF', 1)) & True: # There are two linkages which have opposite sign, cancelling each other
						dri_dict_monthly[mon][var][(dri, lag)] = [dri_dict_monthly[mon][var][(dri, lag)][1]]
					dri_dict_monthly[mon][var][(dri, lag)] = np.mean(dri_dict_monthly[mon][var][(dri, lag)])

		return dri_dict_monthly

	org_dri_dict_monthly = submonthly_to_monthly_step4_dict(org_dri_dict)
	add_dri_dict_monthly = submonthly_to_monthly_step4_dict(add_dri_dict)

	return org_dri_dict_monthly, add_dri_dict_monthly	


if __name__ == "__main__":
	main()
	#dynmiate_CEN()
	#mingfang_CEN()
	#mingfang_CEN_new()
	#main_bs()
