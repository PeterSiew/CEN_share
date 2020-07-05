import xarray as xr
import numpy as np
import datetime as dt
import pandas as pd
import ipdb
import matplotlib.pyplot as plt
import matplotlib

import utilities as ut
import arrows_plotting_temp as apt
import utilities as ut
import CEN as CEN

def fig1_revise():

	### Initiate the plot
	plt.close()
	fig, ((ax1,ax2), (ax3,ax4), (ax5,ax6)) = plt.subplots(3,2, figsize=(24,27))
	fig, (ax1,ax2) = plt.subplots(1,2, figsize=(24,12)) # Fig.3
	fig, ax3 = plt.subplots(1,1, figsize=(12,7)) # Fig.5 Simultaneous relationship

	### Universal parameters: 
	vars = ['ICE', 'THF', 'URALS', 'V*T*','SPV', 'NAO']	
	scheume_autoc = False # Not Allowing to calculate the correlation of autocorrelation in step1
	f_add_n = 'None'
	fname = 'fig3_cen_hl'
	step12_sigs = [0.05,0.1,0.2,0.3,0.4,0.5]
	step12_sigs = [0.05,0.1,0.2]
	step34_sig = 0.05
	pathway_highlight = True

	### Monthly
	CEN_lag = 2; lags = list(range(0, CEN_lag + 1))
	timescales = 'monthly' # monthly / submonthly / weekly / weekly_p

	mon_var_lag_ts = CEN.CEN_mon_var_lag_initiate(timescales, vars, lags)
	org_dri_beta_dict, add_dri_beta_dict, logs = CEN.CEN_step1to4_exe(mon_var_lag_ts, step12_sigs, step34_sig, no_autoc=scheume_autoc, skip_step3=False, step3_fdrc=False)
	apt.causal_arrows_plotting_mod(org_dri_beta_dict, add_dri_beta_dict, step34_sig, CEN_lag, timescales, f_add_n, pw_hl=pathway_highlight, ax_all=(ax1,ax3,ax5))

	### Half-monthly
	CEN_lag = 4; lags = list(range(0, CEN_lag + 1))
	timescales = 'submonthly' # monthly / submonthly / weekly / weekly_p
	half_to_full_mon = True

	mon_var_lag_ts = CEN.CEN_mon_var_lag_initiate(timescales, vars, lags)
	org_dri_beta_dict, add_dri_beta_dict, logs = CEN.CEN_step1to4_exe(mon_var_lag_ts, step12_sigs, step34_sig, no_autoc=scheume_autoc, skip_step3=False, step3_fdrc=False)

	if half_to_full_mon & (timescales == 'submonthly'):
		org_dri_beta_dict, add_dri_beta_dict = CEN.aggregated_submonthly(org_dri_beta_dict, add_dri_beta_dict)
		timescales = 'monthly'
	apt.causal_arrows_plotting_mod(org_dri_beta_dict, add_dri_beta_dict, step34_sig, CEN_lag, timescales, f_add_n, pw_hl=pathway_highlight, half_to_full_mon=half_to_full_mon, ax_all=(ax2,ax4,ax6))

	### Plotting setups
	ax1.set_title('(a) Monthly CEN', loc='left', size=28)
	ax2.set_title('(b) Half-monthly CEN', loc='left', size=28)
	
	colorbar_vmax = 0.8
	mapcolors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
	cmap= matplotlib.colors.ListedColormap(mapcolors)
	cNorm  = matplotlib.colors.Normalize(vmin = -colorbar_vmax, vmax=colorbar_vmax)
	scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm,cmap=cmap)
	#cbaxes1 = fig.add_axes([0.35, 0.01, 0.3, 0.025])
	cbaxes1 = fig.add_axes([0.2, 0, 0.6, 0.025])
	cb1 = matplotlib.colorbar.ColorbarBase(cbaxes1, cmap=cmap, norm=cNorm, orientation='horizontal',
				ticks=np.linspace(-colorbar_vmax, colorbar_vmax, 9))
	cb1.ax.set_xticklabels([round(i,1) for i in np.linspace(-colorbar_vmax, colorbar_vmax, 9)])
	cb1.set_label('Beta coefficient', fontsize=20)
	cb1.ax.tick_params(labelsize=20)
	
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.35)
	plt.savefig('../graphs/%s_fig3_CEN_%s.png' %(dt.date.today(), fname), bbox_inches='tight')
	plt.savefig('../graphs/%s_fig3_CEN_%s.pdf' %(dt.date.today(), fname), bbox_inches='tight', dpi=600)

def fig1():
	vars = ['ICE', 'THF', 'URALS', 'V*T*','SPV', 'NAO']	
	scheume_autoc = False # Not Allowing to calculate the correlation of autocorrelation in step1
	pathway_highlight= False
	f_add_n = 'None'
	fname = 'fig3_cen_hl'

	plt.close()
	fig, ((ax1,ax2), (ax3,ax4), (ax5,ax6)) = plt.subplots(3,2, figsize=(24,27))
	fig, (ax1,ax2) = plt.subplots(1,2, figsize=(24,9))
	#fig, ax3 = plt.subplots(1,1, figsize=(12,6))
	#fig, ax1 = plt.subplots(1,1, figsize=(12,6))

	### Monthly
	sig_lev = 0.01
	CEN_lag = 2
	timescales = 'monthly' # monthly / submonthly / weekly / weekly_p
	
	months = ['Nov', 'Dec', 'Jan', 'Feb', 'Mar']
	monthly_ts_ds = ut.read_ts_cache(vars, year_control=False, dt_start = '1979-09-01', dt_end = '2018-03-01')
	lags = list(range(0, CEN_lag + 1))
	mon_var_lag_ts, mon_var_lag_dt, var_lag_ds, var_lag_dt = CEN.read_CEN_format_ts(vars, months, monthly_ts_ds, lags, timescales)

	poten_dri_dict, step1_log = CEN.CEN_step1(mon_var_lag_ts, sig_lev=sig_lev)
	dri_dict, step2_log = CEN.CEN_step2(mon_var_lag_ts, poten_dri_dict, sig_lev=sig_lev)
	org_dri_beta_dict, add_dri_beta_dict, step3_log = CEN.CEN_step3(mon_var_lag_ts, dri_dict, no_add_dri=False, sig_lev=sig_lev)

	apt.causal_arrows_plotting_mod(org_dri_beta_dict, add_dri_beta_dict, sig_lev, CEN_lag,
						timescales, f_add_n, pw_hl=pathway_highlight, ax_all=(ax1,ax3,ax5))

	# Half-monthly
	sig_lev = 0.01
	CEN_lag = 4
	timescales = 'submonthly' # monthly / submonthly / weekly / weekly_p
	half_to_full_mon = True

	months = ['Nov-1', 'Nov-2', 'Dec-1', 'Dec-2', 'Jan-1', 'Jan-2', 'Feb-1', 'Feb-2', 'Mar-1', 'Mar-2']
	monthly_ts_ds = ut.read_ts_cache(vars, year_control=False, dt_start = '1979-09-01', dt_end = '2018-03-16', resolution='halfmonthly')
	lags = list(range(0, CEN_lag + 1))
	mon_var_lag_ts, mon_var_lag_dt, var_lag_ds, var_lag_dt = CEN.read_CEN_format_ts(vars, months, monthly_ts_ds, lags, timescales)

	poten_dri_dict, step1_log = CEN.CEN_step1(mon_var_lag_ts, sig_lev=sig_lev)
	dri_dict, step2_log = CEN.CEN_step2(mon_var_lag_ts, poten_dri_dict, sig_lev=sig_lev)
	org_dri_beta_dict, add_dri_beta_dict, step3_log = CEN.CEN_step3(mon_var_lag_ts, dri_dict, no_add_dri=False, sig_lev=sig_lev)

	if half_to_full_mon & (timescales == 'submonthly'):
		org_dri_beta_dict, add_dri_beta_dict = CEN.aggregated_submonthly(org_dri_beta_dict, add_dri_beta_dict)
		timescales = 'monthly'

	# Start plotting
	apt.causal_arrows_plotting_mod(org_dri_beta_dict, add_dri_beta_dict, sig_lev, CEN_lag,
						timescales, f_add_n, pw_hl=pathway_highlight, half_to_full_mon=half_to_full_mon, ax_all=(ax2,ax4,ax6))
	
	ax1.set_title('(a) Monthly CEN', loc='left', size=28)
	ax2.set_title('(b) Half-monthly CEN', loc='left', size=28)
	
	colorbar_vmax = 0.8
	mapcolors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
	cmap= matplotlib.colors.ListedColormap(mapcolors)
	cNorm  = matplotlib.colors.Normalize(vmin = -colorbar_vmax, vmax=colorbar_vmax)
	scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm,cmap=cmap)
	#cbaxes1 = fig.add_axes([0.35, 0.01, 0.3, 0.025])
	cbaxes1 = fig.add_axes([0.2, -0.03, 0.6, 0.025])
	cb1 = matplotlib.colorbar.ColorbarBase(cbaxes1, cmap=cmap, norm=cNorm, orientation='horizontal',
				ticks=np.linspace(-colorbar_vmax, colorbar_vmax, 9))
	cb1.ax.set_xticklabels([round(i,1) for i in np.linspace(-colorbar_vmax, colorbar_vmax, 9)])
	cb1.set_label('Beta coefficient', fontsize=20)
	cb1.ax.tick_params(labelsize=20)
	
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.35)
	plt.savefig('../graphs/fig3_CEN_%s.png' %fname, bbox_inches='tight')
	plt.savefig('../graphs/fig3_CENs_%s.pdf' %fname, bbox_inches='tight', dpi=600)

if __name__ == "__main__":
	fig1_revise()
	#fig1()