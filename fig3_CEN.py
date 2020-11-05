import numpy as np
import datetime as dt
import pandas as pd
import ipdb
import matplotlib.pyplot as plt
import matplotlib

import utilities as ut
import arrows_plotting_temp as apt
import CEN as CEN


# Author: Peter Siew (yu.siew@uib.no)
# This script generates the Fig.3 and Fig. 6 in Siew et al. WCD

if __name__ == "__main__":

	### Universal parameters: 
	vars = ['ICE', 'THF', 'URALS', 'V*T*','SPV', 'NAO'] # input indicies. They store in ts_anomaly
	scheume_autoc = False # Not Allowing to calculate the correlation of autocorrelation in step1
	f_add_n = 'None'  # Not useful anymore
	step12_sigs = [0.05,0.1,0.2] # The significane level for the 1st and 2nd steps
	step34_sig = 0.05  # The singificnae level for the 3th and 4th steps
	pathway_highlight = True # Turn False - not to highligh the pathway we are interested in

	# The 1st and 2nd step(PC algorithm - actually only the first step in Siew et al. 2020)  - to identify the the possible casual driver of each index (condition selection)
	# The 3rd step (MCI test - actually the 2nd step in Siew et al. 2020) - to identify the "true" causal drivers based on condition in PC-algorithm
	# The 4th step (actually the 3rd step in Siew et al. 2020) - to identify the strength of causal pathway, as well as simultaneous relationship

	### Initiate the plot
	# ax1, ax2 hold the monthy (half-monthly) CEN causal arrows plot
	# ax3, ax5  hold the monthly (half-monthly) simulataneous relationship
	# Neglect ax4 and ax6 - they are not useful anymore with the new CEN algorithm. They should be removed in the future.
	plt.close()
	fig, ((ax1,ax2), (ax3,ax4), (ax5,ax6)) = plt.subplots(3,2, figsize=(24,27))

	# Decide generateing Fig. 3 or Fig. 6 in Siew et al. 2020
	if True: # Fig.3 in Siew et al. 2020 WCD
		fig, (ax1,ax2) = plt.subplots(1,2, figsize=(24,12)) 
		fname = 'fig3_monthly_and_halfmonthly_CEN'
	elif False: # Fig.6 in Siew et al. 2020 WCD
		fig, ax3 = plt.subplots(1,1, figsize=(12,7)) 
		fname = 'fig6_simultaneous_relationship_CEN'
	else: # This shows all the six panels (not useful)
		fname = 'everything_notuseful'
		pass

	### Do the monthly CEN results
	# Intitiate the parameters
	CEN_lag = 2; lags = list(range(0, CEN_lag + 1))
	timescales = 'monthly' # monthly / submonthly / pendtad
	# Read the monthly timeseries
	mon_var_lag_ts = CEN.CEN_mon_var_lag_initiate(timescales, vars, lags)
	# Do the CEN algorithm step 1 - 4
	org_dri_beta_dict, add_dri_beta_dict, logs = CEN.CEN_step1to4_exe(mon_var_lag_ts, step12_sigs, step34_sig, no_autoc=scheume_autoc, skip_step3=False, step3_fdrc=False)
	# Do the monthly plotting
	apt.causal_arrows_plotting_mod(org_dri_beta_dict, add_dri_beta_dict, step34_sig, CEN_lag, timescales, f_add_n, pw_hl=pathway_highlight, ax_all=(ax1,ax3,ax5))

	### Do the half-monthly CEN results
	# Intitiate the parameters
	CEN_lag = 4; lags = list(range(0, CEN_lag + 1))
	timescales = 'submonthly' # monthly / submonthly / pentead
	half_to_full_mon = True
	# Read the half-monthly timeseries
	mon_var_lag_ts = CEN.CEN_mon_var_lag_initiate(timescales, vars, lags)
	# Do the CEN algorithm step 1 - 4
	org_dri_beta_dict, add_dri_beta_dict, logs = CEN.CEN_step1to4_exe(mon_var_lag_ts, step12_sigs, step34_sig, no_autoc=scheume_autoc, skip_step3=False, step3_fdrc=False)
	# Aggregate the half-monthly result into monthly
	if half_to_full_mon & (timescales == 'submonthly'):
		org_dri_beta_dict, add_dri_beta_dict = CEN.aggregated_submonthly(org_dri_beta_dict, add_dri_beta_dict)
		timescales = 'monthly'
	# Do the half-monthly plotting
	apt.causal_arrows_plotting_mod(org_dri_beta_dict, add_dri_beta_dict, step34_sig, CEN_lag, timescales, f_add_n, pw_hl=pathway_highlight, half_to_full_mon=half_to_full_mon, ax_all=(ax2,ax4,ax6))

	### Plotting the miscancellous - e.g., colobar
	ax1.set_title('(a) Monthly CEN', loc='left', size=28)
	ax2.set_title('(b) Half-monthly CEN', loc='left', size=28)
	
	# Plot the colorbar
	colorbar_vmax = 0.8
	mapcolors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
	cmap= matplotlib.colors.ListedColormap(mapcolors)
	cNorm  = matplotlib.colors.Normalize(vmin = -colorbar_vmax, vmax=colorbar_vmax)
	scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm,cmap=cmap)
	cbaxes1 = fig.add_axes([0.2, 0, 0.6, 0.025])
	cb1 = matplotlib.colorbar.ColorbarBase(cbaxes1, cmap=cmap, norm=cNorm, orientation='horizontal',
				ticks=np.linspace(-colorbar_vmax, colorbar_vmax, 9))
	cb1.ax.set_xticklabels([round(i,1) for i in np.linspace(-colorbar_vmax, colorbar_vmax, 9)])
	cb1.set_label('Beta coefficient', fontsize=20)
	cb1.ax.tick_params(labelsize=20)
	
	plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.35)
	plt.savefig('../graphs/%s_%s.png' %(dt.date.today(), fname), bbox_inches='tight')

