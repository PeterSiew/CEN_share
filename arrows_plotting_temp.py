import numpy as np
import datetime as dt
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
import ipdb



def causal_arrows_plotting_mod(org_dri_dict, add_dri_dict, sig_lev, lag_num, timescales, f_add_n, half_to_full_mon=False, pw_hl=False, ax_all=None, only_keep_ax1=False, **aargus):

	# c_mon here is for label. So it has two more months than dri_dict
	if timescales == 'monthly':
		c_mons = ['Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb', 'Mar']
		x_adj_m = {'Nov':-0.32, 'Dec':-0.16, 'Jan':0, 'Feb':0.16, 'Mar':0.32}
	elif timescales == 'submonthly':
		c_mons = ['Sep-1', 'Sep-2', 'Oct-1', 'Oct-2', 'Nov-1', 'Nov-2', 
				'Dec-1', 'Dec-2', 'Jan-1', 'Jan-2', 'Feb-1', 'Feb-2', 'Mar-1', 'Mar-2']
		x_adj_m = {'Nov-1':-0.38, 'Nov-2':-0.30, 'Dec-1':-0.19, 'Dec-2':-0.10,
					'Jan-1':0, 'Jan-2':0.1, 'Feb-1':0.19, 'Feb-2':0.25, 
					'Mar-1':0.30, 'Mar-2':0.38}
	elif timescales == 'pentad':
		c_mons = ['Sep-1', 'Sep-2', 'Sep-3', 'Sep-4', 'Sep-5', 'Sep-6', 'Oct-1', 'Oct-2', 'Oct-3', 'Oct-4', 'Oct-5', 'Oct-6',
				  'Nov-1', 'Nov-2', 'Nov-3', 'Nov-4', 'Nov-5', 'Nov-6', 'Dec-1', 'Dec-2', 'Dec-3', 'Dec-4', 'Dec-5', 'Dec-6',
		          'Jan-1', 'Jan-2', 'Jan-3', 'Jan-4', 'Jan-5', 'Jan-6', 'Feb-1', 'Feb-2', 'Feb-3', 'Feb-4', 'Feb-5', 'Feb-6',
				  'Mar-1', 'Mar-2', 'Mar-3', 'Mar-4', 'Mar-5', 'Mar-6']
		x_adj_m = {i:0 for i in c_mons}
	elif timescales == 'DJF_pentad':
		c_mons = ['DJF-2p', 'DJF-1g', 'DJF']
		x_adj_m = {i:0 for i in c_mons}
	else:
		raise ValueError('No such timescales')
	
	colorbar_vmax = 0.8
	mons = list(org_dri_dict.keys())
	vars = list(org_dri_dict[mons[-1]].keys())
		
	x_cor_label = ['IR', 'Snow', 'ICE', 'THF', 'THFnew', 'SSLP', 'SLP', 'URALS', 'SLPA', 'SLPL', 'Z500','z500u', 'Z500bks', 'wave1', 'VTstar', 'V*T*1', 'V*T*2', 'V*T*123', 
					'VT0to120', 'VT0to60', 'VT0to180', 'V*T*', 'Z50', 'SPV', 'Z10e1', 'Z10e2', 'Z10e3',
					'U10', 'Z10', 'NAO', 'NAO_SLPbased', 'NAOM', 'VQ70', 'AO', 'OHT', 'R1', 'R2', 'R3',
					'PSLP', 'Aair', 'GBI',
					'Heatwave', 'Prepcipitation', 'Nino34', 'Z200']

	x_cor_label = [i for i in x_cor_label if i in vars]
	x_cor = {var: i+1 for i, var in enumerate(x_cor_label)}
	x_cor_num = len(x_cor_label) + 1
	x = np.arange(1, x_cor_num)

	y_cor = {var:i+1 for i, var in enumerate(c_mons)} #also work for pentads timescales. c_mon is used. Not mov_n
	y_cor_label = c_mons
	y_cor_num = len(c_mons) + 1
	y = np.arange(1, len(c_mons) + 1)

	if ax_all == None:
		plt.close()
		if (timescales == 'monthly'):
			fig, (ax1,ax2,ax3) = plt.subplots(3,1, figsize=(12,27))
		elif timescales in ['submonthly', 'pentad']: 
			fig, (ax1,ax2,ax3) = plt.subplots(3,1, figsize=(12,37))
		else: 
			fig, (ax1,ax2,ax3) = plt.subplots(3,1, figsize=(12,27))
	else: # Design for fig3_CEN.py so that we can make two CEN plots together
		ax1, ax2, ax3 = ax_all
	
	if only_keep_ax1==True:
		fig, ax1 = plt.subplots(1,1, figsize=(12,10)) # Just to crop out the first ax (ax1)
		#fig, ax1 = plt.subplots(1,1, figsize=(12,15)) # For the raw half-monthly

	yticks_fs = {'monthly':30, 'submonthly':30, 'pentad':15}
	yticks_fs['DJF_pentad'] = 30 
	for ax in [ax1, ax2, ax3]:
		ax.set_xticks(x)
		ax.set_xticklabels(x_cor_label, fontsize=30)
		ax.set_yticks(y)
		ax.set_yticklabels(y_cor_label, fontsize=yticks_fs[timescales])
		ax.set_xlim(0.5, x_cor_num-0.5)
		ax.set_ylim(0.5, y_cor_num)
		#ax.set_ylim(0.5, y_cor_num+3)
		ax.set_axisbelow(True)
		ax.spines['right'].set_visible(False)
		ax.spines['top'].set_visible(False)
		ax.spines['left'].set_visible(False)
		ax.spines['bottom'].set_visible(False)
		ax.tick_params(axis='both', which='both',length=0)
	
	ax2.set_ylim(0.5, x_cor_num-0.5)
	ax2.set_xticks(np.insert(x,0,0)+0.5, minor=True);
	ax2.set_yticks(np.insert(x,0,0)+0.5, minor=True);
	ax2.set_yticklabels([i + ' ' for i in x_cor_label], fontsize=30)
	ax2.set_xticklabels(x_cor_label, fontsize=30)
	ax2.grid(which='minor', color='grey', linewidth=0.5)
	ax2.fill_between(np.arange(0, x_cor_num+1), np.repeat(0, x_cor_num+1), np.arange(0, x_cor_num+1), facecolor='#D8D8D8')

	dot_s = {'monthly':150, 'submonthly':150, 'pentad':50}
	dot_s['DJF_pentad']=150
	# the dot in both monthly and weekly timescales
	if timescales in ['monthly', 'submonthly', 'pentad', 'DJF_pentad']:
		for x_pos in range(1, x_cor_num):
			for y_pos in range(1, len(y)+1):
				for ax in [ax1, ax3]:
					ax.scatter(x_pos, y_pos, color='#F5F5F5', alpha = 1, s=dot_s[timescales])

	mapcolors = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
	cmap= matplotlib.colors.ListedColormap(mapcolors)
	cNorm  = matplotlib.colors.Normalize(vmin = -colorbar_vmax, vmax=colorbar_vmax)
	scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm,cmap=cmap)
	curvature = .25
	
	if timescales in ['monthly']:
		astyle="simple,tail_width=5,head_width=25,head_length=25"; akw = dict(arrowstyle=astyle)	
	elif (timescales=='monthly') & half_to_full_mon:
		astyle="simple, tail_width=5,head_width=20,head_length=25"; akw = dict(arrowstyle=astyle)	
	elif timescales in ['submonthly']:
		astyle="simple,tail_width=5,head_width=20,head_length=20"; akw = dict(arrowstyle=astyle)	
	elif timescales =='pentad': # compensate the difference between two different timescales plot
		astyle="simple, tail_width=1,head_width=7,head_length=7"; akw = dict(arrowstyle=astyle)

	##### Plotting original parents. Additonal parents dict has another format
	for mon in mons:		
		for var in vars:

			### Plotting ordinary variables
			dris_lags = list(org_dri_dict[mon][var].keys())
			
			if len(dris_lags) == 0: # It has not driver
				continue
		
			for dri, lag in dris_lags:
				beta = org_dri_dict[mon][var][(dri, lag)]
				colorVal = scalarMap.to_rgba(beta)
				a_alpha = 1
				zorder = 10
				text_color = 'black'
				t_size = 30; t_alpha= 1
				edge_c = colorVal

				if pw_hl & ((timescales == 'monthly') | (timescales == 'submonthly')):
					pathways_selected = pathways_grey(mon, var, dri, lag, timescales, half_to_full_mon)
					if not pathways_selected: # It returns True if it is the seleted pathway
						colorVal, a_alpha, zorder= '#F5F5F5', 1, 0
						edge_c = '#F5F5F5'
						text_color = '#F5F5F5'
						t_size=2
					else: # For the selected pathway
						pass
	
				# For non-autocorrelation. Plotting arrows
				if dri != var:
					a1 = matplotlib.patches.FancyArrowPatch((x_cor[dri], y_cor[mon]-lag), (x_cor[var], y_cor[mon]), 
							connectionstyle="arc3,rad=%f" %curvature, color=colorVal, alpha=a_alpha,zorder=zorder, ec=edge_c, **akw) # (x1, y1), (x2, y2)
					ax1.add_patch(a1)
					
					if 'org_bs_betas' in aargus:
						mid_point = a1.get_path().to_polygons(None)[0][1]
						percentage = len(aargus['org_bs_betas'][mon][var][(dri, lag)]) * 1.0 / (aargus['bs_ss'] * 1.0) * 100
						ax1.text(mid_point[0], mid_point[1], str(int(percentage)) + '%', fontsize=t_size, alpha=t_alpha, zorder = 15, color=text_color)

				else: # Forr auto-coorelation. Ther vertical bar
					if lag == 0: # It only happens for submonthly-to-monthly plots
						continue
					ax1.plot([x_cor[dri], x_cor[var]], [y_cor[mon]-lag + 0.2, y_cor[mon] - 0.2], lw=7, linestyle = '-', 
							color=colorVal, alpha=a_alpha, zorder = zorder)
					
					if 'org_bs_betas' in aargus:
						percentage = len(aargus['org_bs_betas'][mon][var][(dri, lag)]) * 1.0 / (aargus['bs_ss'] * 1.0) * 100
						ax1.text(x_cor[dri], (2*y_cor[mon]-lag)/2.0, str(int(round(percentage))) + '%', fontsize=t_size, 
										alpha=t_alpha, zorder = 15, color=text_color)
	

	##### Plotting additional parents. Additonal parents dict has another format
	# separate from the original one because there is a continue which will also skil this loop if there is no ordinary casual driver
	y_adj_m = -0.2
	for mon in mons:		
		for var in vars:

			dris_lags = list(add_dri_dict[mon][var].keys())
			if len(dris_lags) == 0: # It has not driver
				continue
	
			for dri, lag in dris_lags:
				beta = add_dri_dict[mon][var][(dri, lag)]
				colorVal = scalarMap.to_rgba(beta)
				
				# It shows nothing because when URALS is var, V*T* is its driver. This has stronger beta coeff.
				# This only helps the plottign in my figure - should delete in the future. / or brings it back to the upper function
				if (timescales == 'monthly') & (pw_hl == True):
					if (mon == 'Dec') & (var == 'V*T*') & (dri == 'URALS') & (lag == 0):
						ax1.plot([x_cor[dri], x_cor[var]], [y_cor[mon], y_cor[mon]], c=colorVal, lw=7, zorder=10)
						if 'add_bs_betas' in aargus:
							percentage = len(aargus['add_bs_betas'][mon][var][(dri, lag)]) * 1.0 / (aargus['bs_ss'] * 1.0) * 100
							ax1.text((x_cor[dri]+x_cor[var])/2.0, (y_cor[mon]+y_cor[mon])/2.0, str(int(percentage)) + '%', fontsize=30, alpha=1, zorder = 15, color='black')
				
				# Plotting simultaneous relationship
				if lag == 0:
					if (var, 0) in add_dri_dict[mon][dri]:
						if abs(beta) < abs(add_dri_dict[mon][dri][(var,0)]):
							continue # If that simultaneous relationship is weaker than looking from the other way round, skip this loop
					ax2_x = x_cor[var]
					ax2_y = x_cor[dri]
					if ax2_x < ax2_y: # if y > x, normal plot
						pass
					else:
						ax2_x, ax2_y = ax2_y, ax2_x
					#ax2.scatter(ax2_x + x_adj_m[mon], ax2_y + y_adj_m, color = colorVal, marker=r"$ {} $".format(str(mon)[0]), s=350)
					if timescales == 'monthly':
						ax2.text(ax2_x + x_adj_m[mon] -0.08, ax2_y + y_adj_m, r'$%s$' %mon[0], color = colorVal, fontsize=25, fontweight='bold')
					elif timescales in ['submonthly', 'pentad']:
						ax2.text(ax2_x + x_adj_m[mon] -0.08, ax2_y + y_adj_m, r'$%s_%s$' %(mon[0], mon[4]), color = colorVal, fontsize=10, fontweight='bold')

				else: # if lag is not 0. Additional driver for panel 3
					a1 = matplotlib.patches.FancyArrowPatch((x_cor[dri], y_cor[mon]-lag), (x_cor[var], y_cor[mon]), 
										connectionstyle="arc3,rad=%f" %curvature, color=colorVal, alpha = 1, **akw)
					ax3.add_patch(a1)
		
	if ax_all == None: # It generate the plot by itself if not passing the plotting axis
		if only_keep_ax1:
			cbaxes1 = fig.add_axes([0.2, 0.02, 0.6, 0.025])
			cbaxes2 = fig.add_axes([0.2, 0.02, 0.6, 0.025])
			cbaxes3 = fig.add_axes([0.2, 0.02, 0.6, 0.025])
		else: # Keep the normal 3 ax
			cbaxes1 = fig.add_axes([0.10, 0.64, 0.75, 0.01])
			cbaxes2 = fig.add_axes([0.10, 0.35, 0.75, 0.01])
			cbaxes3 = fig.add_axes([0.10, 0.07, 0.75, 0.01])
			
		for cba in [cbaxes1, cbaxes2, cbaxes3]:
			cb1 = matplotlib.colorbar.ColorbarBase(cba, cmap=cmap, norm=cNorm, orientation='horizontal',
						ticks=np.linspace(-colorbar_vmax, colorbar_vmax, 9))
			cb1.ax.set_xticklabels([round(i,1) for i in np.linspace(-colorbar_vmax, colorbar_vmax, 9)])
			cb1.set_label('Beta coefficient', fontsize=20)
			cb1.ax.tick_params(labelsize=20)

		plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.35)
		plt.savefig('../graphs/%s_%s_%s_%s_%s.png' %(dt.date.today(), f_add_n, sig_lev, timescales, lag_num),bbox_inches='tight', dpi=100)
		plt.savefig('../graphs/%s_%s_%s_%s_%s.pdf' %(dt.date.today(), f_add_n, sig_lev, timescales, lag_num),bbox_inches='tight', dpi=600)


def pathways_grey(mon, var, dri, lag, timescales, half_to_full_mon=False):

	if (timescales == 'monthly') & (half_to_full_mon==False):
		if ((mon == 'Dec' and  var == 'URALS' and  dri == 'ICE') or # stratospheric pathway
			(mon == 'Jan' and var == 'SPV' and dri == 'V*T*') or 
			(mon == 'Feb' and var == 'NAO' and dri == 'SPV')

			#(mon == 'Jan' and var == 'SPV' and dri == 'URALS' and lag==1) or 
			#(mon == 'Jan' and var == 'V*T*' and dri == 'URALS') or
			#(mon == 'Feb' and var == 'SPV' and dri == 'V*T*') or 
			#(mon == 'Mar' and var == 'NAO' and dri == 'SPV')
			):
			return True
		else:
			return False

	
	if timescales == 'submonthly': # For the halfmonthly CEN without aggrgation
		if ((mon == 'Dec-2' and var == 'V*T*' and dri == 'THF' and lag==3) or
			(mon == 'Dec-2' and var == 'URALS' and dri == 'ICE') or 
			(mon == 'Nov-1' and var == 'THF' and dri == 'ICE') or 
			(mon == 'Feb-1' and var == 'NAO' and dri == 'SPV') or 
			(mon == 'Jan-2' and var == 'SPV' and dri == 'V*T*') or 
			(mon == 'Jan-1' and var == 'SPV' and dri == 'V*T*') or 
			(mon == 'Jan-1' and var == 'V*T*' and dri == 'URALS') or 
			(mon == 'Jan-2' and var == 'SPV' and dri =='SPV') or 
			(mon == 'Jan-2' and var == 'SPV' and dri =='SPV') or 

			(mon == 'Jan-2' and var == 'SPV' and dri =='URALS') or 
			(mon == 'Jan-1' and var == 'V*T*' and dri =='ICE') or 
			(mon == 'Feb-1' and var == 'NAO' and dri =='V*T*') or
			(mon == 'Dec-2' and var == 'SPV' and dri =='THF') or
			(mon == 'Jan-1' and var == 'SPV' and dri =='SPV')
			):
			return True
		else:
			return False

	if (timescales == 'monthly') and (half_to_full_mon== True): # For the halfmonthly aggregated CEN
		if ((mon == 'Dec' and var == 'V*T*' and dri == 'THF') or
			(mon == 'Dec' and var == 'URALS' and dri == 'ICE') or 
			(mon == 'Nov' and var == 'THF' and dri == 'ICE' and lag == 1) or 
			(mon == 'Feb' and var == 'NAO' and dri == 'SPV') or 
			(mon == 'Jan' and var == 'SPV' and dri == 'V*T*') or 
			(mon == 'Jan' and var == 'V*T*' and dri == 'ICE') or 
			(mon == 'Jan' and var == 'SPV' and dri == 'URALS') or 
			(mon == 'Dec' and var == 'SPV' and dri == 'THF') or 
			(mon == 'Jan' and var == 'SPV' and dri == 'SPV') or
			(mon == 'Feb' and var == 'SPV' and dri == 'SPV') or
			(mon == 'Mar' and var == 'NAO' and dri == 'SPV')
			):
			return True
		else:
			return False

	if False:
		if (timescales == 'monthly') and (half_to_full_mon== True): # Decide for shorter period (1979to2010) and sig_lev=0.05
			if ((mon == 'Nov' and var == 'SPV' and dri == 'ICE') or
				(mon == 'Dec' and var == 'V*T*' and dri == 'THF') or
				(mon == 'Dec' and var == 'NAO' and dri == 'SPV') or
				(mon == 'Dec' and var == 'URALS' and dri == 'ICE') or 
				(mon == 'Nov' and var == 'THF' and dri == 'ICE' and lag == 1) or 
				(mon == 'Feb' and var == 'NAO' and dri == 'SPV') or 
				(mon == 'Jan' and var == 'SPV' and dri == 'URALS') or 
				#(mon == 'Jan' and var == 'V*T*' and dri == 'ICE') or 
				(mon == 'Jan' and var == 'SPV' and dri == 'V*T*' and lag==1)):
				return True
			else:
				return False
	if False:
		if (timescales == 'monthly') & (half_to_full_mon==False): # A figure for revision showing all tropospheric and stratospheric pathway
			if ((mon == 'Dec' and  var == 'URALS' and  dri == 'ICE') or 
				(mon == 'Jan' and var == 'V*T*' and dri == 'URALS') or
				(mon == 'Feb' and var == 'SPV' and dri == 'V*T*') or 
				(mon == 'Mar' and var == 'NAO' and dri == 'SPV') or
				(mon == 'Nov' and var == 'THF' and dri == 'ICE') or
				(mon == 'Jan' and var == 'URALS' and dri == 'THF') or 
				(mon == 'Feb' and var == 'URALS' and dri == 'URALS') or 
				(mon == 'Mar' and var == 'NAO' and dri == 'URALS') or
				(mon == 'Feb' and var == 'URALS' and dri == 'ICE') or
				(mon == 'Jan' and var == 'SPV' and dri == 'URALS' and lag==1) or 
				(mon == 'Feb' and var == 'NAO' and dri == 'SPV') or
				(mon == 'Jan' and var == 'SPV' and dri == 'V*T*')
				):
				return True
			else:
				return False

