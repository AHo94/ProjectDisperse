import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec

s_variable = 0.7

def Call_plot(xdata, ydata, xlabel, ylabel, legend, style='-', **kwargs):
	""" 
	Calls plotting, xdata and ydata must be same size 
	Extra parameter logscal determines the plotting method
	normal = plot(), logx = semilogx(), etc.
	"""
	if len(xdata) != len(ydata):
		raise ValueError("xdata and ydata not of same length!")
	figure = plt.figure()
	plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
	for kw in kwargs:
		if kw == 'color':
			colors = kwargs[kw]
			for i in range(len(xdata)):
				plt.plot(xdata[i], ydata[i], style, color=colors[i])
		else:
			for i in range(len(xdata)):
				plt.plot(xdata[i], ydata[i], style)
			
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.legend(legend)
	for kw in kwargs:
		if kw == 'title':
			plt.title(title)
		if kw == 'logscale':
			if kwargs[kw] == 'logx':
				plt.xscale('log')
			elif kwargs[kw] == 'logy':
				plt.yscale('log')
			elif kwargs[kw] == 'loglog':
				plt.xscale('log')
				plt.yscale('log')
			else:
				raise ValueError("Log scale parameter inserted wrong!")
	return figure

def Call_plot_sameX_OLD(xdata, ydata, xlabel, ylabel, legend, style='-', **kwargs):
	""" 
	Calls plotting, xdata to be common lengths/bins for all ydata
	Extra parameter logscal determines the plotting method
	normal = plot(), logx = semilogx(), etc.
	"""
	#if len(xdata) != len(ydata[0]):
	#	raise ValueError("xdata and ydata not of same length!")
	figure = plt.figure()
	plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
	Do_error = 0
	for kw in kwargs:
		if kw == 'error':
			Errors = kwargs[kw]
			Do_error = 1
		break
	for kw in kwargs:
		if kw == 'color':
			colors = kwargs[kw]
			if Do_error:
				for i in range(len(ydata)):
					plt.plot(xdata, ydata[i], style, color=colors[i])
					plt.fill_between(xdata, ydata[i]-Errors[i], ydata[i]+Errors[i], alpha=0.3, facecolor=colors[i])
			else:
				for i in range(len(ydata)):
					plt.plot(xdata, ydata[i], style, color=colors[i])
			break
		else:
			for i in range(len(ydata)):
				plt.plot(xdata, ydata[i], style)
			break
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.legend(legend)
	if np.min(xdata) > 1e12:
		xfmt = plt.ScalarFormatter()
		xfmt.set_powerlimits((0,0))
		plt.gca().xaxis.set_major_formatter(xfmt)

	for kw in kwargs:
		if kw == 'title':
			plt.title(title)
		if kw == 'logscale':
			if kwargs[kw] == 'logx':
				plt.xscale('log')
			elif kwargs[kw] == 'logy':
				plt.yscale('log')
			elif kwargs[kw] == 'loglog':
				plt.xscale('log')
				plt.yscale('log')
			else:
				raise ValueError("Log scale parameter inserted wrong!")
	return figure

def Call_plot_sameX(xdata, ydata, xlabel, ylabel, legend, colors, **kwargs):
	""" 
	Calls plotting, xdata to be common lengths/bins for all ydata
	Extra parameter logscal determines the plotting method
	normal = plot(), logx = semilogx(), etc.
	"""
	do_fill_between = False
	set_ylimits = False
	set_xlimits = False
	New_figure_size = False
	Remove_y_ticks = True
	Nrows = 1
	Ncols = len(ydata)
	fb_alpha = 0.3
	anchor_legend = True
	titles = []
	Change_xscales = False
	Change_yscales = False
	Plot_LCDMDiff = False
	linestyles = ['-']*len(ydata[0])
	error = np.array([])
	for kw in kwargs:
		if kw == 'fillbetween':
			do_fill_between = kwargs[kw]
		elif kw == 'xlim':		# Sets xlimit if called
			set_xlimits = True
			xlims = kwargs[kw]
			if type(xlims) != tuple:
				raise ValueError("Keyword argument xlims must be a tuple!")
			if not xlims[0] and not xlims[1]:
				set_xlimits = False
		elif kw == 'ylim':		# Sets ylimit if called
			set_ylimits = True
			ylims = kwargs[kw]
			if type(ylims) != tuple:
				raise ValueError("Keyword argument ylims must be a tuple!")
			if not ylims[0] and not ylims[1]:
				set_ylimits = False
		elif kw == 'figsize':   	# Figure size of plot
			figure_size = kwargs[kw]
			New_figure_size = True
		elif kw == 'NoYticks':   # Removes Y-ticks of the latter plots, removes clutted.
			Remove_y_ticks = kwargs[kw]
		elif kw == 'fb_alpha':   # Alpha used for fill_between
			fb_alpha = kwargs[kw]
			if type(fb_alpha) != float:
				raise ValueError("Keyword argument fb_alpha must be a float number!")
		elif kw == 'legend_anchor':   # Anchors the legend outside the figure if True
			anchor_legend = kwargs[kw]
		elif kw == 'title':    	# Titles above subplots
			titles = kwargs[kw]
			if type(titles) != list:
				raise ValueError("Keyword argument titles must be a list!")
		elif kw == 'xscale':	# Changes x plot scales based on input
			logXscale_name = kwargs[kw]
			Change_xscales = True
		elif kw == 'yscale':	# Changes y plot scales based on input
			logYscale_name = kwargs[kw]
			Change_yscales = True
		elif kw == 'linestyles': 	# Choose linestyle of plotting
			linestyles = kwargs[kw]
		elif kw == 'reldiff':	# Plots a line that distinguishes the LCDM model as a zero line
			Plot_LCDMDiff = kwargs[kw] 
		elif kw == 'error':
			error = kwargs[kw]
		else:
			raise ValueError("Keyword argument " + kw + " not recognized!")
	if type(error) == list:
		if do_fill_between and not error:
			do_fill_between = False
	else:        
		if do_fill_between and not error.any():
			do_fill_between = False
		#print "Warning: fill_between enabled, but no error!"
	if New_figure_size:
		figure = plt.figure(figsize=figure_size)
		ax = plt.axes()
	else:
		figure = plt.figure()
		ax = plt.axes()
	plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
	if do_fill_between:
		if Plot_LCDMDiff:
			plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM')
		for i in range(len(ydata)):
			plt.plot(xdata, ydata[i], label=legend[i], color=colors[i], linestyle=linestyles[i])
			plt.fill_between(xdata, ydata[i]-error[i], ydata[i]+error[i], alpha=fb_alpha, facecolor=colors[i])
	else:
		if Plot_LCDMDiff:
			plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM')
		for i in range(len(ydata)):
			plt.plot(xdata, ydata[i], label=legend[i], color=colors[i], linestyle=linestyles[i])
	if titles:
		plt.title(titles, fontsize=10)
	if set_xlimits:
		plt.xlim(xlims)
	if set_ylimits:
		plt.ylim(ylims)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	if np.min(xdata) > 1e12:
		xfmt = plt.ScalarFormatter()
		xfmt.set_powerlimits((0,0))
		plt.gca().xaxis.set_major_formatter(xfmt)
	if anchor_legend:
		ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
	else:
		plt.legend(legend)
	if Change_xscales:
		plt.xscale(logXscale_name)
	if Change_yscales:
		plt.yscale(logYscale_name)
	plt.tight_layout()
	return figure

def Plot_differences(self, xdata, ydata, xlabel, ylabel, legend, logscale='normal', style='-', title='None', diff='rel'):
	""" Plots relative or absolute differences. Base data assumed to be the first element of xdata and ydata """
	if len(xdata) != len(ydata):
		raise ValueError("xdata and ydata not of same length!")
		
	deltas = []
	if diff == 'rel':
		for i in range(1, len(xdata)):
			deltas.append(self.relative_deviation(ydata, i))
	elif diff == 'abs':
		for i in range(1, len(xdata)):
			deltas.append(ydata[i] - ydata[0])
	else:
		raise ValueError("Argument diff not properly set! Use either diff='rel' or diff='abs'")

	figure = plt.figure()
	plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
	if logscale == 'normal':
		plt.plot(xdata[0], np.zeros(len(xdata[0])), style)
		for i in range(1, len(xdata)):
			plt.plot(xdata[i], deltas[i-1], style)
	elif logscale == 'logx':
		plt.semilogx(xdata[0], np.zeros(len(xdata[0])), style)
		for i in range(1, len(xdata)):
			plt.semilogx(xdata[i], deltas[i-1], style)
	elif logscale == 'logy':
		plt.semilogy(xdata[0], np.zeros(len(xdata[0])), style)
		for i in range(1, len(xdata)):
			plt.semilogy(xdata[i], deltas[i-1], style)
	elif logscale == 'loglog':
		plt.loglog(xdata[0], np.zeros(len(xdata[0])), style)
		for i in range(1, len(xdata)):
			plt.loglog(xdata[i], deltas[i-1], style)
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.legend(legend)
	plt.title('') if title == 'None' else plt.title(title)
	return figure

def Plot_errobar_sameX(self, xdata, ydata, error, xlabel, ylabel, legend, logscale='normal', fill_between=False, diff=False, limit_yaxis=True, Lengths=True):
	if len(ydata) != len(error):
		raise ValueError("ydata and error data not of same length!")
	self.Check_ok_parameters(logscale)

	Max_yval = np.array([np.max(ydata[i]) for i in range(len(ydata))])
	Min_yval = np.array([np.min(ydata[i]) for i in range(len(ydata))])
	figure = plt.figure()
	plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
	ax = plt.axes()
	if diff:
		plt.plot(xdata, np.zeros(len(xdata)), label='$\Lambda$CDM')
		if fill_between:
			plt.fill_between(xdata, np.zeros(len(xdata)), np.zeros(len(xdata)), alpha=0)
	if fill_between:
		for i in range(len(ydata)):
			plt.plot(xdata, ydata[i], label=legend[i])
			plt.fill_between(xdata, ydata[i]-error[i], ydata[i]+error[i], alpha=0.3)
		if limit_yaxis:
			plt.ylim((-1.0, 1))
	else:
		for i in range(len(ydata)):
			plt.errorbar(xdata, ydata[i], error[i], label=legend[i])
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)

	if Lengths:
		plt.xlim((1, np.max(xdata)))
		
	ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
	if logscale == 'logx':
		plt.xscale('log')
	elif logscale == 'logy':
		plt.yscale('log')
	elif logscale == 'loglog':
		plt.xscale('log')
		plt.yscale('log')
	#plt.legend(legend)
	return figure

def Do_subplots_sameX(xdata, ydata, xlabel, ylabel, legend, colors, error=[], **kwargs):
	""" Plots subplots for multiple ydata sets but all using the same xdata """
	# Default arguments
	do_fill_between = False
	set_ylimits = False
	set_xlimits = False
	New_figure_size = False
	Remove_y_ticks = True
	Nrows = 1
	Ncols = len(ydata)
	fb_alpha = 0.3
	anchor_legend = True
	titles = []
	Change_xscales = False
	Change_yscales = False
	Plot_LCDMDiff = False
	linestyles = ['-']*len(ydata[0])
	# Iterate through keyword arguments to check if anything special will happen. Changes some default arguments
	for kw in kwargs:
		if kw == 'fillbetween':
			do_fill_between = kwargs[kw]
		elif kw == 'xlim':		# Sets xlimit if called
			set_xlimits = True
			xlims = kwargs[kw]
			if type(xlims) != tuple:
				raise ValueError("Keyword argument xlims must be a tuple!")
			if not xlims[0] and not xlims[1]:
				set_xlimits = False
		elif kw == 'ylim':		# Sets ylimit if called
			set_ylimits = True
			ylims = kwargs[kw]
			if type(ylims) != tuple:
				raise ValueError("Keyword argument ylims must be a tuple!")
			if not ylims[0] and not ylims[1]:
				set_ylimits = False
		elif kw == 'figsize':   	# Figure size of plot
			figure_size = kwargs[kw]
			New_figure_size = True
		elif kw == 'rowcol':   	# Number of rows and columns for subplots.
			row_and_col = kwargs[kw]
			if type(row_and_col) != list:
				raise ValueError("Keyword argument rowcol must be a list!")
			Nrows = row_and_col[0]
			Ncols = row_and_col[1]
		elif kw == 'NoYticks':   # Removes Y-ticks of the latter plots, removes clutted.
			Remove_y_ticks = kwargs[kw]
		elif kw == 'fb_alpha':   # Alpha used for fill_between
			fb_alpha = kwargs[kw]
			if type(fb_alpha) != float:
				raise ValueError("Keyword argument fb_alpha must be a float number!")
		elif kw == 'legend_anchor':   # Anchors the legend outside the figure if True
			anchor_legend = kwargs[kw]
		elif kw == 'title':    	# Titles above subplots
			titles = kwargs[kw]
			if type(titles) != list:
				raise ValueError("Keyword argument titles must be a list!")
		elif kw == 'xscale':	# Changes x plot scales based on input
			logXscale_name = kwargs[kw]
			Change_xscales = True
		elif kw == 'yscale':	# Changes y plot scales based on input
			logYscale_name = kwargs[kw]
			Change_yscales = True
		elif kw == 'linestyles': 	# Choose linestyle of plotting
			linestyles = kwargs[kw]
		elif kw == 'reldiff':	# Plots a line that distinguishes the LCDM model as a zero line
			Plot_LCDMDiff = kwargs[kw]
		else:
			raise ValueError("Keyword argument " + kw + " not recognized!")

	# Quick checks of error data vs ydata and titles
	if not error and do_fill_between:
		print 'Warning: fillbetween is True but no error data found'
		do_fill_between = False
	elif error and do_fill_between:
		if len(error) != len(ydata):
			raise ValueError("Error data not the same size as ydata!")
	if titles:
		if len(titles) != len(ydata):
			raise ValueError("title list not the same length as ydata!")
	if len(legend) != len(ydata[0]):
		raise ValueError("Number of legends not the same as number of ydata!")

	if New_figure_size:
		figure = plt.figure(figsize=figure_size)
	else:
		figure = plt.figure()
	plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
	ax = plt.subplot(Nrows, Ncols, 1)
	if do_fill_between:
		for j in range(len(ydata)):
			if j > 0:
				ax = plt.subplot(Nrows,Ncols, j+1, sharey=ax, sharex=ax)
				plt.setp(ax.get_yticklabels(), visible=False) if Remove_y_ticks else plt.setp(ax.get_yticklabels(), visible=True)
			if Plot_LCDMDiff:
				#plt.plot(xdata, np.zeros(len(xdata)), color='b', label='$\Lambda$CDM', linestyle=(0, (3, 10, 1, 10, 1, 10)))
				plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM', linestyle='-')
			for i in range(len(ydata[j])):
				plt.plot(xdata, ydata[j][i], label=legend[i], color=colors[i], linestyle=linestyles[i])
				plt.fill_between(xdata, ydata[j][i]-error[j][i], ydata[j][i]+error[j][i], alpha=fb_alpha, facecolor=colors[i])
			if titles:
				plt.title(titles[j], fontsize=10)
	else:
		for j in range(len(ydata)):
			if j > 0:
				ax = plt.subplot(Nrows,Ncols, j+1, sharey=ax, sharex=ax)
				plt.setp(ax.get_yticklabels(), visible=False) if Remove_y_ticks else plt.setp(ax.get_yticklabels(), visible=True)
			if Plot_LCDMDiff:
				plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM')
			for i in range(len(ydata[j])):
				plt.plot(xdata, ydata[j][i], label=legend[i], color=colors[i], linestyle=linestyles[i])
			if titles:
				plt.title(titles[j], fontsize=10)
	if set_xlimits:
		plt.xlim(xlims)
	if set_ylimits:
		plt.ylim(ylims)

	if anchor_legend:
		ax.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5), ncol=1, fancybox=True)
	figure.text(0.5, 0, xlabel, ha='center', fontsize=10)
	figure.text(0, 0.5, ylabel, ha='center', va='center', rotation='vertical', fontsize=10)
	if Change_xscales:
		plt.xscale(logXscale_name)
	if Change_yscales:
		plt.yscale(logYscale_name)
	plt.tight_layout()
	return figure

def Do_gridspec_sameX(xdata, primaryY, secondaryY, xlabel, ylabel1, ylabel2, legend, colors, Primerror=[], Secerror=[], **kwargs):
	do_fill_between = False
	set_ylimits = False
	set_xlimits = False
	set_ylimits_diff = False
	set_xlimits_diff = False
	New_figure_size = False
	Remove_y_ticks = True
	Nrows = 3
	Ncols = len(primaryY)
	fb_alpha = 0.3
	anchor_legend = True
	titles = []
	Change_xscales = False
	Change_yscales = False
	Plot_LCDMDiff = False
	linestyles = ['-']*10
	# Iterate through keyword arguments to check if anything special will happen. Changes some default arguments
	for kw in kwargs:
		if kw == 'fillbetween':
			do_fill_between = kwargs[kw]
		elif kw == 'xlim':		# Sets xlimit if called
			set_xlimits = True
			xlims = kwargs[kw]
			if type(xlims) != tuple:
				raise ValueError("Keyword argument xlims must be a tuple!")
			if not xlims[0] and not xlims[1]:
				set_xlimits = False
		elif kw == 'ylim':		# Sets ylimit if called
			set_ylimits = True
			ylims = kwargs[kw]
			if type(ylims) != tuple:
				raise ValueError("Keyword argument ylims must be a tuple!")
			if not ylims[0] and not ylims[1]:
				set_ylimits = False
		elif kw == 'xlim_diff':		# Sets xlimit of reldiff plot if called
			set_xlimits_diff = True
			xlims_d = kwargs[kw]
			if type(xlims_d) != tuple:
				raise ValueError("Keyword argument xlims must be a tuple!")
			if not xlims_d[0] and not xlims_d[1]:
				set_xlimits_diff = False
		elif kw == 'ylim_diff':		# Sets ylimit of reldiff plot if called
			set_ylimits_diff = True
			ylims_d = kwargs[kw]
			if type(ylims_d) != tuple:
				raise ValueError("Keyword argument ylims must be a tuple!")
			if not ylims_d[0] and not ylims_d[1]:
				set_ylimits_diff = False
		elif kw == 'figsize':   	# Figure size of plot
			figure_size = kwargs[kw]
			New_figure_size = True
		elif kw == 'rowcol':   	# Number of rows and columns for subplots.
			row_and_col = kwargs[kw]
			if type(row_and_col) != list:
				raise ValueError("Keyword argument rowcol must be a list!")
			Nrows = row_and_col[0]
			Ncols = row_and_col[1]
		elif kw == 'NoYticks':   # Removes Y-ticks of the latter plots, removes clutted.
			Remove_y_ticks = kwargs[kw]
		elif kw == 'fb_alpha':   # Alpha used for fill_between
			fb_alpha = kwargs[kw]
			if type(fb_alpha) != float:
				raise ValueError("Keyword argument fb_alpha must be a float number!")
		elif kw == 'legend_anchor':   # Anchors the legend outside the figure if True
			anchor_legend = kwargs[kw]
		elif kw == 'title':    	# Titles above subplots
			titles = kwargs[kw]
			if type(titles) != list:
				raise ValueError("Keyword argument titles must be a list!")
		elif kw == 'xscale':	# Changes x plot scales based on input
			logXscale_name = kwargs[kw]
			Change_xscales = True
		elif kw == 'yscale':	# Changes y plot scales based on input
			logYscale_name = kwargs[kw]
			Change_yscales = True
		elif kw == 'linestyles': 	# Choose linestyle of plotting
			linestyles = kwargs[kw]
		elif kw == 'reldiff':	# Plots a line that distinguishes the LCDM model as a zero line
			Plot_LCDMDiff = kwargs[kw]
		else:
			raise ValueError("Keyword argument " + kw + " not recognized!")

	# Quick checks of error data vs ydata and titles
	#if len(primaryY) != secondaryY:
	#	raise ValueError("Primary and secondary Ydata not of same size!")
	if type(Primerror) == list and type(Secerror) == list:
		if (not Primerror and do_fill_between) and (not Secerror and do_fill_between):
			print 'Warning: fillbetween is True but no error data found'
			do_fill_between = False
	elif type(Primerror) != type(Secerror):
		raise ValueError("Error inputs not of same type!")
	else:
		if (not Primerror.any() and do_fill_between) and (not Secerror.any() and do_fill_between):
			print 'Warning: fillbetween is True but no error data found'
			do_fill_between = False

	if Primerror and do_fill_between:
		if len(Primerror) != len(primaryY):
			raise ValueError("Primary error data not the same size as primary ydata!")
	if Secerror and do_fill_between:
		if len(Secerror) != len(secondaryY):
			raise ValueError("Secondary error data not the same size as secondary ydata!")
	if titles:
		if len(titles) != len(primaryY):
			raise ValueError("title list not the same length as ydata!")
	if len(legend) != len(primaryY[0]):
		raise ValueError("Number of legends not the same as number of ydata!")

	if New_figure_size:
		figure = plt.figure(figsize=figure_size)
	else:
		figure = plt.figure(figsize=(8,6))
	plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
	gs = gridspec.GridSpec(Nrows, Ncols) 
	if do_fill_between:
		subfactor = 0
		for j in range(len(primaryY)):
			ax0 = plt.subplot(gs[0:Nrows-1, j]) if j == 0 else plt.subplot(gs[0:Nrows-1, j], sharey=ax0)
			plt.setp(ax0.get_xticklabels(), visible=False)
			plt.setp(ax0.get_yticklabels(), visible=False) if j > 0 else plt.setp(ax0.get_yticklabels(), visible=True)
			for i in range(len(primaryY[j])):
				plt.plot(xdata, primaryY[j][i], label=legend[i], color=colors[i], linestyle=linestyles[i])
			if Primerror:
				plt.fill_between(xdata, primaryY[j][i]-Primerror[j][i], primaryY[j][i]+Primerror[j][i], alpha=fb_alpha, facecolor=colors[i])
			if titles:
				plt.title(titles[j], fontsize=7)
			plt.ylabel(ylabel1) if j == 0 else plt.ylabel('')
			ax1 = plt.subplot(gs[Nrows-1, j], sharex=ax0) if j == 0 else plt.subplot(gs[Nrows-1, j], sharex=ax0, sharey=ax1)
			plt.setp(ax1.get_yticklabels(), visible=False) if j > 0 else plt.setp(ax1.get_xticklabels(), visible=True)
			if Plot_LCDMDiff:
				plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM', linestyle='-')
				subfactor = 1
			for i in range(len(secondaryY[j])):
				plt.plot(xdata, secondaryY[j][i], label=legend[i+subfactor], color=colors[i+subfactor], linestyle=linestyles[i+subfactor])
				plt.fill_between(xdata, secondaryY[j][i]-Secerror[j][i], secondaryY[j][i]+Secerror[j][i], alpha=fb_alpha, facecolor=colors[i+subfactor])
			#plt.xlabel(xlabel)
			plt.ylabel(ylabel2) if j == 0 else plt.ylabel('')
			if set_xlimits:
				ax0.set_xlim(xlims)
			if set_ylimits:
				ax0.set_ylim(ylims)
			if set_xlimits_diff:
				ax1.set_xlim(xlims_d)
			if set_ylimits_diff:
				ax1.set_ylim(ylims_d)
			if Change_xscales:
				plt.xscale(logXscale_name)
			if Change_yscales:
				plt.yscale(logYscale_name)
	else:
		subfactor = 0
		for j in range(len(primaryY)):
			ax0 = plt.subplot(gs[0:Nrows-1, j])
			plt.setp(ax0.get_xticklabels(), visible=False)
			plt.setp(ax0.get_yticklabels(), visible=False) if j > 0 else plt.setp(ax0.get_yticklabels(), visible=True)
			for i in range(len(primaryY[j])):
				plt.plot(xdata, primaryY[j][i], label=legend[i], color=colors[i], linestyle=linestyles[i])
			if titles:
				plt.title(titles[j], fontsize=7)
			plt.ylabel(ylabel1) if j == 0 else plt.ylabel('')
			ax1 = plt.subplot(gs[Nrows-1, j], sharex=ax0) if j == 0 else plt.subplot(gs[Nrows-1, j], sharex=ax0, sharey=ax1)
			plt.setp(ax1.get_yticklabels(), visible=False) if j > 0 else plt.setp(ax1.get_xticklabels(), visible=True)
			if Plot_LCDMDiff:
				plt.plot(xdata, np.zeros(len(xdata)), color='k', label='$\Lambda$CDM', linestyle='-')
				subfactor = 1
			for i in range(len(secondaryY[j])):
				plt.plot(xdata, secondaryY[j][i], label=legend[i+subfactor], color=colors[i+subfactor], linestyle=linestyles[i+subfactor])
			plt.ylabel(ylabel2) if j == 0 else plt.ylabel('')
			if set_xlimits:
				plt.xlim(xlims)
			if set_ylimits:
				plt.ylim(ylims)
			if Change_xscales:
				plt.xscale(logXscale_name)
			if Change_yscales:
				plt.yscale(logYscale_name)

	if anchor_legend:
		ax0.legend(loc = 'lower left', bbox_to_anchor=(1.0,0.5*(Nrows-1)/Nrows), ncol=1, fancybox=True)
	else:
		h,l=ax0.get_legend_handles_labels()
		ax0.legend(h,l)
	figure.text(0.5, 0, xlabel, ha='center', fontsize=10)
	#figure.text(0, 0.5, ylabel, ha='center', va='center', rotation='vertical', fontsize=10)
	plt.tight_layout()
	return figure