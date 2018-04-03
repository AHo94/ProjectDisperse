import matplotlib.pyplot as plt
import numpy as np

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

def Call_plot_sameX(xdata, ydata, xlabel, ylabel, legend, style='-', **kwargs):
	""" 
	Calls plotting, xdata to be common lengths/bins for all ydata
	Extra parameter logscal determines the plotting method
	normal = plot(), logx = semilogx(), etc.
	"""
	#if len(xdata) != len(ydata[0]):
	#	raise ValueError("xdata and ydata not of same length!")
	figure = plt.figure()
	plt.gcf().set_size_inches((8*s_variable, 6*s_variable))
	for kw in kwargs:
		if kw == 'color':
			colors = kwargs[kw]
			for i in range(len(ydata)):
				plt.plot(xdata, ydata[i], style, color=colors[i])
		else:
			for i in range(len(ydata)):
				plt.plot(xdata, ydata[i], style)
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
