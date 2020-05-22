import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits import mplot3d

import matplotlib.gridspec as gridspec
import seaborn as sns

def plot(plot_type,x,y=None,yerr=None,xerr=None,ax=None,x_lab='',y_lab='',fontsize=18,figsize=(12,12),
		 x_name=None,y_name=None,label_name=None,**kwargs):
	if ax is None and plot_type != 'joint':
		fig=plt.figure(figsize=figsize)
		ax=fig.gca()

	if plot_type=='scatter':
		ax.scatter(x,y,**kwargs)
	elif plot_type=='plot':
		ax.plot(x,y,**kwargs)
	elif plot_type=='errorbar':
		ax.errorbar(x,y,xerr=xerr,yerr=yerr,**kwargs)
	elif plot_type=='hist':
		ax.hist(x,**kwargs)
	elif plot_type=='joint':
		g=multivariateGrid(x_name, y_name, label_name, df=x)
		fig=g.ax_joint.__dict__['figure']
		ax=fig.gca()
		fig.set_size_inches(figsize[0],figsize[1])
	else:
		raise RuntimeError('What plot are you trying to do.')
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(16)
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(16)
	ax.set_xlabel(x_lab,fontsize=fontsize)
	ax.set_ylabel(y_lab,fontsize=fontsize)
	return(ax)

def plot3d(x,y,z,ax=None,x_lab='',y_lab='',fontsize=18,figsize=(12,12),
		 x_name=None,y_name=None,label_name=None,cmap='viridis',edgecolor='none',**kwargs):
	fig=plt.figure(figsize=figsize)
	ax=fig.gca()

	xsurf = np.outer(x, np.ones(len(y)))
	ysurf = np.outer(y,np.ones(len(x))).T
	ax.plot_surface(xsurf, ysurf, z,cmap=cmap, edgecolor=edgecolor,**kwargs)
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize(16)
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize(16)
	ax.set_xlabel(x_lab,fontsize=fontsize)
	ax.set_ylabel(y_lab,fontsize=fontsize)
	return(ax)

def split_plot(ax,plot_type,x,y=None,yerr=None,xerr=None,x_lab='',y_lab='',xticks=False,fontsize=18,split_size='50%',**kwargs):
	ax_divider = make_axes_locatable(ax)
	ax_ml = ax_divider.append_axes("bottom", size=split_size, pad=.2)
	ticks=[]
	for tick in ax_ml.xaxis.get_major_ticks():
		ticks.append('')
		tick.label.set_fontsize(16)
	if not xticks:
		ax.set_xticklabels(ticks)
	for tick in ax_ml.yaxis.get_major_ticks():
		tick.label.set_fontsize(16)
	if plot_type=='scatter':
		ax_ml.scatter(x,y,**kwargs)
	elif plot_type=='plot':
		ax_ml.plot(x,y,**kwargs)
	elif plot_type=='errorbar':
		ax_ml.errorbar(x,y,xerr=xerr,yerr=yerr,**kwargs)
	elif plot_type=='hist':
		ax_ml.hist(x,**kwargs)
	else:
		raise RuntimeError('What plot are you trying to do.')
	for tick in ax_ml.xaxis.get_major_ticks():
		tick.label.set_fontsize(16)
	for tick in ax_ml.yaxis.get_major_ticks():
		tick.label.set_fontsize(16)
	ax_ml.set_xlabel(x_lab,fontsize=fontsize)
	ax_ml.set_ylabel(y_lab,fontsize=fontsize)
	return(ax,ax_ml)

def grid_plot(grid_x,grid_y,figsize=(12,12)):
	fig=plt.figure(figsize=figsize)
	gs = gridspec.GridSpec(grid_x, grid_y)
	return(fig,gs)


def multivariateGrid(col_x, col_y, col_k, df, k_is_color=False, 
					 scatter_alpha=.5,global_hist=False):
	def colored_scatter(x, y, c=None):
		def scatter(*args, **kwargs):
			args = (x, y)
			if c is not None:
				kwargs['c'] = c
			kwargs['alpha'] = scatter_alpha
			plt.scatter(*args, **kwargs)
			ax=plt.gca()

		return scatter
	g = sns.JointGrid(
		x=col_x,
		y=col_y,
		data=df
	)

	color = None
	legends=[]
	for name, df_group in df.groupby(col_k):
		legends.append(name)
		if k_is_color:
			color=name

		g.plot_joint(
			colored_scatter(df_group[col_x],df_group[col_y],color),
		)
		sns.distplot(
			df_group[col_x].values,
			ax=g.ax_marg_x,
			color=color,
		)
		sns.distplot(
			df_group[col_y].values,
			ax=g.ax_marg_y,
			color=color,            
			vertical=True
		)
	# Do also global Hist:
	if global_hist:
		sns.distplot(
			df[col_x].values,
			ax=g.ax_marg_x,
			color='grey'
		)
		sns.distplot(
			df[col_y].values.ravel(),
			ax=g.ax_marg_y,
			color='grey',
			vertical=True
		)
	plt.legend(legends,fontsize=20)
	return(g)

