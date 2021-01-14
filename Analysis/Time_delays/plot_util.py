import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits import mplot3d
import numpy as np
import pickle
import matplotlib.lines as mlines

from astropy.table import Table
import matplotlib.gridspec as gridspec

def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

def plot(plot_type,x,y=None,yerr=None,xerr=None,ax=None,x_lab='',y_lab='',fontsize=18,figsize=(12,12),
         x_name=None,y_name=None,label_name=None,ticklabelsize=16,
         majorticklength=7,minorticklength=3,tickwidth=2,**kwargs):
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
        tick.label.set_fontsize(ticklabelsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(ticklabelsize)
    ax.set_xlabel(x_lab,fontsize=fontsize)
    ax.set_ylabel(y_lab,fontsize=fontsize)
    ax.tick_params('both',length=majorticklength,width=tickwidth)
    ax.tick_params('both',length=minorticklength,width=tickwidth,which='minor')
    return(ax)



def plot3d(x,y,z,fig=None,x_lab='',y_lab='',fontsize=18,figsize=(12,12),line_dict=dict(color='black',width=8),
         x_name='x',y_name='y',z_name='z',label_name='z',fig_title='My Figure',hovertext=None,**kwargs):
    import chart_studio.plotly as py
    import plotly.graph_objects as go
    import plotly.tools as tls
    import plotly.express as px

    if hovertext is None:
        hovertext=['%s: %.2f<br>%s: %.2f<br>%s: %.2f'%(x_name,x[j],y_name,y[j],z_name,z[j]) for j in range(len(x))]
    if fig is None:
        fig=go.FigureWidget(go.Scatter3d(x=x, y=y, z=z,
            hoverinfo='text',hovertext=hovertext,**kwargs))
        layout = dict(
        width=800,
        height=800,
        autosize=False,
        title=fig_title,
        scene=dict(
                xaxis=dict(
                    gridcolor='rgb(255, 255, 255)',
                    zerolinecolor='rgb(255, 255, 255)',
                    showbackground=True,
                    backgroundcolor='rgb(230, 230, 230)',
                    title=x_name,
                    autorange='reversed',
                    zerolinewidth=5,
                    tick0=0

                    #titlefont=dict(size=18,color='rgb(255,255,255)'),
                ),
                yaxis=dict(
                    gridcolor='rgb(255, 255, 255)',
                    zerolinecolor='rgb(255, 255, 255)',
                    showbackground=True,
                    backgroundcolor='rgb(230, 230, 230)',
                    title=y_name,
                    autorange='reversed'
                    #titlefont=dict(size=18,color='rgb(255,255,255)'),

                ),
                zaxis=dict(
                    gridcolor='rgb(255, 255, 255)',
                    zerolinecolor='rgb(255, 255, 255)',
                    showbackground=True,
                    backgroundcolor='rgb(230, 230, 230)',
                    title=z_name,
                    #titlefont=dict(size=18,color='rgb(255,255,255)'),

                ),
                camera=dict(
                    up=dict(
                        x=0,
                        y=0,
                        z=1
                    ),
                    eye=dict(
                       x=-2.5,#-1.7428,
                       y=0,#1.0707,
                       z=0#0.7100,
                    ),
                    center=dict(x=0,y=0,z=0),
                    projection=dict(type='orthographic')
                ),
                #aspectratio = dict( x=1, y=1, z=1 ),
                #aspectmode = 'manual'
            ),
        )

        fig['layout']=layout
    else:
        fig.add_trace(go.Scatter3d(x=x,y=y,z=z,hoverinfo='text',hovertext=hovertext,**kwargs))
    return(fig)



def plot3d_Volume():
    import chart_studio.plotly as py
    import plotly.graph_objects as go
    import plotly.tools as tls
    import plotly.express as px

    X, Y, Z = np.mgrid[-5:5:40j, -5:5:40j, -5:5:40j]

    # ellipsoid
    values = X * X * 0.75 + Y * Y + Z * Z * 3
    values/=np.max(values)
    
    fig = go.Figure(data=go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=values.flatten(),
        isomin=.05,
        showscale=False,
        colorscale='blues',
        isomax=.2,
        surface_count=1,
        caps=dict(x_show=False, y_show=False)
        ))
    X2, Y2, Z2 = np.mgrid[-5:5:40j, -5:5:40j, -5:5:40j]
    values2 = X2 * X2 * 0.75 + Y2 * Y2 + Z2 * Z2 * 3
    values2/=np.max(values2)
    
    fig.add_trace(go.Isosurface(x=X2.flatten()-3,
        y=Y2.flatten()-9,
        z=Z2.flatten(),
        value=values2.flatten(),
        isomin=.05,
        colorscale='reds',
        showscale=False,
        isomax=.15,
        surface_count=1,
        caps=dict(x_show=False, y_show=False)))

    fig.add_trace(go.Scatter3d(x=[0],y=[20],z=[0],mode='markers'))
    x=np.linspace(-2,10,100)
    y=np.linspace(-9,0,100)
    z=np.zeros(len(y))
    
    fig.add_trace(go.Scatter3d(x=x,y=y,z=z,mode='lines',line=dict(color='white')))
    x=np.linspace(-2,-12,100)
    fig.add_trace(go.Scatter3d(x=x,y=y,z=z,mode='lines',line=dict(color='white')))
    x=np.linspace(10,0,100)
    y=np.linspace(0,20,100)
    z=np.zeros(len(y))
    fig.add_trace(go.Scatter3d(x=x,y=y,z=z,mode='lines',line=dict(color='white')))
    x=np.linspace(-12,0,100)
    fig.add_trace(go.Scatter3d(x=x,y=y,z=z,mode='lines',line=dict(color='white')))

    layout = dict(
        width=800,
        height=700,
        autosize=False,
        title='SN Requiem',
        #plot_bgcolor='#000000',
        scene=dict(
            xaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(0, 0, 0)',
                #title='Age (Observer Days)',
                #titlefont=dict(size=18,color='rgb(255,255,255)'),
                #autorange='reversed'
            ),
            yaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(0, 0, 0)',

                #title='F105W-F160W Color',
                #titlefont=dict(size=18,color='rgb(255,255,255)'),
                autorange='reversed'

            ),
            zaxis=dict(
                gridcolor='rgb(255, 255, 255)',
                zerolinecolor='rgb(255, 255, 255)',
                showbackground=True,
                backgroundcolor='rgb(0, 0, 0)',
                #title='F160W AB Magnitude',
                #titlefont=dict(size=18,color='rgb(255,255,255)'),
                autorange='reversed'

            ),
            camera=dict(
                up=dict(
                    x=0,
                    y=0,
                    z=1
                ),
                #eye=dict(
                #    x=-1.7428,
                #    y=1.0707,
                #    z=0.7100,
                #)
            ),
            aspectratio = dict( x=1, y=1, z=0.7 ),
            aspectmode = 'manual'
            ),
        )

    fig['layout']=layout

    return(fig)


def split_plot(ax,plot_type,x,y=None,yerr=None,xerr=None,x_lab='',y_lab='',xticks=False,fontsize=18,
                    split_size='50%',ticklabelsize=16,**kwargs):
    ax_divider = make_axes_locatable(ax)
    ax_ml = ax_divider.append_axes("bottom", size=split_size, pad=.2)
    ticks=[]
    for tick in ax_ml.xaxis.get_major_ticks():
        ticks.append('')
        tick.label.set_fontsize(ticklabelsize)
    if not xticks:
        ax.set_xticklabels(ticks)
    for tick in ax_ml.yaxis.get_major_ticks():
        tick.label.set_fontsize(ticklabelsize)
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
        tick.label.set_fontsize(ticklabelsize)
    for tick in ax_ml.yaxis.get_major_ticks():
        tick.label.set_fontsize(ticklabelsize)
    ax_ml.set_xlabel(x_lab,fontsize=fontsize)
    ax_ml.set_ylabel(y_lab,fontsize=fontsize)
    return(ax,ax_ml)



def grid_plot(grid_x,grid_y,figsize=(12,12),sharex=False,sharey=False):
    fig=plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(grid_x, grid_y)
    
    axes=[]
    if grid_x==1:
        for j in range(grid_y):
            axes.append(fig.add_subplot(gs[0,j]))
            if sharex and j>0:
                axes[-1].get_shared_x_axes().join(axes[0], axes[-1])
            if sharey and j>0:
                axes[-1].get_shared_y_axes().join(axes[0], axes[-1])
    elif grid_y==1:
        for i in range(grid_x):
            axes.append(fig.add_subplot(gs[i,0]))
            if sharex and i>0:
                axes[-1].get_shared_x_axes().join(axes[0], axes[-1])
            if sharey and i>0:
                axes[-1].get_shared_y_axes().join(axes[0], axes[-1])
    else:
        for i in range(grid_x):
            tempax=[]
            for j in range(grid_y):
                tempax.append(fig.add_subplot(gs[i,j]))
            axes.append(tempax)
        for i in range(grid_x):
            for j in range(grid_y):
                if sharex and i>0:
                    axes[i][j].get_shared_x_axes().join(axes[0][j], axes[i][j])
                if sharey and j>0:
                    axes[i][j].get_shared_y_axes().join(axes[i][0], axes[i][j])
    
    
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    return(fig,axes,gs)


def multivariateGrid(col_x, col_y, col_k, df, k_is_color=False, 
                     scatter_alpha=.5,global_hist=False,kind='scatter',dist_kde=True,dist_hist=False,
                     replacements=None,xlim=None,ylim=None,bins_x=None,bins_y=None,fontsize=16,**kwargs):
    import seaborn as sns

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
        data=df,
        xlim=xlim,
        ylim=ylim,
        **kwargs
    )

    color = None
    legends=[]
    handles=[]
    for name, df_group in df.groupby(col_k):
        legends.append(str(name))
        if k_is_color:
            color=name
        if kind=='scatter':
            g.plot_joint(
                colored_scatter(df_group[col_x],df_group[col_y],c=color),
            )
        elif kind=='kde':
            g.x=df_group[col_x]
            g.y=df_group[col_y]
            handles.append(g.plot_joint(sns.kdeplot,label=str(name),**kwargs))
    
    for name, df_group in df.groupby(col_k):
        sns.distplot(
            df_group[col_x].values,
            ax=g.ax_marg_x,
            color=color,
            kde=dist_kde,
            hist=dist_hist,
            bins=bins_x
        )
        sns.distplot(
            df_group[col_y].values,
            ax=g.ax_marg_y,
            color=color,            
            vertical=True,
            kde=dist_kde,
            hist=dist_hist,
            bins=bins_y
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
    if replacements is not None:
        xlabel = g.ax_joint.get_xlabel()
        ylabel = g.ax_joint.get_ylabel()
        if xlabel in replacements.keys():
            g.ax_joint.set_xlabel(replacements[xlabel],fontsize=fontsize)
        if ylabel in replacements.keys():
            g.ax_joint.set_ylabel(replacements[ylabel],fontsize=fontsize)
        g.ax_joint.legend([x if x not in replacements.keys() else replacements[x] for x in legends])
    return(g)

