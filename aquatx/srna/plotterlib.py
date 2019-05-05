""" Plotting functions for small RNA data. 

This module contains functions to create relevant plots for small RNA data for use
with the AQuATx pipeline. The plots are built off of matplotlib, but updated to
use the plot style of this tool. Other color schemes and built-in matplotlib styles
can be used. 
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings(action='once')

def set_aquatx_style(**kwargs):
    """Set parameters for the stylesheet for AQuATx if parameters not set by user."""

    # Define default aquatx settings
    srna_colors =  ['F1605D', '51B9CF', 'FDC010', 'A5D38E', 'ED2891', '989898']
    aq_kwargs = {'figure': {'figsize': (6,6), 
                            'dpi': 100, 
                            'facecolor': 'white', 
                            'edgecolor': '0.5'},
                 'xtick': {'color': '333333',
                           'direction': 'out',
                           'major.size': 3.0,
                           'minor.size': 1.5,
                           'major.width': 0.8,
                           'minor.width': 0.6,
                           'major.pad': 3.5,
                           'labelsize': 16},
                 'ytick': {'color': '333333',
                           'direction': 'out',
                           'major.size': 3.0,
                           'minor.size': 1.5,
                           'major.width': 0.8,
                           'minor.width': 0.6,
                           'major.pad': 3.5,
                           'labelsize': 16},
                 'axes': {'labelsize': 16,
                          'titlesize': 22,
                          'facecolor': 'white',
                          'edgecolor': '333333',
                          'linewidth': 1,
                          'grid': True,
                          'labelcolor': '333333',
                          'labelpad': 4.0,
                          'axisbelow': True,
                          'autolimit_mode': 'round_numbers',
                          'prop_cycle': mpl.cycler('color', srna_colors),
                          'xmargin': 0.5,
                          'ymargin': 0.5},
                 'legend': {'fontsize': 16,
                            'framealpha': 0.8,
                            'edgecolor': '0.8',
                            'markerscale': 1.0,
                            'borderpad': 0.4,
                            'labelspacing': 0.5,
                            'handletextpad': 0.8,
                            'borderaxespad': 0.5,
                            'columnspacing': 2.0},
                 'pdf': {'fonttype': 42},
                 'ps': {'fonttype': 42},
                 'savefig': {'pad_inches': 0.1},
                 'font': {'family': 'sans-serif',
                          'sans-serif': 'Arial',
                          'size': 10},
                 'lines': {'linewidth': 3.0,
                           'markersize': 11.0,
                           'markeredgewidth': 0},
                 'grid': {'color': '333333',
                          'linestyle': '--',
                          'linewidth': 0.8,
                          'alpha': 0.2},
                }

    # Add default values if user has not set
    for key, val in aq_kwargs.items():
        kwargs.setdefault(key, val)

    # Set all values in rcParams, then remove from kwargs
    kwargskeys = list(kwargs.keys())
    for key in kwargskeys:
        try:
            mpl.rc(key, **kwargs[key])
            _ = kwargs.pop(key)
        except (KeyError, TypeError):
            pass
        
        try:
            mpl.rcParams[key] = kwargs[key]
            _ = kwargs.pop(key)
        except KeyError:
            pass
    
    return kwargs

def size_dist_bar(size_df, **kwargs):
    """Creates a size distribution plot.

    Args:
        size_df: A pandas dataframe containing the size x 5'nt raw counts

    Returns:
        sizeb: A stacked barplot of size + 5'nt data
    """
    # Set the plot style
    kwargs = set_aquatx_style(figure={'figsize': (6,4)}, **kwargs)

    # Convert reads to proportion
    size_prop = size_df / size_df.sum().sum()

    # Create the plot
    sizeb = size_prop.plot(kind='bar', stacked=True, **kwargs)
    sizeb.set_title('Distribution of aligned reads')
    sizeb.set_ylim(0,np.max(np.sum(size_prop, axis=1))+0.025)
    sizeb.set_ylabel('Proportion of Reads')
    sizeb.set_xlabel('Length of Sequence')
    sizeb.set_xticklabels(sizeb.get_xticklabels(), rotation=0)
    
    return sizeb

def class_pie(class_df, **kwargs):
    """Creates a pie chart of sRNA classes.
  
    Args:
        class_df: A pandas dataframe containing counts per class

    Returns:
        cpie: A pie chart of sRNA classes
    """
    # Set the plot style
    kwargs = set_aquatx_style(**kwargs)

    # Convert reads to proportion
    class_prop = class_df/(class_df.sum())

    # Create the plot
    cpie = class_prop.plot(kind='pie', **kwargs)
    cpie.set_aspect("equal")
    cpie.legend(loc='best', bbox_to_anchor= (1,0.5), fontsize=10, labels=class_prop.index)
    cpie.set_ylabel('')
    cpie.set_xlabel('')
    
    return cpie

def class_barh(class_df, **kwargs):
    """Creates a horizontal bar chart of sRNA classes.

    Args:
        class_df: A pandas dataframe containing counts per class

    Returns:
        cbar: A horizontal bar chart of sRNA classes
    """
    # Set the plot style
    kwargs = set_aquatx_style(**kwargs)

    # Convert reads to proportion
    class_prop = class_df/(class_df.sum())

    # Create the plot
    cbar = class_prop.barh(color=plt.rcParams['axes.prop_cycle'].by_key()['color'])
    cbar.set_xlabel('Proportion of reads')
    cbar.set_xlim(0,1)
    
    return cbar

def scatter_range(df):
    """ Find an appropriate range for x,y limits of a scatter plot.

    Args:
        df: A dataframe being plotted

    Returns:
        lim_min: The minimum value to set axis limits
        lim_max: The maximum value to set axis limits
    """


    if np.min(np.min(df)) == -np.inf:
        df_min = 0
    else:
        df_min = np.min(np.min(df))

    if np.max(np.max(df)) == np.inf:
        df_max = np.max(np.max(df[~(df == np.inf)]))
    else:
        df_max = np.max(np.max(df))

    intv = (df_max - df_min) / 12
    
    lim_min = df_min - intv
    lim_max = df_max + intv

    return lim_min, lim_max

def scatter_simple(count_df, **kwargs):
    """Creates a simple scatter plot of counts.

    Args:
        count_df: A pandas dataframe containing the counts per feature

    Returns:
        sscat: A simple scatter plot of counts
    """
    # Set the plot style
    kwargs = set_aquatx_style(figure={'figsize': (8,8)})
    
    # Create the plot
    sscat = count_df.plot(kind='scatter', **kwargs)
    sscat.set_xticklabels(sscat.get_xticklabels(), rotation=0)
    sscat_lims = scatter_range(count_df)
    sscat.set_xlim(sscat_lims)
    sscat.set_ylim(sscat_lims)

    return sscat
