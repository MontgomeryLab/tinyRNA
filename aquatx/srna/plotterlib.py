""" Plotting functions for small RNA data. 

This module contains functions to create relevant plots for small RNA data for use
with the AQuATx pipeline. The plots are built off of matplotlib, but updated to
use the plot style of this tool. Other color schemes and built-in matplotlib styles
can be used. 
"""
import itertools
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import warnings; warnings.filterwarnings(action='once')

class plotterlib:
    def set_aquatx_style(self, **kwargs):
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

    def size_dist_bar(self, size_df: pd.DataFrame, **kwargs):
        """Creates a size distribution plot.

        Args:
            size_df: A pandas dataframe containing the size x 5'nt raw counts

        Returns:
            sizeb: A stacked barplot of size + 5'nt data
        """
        # Set the plot style
        kwargs = self.set_aquatx_style(figure={'figsize': (6,4)}, **kwargs)

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

    def class_pie(self, class_df, **kwargs):
        """Creates a pie chart of sRNA classes.

        Args:
            class_df: A pandas dataframe containing counts per class

        Returns:
            cpie: A pie chart of sRNA classes
        """
        # Set the plot style
        kwargs = self.set_aquatx_style(**kwargs)

        # Convert reads to proportion
        class_prop = class_df/(class_df.sum())

        # Create the plot
        cpie = class_prop.plot(kind='pie', **kwargs)
        cpie.set_aspect("equal")
        cpie.legend(loc='best', bbox_to_anchor= (1,0.5), fontsize=10, labels=class_prop.index)
        cpie.set_ylabel('')
        cpie.set_xlabel('')

        return cpie

    def class_barh(self, class_df, **kwargs):
        """Creates a horizontal bar chart of sRNA classes.

        Args:
            class_df: A pandas dataframe containing counts per class

        Returns:
            cbar: A horizontal bar chart of sRNA classes
        """
        # Set the plot style
        kwargs = self.set_aquatx_style(**kwargs)

        # Convert reads to proportion
        class_prop = class_df/(class_df.sum())

        # Create the plot
        colors = kwargs.get('colors', plt.rcParams['axes.prop_cycle'].by_key()['color'])
        cbar = class_prop.plot(kind='barh', color=colors)
        cbar.set_xlabel('Proportion of reads')
        cbar.set_xlim(0,1)

        return cbar

    def class_pie_barh(self, class_df, **kwargs):
        """Creates both a pie & bar chart in the same figure

        Args:
            class_df: A pandas dataframe containing counts per class

        Returns:
            cplots: A pie chart & horizontal bar chart of sRNA classes
        """
        # Set the plot style
        kwargs = self.set_aquatx_style(**kwargs, figure={'figsize':(8,4)})

        # Convert reads to proportion
        class_prop = class_df/(class_df.sum())

        # Create the plots
        fig, ax = plt.subplots(nrows=1, ncols=2)

        # pie
        class_prop.plot(kind='pie', subplots=True, ax=ax[0], labels=None, **kwargs)
        ax[0].set_aspect("equal")
        ax[0].legend(loc='best', bbox_to_anchor= (1,0.5), fontsize=10, labels=class_prop.index)
        ax[0].set_ylabel('')
        ax[0].set_xlabel('')

        # bar
        colors = kwargs.get('colors', plt.rcParams['axes.prop_cycle'].by_key()['color'])
        class_prop.plot(kind='barh', subplots=True, ax=ax[1], color=colors)
        ax[1].legend().set_visible(False)
        ax[1].set_title('')
        ax[1].set_xlabel('Proportion of reads')
        ax[1].set_ylabel('')
        ax[1].set_xlim(0,1)

        # finalize & save figure
        fig.suptitle("Proportion of classes of small RNAs", fontsize=22)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        fig.subplots_adjust(top=0.85)

        return fig, ax

    def scatter_range(self, df):
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

    def scatter_simple(self, count_x, count_y, log_norm=False, **kwargs):
        """Creates a simple scatter plot of counts.

        Args:
            count_x: A pandas dataframe/series of counts per feature (X axis)
            count_y: A pandas dataframe/series of counts per feature (Y axis)
            log_norm: Apply log2 normalization to the data
            kwargs: Additional arguments to pass to matplotlib rc or plot

        Returns:
            sscat: A simple scatter plot of counts
        """
        # Set the plot style
        kwargs = self.set_aquatx_style(figure={'figsize': (8,8)}, **kwargs)

        fig, sscat = plt.subplots()

        # log2 normalize data if requested
        if log_norm:
            count_x = count_x.apply(np.log2)
            count_y = count_y.apply(np.log2)
            sscat_lims = self.scatter_range(pd.concat([count_x, count_y]))
            sscat.set_xlim(sscat_lims)
            sscat.set_ylim(sscat_lims)
            sscat.scatter(count_x, count_y, **kwargs)

            oldticks = sscat.get_xticks()
            newticks = np.empty([len(oldticks)-1, 8])

            for i in range(1,len(oldticks)-1):
                newticks[i,:] = np.arange(2**oldticks[i-1], 2**oldticks[i], (2**oldticks[i] - 2**oldticks[i-1])/8)

            newticks = np.sort(newticks[2:,:].flatten())
            sscat.set_xticks(np.log2(newticks), minor=True)
            sscat.set_xticklabels(np.round(2**oldticks))
            sscat.set_yticks(np.log2(newticks), minor=True)
            sscat.set_yticklabels(np.round(2**oldticks))

        else:
            sscat.scatter(count_x, count_y, **kwargs)
            sscat_lims = self.scatter_range(pd.concat([count_x, count_y]))
            sscat.set_xlim(sscat_lims)
            sscat.set_ylim(sscat_lims)

        return sscat

    def scatter_grouped(self, count_x, count_y, *args, log_norm=False, labels=None, **kwargs):
        """Creates a scatter plot with different groups highlighted.

        Args:
            count_x: A pandas dataframe/series of counts per feature (X axis)
            count_y: A pandas dataframe/series of counts per feature (Y axis)
            args: A list of features to highlight, can pass multiple lists
            log_norm: whether or not the data should be log-normalized
            kwargs: Additional arguments to pass to matplotlib rc or plot

        Returns:
            gscat: A scatter plot containing groups highlighted different colors
        """
        kwargs = self.set_aquatx_style(figure={'figsize': (8,8)}, **kwargs)

        # Subset the base points to avoid overplotting
        count_x_base = count_x.drop(list(itertools.chain(*args)))
        count_y_base = count_y.drop(list(itertools.chain(*args)))

        # Create a base plot using counts not in *args
        colors = iter(kwargs.get('colors', plt.rcParams['axes.prop_cycle'].by_key()['color']))
        gscat = self.scatter_simple(count_x_base, count_y_base, log_norm=log_norm, color='#888888', marker='s', alpha=0.3, s=50, edgecolors='none', **kwargs)

        if log_norm:
            count_x = count_x.apply(np.log2).replace(-np.inf, 0)
            count_y = count_y.apply(np.log2).replace(-np.inf, 0)

        if labels is None:
            labels = list(range(len(args)))

        # Add points for each *args
        for group in args:
            gscat.scatter(count_x.loc[group], count_y.loc[group], color=next(colors), marker='s', alpha=0.9, s=50, edgecolors='none', **kwargs)

        gscat.legend(labels=labels)

        return gscat
