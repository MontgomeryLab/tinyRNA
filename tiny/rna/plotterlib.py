""" Plotting functions for small RNA data. 

This module contains functions to create relevant plots for small RNA data for use
with the tinyRNA pipeline. The plots are built using matplotlib and our style sheet.
You may override these styles by obtaining a copy of the style sheet (tiny get-template),
modifying it, and passing it to tiny-plot via the -s/--style-sheet argument. If
using this module directly, it may be passed at construction time.
"""
import pandas as pd
import numpy as np
import itertools
import locale
import os

# cwltool appears to unset all environment variables including those related to locale
# This leads to warnings from plt's FontConfig manager, but only for pipeline/cwl runs
curr_locale = locale.getlocale()
if curr_locale[0] is None:
    # Default locale otherwise unset
    os.environ['LC_CTYPE'] = 'en_US.UTF-8'

import matplotlib as mpl
import matplotlib.ticker as tix
import matplotlib.pyplot as plt

from typing import Union


class plotterlib:

    def __init__(self, user_style_sheet):

        # Set global plot style once
        plt.style.use(user_style_sheet)

        # Improves performance (sacrifices interactivity)
        mpl.use("PDF")

        # Create one subplot per plot type to reuse between calls
        fig_args = {
            'class_pie_barh': {'figsize': (8, 4), 'nrows': 1, 'ncols': 2, 'tight_layout': True},
            'len_dist_bar': {'figsize': (6, 4)},
            'scatter_simple': {'figsize': (8, 8)}
        }

        self.subplots = {}
        for plot in fig_args:
            fig, ax = plt.subplots(**fig_args[plot])
            self.subplots[plot] = {'fig': fig, 'ax': ax}

    def len_dist_bar(self, size_df: pd.DataFrame, **kwargs) -> plt.Axes:
        """Creates a stacked barplot of 5' end nucleotides by read length

        Args:
            size_df: A dataframe containing the size x 5'nt raw counts
            kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()

        Returns:
            sizeb: A stacked barplot of size + 5'nt data
        """

        # Retrieve axis and styles for this plot type
        fig, ax = self.reuse_subplot("len_dist_bar")

        # Convert reads to proportion
        size_prop = size_df / size_df.sum().sum()

        # Override default colors. User may override with kwargs. (Orange, Yellow-green, Blue, Pink)
        colors = {'axes.prop_cycle': mpl.cycler(color=['#F78E2D', '#CBDC3F', '#4D8AC8', '#E06EAA'])}

        # The style context allows us to use temporary styles
        with plt.style.context(colors):
            plt.sca(ax)
            sizeb = size_prop.plot(kind='bar', stacked=True, reuse_plot=True, **kwargs)
            sizeb.set_title('Distribution of aligned reads')
            sizeb.set_ylim(0,np.max(np.sum(size_prop, axis=1))+0.025)
            sizeb.set_ylabel('Proportion of Reads')
            sizeb.set_xlabel('Length of Sequence')
            sizeb.set_xticklabels(sizeb.get_xticklabels(), rotation=0)

        return sizeb

    def class_pie(self, class_s: pd.Series, **kwargs) -> plt.Axes:
        """Creates a pie chart of sRNA classes.

        Args:
            class_s: A pandas Series containing counts per class
            kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()

        Returns:
            cpie: A pie chart of sRNA classes
        """

        # Convert reads to proportion
        class_prop = class_s / class_s.sum()

        # Create the plot
        cpie = class_prop.plot(kind='pie', normalize=True, **kwargs)
        cpie.legend(loc='best', bbox_to_anchor=(1, 0.5), fontsize=10, labels=class_prop.index)
        cpie.set_aspect("equal")
        cpie.set_ylabel('')
        cpie.set_xlabel('')

        return cpie

    def class_barh(self, class_s: pd.Series, **kwargs) -> plt.Axes:
        """Creates a horizontal bar chart of sRNA classes.

        Args:
            class_s: A pandas Series containing a single library's counts per class
            kwargs: Additional keyword arguments to pass to pandas.Series.plot()

        Returns:
            cbar: A horizontal bar chart of sRNA classes
        """

        # Convert reads to proportion
        class_prop = class_s / class_s.sum()

        # df.plot(kind=barh) ignores axes.prop_cycle... (ugh)
        colors = kwargs.get('colors', plt.rcParams['axes.prop_cycle'].by_key()['color'])

        # Create the plot
        cbar = class_prop.plot(kind='barh', color=colors, **kwargs)
        cbar.set_xlabel('Proportion of reads')
        cbar.legend().set_visible(False)
        cbar.set_title('')
        cbar.set_ylabel('')
        cbar.set_xlim(0,1)

        return cbar

    def class_pie_barh(self, class_s: pd.Series, **kwargs) -> plt.Figure:
        """Creates both a pie & bar chart in the same figure

        Args:
            class_s: A pandas Series containing counts per class
            kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()

        Returns:
            cplots: A pie chart & horizontal bar chart of sRNA classes
        """

        # Retrieve axis and styles for this plot type
        fig, ax = self.reuse_subplot("class_pie_barh")

        # Plot pie and barh on separate axes
        self.class_pie(class_s, ax=ax[0], labels=None, **kwargs)
        self.class_barh(class_s, ax=ax[1], legend=None, title=None, ylabel=None, **kwargs)

        # finalize & save figure
        fig.suptitle("Proportion of classes of small RNAs", fontsize=22)
        fig.subplots_adjust(top=0.85)

        return fig

    def scatter_range(self, df: pd.DataFrame) -> (int, int):
        """Find an appropriate range for x,y limits of a scatter plot.

        Args:
            df: A dataframe being plotted

        Returns:
            lim_min: The minimum value to set axis limits
            lim_max: The maximum value to set axis limits
        """

        if np.min(df) == -np.inf:
            df_min = 0
        else:
            df_min = np.min(df)

        if np.max(df) == np.inf:
            df_max = np.max(df[~(df == np.inf)])
        else:
            df_max = np.max(df)

        intv = (df_max - df_min) / 12

        lim_min = df_min - intv
        lim_max = df_max + intv

        return lim_min, lim_max

    def scatter_simple(self, count_x: pd.DataFrame, count_y: pd.DataFrame, log_norm=False, lim=None, **kwargs) -> plt.Axes:
        """Creates a simple scatter plot of counts.

        Args:
            count_x: A pandas dataframe/series of counts per feature (X axis)
            count_y: A pandas dataframe/series of counts per feature (Y axis)
            log_norm: Apply log2 normalization to the data
            kwargs: Additional keyword arguments to pass to pyplot.Axes.scatter()

        Returns:
            ax: A simple scatter plot of counts
        """

        # Retrieve axis and styles for this plot type
        fig, ax = self.reuse_subplot("scatter_simple")

        # log2 normalize data if requested
        if log_norm:
            count_x = count_x.apply(np.log2)
            count_y = count_y.apply(np.log2)
            # Shouldn't this be done for the entire plotted dataset and not for the initial omitted subset?
            sscat_lims = lim if lim is not None else self.scatter_range(pd.concat([count_x, count_y]))
            if not np.isnan(sscat_lims).any():
                ax.set_xlim(sscat_lims)
                ax.set_ylim(sscat_lims)
            ax.scatter(count_x, count_y, **kwargs)

            oldticks = ax.get_xticks()
            newticks = np.empty([len(oldticks)-1, 8])

            for i in range(1,len(oldticks)-1):
                newticks[i,:] = np.arange(2**oldticks[i-1], 2**oldticks[i], (2**oldticks[i] - 2**oldticks[i-1])/8)

            newticks = np.sort(newticks[2:,:].flatten())

            # These lines have been added to address the FixedLocator warning
            ax.xaxis.set_major_locator(tix.FixedLocator(oldticks))
            ax.yaxis.set_major_locator(tix.FixedLocator(oldticks))

            ax.set_xticks(np.log2(newticks), minor=True)
            ax.set_xticklabels(np.round(2**oldticks).astype(int))
            ax.set_yticks(np.log2(newticks), minor=True)
            ax.set_yticklabels(np.round(2**oldticks).astype(int))
        else:
            ax.scatter(count_x, count_y, **kwargs)
            sscat_lims = self.scatter_range(pd.concat([count_x, count_y]))
            ax.set_xlim(sscat_lims)
            ax.set_ylim(sscat_lims)

        ax.axline([0, 1], [1, 2], color='#CCCCCC', label='_nolegend_')
        ax.axline([0, 0], [1, 1], color='#CCCCCC', label='_nolegend_')
        ax.axline([0,-1], [1, 0], color='#CCCCCC', label='_nolegend_')
        return ax

    def scatter_grouped(self, count_x: pd.DataFrame, count_y: pd.DataFrame, *args, log_norm=False, labels=None, **kwargs):
        """Creates a scatter plot with different groups highlighted.

        Args:
            count_x: A pandas dataframe/series of counts per feature (X axis)
            count_y: A pandas dataframe/series of counts per feature (Y axis)
            args: A list of features to highlight, can pass multiple lists
            log_norm: whether or not the data should be log-normalized
            kwargs: Additional arguments to pass to pyplot.Axes.scatter()

        Returns:
            gscat: A scatter plot containing groups highlighted with different colors
        """

        # Subset counts not in *args (for example, points with p-val above threshold)
        count_x_base = count_x.drop(list(itertools.chain(*args)))
        count_y_base = count_y.drop(list(itertools.chain(*args)))

        if log_norm:
            count_x = count_x.apply(np.log2).replace(-np.inf, 0)
            count_y = count_y.apply(np.log2).replace(-np.inf, 0)

        if labels is None:
            labels = list(range(len(args)))

        colors = iter(kwargs.get('colors', plt.rcParams['axes.prop_cycle'].by_key()['color']))
        ax_lim = self.scatter_range(pd.concat([count_x, count_y]))
        argsit = iter(args)

        if any([len(base_axis) == 0 for base_axis in [count_x_base, count_y_base]]):
            group = next(argsit)
            gscat = self.scatter_simple(count_x.loc[group], count_y.loc[group],
                log_norm=log_norm, color=next(colors), edgecolors='none', lim=ax_lim, **kwargs)
        else:
            # Create a base plot using counts not in *args
            gscat = self.scatter_simple(count_x_base,count_y_base,
                log_norm=log_norm, color='#B3B3B3', edgecolors='none', lim=ax_lim, **kwargs)

        # Add points for each *args
        for group in argsit:
            gscat.scatter(count_x.loc[group], count_y.loc[group],
                          color=next(colors), edgecolors='none', **kwargs)

        gscat.legend(labels=labels)

        return gscat

    def reuse_subplot(self, plot_type: str) -> (plt.Figure, Union[plt.Axes, np.ndarray]):
        """Retrieves the reusable subplot for this plot type

        Args:
            plot_type: The reusable plot name

        Returns:
            fig: The subplot's pyplot.Figure
            ax: The subplot's pyplot.Axes, or an array of axes if subplot's nrows or ncols is >1
        """

        fig, ax = self.subplots[plot_type].values() # Each plot type has a dedicated figure and axis
        if type(ax) == np.ndarray:
            for subax in ax: subax.clear()  # Remove previous plot data
        else:
            ax.clear()
        return fig, ax
