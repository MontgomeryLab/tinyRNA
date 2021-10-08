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
import math
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

    def __init__(self, user_style_sheet, debug=False):

        # Set global plot style once
        plt.style.use(user_style_sheet)

        if debug:
            mpl.use('TkAgg')
            mpl.rcParams['savefig.dpi'] = 100
        else:
            # Slightly better performance
            mpl.use('PDF')

        # Create one subplot per plot type to reuse between calls
        fig_args = {
            'class_pie_barh': {'figsize': (8, 4), 'nrows': 1, 'ncols': 2, 'tight_layout': True},
            'len_dist_bar': {'figsize': (7, 4)},
            'scatter_simple': {'figsize': (8, 8)}
        }

        self.subplots = {}
        for plot in fig_args:
            fig, ax = plt.subplots(**fig_args[plot])
            self.subplots[plot] = {'fig': fig, 'ax': ax}

        self.sample_lims = {}
        self.axis_cache = {}

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

    def scatter_simple(self, count_x: pd.Series, count_y: pd.Series, log_norm=False, **kwargs) -> plt.Axes:
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
        ax: plt.Axes

        # log2 normalize data if requested
        if log_norm:
            # Set log2 scale
            ax.set_xscale('log', base=2)
            ax.set_yscale('log', base=2)

            # Unset default locators and formatters
            for axis in [ax.xaxis, ax.yaxis]:
                axis.set_major_locator(tix.NullLocator())
                axis.set_minor_locator(tix.NullLocator())
                axis.set_major_formatter(tix.LogFormatter(base=2))
                axis.set_minor_formatter(tix.NullFormatter())

            ax.scatter(count_x, count_y, edgecolor='none', **kwargs)

            # Draw y = x, x +/- 1 lines (log2 scale) using point pairs for the 3 lines
            for p1,p2 in [((1, 2), (2, 4)), ((1, 1), (2, 2)), ((1, 0.5), (2, 1))]:
                ax.axline(p1, p2, color='#CCCCCC', label='_nolegend_')
        else:
            # Set linear scale
            ax.set_xscale('linear')
            ax.set_yscale('linear')

            # Set up axis ticks and labels
            for axis in [ax.xaxis, ax.yaxis]:
                axis.set_major_locator(tix.MaxNLocator(6))
                axis.set_major_formatter(tix.ScalarFormatter())

            ax.scatter(count_x, count_y, **kwargs)

            # Draw y = x, x +/- 1 lines using point pairs for the 3 lines
            for p1,p2 in [((0, 1), (1, 2)), ((0, 0), (1, 1)), ((0,-1), (1, 0))]:
                ax.axline(p1, p2, color='#CCCCCC', label='_nolegend_')

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
        count_x_base = count_x.drop(itertools.chain(*args))
        count_y_base = count_y.drop(itertools.chain(*args))

        if labels is None:
            labels = list(range(len(args)))

        colors = iter(kwargs.get('colors', plt.rcParams['axes.prop_cycle'].by_key()['color']))
        argsit = iter(args)

        if any([len(outgroup) == 0 for outgroup in [count_x_base, count_y_base]]):
            # There is no outgroup, plot the first group with scatter_simple() to set scale and lines
            group = next(argsit)
            gscat = self.scatter_simple(count_x.loc[group], count_y.loc[group],
                                        log_norm=log_norm, color=next(colors), **kwargs)
        else:
            # Plot the outgroup in light grey (these are counts not in *args)
            gscat = self.scatter_simple(count_x_base,count_y_base, log_norm=log_norm, color='#B3B3B3', **kwargs)

        # Add each group to plot with a different color
        for group in argsit:
            gscat.scatter(count_x.loc[group], count_y.loc[group], color=next(colors), edgecolor='none', **kwargs)

        self.set_scatter_lims(gscat)
        gscat.legend(labels=labels)
        return gscat

    def set_scatter_lims(self, ax: plt.Axes):
        """Scatter plot will be centered about diagonal, without lower left tick label

        Args:
            ax: A scatter plot Axes object
        """

        # Let autoscale do most of the work
        ax.autoscale(True)

        # Get maximum and minimum limits of all axes
        lim = ax.viewLim.bounds
        ax_min, ax_max = min(lim), max(lim)
        tick_locs = [2 ** x for x in range(math.floor(np.log2(ax_min)), math.ceil(np.log2(ax_max)))]
        min_tick, max_tick = tick_locs[0], tick_locs[-1]

        ax.xaxis.set_major_locator(tix.FixedLocator(tick_locs))
        ax.yaxis.set_major_locator(tix.FixedLocator(tick_locs))

        # Square up limits to center about diagonal
        ax.set_xlim(left=ax_min, right=ax_max)
        ax.set_ylim(bottom=ax_min, top=ax_max)

        every_nth_label = 3

        # Hide ticks near origin and set minor tick parameters
        for axis, spine in [(ax.xaxis, ax.spines.bottom), (ax.yaxis, ax.spines.left)]:
            major_ticks = axis.get_major_ticks()
            major_ticks[0].set_visible(False)
            for i, tick in enumerate(major_ticks):
                if i % every_nth_label != 0:
                    tick.label1.set_visible(False)

            axis.set_minor_locator(tix.LogLocator(
                base=2.0,
                numticks=self.get_min_LogLocator_numticks(axis),
                subs=np.log2(np.linspace(2 ** 2, 2 ** 4, 10))[:-1]))

            # THUNDEROUS PRAYER HANDS, FREE US FROM OUR EARTHLY BONDS [[[thunder]]]
            ax.figure.canvas.draw()

            for tick in major_ticks:
                line = tick.tick1line
                if line._xy is None or line._xy.max() < min_tick or line._xy.max() > max_tick:
                    line.set_visible(False)

            for tick in axis.get_minor_ticks():
                line = tick.tick1line
                if line._xy is None or line._xy.max() < 0.25 or line._xy.max() > max_tick:
                    line.set_visible(False)

            pass
        ## NOTE TO FUTURE DEVELOPERS
        # You may have encountered a pairwise comparison which differs in view limits between
        # scatter_by_dge and scatter_by_dge_class. This is going to be a painful problem to solve
        # due to the amount of explicit and implicit redundancy in mpl (also: keep an eye on those
        # @property decorators...). You'll want to start at _AxesBase.autoscale_view(). Good luck.

    @staticmethod
    def get_min_LogLocator_numticks(axis: plt.Axes):
        vmin, vmax = axis.get_view_interval()
        log_vmin = math.log(vmin) / math.log(2)
        log_vmax = math.log(vmax) / math.log(2)

        numdec = math.floor(log_vmax) - math.ceil(log_vmin)
        return numdec + 2  # Want: [ (numdec + 1) // nticks + 1 ] == 1

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
