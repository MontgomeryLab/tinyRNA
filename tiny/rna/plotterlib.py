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
import matplotlib.pyplot as plt
import matplotlib.ticker as tix
import matplotlib.axis
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox
from matplotlib.scale import LogTransform

from typing import Union, Tuple, List, Optional


class plotterlib:

    def __init__(self, user_style_sheet, debug=False):

        if debug:
            mpl.use("TkAgg", force=True)
            mpl.rcParams['savefig.dpi'] = 100
        else:
            mpl.use("PDF", force=True)

        # Must occur after mpl.use()
        import matplotlib.pyplot as plt

        # Set global plot style once
        plt.style.use(user_style_sheet)

        self.subplots = self.init_subplots(debug)
        self.dge_scatter_tick_cache = {}

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
            log_norm: Plot on log scale rather than linear scale
            kwargs: Additional keyword arguments to pass to pyplot.Axes.scatter()

        Returns:
            ax: A simple scatter plot of counts
        """

        # Retrieve axis and styles for this plot type
        fig, ax = self.reuse_subplot("scatter")
        ax: plt.Axes

        # log2 normalize data if requested
        if log_norm:
            # Set log2 scale
            ax.set_xscale('log', base=2)
            ax.set_yscale('log', base=2)

            # Unset default formatters
            for axis in [ax.xaxis, ax.yaxis]:
                axis.set_major_formatter(tix.LogFormatterExponent(base=2))
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

    def scatter_grouped(self, count_x: pd.DataFrame, count_y: pd.DataFrame,
                        view_lims: Tuple[float, float] = None, *args, log_norm=False, labels=None, **kwargs):
        """Creates a scatter plot with different groups highlighted.

        Args:
            count_x: A pandas dataframe/series of counts per feature (X axis)
            count_y: A pandas dataframe/series of counts per feature (Y axis)
            view_lims: Optional plot view limits as tuple(min, max)
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

        self.set_square_scatter_view_lims(gscat, view_lims)
        self.set_scatter_ticks(gscat)
        gscat.legend(labels=labels)
        return gscat

    @staticmethod
    def get_scatter_view_lims(counts_df: pd.DataFrame) -> Tuple[float, float]:
        """Calculates scatter view limits for the counts dataframe"""

        x0 = counts_df.min(axis='columns').where(lambda x: x != 0).dropna().min()
        x1 = np.max(counts_df).max()
        minpos = 1e-300

        if not np.isfinite([x0, x1]).all() or not isinstance(x0, np.float) or x1 <= 0:
            print("The provided dataset contains invalid values.")
            return (minpos, minpos)

        x0, x1 = (minpos if x0 <= 0 else x0,
                  minpos if x1 <= 0 else x1)

        transform = LogTransform(base=2)
        inverse_trans = transform.inverted()

        x0t, x1t = transform.transform([x0, x1])
        delta = (x1t - x0t) * mpl.rcParams.get('axes.xmargin', 0)
        if not np.isfinite(delta): delta = 0

        return inverse_trans.transform([x0t - delta, x1t + delta])

    @staticmethod
    def set_square_scatter_view_lims(ax: plt.Axes, min_max=None):
        """Adjusts the scatter plot to display the same interval on both axes"""

        if min_max is not None:
            ax_min, ax_max = min_max
        else:
            lim = ax.viewLim.bounds
            ax_min, ax_max = min(lim), max(lim)

        # Square up limits to center about diagonal
        ax.set_xlim(left=ax_min, right=ax_max)
        ax.set_ylim(bottom=ax_min, top=ax_max)

    @staticmethod
    def get_min_LogLocator_numticks(axis: plt.Axes) -> int:
        """Calculates an axis' threshold numticks value for LogLocator to calculate (display) minor ticks

        Matplotlib's LogLocator will not locate ticks if its `numticks` parameter
        is below threshold for the view interval. Providing a `numticks` value below
        threshold results in minor ticks not being drawn.

        Args:
            axis: an x or y axis

        Returns: minimum `numticks` parameter value
        """

        vmin, vmax = axis.get_view_interval()
        log_vmin = math.log(vmin) / math.log(2)
        log_vmax = math.log(vmax) / math.log(2)

        numdec = math.floor(log_vmax) - math.ceil(log_vmin)
        return numdec + 2  # Want: [ (numdec + 1) // nticks + 1 ] == 1

    @staticmethod
    def get_fixed_majorticklocs(view_lims: Tuple[float, float, float, float]) -> Tuple[List[int], float, float]:
        """Produces a list of locations for major ticks for the given view limit"""

        ax_min, ax_max = min(view_lims), max(view_lims)
        floor, ceil, log2 = math.floor, math.ceil, np.log2
        locs = [2 ** x for x in range(floor(log2(ax_min)), ceil(log2(ax_max)))]
        return locs, ax_min, ax_max

    def set_scatter_ticks(self, ax: plt.Axes, minor_ticks=False):
        """Intelligently creates major and minor ticks for a square scatter plot while avoiding crowding"""

        # Get tick locations corresponding to the current view limits
        major_locs, ax_min, ax_max = self.get_fixed_majorticklocs(ax.viewLim.bounds)

        ax.xaxis.set_major_locator(tix.FixedLocator(major_locs))
        ax.yaxis.set_major_locator(tix.FixedLocator(major_locs))

        if len(self.dge_scatter_tick_cache):
            self.restore_ticks(ax, "both")
            return

        for axis in [ax.xaxis, ax.yaxis]:
            # Only display every nth major tick label
            ticks_displayed, last_idx = self.every_nth_label(axis, 3)

            if minor_ticks:
                axis.set_minor_locator(tix.LogLocator(
                    base=2.0,
                    numticks=self.get_min_LogLocator_numticks(axis),
                    subs=np.log2(np.linspace(2 ** 2, 2 ** 4, 10))[:-1]))

            min_tick = 2 ** (np.log2(ax_min)+1)
            max_tick = major_locs[last_idx]
            self.set_tick_bounds(axis, min_tick=min_tick, max_tick=max_tick, minor=minor_ticks)
            self.cache_ticks(axis, axis.__name__)

    def every_nth_label(self, axis: mpl.axis.Axis, n: int) -> Tuple[List[mpl.axis.Tick], int]:
        """Major ticks: hides all labels except every nth and mitigates crowding on lower and upper axis"""

        major_ticks = axis.get_major_ticks()
        ticks_displayed = []
        last_idx = 0

        for i, tick in enumerate(major_ticks):
            if i % n != 0:
                tick.label1.set_visible(False)
            else:
                ticks_displayed.append(tick)
                last_idx = i

        # Hide tick labels in the lower left corner, regardless
        major_ticks[0].label1.set_visible(False)

        # If the last tick label on the x-axis will extend past the plot space,
        # then hide it and its corresponding tick on the y-axis
        if axis.__name__ is "xaxis" and axis.get_tick_space() == len(ticks_displayed):
            major_ticks[last_idx].label1.set_visible(False)
            yaxis = axis.axes.yaxis
            yaxis.get_major_ticks()[last_idx].label1.set_visible(False)

        return ticks_displayed, last_idx

    def set_tick_bounds(self, axis: mpl.axis.Axis, min_tick: float, max_tick: float, minor=False):
        """Hide major and minor ticks that lie outside of the bounds defined

        For both major ticks and minor ticks, we need to first call their getters
        to build a list of partially constructed tick objects. Hiding the tick
        object excludes both the label and tickline subcomponents from the axis's
        draw() function at render time.
        """

        axis.get_major_ticks()
        for i, loc in enumerate(axis.get_majorticklocs()):
            if loc < min_tick or loc > max_tick:
                axis.majorTicks[i].label1.set_visible(False)
                axis.majorTicks[i].tick1line.set_visible(False)
                axis.majorTicks[i].gridline.set_visible(True)

        if minor:
            axis.get_minor_ticks()
            for i, loc in enumerate(axis.get_minorticklocs()):
                if loc < min_tick or loc > max_tick:
                    axis.minorTicks[i].set_visible(False)

    def cache_ticks(self, axis: mpl.axis.Axis, name: str):
        """Cache major and minor tick objects, which contain expensive data"""

        for type in ["major", "minor"]:
            self.dge_scatter_tick_cache[f"{name}_{type}_loc"] = getattr(axis, type).locator
            self.dge_scatter_tick_cache[f"{name}_{type}_tix"] = getattr(axis, f"{type}Ticks")

    def restore_ticks(self, ax: plt.Axes, axis: str):
        """Restore tick objects from previous render"""

        axes = [ax.xaxis, ax.yaxis] if axis is "both" else [getattr(ax, axis)]
        for axis in axes:
            name = axis.__name__
            for type in ["major", "minor"]:
                getattr(axis, type).locator = self.dge_scatter_tick_cache[f'{name}_{type}_loc']
                setattr(axis, f"{type}Ticks", self.dge_scatter_tick_cache[f'{name}_{type}_tix'])

    def cache_rendered_axis(self, ax: plt.Axes, axis: str):
        """Caches the current axis as rendered (spine, ticks, and labels) as raw 2D rgb data

        Note: this is a bitmap copy so artist information will be lost when the cached copy
        is restored. Regardless of backend one selects, they will all ultimately use the same
        pdf backend to .savefig(). Rerendering the entire canvas, including expensive minor
        ticks, is therefore unavoidable. For now it makes more sense to reuse locators and
        their associated major and minor tick objects.
        """

        assert "agg" in mpl.get_backend().lower(), "Render caching requires an Agg backend."
        canvas = ax.figure.canvas
        rr = canvas.renderer
        canvas.draw()

        if axis in ['x', 'both']:
            box = ax.spines['bottom'].get_tightbbox(rr)
            box = Bbox.from_bounds(x0=0, y0=0, width=ax.figure.bbox.width, height=box.y1)
            self.dge_scatter_tick_cache['x'] = (box, canvas.copy_from_bbox(box))
        if axis in ['y', 'both']:
            box = ax.spines['left'].get_tightbbox(rr)
            box = Bbox.from_bounds(x0=0, y0=0, width=box.x1, height=ax.figure.bbox.height)
            self.dge_scatter_tick_cache['y'] = (box, canvas.copy_from_bbox(box))

    def restore_rendered_axis(self, ax: plt.Axes, axis: str):
        """Restores the cached 2D rgb region from a previous plot render"""

        assert "agg" in mpl.get_backend().lower(), "Render caching requires an Agg backend."

        canvas = ax.figure.canvas
        if axis in ['x', 'both']:
            box, rendered_axis = self.dge_scatter_tick_cache['x']
            canvas.restore_region(rendered_axis)
            canvas.blit(box)
        if axis in ['y', 'both']:
            box, rendered_axis = self.dge_scatter_tick_cache['y']
            canvas.restore_region(rendered_axis)
            canvas.blit(box)

    def box_the_artist(self, artist):  # FIGHT!
        ax = artist.axes
        box = self.get_artist_bbox(artist)
        self.draw_bbox_rectangle(ax, box)

    def get_artist_bbox(self, artist: plt.Artist) -> Optional[Bbox]:
        """Attempts to obtain the rectangular coordinates from an artist"""

        ax = artist.axes
        rr = ax.figure.canvas.renderer

        if hasattr(artist, "get_window_extent") \
                and np.any(artist.get_window_extent(rr)) \
                and np.isfinite(artist.get_window_extent(rr)).all():
            box = artist.get_window_extent(rr)
        elif hasattr(artist, "get_tightbbox") and artist.get_tightbbox(rr) is not None:
            box = artist.get_tightbbox(rr)
        elif hasattr(artist, "clipbox") and artist.clipbox is not None:
            box = artist.clipbox
        else:
            print("Couldn't obtain a Bbox for this artist.")
            return

        return box

    @staticmethod
    def draw_bbox_rectangle(ax: plt.Axes, box: Bbox):
        """Draws a green rectangle around the specified Bbox"""

        rect = Rectangle(xy=(box.x0, box.y0), width=box.width, height=box.height,
                                  transform=ax.get_transform(), clip_on=False, color="green", fill=False)
        ax.add_patch(rect)

    def init_subplots(self, debug):
        """Create one subplot per plot type to reuse between calls"""

        fig_args = {
            'class_pie_barh': {'figsize': (8, 4), 'nrows': 1, 'ncols': 2},
            'len_dist_bar': {'figsize': (7, 4)},
            'scatter': {'figsize': (8, 8), 'tight_layout': False}
        }

        subplots = {}
        for plot in fig_args:
            fig, ax = plt.subplots(**fig_args[plot])
            subplots[plot] = {'fig': fig, 'ax': ax}

        if debug: plt.show(block=False)
        return subplots

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
            # Figure for class_charts has 2 subaxes
            for subax in ax: subax.clear()
        else:
            # Clear only the points for scatter plots
            if len(ax.collections): ax.collections.clear()
            else: ax.clear()
        return fig, ax
