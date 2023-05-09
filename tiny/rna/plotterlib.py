""" Plotting functions for small RNA data. 

This module contains functions to create relevant plots for small RNA data for use
with the tinyRNA pipeline. The plots are built using matplotlib and our style sheet.
You may override these styles by obtaining a copy of the style sheet (tiny get-templates),
modifying it, and passing it to tiny-plot via the -s/--style-sheet argument. If
using this module directly, it may be passed at construction time.
"""
import pandas as pd
import numpy as np
import itertools
import locale
import math
import sys
import os
import re

# This has to be done before importing matplotlib.pyplot
# cwltool appears to unset all environment variables including those related to locale
# This leads to warnings from plt's FontConfig manager, but only for pipeline/cwl runs
curr_locale = locale.getlocale()
if curr_locale[0] is None:
    # Default locale otherwise unset
    os.environ['LC_CTYPE'] = 'en_US.UTF-8'

import matplotlib as mpl; mpl.use("PDF")
import matplotlib.pyplot as plt
import matplotlib.ticker as tix
from matplotlib.patches import Rectangle
from matplotlib.transforms import Bbox
from matplotlib.scale import LogTransform

from typing import Union, Tuple, List, Optional
from abc import ABC, abstractmethod

from tiny.rna.util import sorted_natural


class plotterlib:

    def __init__(self, user_style_sheet, **prefs):

        self.debug = self.is_debug_mode()
        if self.debug:
            mpl.use("TkAgg", force=True)
            mpl.rcParams['savefig.dpi'] = 100

        # Set global plot style once
        plt.style.use(user_style_sheet)

        self.prefs = prefs
        self.subplot_cache = {}
        self.dge_scatter_tick_cache = {}

    def len_dist_bar(self, size_prop: pd.DataFrame, subtype: str, **kwargs) -> plt.Axes:
        """Creates a stacked barplot of 5' end nucleotides by read length

        Args:
            size_prop: A dataframe of size x 5'nt read count proportions
            subtype: The subtype of this len_dist plot so the title can be properly set
            kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()

        Returns:
            sizeb: A stacked barplot of size + 5'nt data
        """

        # Retrieve axis and styles for this plot type
        fig, ax = self.reuse_subplot("len_dist")

        # Ensure xaxis tick labels won't be too crowded
        font_size = self.get_xtick_labelsize_for_axis(ax.xaxis, size_prop.index)

        # Override default colors. User may override with kwargs. (Orange, Yellow-green, Blue, Pink, Gray)
        colors = {'axes.prop_cycle': mpl.cycler(color=['#F78E2D', '#CBDC3F', '#4D8AC8', '#E06EAA', '#B3B3B3'])}

        # The style context allows us to use temporary styles
        with plt.style.context(colors):
            plt.sca(ax)
            sizeb = size_prop.plot(kind='bar', stacked=True, reuse_plot=True, **kwargs)
            sizeb.tick_params(axis='x', labelsize=font_size, rotation=0)
            sizeb.set_ylim(0, np.max(np.sum(size_prop, axis=1)) + 0.025)
            sizeb.set_title(f'Distribution of {subtype} reads')
            sizeb.set_ylabel('Proportion of reads')
            sizeb.set_xlabel('Length of sequence')

        return sizeb

    def pie_proportion(self, prop_ds: pd.Series, ax: plt.Axes, scale=2, **kwargs) -> plt.Axes:
        """Creates a pie chart of read proportions

        Args:
            prop_ds: A pandas Series with proportions to plot
            ax: The subplot axes to use for the pie chart
            scale: The desired number of *percentage* decimal places for labels
            kwargs: Additional keyword arguments to pass to pandas.DataFrame.plot()

        Returns:
            cpie: A pie chart axes object
        """

        # Create the plot
        cpie = prop_ds.plot(kind='pie', normalize=True, ax=ax, **kwargs)
        cpie.legend(
            loc='lower center',
            bbox_to_anchor=(0.5, -0.3),
            fontsize=10,
            ncol=math.ceil(len(prop_ds) / 5),
            labels=[f"{c} ({v:.{scale}f}%)" for c,v in prop_ds.to_dict().items()])
        cpie.set_aspect("equal")
        cpie.set_ylabel('')
        cpie.set_xlabel('')

        return cpie

    def barh_proportion(self, prop_ds: pd.Series, max_prop=1.0, scale=2, **kwargs) -> plt.Axes:
        """Creates a horizontal bar chart of read proportions (as percentages) for a library

        Args:
            prop_ds: A pandas Series with proportions to chart
            max_prop: The xaxis upper bound as a float between 0 and 1
            scale: The desired number of *percentage* decimal places for labels
            kwargs: Additional keyword arguments to pass to pandas.Series.plot()

        Returns:
            cbar: A horizontal bar chart axes object
        """

        # Retrieve axis and styles for this plot type
        fig, ax = self.reuse_subplot("proportions")
        prop_ds = prop_ds.dropna()[::-1]

        # Set bar and text colormaps
        bar_colors = kwargs.get('colors', plt.get_cmap("viridis")(prop_ds))
        text_colors = plt.get_cmap("Greys")

        # Create the plot and set plot attributes
        cbar = (prop_ds * 100).plot(kind='barh', ax=ax, color=bar_colors, **kwargs)
        cbar.set_xlabel('Percentage of reads')
        cbar.set_xlim(0, min([(max_prop * 100) + 10, 100]))

        # Remove irrelevant plot attributes
        cbar.legend().set_visible(False)
        cbar.grid(False, axis='y')

        # For converting data coordinates to axes fraction coordinates
        frac_tf = cbar.transData + cbar.transAxes.inverted()

        # Place percentage annotations
        for bar in cbar.patches:
            width = bar.get_width()
            height = bar.get_height()
            axes_frac = frac_tf.transform((width, 0))[0]
            fontsize = min(self.data_val_to_points(cbar, height), plt.rcParams['axes.labelsize'])
            x_offset = 0.8 * fontsize if axes_frac < 0.8 else -4 * fontsize
            cbar.annotate(
                f"{width:.{scale}f}%",
                (1.0, 0.5), xycoords=bar,  # anchor text to bar
                xytext=(x_offset, 0), textcoords='offset points',
                color=text_colors(width / 100 + 0.4),
                fontproperties={'size': fontsize},
                va='center_baseline',
            )

        return cbar

    def scatter_simple(self, count_x: pd.Series, count_y: pd.Series, log_scale=True, **kwargs) -> plt.Axes:
        """Creates a simple scatter plot of counts.

        Args:
            count_x: A pandas dataframe/series of counts per feature (X axis)
            count_y: A pandas dataframe/series of counts per feature (Y axis)
            log_scale: Plot on log scale rather than linear scale
            kwargs: Additional keyword arguments to pass to pyplot.Axes.scatter()

        Returns:
            ax: A simple scatter plot of counts
        """

        # Retrieve axis and styles for this plot type
        fig, ax = self.reuse_subplot("scatter")
        ax: plt.Axes

        # log2 normalize data if requested
        if log_scale:
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

    def scatter_grouped(self, count_x: pd.DataFrame, count_y: pd.DataFrame, *groups, colors: dict, pval=0.05,
                        view_lims: Tuple[float, float] = None, labels=None, **kwargs):
        """Creates a scatter plot with different groups highlighted.

        Args:
            count_x: A pandas dataframe/series of counts per feature (X axis)
            count_y: A pandas dataframe/series of counts per feature (Y axis)
            groups: An iterable of lists each representing a group
            labels: An iterable of labels corresponding to each group
            colors: A dictionary of label -> color for each group
            view_lims: Optional plot view limits as tuple(min, max)
            pval: The p-value to use for the outgroup label

        Keyword Args:
            log_scale: Data is plotted on log scale if true (default: true)
            kwargs: Additional arguments to pass to pyplot.Axes.scatter()

        Returns:
            gscat: A scatter plot containing groups highlighted with different colors
        """

        # Subset counts not in *groups (e.g., p-val above threshold)
        count_x_out = count_x.drop(itertools.chain(*groups))
        count_y_out = count_y.drop(itertools.chain(*groups))
        has_outgroup = all(co.replace(0, pd.NA).dropna().any()
                           for co in (count_x_out, count_y_out))

        # Make all counts log-compatible, or drop the group if it is zero in both conditions
        plottable_groups = self.nonzero_group_indexes(count_x, count_y, groups, view_lims)
        plot_labels = [labels[i] for i in plottable_groups]
        plot_groups = [groups[i] for i in plottable_groups]
        group_it = iter(plot_groups)

        if has_outgroup:
            x, y = count_x_out, count_y_out
            gscat = self.scatter_simple(x, y, color='#B3B3B3', **kwargs)
        elif plottable_groups:
            group = next(group_it)
            x, y = count_x.loc[group], count_y.loc[group]
            gscat = self.scatter_simple(x, y, **kwargs)
        else:
            has_outgroup = None
            x = y = pd.Series(dtype='float64')
            gscat = self.scatter_simple(x, y, **kwargs)

        # Add remaining groups
        for group in group_it:
            x, y = count_x.loc[group], count_y.loc[group]
            gscat.scatter(x, y, edgecolor='none', **kwargs)

        self.sort_point_groups_and_label(gscat, plot_groups, plot_labels, colors, has_outgroup, pval)
        self.set_square_scatter_view_lims(gscat, view_lims)
        self.set_scatter_ticks(gscat)

        return gscat

    @staticmethod
    def nonzero_group_indexes(count_x, count_y, groups, view_lims):
        """When scatter plotting features on a log scale, if the feature has a count of
        zero in either condition then it is omitted from the plot by default due to the
        singularity at 0, but we want to represent them nonetheless. So, we set their count
        in the zero condition to the lower plot limit so that they are plotted on the very
        edge of the plot space. Otherwise, approximating zero for these features would
        shrink the plot.

        We still want to omit features that are zero in both conditions, and if an entire
        group consists of zero counts, its label should be omitted from the legend. This
        is accomplished by omitting its index from the returned list."""

        non_zero_groups = []
        minpos = min(view_lims)
        for i, group in enumerate(groups):
            x, y = count_x.loc[group], count_y.loc[group]
            x_is_zeros = x.replace(0, pd.NA).dropna().empty
            y_is_zeros = y.replace(0, pd.NA).dropna().empty

            if not (x_is_zeros and y_is_zeros):
                non_zero_groups.append(i)

                # Replace counts that are zero in only one condition
                x.loc[(x == 0) & (y != 0)] = minpos
                y.loc[(y == 0) & (x != 0)] = minpos
                count_x.loc[group] = x
                count_y.loc[group] = y

        return non_zero_groups

    @staticmethod
    def sort_point_groups_and_label(axes: plt.Axes, groups, labels, colors, outgroup: Optional[bool], pval):
        """Sorts scatter groups so that those with fewer points are rendered on top of the stack.
        After sorting, group colors and labels are assigned, and the legend is created. Labels
        in the legend are sorted by natural order with the outgroup always listed last.
            Args:
                axes: The scatter plot Axes object
                groups: A list of DataFrames that were able to be plotted
                labels: A list of names, one for each group, for the corresponding index in `groups`
                colors: A dictionary of group labels and their assigned colors
                outgroup: True if an out group was plotted, None if empty plot (no groups or out groups)
        """

        lorder = np.argsort([len(grp) for grp in groups if len(grp)])[::-1]   # Index of groups by size
        offset = int(bool(outgroup))
        layers = axes.collections

        if outgroup:
            layers[0].set_label('p ≥ %g' % pval)
        if labels is None:
            labels = list(range(len(groups)))
        if outgroup is None:
            return

        groupsize_sorted = [(labels[i], layers[i + offset]) for i in lorder]
        for z, (label, layer) in enumerate(groupsize_sorted, start=offset+1):
            layer.set_label(re.sub(r'^_', ' _', label))                       # To allow labels that start with _
            layer.set_facecolor(colors[label])
            layer.set_zorder(z)                                               # Plot in order of group size

        # Ensure lines remain on top of points
        for line in axes.lines:
            line.set_zorder(len(layers) + 1)

        # Sort the legend with outgroup last while retaining layer order
        handles = sorted_natural(layers[offset:], key=lambda x: x.get_label())
        if outgroup: handles.append(layers[0])
        axes.legend(handles=handles)

    @staticmethod
    def assign_class_colors(classes):
        """Assigns a color to each class for consistency across samples"""

        stylesheet_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        if len(classes) <= len(stylesheet_colors):
            colors = iter(stylesheet_colors)
        else:
            colors = iter(plt.get_cmap("tab20"))

        return {cls: next(colors) for cls in sorted_natural(classes)}

    def set_dge_class_legend_style(self):
        """Widens the "scatter" figure and moves plot to the left to accommodate legend"""

        expand_width_inches = 3

        fig, scatter = self.reuse_subplot("scatter")
        transFigure = fig.transFigure
        orig_axes_pos = scatter.get_position().transformed(transFigure)
        orig_fig_size = fig.get_size_inches()

        # Expand the figure and move the plot back to its original position at left
        fig.set_size_inches(orig_fig_size[0] + expand_width_inches, orig_fig_size[1])
        scatter.set_position(orig_axes_pos.transformed(transFigure.inverted()))

    @staticmethod
    def get_scatter_view_lims(counts_df: pd.DataFrame, vmin: int = None, vmax: int = None) -> Tuple[float, float]:
        """Calculates scatter view limits for the counts dataframe

        Args:
            counts_df: A pandas dataframe of counts per feature
            vmin: Optional log2 minimum view limit
            vmax: Optional log2 maximum view limit
        """

        # For transforming values to/from log2 scale
        transform = LogTransform(base=2)
        inverse_trans = transform.inverted()

        # User-specified min & max, no calculation necessary
        if (vmin, vmax) != (None, None):
            return inverse_trans.transform([vmin, vmax])

        # Obtain the minimum and maximum counts from the counts dataframe
        x0 = counts_df.replace(0, pd.NA).min(axis="columns").dropna().min()
        x1 = counts_df.max().max()
        minpos = 1e-300

        if not np.isfinite([x0, x1]).all() or not isinstance(x0, np.float) or x1 <= 0:
            print("The provided dataset contains invalid values.", file=sys.stderr)
            return (minpos, minpos)

        # Avoid log2(0) errors
        x0, x1 = (minpos if x0 <= 0 else x0,
                  minpos if x1 <= 0 else x1)

        # Get axes margin preferences from style sheet
        rc_mar = {mpl.rcParams.get(f"axes.{m}", 0)
                  for m in ('xmargin', 'ymargin')}

        margin = max(rc_mar)
        if len(rc_mar) != 1:
            print("Style sheet values for axes.xmargin and axes.ymargin differ. "
                  "The larger value will be chosen for the scatter plot margin.",
                  file=sys.stderr)

        # Calculate plot margins
        x0t, x1t = transform.transform([x0, x1])
        delta = (x1t - x0t) * margin
        if not np.isfinite(delta): delta = 0

        if vmin is None: vmin = x0t - delta
        if vmax is None: vmax = x1t + delta
        return inverse_trans.transform([vmin, vmax])

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
        ceil, log2 = math.ceil, np.log2
        locs = [2 ** x for x in range(ceil(log2(ax_min)), ceil(log2(ax_max)))]
        return locs, ax_min, ax_max

    def set_scatter_ticks(self, ax: plt.Axes, minor_ticks=False):
        """Intelligently creates major and minor ticks for a square scatter plot while avoiding crowding"""

        # Get tick locations corresponding to the current view limits
        major_locs, ax_min, ax_max = self.get_fixed_majorticklocs(ax.viewLim.extents)

        ax.xaxis.set_major_locator(tix.FixedLocator(major_locs))
        ax.yaxis.set_major_locator(tix.FixedLocator(major_locs))

        if len(self.dge_scatter_tick_cache):
            self.restore_ticks(ax, "both")
            return

        for axis in [ax.xaxis, ax.yaxis]:
            # Only display every nth major tick label
            n = int(np.log2(len(major_locs)) - 1) or 1
            ticks_displayed, last_idx = self.every_nth_label(axis, n)

            if minor_ticks:
                axis.set_minor_locator(tix.LogLocator(
                    base=2.0,
                    numticks=self.get_min_LogLocator_numticks(axis),
                    subs=np.log2(np.linspace(2 ** 2, 2 ** 4, 10))[:-1]))

            min_tick = ax_min
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

        if self.prefs.get('cache_scatter_ticks', True):
            for type in ["major", "minor"]:
                self.dge_scatter_tick_cache[f"{name}_{type}_loc"] = getattr(axis, type).locator
                self.dge_scatter_tick_cache[f"{name}_{type}_tix"] = getattr(axis, f"{type}Ticks")

    def restore_ticks(self, ax: plt.Axes, axis: str):
        """Restore tick objects from previous render"""

        axes = [ax.xaxis, ax.yaxis] if axis == "both" else [getattr(ax, axis)]
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
    def data_val_to_points(ax, val):
        """Converts a length from data coordinates (relative coords) to points (absolute coords)"""

        # Map a (1, 1) point from data coordinates to display coordinates
        dd = ax.transData.transform([[0, 0], [1, 1]])  # Units: px

        # Calculate the data/points coefficient
        co = ax.get_figure().get_dpi() / round(dd[1, 1] - dd[0, 1]) / 72  # Units: data/pts

        return val / co  # Units: data * pts/data = pts

    @staticmethod
    def get_xtick_labelsize_for_axis(axis: mpl.axis.XAxis, index: pd.Index, min_size=5) -> int:
        """Calculates a new labelsize for the xaxis if the default size will cause crowding

        Args:
            axis: The xaxis to fit
            index: The index representing plotted xaxis values
            min_size: The minimum acceptable font size

        Returns: an adjusted fontsize if necessary, otherwise the default font size
        """

        # Get width of axis in points
        ax_length = Bbox.from_bounds(0, 0, 1, 1) \
                        .transformed(axis.axes.transAxes - axis.figure.dpi_scale_trans) \
                        .width * 72

        # Get current xtick labelsize in points
        default_size = int(axis._get_tick_label_size('x'))

        # Assume text for values <100 will have aspect ratio of 1.1:1
        # Solve for size: len(index) = ax_length / (size * 1.1)
        size = int(np.floor(ax_length / (len(index) * 1.1)))

        if size > default_size:
            size = default_size
        elif size < min_size:
            size = min_size
            print(f"WARNING: minimum font size ({min_size}) reached while attempting "
                  "to reduce xaxis tick label crowding.", file=sys.stderr)

        return size

    @staticmethod
    def draw_bbox_rectangle(ax: plt.Axes, box: Bbox):
        """Draws a green rectangle around the specified Bbox"""

        rect = Rectangle(xy=(box.x0, box.y0), width=box.width, height=box.height,
                                  transform=ax.get_transform(), clip_on=False, color="green", fill=False)
        rect.set_in_layout(False)
        ax.add_patch(rect)

    @staticmethod
    def is_debug_mode():
        return getattr(sys, 'gettrace', lambda: None)() is not None

    def reuse_subplot(self, plot_type: str) -> Tuple[plt.Figure, Union[plt.Axes, List[plt.Axes]]]:
        """Retrieves the reusable subplot for this plot type

        Args:
            plot_type: The reusable plot name

        Returns:
            fig: The subplot's pyplot.Figure
            ax: The subplot's pyplot.Axes, or an array of axes if subplot's nrows or ncols is >1
        """

        if plot_type in self.subplot_cache:
            cache = self.subplot_cache[plot_type]
        else:
            self.subplot_cache[plot_type] = cache = {
                "proportions": ClassChartCache,
                "len_dist": LenDistCache,
                "scatter": ScatterCache
            }[plot_type]()

        if self.debug: plt.show(block=False)
        return cache.get()


class CacheBase(ABC):
    @abstractmethod
    def get(self): pass

    def __del__(self):
        plt.close(self.fig)


class ClassChartCache(CacheBase):
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(8, 6))

    def get(self) -> Tuple[plt.Figure, plt.Axes]:
        self.ax.clear()
        return self.fig, self.ax


class LenDistCache(CacheBase):
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(8, 4))

    def get(self) -> Tuple[plt.Figure, plt.Axes]:
        self.ax.clear()
        return self.fig, self.ax


class ScatterCache(CacheBase):
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(8, 8), tight_layout=False)
        self.ax.set_aspect('equal')

    def get(self) -> Tuple[plt.Figure, plt.Axes]:
        if len(self.ax.collections): self.ax.collections.clear()
        return self.fig, self.ax