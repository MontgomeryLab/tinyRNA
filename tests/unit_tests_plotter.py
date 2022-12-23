import sys
import unittest
import numpy as np
import pandas as pd

from unittest.mock import patch
from pandas.testing import assert_frame_equal
from pkg_resources import resource_filename

import tiny.rna.plotter as plotter
import tiny.rna.plotterlib as lib

class PlotterTests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.stylesheet = resource_filename('tiny', 'templates/tinyrna-light.mplstyle')

    #====== HELPER METHODS ===================================================

    def get_label_width_pairs_from_annotations_mock(self, annotations):
        bar_widths = [i[1]['xycoords'].get_width() for i in annotations.call_args_list]
        bar_labels = [i[0][0] for i in annotations.call_args_list]
        return list(zip(bar_labels, bar_widths))

    def aqplt_mock(self):
        return patch(
            'tiny.rna.plotter.aqplt',
            lib.plotterlib(
                resource_filename('tiny', 'templates/tinyrna-light.mplstyle'),
                **{'cache_scatter_ticks': False}
            )
        )

    def get_empty_scatter_dge_dataframes(self):
        counts = pd.DataFrame(
            columns=['Feature ID', 'Classifier', 'ConditionA', 'ConditionB']
        ).set_index(['Feature ID', 'Classifier'])

        dge = pd.DataFrame(
            columns=['Feature ID', 'Classifier', 'ConditionA_vs_ConditionB']
        ).set_index(['Feature ID', 'Classifier'])

        return counts, dge

    #====== TESTS =============================================================

    """Are class counts properly calculated?"""

    def test_class_counts(self):
        """If a feature has multiple associated classes, its counts should
        be divided by the number of classes before being summed."""

        # Each feature contributes a single count to its listed classes
        raw_counts_df = pd.DataFrame.from_dict(
            {('feat1', 'wago'):  ['', 1.0, 1.0, 1.0],
             ('feat2', 'csr'):   ['', 1.0, 1.0, 1.0],
             ('feat2', 'wago'):  ['', 1.0, 1.0, 1.0],
             ('feat3', 'wago'):  ['', 1.0, 1.0, 1.0],
             ('feat3', 'csr'):   ['', 1.0, 1.0, 1.0],
             ('feat3', 'other'): ['', 1.0, 1.0, 1.0]},
            orient='index',
            columns=['Feature Name', 'lib1', 'lib2', 'lib3'])
        raw_counts_df.index = pd.MultiIndex.from_tuples(raw_counts_df.index)

        actual = plotter.get_class_counts(raw_counts_df)
        expected = pd.DataFrame.from_dict(
                    {'other': [1.0, 1.0, 1.0],
                     'csr':   [2.0, 2.0, 2.0],
                     'wago':  [3.0, 3.0, 3.0]},
                    orient='index', columns=['lib1', 'lib2', 'lib3'])

        assert_frame_equal(actual, expected, check_like=True)

    """Are proportions properly labeled as percentages?"""

    def test_proportion_chart_percentage_labels(self):
        group_props = pd.Series({'group1': 0.75, 'group2': 0.50, 'group3': 0.25})

        with patch.object(lib.plt.Axes, 'annotate') as annotations:
            plib = lib.plotterlib(self.stylesheet)
            plib.barh_proportion(group_props)

        actual = self.get_label_width_pairs_from_annotations_mock(annotations)
        expected = [
            ('25.00%', 25.0),
            ('50.00%', 50.0),
            ('75.00%', 75.0)
        ]

        self.assertListEqual(actual, expected)

    """Do plotted proportions have the correct percentage label scale?"""

    def test_proportion_chart_percentage_label_scale(self):
        prop_df = plotter.get_proportions_df(
            pd.DataFrame.from_dict({'group1': [1]}, orient='index', columns=['lib1']),
            mapped_totals=pd.Series({'lib1': 3}),
            un="Unassigned"
        )

        with patch.object(lib.plt.Axes, 'annotate') as annotations:
            plib = lib.plotterlib(self.stylesheet)
            plib.barh_proportion(prop_df['lib1'])

        actual = self.get_label_width_pairs_from_annotations_mock(annotations)
        expected = [
            ('33.33%', 33.33),  # group1
            ('66.67%', 66.67)   # unassigned
        ]

        self.assertListEqual(actual, expected)

    """Do proportion charts display "Unassigned" at the correct threshold?"""

    def test_proportion_chart_percentage_unassigned(self):
        """This problem becomes very interesting when you consider floating point precision.
        Percentage values are rounded according to the scale for simplicity, which we assume
        to be 2 for percentages (4 for decimal). With this in mind the threshold value for
        "Unassigned" is the rounding threshold at scale + 1"""

        class_count = 1.0
        df_kwargs = {'orient': 'index', 'columns': ['lib1']}

        # Value refers to the decimal value of the "Unassigned" category
        # Total refers to the total Mapped Reads required to push class_count proportion above/below threshold
        above_thresh_value = 0.00005
        above_thresh_total = class_count / (1 - above_thresh_value)
        above_thresh_series = pd.Series({'lib1': above_thresh_total})

        below_thresh_value = 0.00004
        below_thresh_total = class_count / (1 - below_thresh_value)
        below_thresh_series = pd.Series({'lib1': below_thresh_total})

        class_df = pd.DataFrame.from_dict({'group1': [class_count]}, **df_kwargs)

        # ===== ABOVE THRESHOLD ======================================================================
        actual_above_thresh = plotter.get_proportions_df(class_df, mapped_totals=above_thresh_series, un="Unassigned")
        expected_above_thresh = pd.DataFrame.from_dict({
            'group1':     0.9999,
            'Unassigned': 0.0001
        }, **df_kwargs)

        assert_frame_equal(actual_above_thresh, expected_above_thresh, check_like=True)

        # ===== BELOW THRESHOLD ======================================================================
        actual_below_thresh = plotter.get_proportions_df(class_df, mapped_totals=below_thresh_series, un="Unassigned")
        expected_below_thresh = pd.DataFrame.from_dict({
            'group1':     1.0000,
            'Unassigned': pd.NA  # Ok to have NA values here because they are dropped in barh_proportion()
        }, **df_kwargs)

        assert_frame_equal(actual_below_thresh, expected_below_thresh, check_dtype=False, check_like=True)

    """Does scatter_by_dge do the right thing when DataFrame inputs are empty?"""

    def test_scatter_by_dge_empty_dataframes(self):
        counts, dge = self.get_empty_scatter_dge_dataframes()

        with patch('tiny.rna.plotter.save_plot') as save_plot, self.aqplt_mock():
            plotter.scatter_by_dge(counts, dge, 'dummy_prefix', (0, 0))

        save_plot.assert_not_called()

    """Does scatter_by_dge_class do the right thing when DataFrame inputs are empty?"""

    def test_scatter_by_dge_class_empty_dataframes(self):
        counts, dge = self.get_empty_scatter_dge_dataframes()

        with patch('tiny.rna.plotter.save_plot') as save_plot, self.aqplt_mock():
            plotter.scatter_by_dge_class(counts, dge, 'dummy_prefix', (0, 0))

        save_plot.assert_not_called()

    """Does scatter_by_dge_class properly handle empty inclusive filter lists?"""

    def test_scatter_dge_class_empty_inclusive_filter(self):
        counts, dge = self.get_empty_scatter_dge_dataframes()

        with patch('tiny.rna.plotter.plotterlib.scatter_grouped') as scatter, self.aqplt_mock():
            plotter.scatter_by_dge_class(counts, dge, 'dummy_prefix', (0, 0), include=[])
            scatter.assert_not_called()
            plotter.scatter_by_dge_class(counts, dge, 'dummy_prefix', (0, 0), exclude=[])
            scatter.assert_not_called()

    """Do scatter plots show the appropriate major ticks through a range of view limits?"""

    @unittest.skip("Long-running test, execute manually if needed")
    def test_scatter_major_ticks(self):
        counts, dge = self.get_empty_scatter_dge_dataframes()
        min_non_zero = 1 / sys.maxsize  # avoid zero on log scale
        fps = 3

        counts.loc[('featA', 'featClassA'), 'ConditionA'] = min_non_zero
        counts.loc[('featA', 'featClassA'), 'ConditionB'] = min_non_zero
        dge.loc[('featA', 'featClassA'), 'ConditionA_vs_ConditionB'] = 0

        for x in range(0, 121):
            with self.aqplt_mock():
                x /= fps  # Range only produces integer values, but we want fractional powers of 2 in the demo

                # Rolling view limit window
                lo_bound = 2**(-6 + (x/2))  # Walk lower bound slowly forward from 2 decimal log2 minimum
                hi_bound = 2**x + x         # Walk upper bound forward much faster
                vlim = np.array((lo_bound, hi_bound))

                # title = f"Range: 2^{int(np.log2(view_lims[0]))} .. 2^{np.log2(view_lims[1]):.1f}"
                # ^ must be set within scatter_* functions in plotter.py, not worth refactoring to support
                plotter.scatter_by_dge(counts, dge, f'lim_{x:.2f}', vlim)



if __name__ == '__main__':
    unittest.main()
