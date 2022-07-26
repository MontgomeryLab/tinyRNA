import sys
import unittest
from unittest.mock import patch, call

import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from pkg_resources import resource_filename

import tiny.rna.plotter as plotter
import tiny.rna.plotterlib as lib

class MyTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.stylesheet = resource_filename('tiny', 'templates/tinyrna-light.mplstyle')

    def get_label_width_pairs_from_annotations_mock(self, annotations):
        bar_widths = [i[1]['xycoords'].get_width() for i in annotations.call_args_list]
        bar_labels = [i[0][0] for i in annotations.call_args_list]
        return list(zip(bar_labels, bar_widths))

    """Are class counts properly calculated?"""

    def test_class_counts(self):
        """If a feature has multiple associated classes, its counts should
        be divided by the number of classes before being summed."""

        # Each feature contributes a single count to its listed classes
        raw_counts_df = plotter.tokenize_feature_classes(
            pd.DataFrame.from_dict(
                {('feat1', pd.NA): ['', 'wago', 1, 1, 1],
                 ('feat2', pd.NA): ['', 'csr,wago', 2, 2, 2],
                 ('feat3', pd.NA): ['', 'wago,csr,other', 3, 3, 3]},
                orient='index',
                columns=['Feature Name', 'Feature Class', 'lib1', 'lib2', 'lib3'])
        )

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


if __name__ == '__main__':
    unittest.main()
