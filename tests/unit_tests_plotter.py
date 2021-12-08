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

    """Are class counts properly calculated?"""

    def test_class_counts(self):
        """If a feature has multiple associated classes, its counts should
        be divided by the number of classes before being summed."""

        # Each feature contributes a single count to its listed classes
        counts = {'feat1': ['wago', 1, 1, 1],
                  'feat2': ['csr,wago', 2, 2, 2],
                  'feat3': ['wago,csr,other', 3, 3, 3]}
        raw_counts_df = pd.DataFrame.from_dict(counts, orient='index',
                        columns=['Feature Class', 'lib1', 'lib2', 'lib3'])

        actual = plotter.get_class_counts(raw_counts_df)
        expected = pd.DataFrame.from_dict(
                    {'other': [1.0, 1.0, 1.0],
                     'csr':   [2.0, 2.0, 2.0],
                     'wago':  [3.0, 3.0, 3.0]},
                    orient='index', columns=['lib1', 'lib2', 'lib3'])

        assert_frame_equal(actual, expected, check_like=True)

    """Are class proportion percentages calculated properly?"""

    def test_class_chart_table_percentage(self):
        class_s = pd.Series({'csr': 75, 'wago': 50, 'other': 25})

        with patch.object(lib.plt.Axes, 'table') as table:
            plib = lib.plotterlib(self.stylesheet)
            plib.class_pie_barh(class_s, 100)

        expected = np.array([
            ['csr', '75.00%'],
            ['wago', '50.00%'],
            ['other', '25.00%']
        ], dtype=object)

        np.testing.assert_array_equal(table.call_args[1]['cellText'], expected)

    """Do class proportion percentages have the correct scale?"""

    def test_class_chart_table_percentage_scale(self):
        class_s = pd.Series({'csr': 1})

        with patch.object(lib.plt.Axes, 'table') as table:
            plib = lib.plotterlib(self.stylesheet)
            plib.class_pie_barh(class_s, 3)

        expected = np.array([
            ['csr', '33.33%'],
            ['Unassigned', '66.67%']
        ], dtype=object)

        np.testing.assert_array_equal(table.call_args[1]['cellText'], expected)

    """Do class proportion percentages display "Unassigned" at the correct threshold?"""

    def test_class_chart_table_percentage_unassigned(self):
        """This problem becomes very interesting when you consider floating point precision.
        Percentage values are rounded according to the scale for simplicity, which we assume
        to be 2 for percentages (4 for decimal). With this in mind the threshold value for
        "Unassigned" is the rounding threshold at scale + 1"""

        class_count = 1
        class_s = pd.Series({'csr': class_count})

        # Value here refers to the decimal value of the "Unassigned" category
        above_thresh_value = 0.00005
        above_thresh_total = class_count / (1 - above_thresh_value)
        below_thresh_value = 0.00004
        below_thresh_total = class_count / (1 - below_thresh_value)

        with patch.object(lib.plt.Axes, 'table') as table:
            plib = lib.plotterlib(self.stylesheet)
            plib.class_pie_barh(class_s, above_thresh_total)
            plib.class_pie_barh(class_s, below_thresh_total)

        expected_above_thresh = np.array([
            ['csr', '99.99%'],
            ['Unassigned', '0.01%']
        ])

        expected_below_thresh = np.array([
            ['csr', '100.00%']
        ])

        np.testing.assert_array_equal(table.call_args_list[0][1]['cellText'], expected_above_thresh)
        np.testing.assert_array_equal(table.call_args_list[1][1]['cellText'], expected_below_thresh)

if __name__ == '__main__':
    unittest.main()
