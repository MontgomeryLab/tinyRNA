#!/usr/bin/env python
""" Testing module for aquatx. """

import unittest
import aquatx

class test_placeholder(unittest.TestCase):
    """ 
    Currently creating a placeholder to have TravisCI set up while
    I create the tests. Old tests had to be replaced as the scripts
    were revamped quite a bit so I decided to wait until things were
    more complete and cleaned up before writing a bunch of tests only
    to have to delete them all and start again. 
    
    So all 'tests' will pass for now, see unit_tests_counter.py for old
    tests that are no longer run, but left for reference when I make new
    ones. This is purely for initial public repository set up and to
    have Travis set up for eventual testing.

    """
    def test_math(self):
        self.assertEqual(2+2, 4)

if __name__ == '__main__':
    unittest.main()
