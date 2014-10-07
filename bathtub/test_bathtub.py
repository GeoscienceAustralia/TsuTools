#!/usr/bin/env python

import unittest

class Test_bathtub(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_example_usage(self):
        """
            Check we can run example_usage without errors
        """
        cmd = 'python example_usage.py'
        import os
        os.system(cmd)
        return


if __name__ == "__main__":
    suite = unittest.makeSuite(Test_bathtub, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
