import unittest

from aquatx.srna.Configuration import Configuration


class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.file = "./testdata/run_config_template.yml"

    def test_init(self):
        result = Configuration(self.file)
        result.write_processed_config()

    def test_rundir(self):
        from aquatx.aquatx import run
        import os

        if os.path.isdir('run_directory'): os.removedirs('run_directory')
        run("../aquatx/cwl", "./testdata/run_config_template.yml")

if __name__ == '__main__':
    unittest.main()
