import unittest

from aquatx.srna.Configuration import Configuration


class MyTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.file = "./testdata/run_config_template.yml"

    def test_init(self):
        result = Configuration(self.file)
        result.write_config_file()

if __name__ == '__main__':
    unittest.main()
