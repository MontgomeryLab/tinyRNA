import unittest
from unittest.mock import patch, mock_open

from rna.counter.validation import GFFValidator, ReportFormatter


class ValidationTests(unittest.TestCase):

    """======== Helper functions =========================== """

    def make_gff_row(self, **cols):
        template = {
            'seqid': "I",
            'source': "Wormbase",
            'type': "gene",
            'start': 3747,
            'end': 3909,
            'score': ".",
            'strand': "+",
            'phase': ".",
            'attrs':  {
                'ID': 'featid',
                'gene_id': 'featid',
                'Parent': 'parentid'
            }
        }

        new_gff = dict(template, **cols)
        new_gff['attrs'] = ';'.join([f"{k}={v}" for k, v in new_gff['attrs'].items()])
        return '\t'.join(map(str, new_gff.values()))

    def make_gff_validator(self) -> GFFValidator:
        return GFFValidator({}, {})

    def mock_gff_open(self, contents):
        mo = mock_open(read_data=contents)
        return patch('tiny.rna.counter.hts_parsing.HTSeq.utils.open', mo)

    """Does GFFValidator correctly validate strand information?"""

    def test_gff_strand_validation(self):
        validator = self.make_gff_validator()
        mock_filename = '/dev/null'
        mock_gff = '\n'.join([
            self.make_gff_row(strand="+"),                      # Valid
            self.make_gff_row(strand="-"),                      # Valid
            self.make_gff_row(strand=".", type="chromosome"),   # Valid
            self.make_gff_row(strand=".")                       # Invalid
        ])

        expected = '\n'.join([
            "The following issues were found in the GFF files provided:",
            "\t" + f"{mock_filename}: ",
            "\t\t" + f"{validator.targets['strand']}: 1"
        ])

        with self.mock_gff_open(mock_gff) as p:
            validator.parse_and_validate_gffs({mock_filename: []})

        self.assertListEqual(validator.report.errors, [expected])
        self.assertListEqual(validator.report.warnings, [])

    """Does GFFValidator correctly validate ID attributes?"""

    def test_gff_id_validation(self):
        validator = self.make_gff_validator()
        mock_filename = '/dev/null'
        mock_gff = '\n'.join([
            self.make_gff_row(attrs={'ID':      'feat1'}),  # Valid
            self.make_gff_row(attrs={'gene_id': 'feat3'}),  # Valid
            self.make_gff_row(attrs={'Parent':  'feat5'}),  # Valid
            self.make_gff_row(attrs={'Gene id': 'feat6'}),  # Invalid
            self.make_gff_row(attrs={'other':   'feat7'}),  # Invalid
            self.make_gff_row(attrs={}),                    # Invalid
        ])

        expected = '\n'.join([
            "The following issues were found in the GFF files provided:",
            "\t" + f"{mock_filename}: ",
            "\t\t" + f"{validator.targets['ID attribute']}: 3"
        ])

        with self.mock_gff_open(mock_gff) as p:
            validator.parse_and_validate_gffs({mock_filename: []})

        self.assertListEqual(validator.report.errors, [expected])
        self.assertListEqual(validator.report.warnings, [])

    """Does GFFValidator.chroms_shared_with_ebwt() correctly identify ebwt chromosomes?"""

    def test_ebwt_chroms(self):
        validator = self.make_gff_validator()
        ebwt_prefix = "./testdata/counter/validation/ebwt/ram1"

        # Chroms are shared
        validator.chrom_set = {'ram1'}
        shared, ebwt_chroms = validator.chroms_shared_with_ebwt(ebwt_prefix)
        self.assertSetEqual(shared, validator.chrom_set)
        self.assertSetEqual(shared, ebwt_chroms)

        # Chroms aren't shared
        validator.chrom_set = {'chr1', 'chr2', 'chr3'}
        shared, ebwt_chroms = validator.chroms_shared_with_ebwt(ebwt_prefix)
        self.assertSetEqual(shared, set())
        self.assertSetEqual(ebwt_chroms, {'ram1'})

    """Does GFFValidator.chroms_shared_with_genomes() correctly identify genome chromosomes?"""

    def test_genome_chroms(self):
        validator = self.make_gff_validator()
        fasta_file = "./testdata/counter/validation/genome/genome.fasta"

        # Chroms are shared
        validator.chrom_set = {'chr1', 'chr2', 'chr3'}
        shared, genome_chroms = validator.chroms_shared_with_genomes([fasta_file])
        self.assertSetEqual(shared, validator.chrom_set)
        self.assertSetEqual(shared, genome_chroms)

        # Chroms aren't shared
        validator.chrom_set = {'ram1'}
        shared, genome_chroms = validator.chroms_shared_with_genomes([fasta_file])
        self.assertSetEqual(shared, set())
        self.assertSetEqual(genome_chroms, {'chr1', 'chr2', 'chr3'})

    """Does GFFValidator.generate_gff_report() correctly process an infractions report?"""

    def test_gff_report_output(self):
        validator = self.make_gff_validator()
        infractions = {
            "gff1": {'ID attribute': 10, 'strand': 5},
            "gff2": {'ID attribute': 1},
            "gff3": {'strand': 1},
            "gff4": {}
        }

        expected = '\n'.join([
            "The following issues were found in the GFF files provided:",
            "\tgff1: ",
            "\t\t" + f"{validator.targets['ID attribute']}: 10",
            "\t\t" + f"{validator.targets['strand']}: 5",
            "\tgff2: ",
            "\t\t" + f"{validator.targets['ID attribute']}: 1",
            "\tgff3: ",
            "\t\t" + f"{validator.targets['strand']}: 1"
        ])

        validator.generate_gff_report(infractions)
        self.assertListEqual(validator.report.errors, [expected])

    def test_chrom_report_output(self):
        validator = self.make_gff_validator()
        validator.chrom_set = {'chr1', 'chr2'}
        seq_chroms = {'chr3', 'chr4'}
        shared_chroms = validator.chrom_set & seq_chroms

        exp_errors = '\n'.join([
            "GFF files and sequence files don't share any chromosome identifiers.",
            "\tChromosomes are present in GFF files: ",
            "\t\tchr1",
            "\t\tchr2",
            "\tThe following chromosomes are present in sequence files: ",
            "\t\tchr3",
            "\t\tchr4"
        ])

        validator.generate_chrom_report(shared_chroms, seq_chroms)
        # self.assertEqual(validator.report.errors[0], exp_errors)
        validator.report.errors *= 3
        validator.report.print_report()

if __name__ == '__main__':
    unittest.main()
