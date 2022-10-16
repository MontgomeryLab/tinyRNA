import contextlib
import unittest
import time
import io

from glob import glob
from unittest.mock import patch, mock_open

from rna.counter.validation import GFFValidator, ReportFormatter


class ValidationTests(unittest.TestCase):

    """======== Helper functions =========================== """

    @classmethod
    def setUpClass(self):
        self.strand_ext_header = \
            'Unstranded features are allowed, but they can lead to potentially unexpected results.\n' \
            'These features will match "sense", "antisense", and "both" strand selectors. 5\'/3\' anchored\n' \
            "overlap selectors for these features will evaluate for termini shared with the alignment,\n" \
            "but will not distinguish between the alignment's 5' and 3' ends."

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
            "The following issues were found in the GFF files provided. ",
            self.strand_ext_header,
            "\t" + f"{validator.targets['strand']}: ",
            "\t\t" + f"1 missing in {mock_filename}"
        ])

        with self.mock_gff_open(mock_gff) as p:
            validator.parse_and_validate_gffs({mock_filename: []})

        self.assertListEqual(validator.report.warnings, [expected])
        self.assertListEqual(validator.report.errors, [])

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
            "The following issues were found in the GFF files provided. ",
            "\t" + f"{validator.targets['ID attribute']}: ",
            "\t\t" + f"3 missing in {mock_filename}"
        ])

        with patch('tiny.rna.counter.hts_parsing.HTSeq.utils.open', mock_open(read_data=mock_gff)) as p:
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

    """Does GFFValidator's alignments heuristic identify potentially incompatible SAM files?"""

    def test_alignments_heuristic(self):
        validator = self.make_gff_validator()
        sam_files = ['./testdata/counter/identity_choice_test.sam',
                     './testdata/counter/single.sam']

        sam_chroms = {
            sam_files[0]: {'I', 'V', 'MtDNA'},
            sam_files[1]: {'I'}
        }

        # Some chroms are shared
        validator.chrom_set = {'I', 'not_shared', 'also_not_shared'}
        bad_sams = validator.alignment_chroms_mismatch_heuristic(sam_files)
        self.assertDictEqual(bad_sams, {})

        # Some chroms aren't shared
        validator.chrom_set = {'V', 'not_shared', 'also_not_shared'}
        bad_sams = validator.alignment_chroms_mismatch_heuristic(sam_files)
        self.assertDictEqual(bad_sams, {sam_files[1]: sam_chroms[sam_files[1]]})

        # Chroms aren't shared
        validator.chrom_set = {'not_shared', 'also_not_shared'}
        bad_sams = validator.alignment_chroms_mismatch_heuristic(sam_files)
        self.assertDictEqual(bad_sams, sam_chroms)

        # Chroms aren't shared, single SAM file input
        validator.chrom_set = {'V'}
        sam_file = sam_files[1]
        bad_sams = validator.alignment_chroms_mismatch_heuristic([sam_file])
        self.assertDictEqual(bad_sams, {sam_file: sam_chroms[sam_file]})

    """Does GFFValidator.generate_gff_report() correctly process an infractions report?"""

    def test_gff_report_output(self):
        validator = self.make_gff_validator()
        infractions = {
            'ID attribute': {'gff1': 10, 'gff2': 1},
            'strand': {'gff1': 5, 'gff3': 1}
        }

        gff_sections_header = "The following issues were found in the GFF files provided. "
        expected_errors = '\n'.join([
            gff_sections_header,
            "\t" + f"{validator.targets['ID attribute']}: ",
            "\t\t" + "10 missing in gff1",
            "\t\t" + "1 missing in gff2"
        ])

        expected_warnings = '\n'.join([
            '\n'.join([gff_sections_header, self.strand_ext_header]),
            "\t" + f"{validator.targets['strand']}: ",
            "\t\t" + "5 missing in gff1",
            "\t\t" + "1 missing in gff3"
        ])

        validator.generate_gff_report(infractions)
        self.assertEqual(validator.report.errors[0], expected_errors)
        self.assertEqual(validator.report.warnings[0], expected_warnings)

    """Does GFFValidator.generate_chrom_report() correctly process an infractions report?"""

    def test_chrom_report_output(self):
        validator = self.make_gff_validator()
        validator.chrom_set = {'chr1', 'chr2'}
        seq_chroms = {'chr3', 'chr4'}
        shared_chroms = validator.chrom_set & seq_chroms

        exp_errors = '\n'.join([
            "GFF files and sequence files don't share any chromosome identifiers.",
            "\t" + f"{validator.targets['gff chromosomes']}: ",
            "\t\t" + "chr1",
            "\t\t" + "chr2",
            "\t" + f"{validator.targets['seq chromosomes']}: ",
            "\t\t" + "chr3",
            "\t\t" + "chr4"
        ])

        validator.generate_chrom_report(shared_chroms, seq_chroms)
        self.assertEqual(validator.report.errors[0], exp_errors)

    """Does GFFValidator.generate_chrom_heuristics_report() correctly process an infractions report?"""

    def test_chrom_heuristics_report_output(self):
        validator = self.make_gff_validator()
        validator.chrom_set = {'chr1', 'chr2'}
        suspect_files = {
            'sam1': {'chr3', 'chr4'},
            'sam2': {'chr5', 'chr6'}
        }

        exp_warnings = '\n'.join([
            "GFF files and sequence files might not contain the same chromosome identifiers.",
            "This is determined from a subset of each sequence file, so false positives may be reported.",
            "\t" + f"{validator.targets['sam files']}: ",
            "\t\t" + "sam1: ",
            "\t\t\t" + f"Chromosomes sampled: {', '.join(sorted(suspect_files['sam1']))}",
            "\t\t" + "sam2: ",
            "\t\t\t" + f"Chromosomes sampled: {', '.join(sorted(suspect_files['sam2']))}",
            "\t" + f"{validator.targets['gff chromosomes']}: ",
            "\t\t" + "chr1",
            "\t\t" + "chr2"
        ])

        validator.generate_chrom_heuristics_report(suspect_files)
        self.assertEqual(validator.report.warnings[0], exp_warnings)

    """Does ReportFormatter add and print all sections with correct formatting?"""

    def test_report_multi_section(self):
        key_mapper = {'short': "long description"}
        formatter = ReportFormatter(key_mapper)

        formatter.add_warning_section("Header only")
        formatter.add_warning_section("Header 2", {'short': 5, 'short2': [1,2,3]})
        formatter.add_error_section("Header 3", {'short': {'level2a': {4,5,6}, 'level2b': 'msg'}})

        expected_report = '\n'.join([
            ReportFormatter.error_header,
            "Header 3",
            "\t" + "long description: ",
            "\t\t" + "level2a: ",
            *["\t\t\t" + str(i) for i in [4,5,6]],
            "\t\t" + "level2b: msg",
            "",
            ReportFormatter.warning_header,
            "Header only",
            "",
            ReportFormatter.warning_header,
            "Header 2",
            "\t" + "long description: 5",
            "\t" + "short2: ",
            *["\t\t" + str(i) for i in [1,2,3]],
            "",
            ""
        ])

        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            formatter.print_report()

        self.assertEqual(stdout.getvalue(), expected_report)

    """Do chromosome heuristics run in 2 seconds or less for a full-size test dataset?"""

    def test_chrom_heuristics_runtime(self):
        validator = self.make_gff_validator()
        validator.chrom_set = {'none match'}
        files = glob("./testdata/local_only/sam/full/*.sam")

        start = time.time()
        _ = validator.alignment_chroms_mismatch_heuristic(files)
        end = time.time()

        print(f"Chromosome heuristics runtime: {end-start:.2f}s")
        self.assertLessEqual(end-start, 2)



if __name__ == '__main__':
    unittest.main()
