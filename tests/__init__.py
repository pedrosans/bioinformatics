import unittest
from unittest import TestSuite
from tests.amber_ff_test import AmberTestCase
from tests.fasta_test import FastaTestCase
from tests.minimization_test import MinimizationTestCase
from tests.pdb_test import PdbTestCase
from tests.topology_test import TopologyTestCase

test_cases = (AmberTestCase, FastaTestCase,
				MinimizationTestCase, PdbTestCase,
				TopologyTestCase,)


def load_tests(loader, tests, pattern):
	suite = TestSuite()
	for test_class in test_cases:
		tests = loader.loadTestsFromTestCase(test_class)
		suite.addTests(tests)
	return suite
