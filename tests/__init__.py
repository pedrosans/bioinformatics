import unittest
from unittest import TestSuite
from tests.amber_ff_test import AmberTestCase
from tests.fasta_test import FastaTestCase

test_cases = (AmberTestCase, FastaTestCase,)


def load_tests(loader, tests, pattern):
    suite = TestSuite()
    for test_class in test_cases:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    return suite
