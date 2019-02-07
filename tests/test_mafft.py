from unittest import TestCase

from minimock import mock, restore

import alignment

class TestMafft (TestCase): 

	def tearDown(self):
		restore()

	def test_simple(self):
		expected = {"A": "THISISASENTENCE", "B": "---THATSENTENCE"}
		assert alignment.run_mafft({"A": "THISISASENTENCE", "B": "THATSENTENCE"}) == expected

	def test_empty(self):
		with self.assertRaisesRegex(ValueError, "Requires at least two sequences"):
			alignment.run_mafft({})

	def test_single(self):
		with self.assertRaisesRegex(ValueError, "Requires at least two sequences"):
			alignment.run_mafft({"A": "THISISASENTENCE"})

	def test_bad_sequence(self):
		mock("alignment.is_invalid_sequence", returns=True)
		with self.assertRaisesRegex(ValueError, "Invalid sequence"):
			alignment.run_mafft({"A": "1234", "B": "THISSENTENCE"})


