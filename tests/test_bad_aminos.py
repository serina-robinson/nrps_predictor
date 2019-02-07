from unittest import TestCase

import alignment

class TestAminos(TestCase):
	def test_empty(self):
		with self.assertRaisesRegex(ValueError, "No sequence provided"):
			alignment.is_invalid_sequence([])

	def test_bad(self):
		assert alignment.is_invalid_sequence(("AB1234")) == True

	def test_good(self):
		assert alignment.is_invalid_sequence(("THISSENTENCE")) == False