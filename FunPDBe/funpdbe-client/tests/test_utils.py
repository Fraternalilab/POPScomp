#!/usr/bin/env python3

# Copyright 2018 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
# either express or implied. See the License for the specific
# language governing permissions and limitations under the
# License.

from unittest import TestCase
from unittest.mock import patch
from funpdbe_client.utils import set_attribute, map_url
from funpdbe_client.constants import LOCAL_API_URL, DEV_API_URL, PROD_API_URL


class TestUser(TestCase):

    def test_set_attribute(self):
        with patch('builtins.input', side_effect='f'):
            new_value = set_attribute(None, 'user prompt message')
        self.assertEqual('f', new_value)

    def test_map_url(self):
        self.assertEqual(map_url("prod"), PROD_API_URL)
        self.assertEqual(map_url("dev"), DEV_API_URL)
        self.assertEqual(map_url("local"), LOCAL_API_URL)
        self.assertIsNone(map_url(None))
