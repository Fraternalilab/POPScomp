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
from funpdbe_client.plugin import Plugin

mock_data = {
    "data_resource": "foo",
    "pdb_id": "bar"
}


class MockClient(object):

    def __init__(self):
        pass

    def post(self, path, resource, plugin, overwrite):
        return "success"

    def delete_one(self, pdb_id, resource):
        pass


class TestPlugin(TestCase):

    def setUp(self):
        self.plugin = Plugin("user", "pwd", mock_data)

    def test_set_value(self):
        self.assertEqual(self.plugin.set_value("data_resource"), "foo")
        self.assertEqual(self.plugin.set_value("pdb_id"), "bar")
        self.assertIsNone(self.plugin.set_value("invalid_key"))

    def test_init_user(self):
        self.assertIsNotNone(self.plugin.initialize_client())
        self.plugin = Plugin(None, None, mock_data)
        self.assertIsNotNone(self.plugin.initialize_client())

    def test_init_client(self):
        self.assertIsNotNone(self.plugin.initialize_client())

    def test_check_if_variables_are_set(self):
        self.assertTrue(self.plugin.check_if_variables_are_set())
        self.plugin.pdb_id = None
        self.assertFalse(self.plugin.check_if_variables_are_set())
        self.plugin.pdb_id = "foo"
        self.plugin.resource = None
        self.assertFalse(self.plugin.check_if_variables_are_set())
        self.plugin.resource = "foo"
        self.plugin.user = None
        self.assertFalse(self.plugin.check_if_variables_are_set())
        self.plugin.user = "foo"
        self.plugin.pwd = None
        self.assertFalse(self.plugin.check_if_variables_are_set())

    def test_post(self):
        self.plugin.client = MockClient()
        self.assertIsNotNone(self.plugin.post())
        self.assertIsNotNone(self.plugin.post(update=True))
        self.plugin.pdb_id = None
        self.assertFalse(self.plugin.post())
