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
from unittest import mock
from funpdbe_client.schema import Schema


def mocked_requests_get(*args, **kwargs):
    class MockResponse:
        def __init__(self, text, status_code):
            self.text = text
            self.status_code = status_code
    if args[0].endswith("funpdbe_schema.json"):
        return MockResponse('{"foo":"bar"}', 200)
    elif args[0].endswith("bad.json"):
        return MockResponse("asd", 200)
    return MockResponse(None, 404)


class TestSchema(TestCase):

    def setUp(self):
        self.schema = Schema()

    @mock.patch('requests.get', side_effect=mocked_requests_get)
    def test_get_schema(self, mock):
        self.schema.get_schema()
        self.assertIsNotNone(self.schema.json_schema)

    @mock.patch('requests.get', side_effect=mocked_requests_get)
    def test_get_schema_bad_json(self, mock):
        self.schema.json_url = "bad.json"
        self.schema.get_schema()
        self.assertIsNone(self.schema.json_schema)

    @mock.patch('requests.get', side_effect=mocked_requests_get)
    def test_get_schema_bad_data(self, mock):
        self.schema.json_url = "https://www.ebi.ac.uk"
        self.schema.get_schema()
        self.assertIsNone(self.schema.json_schema)

    def test_validate_json_with_missing_schema(self):
        mock_schema = {
            "type": "object",
            "properties": {
                "foo": {
                    "type": "string",
                    "description": "bar"
                }
            }
        }
        self.schema.json_schema = mock_schema
        self.assertFalse(self.schema.validate_json(42))

    def test_validate_json_with_missing_data(self):
        self.schema.json_schema = {"foo": "bar"}
        self.assertFalse(self.schema.validate_json(None))

    def test_validate_json(self):
        self.schema.json_schema = {"foo": "bar"}
        self.assertTrue(self.schema.validate_json({"foo": "bar"}))

    def test_clean_json(self):
        mock_data = {
            "sites": [
                {"source_database": "FOO"},
                {"source_database": "BAR"}
            ],
            "pdb_id": "1ABC"
        }
        expected = {
            "sites": [
                {"source_database": "foo"},
                {"source_database": "bar"}
            ],
            "pdb_id": "1abc"
        }
        self.assertEqual(expected, self.schema.clean_json(mock_data))
