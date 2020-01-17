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

import os
from unittest import TestCase
from funpdbe_client.control import Control


class MockUser(object):

    def __init__(self):
        self.user_name = None
        self.user_pwd = None


class MockObject(object):

    def __init__(self, mock_user):
        self.user = mock_user
        self.schema = None
        self.pwd = None

    @staticmethod
    def get_one(arg1, arg2):
        return MockResponse2()

    @staticmethod
    def get_all(arg1):
        return True

    @staticmethod
    def post(arg1, arg2, arg3):
        return MockResponse()

    @staticmethod
    def set_api_url(arg1):
        return True

    @staticmethod
    def delete_one(arg1, arg2):
        return True

    @staticmethod
    def parse_json(path):
        if path:
            return True
        return False

    @staticmethod
    def validate_json():
        return True


class MockResponse(object):

    def __init__(self):
        self.status_code = 201


class MockResponse2(object):

    def __init__(self):
        self.status_code = 200


class TestControl(TestCase):

    def setUp(self):
        mock_opts = [("--user", "test"), ("--pwd", "test")]
        self.control = Control(mock_opts, MockObject(MockUser()))

    def mock_function(self):
        return True

    def test_run_no_mode(self):
        self.control.debug = True
        self.assertIsNone(self.control.run())
        self.control.debug = False
        self.assertIsNone(self.control.run())

    def test_run_help(self):
        mock_opts = [("--help", "help")]
        self.control = Control(mock_opts, MockObject(MockUser()))
        self.assertIsNone(self.control.run())

    def test_set_help_and_overwrite_flags(self):
        mock_opts = [("--help", "help")]
        self.control = Control(mock_opts, MockObject(MockUser()))
        self.control.set_help_and_overwrite_flags()
        self.assertTrue(self.control.help)
        mock_opts = [("--overwrite", "overwrite")]
        self.control = Control(mock_opts, MockObject(MockUser()))
        self.control.set_help_and_overwrite_flags()
        self.assertTrue(self.control.overwrite)

    def test_check_entry_exists(self):
        self.control = Control([("--mode", "get")], MockObject(MockUser()))
        self.assertTrue(self.control.check_entry_exists())
        self.control.client.get_one = lambda x,y: MockResponse()
        self.assertFalse(self.control.check_entry_exists())

    def test_run(self):
        self.control = Control([("--mode", "get")], MockObject(MockUser()))
        self.control.api_url = "local"
        self.control.get = self.mock_function
        self.assertTrue(self.control.run())
        self.control = Control([("--mode", "post")], MockObject(MockUser()))
        self.control.post = self.mock_function
        self.assertTrue(self.control.run())
        self.control = Control([("--mode", "delete")], MockObject(MockUser()))
        self.control.delete = self.mock_function
        self.assertTrue(self.control.run())
        self.control = Control([("--mode", "foo")], MockObject(MockUser()))
        self.assertIsNone(self.control.run())

    def test_get(self):
        self.control.client = MockObject(MockUser())
        self.assertIsNone(self.control.get())
        self.control.pdb_id = "foo"
        self.assertIsNone(self.control.get())
        self.control.resource = "bar"
        self.assertIsNotNone(self.control.get())
        self.control.pdb_id = None
        self.assertIsNotNone(self.control.get())

    def test_post(self):
        self.control.client = MockObject(MockUser())
        self.assertIsNone(self.control.post())
        self.control.path = ".json"
        self.assertIsNotNone(self.control.post())
        os.system("touch foo.json")
        self.control.path = "./"
        self.assertIsNotNone(self.control.post())
        os.system("rm foo.json")

    def test_batch_post(self):
        self.control.client = MockObject(MockUser())
        self.assertIsNone(self.control.batch_post())

    def test_delete(self):
        self.control.client = MockObject(MockUser())
        self.assertIsNotNone(self.control.delete())

    def test_single_validate(self):
        self.assertFalse(self.control.single_validate(None))
        self.assertTrue(self.control.single_validate("path/to/json"))

    def test_validate(self):
        self.control.path = None
        self.assertFalse(self.control.validate())
        self.control.path = "path/to/file.json"
        self.assertTrue(self.control.validate())
        self.control.path = "tests/"
        self.assertTrue(self.control.validate())
        self.control.single_validate = lambda x: False
        self.assertFalse(self.control.validate())
