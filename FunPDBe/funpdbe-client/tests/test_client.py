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
import os
from funpdbe_client.client import Client


class MockUser(object):

    def __init__(self):
        self.user_name = "foo"
        self.user_pwd = "bar"

    @staticmethod
    def set_user():
        return True

    @staticmethod
    def set_pwd():
        return True


class MockSchema(object):

    def __init__(self):
        self.json_schema = "foo"

    @staticmethod
    def get_schema():
        return True

    @staticmethod
    def clean_json(data):
        return data

    @staticmethod
    def validate_json(data):
        print(data)
        if data == {"foo": "bar"}:
            return True
        return False


def mocked_requests_get(*args, **kwargs):
    class MockResponse:
        def __init__(self, json_data, status_code):
            self.text = "foo"
            self.json_data = json_data
            self.status_code = status_code
    if args[0].endswith("resource/cath-funsites/1abc/"):
        return MockResponse({"resource": "ok"}, 200)
    elif args[0].endswith("resource/cath-funsites/"):
        return MockResponse({"resource": "ok"}, 200)
    elif args[0].endswith("resource/invalid/"):
        return MockResponse({"resource": "bad"}, 400)
    return MockResponse(None, 404)


def mocked_requests_post(*args, **kwargs):
    class MockResponse:
        def __init__(self, json_data, status_code):
            self.text = "foo"
            self.json_data = json_data
            self.status_code = status_code
    if args[0].endswith("cath-funsites/") or args[0].endswith("1abc/"):
        return MockResponse({"bar": "ok"}, 201)
    return MockResponse(None, 404)


def mocked_requests_delete(*args, **kwargs):
    class MockResponse:
        def __init__(self, json_data, status_code):
            self.text = "foo"
            self.json_data = json_data
            self.status_code = status_code

    if args[0].endswith("1abc/"):
        return MockResponse({"foo": "bar"}, 301)
    return MockResponse(None, 404)


class TestClient(TestCase):

    def setUp(self):
        self.mock_user = MockUser()
        self.mock_schema = MockSchema()
        self.client = Client(self.mock_schema, self.mock_user)

    def test_help_text(self):
        self.assertIsNotNone(self.client.__str__())

    def test_no_pdb_id(self):
        self.assertFalse(self.client.check_pdb_id(None))

    def test_no_resource_name(self):
        self.assertFalse(self.client.check_resource(None))

    def test_get_one_no_pdb_id(self):
        self.assertIsNone(self.client.get_one(pdb_id=None, resource="nod"))

    def test_pattern_mismatch(self):
        self.assertIsNone(self.client.get_one(pdb_id="invalid", resource="nod"))

    @mock.patch('requests.get', side_effect=mocked_requests_get)
    def test_get_all(self, mock):
        call = self.client.get_all("cath-funsites")
        self.assertEqual({"resource": "ok"}, call.json_data)

    def test_get_all_bad_resource(self):
        call = self.client.get_all("invalid")
        self.assertIsNone(call)

    @mock.patch('requests.get', side_effect=mocked_requests_get)
    def test_get_one_with_id_not_exist(self, mock):
        call = self.client.get_one("2abc", "nod")
        self.assertIsNone(call.json_data)

    @mock.patch('requests.get', side_effect=mocked_requests_get)
    def test_get_one_with_resource(self, mock):
        call = self.client.get_one("1abc", "cath-funsites")
        self.assertEqual({"resource": "ok"}, call.json_data)

    def test_parse_json_no_path(self):
        self.assertIsNone(self.client.parse_json(None))

    def test_parse_json_no_file(self):
        self.assertFalse(self.client.parse_json("random"))

    def test_parse_json(self):
        with open("tmp.json", "w") as tmp:
            tmp.write('{"foo": "bar"}')
        self.assertTrue(self.client.parse_json("tmp.json"))
        os.system("rm tmp.json")

    def test_parse_bad_json(self):
        with open("tmp.json", "w") as tmp:
            tmp.write("foo")
        self.assertFalse(self.client.parse_json("tmp.json"))
        os.system("rm tmp.json")

    def test_validate_json_no_schema(self):
        self.client.schema.json_schema = None
        self.client.json_data = {"foo": "bar"}
        self.assertTrue(self.client.validate_json())

    def test_validate_json_valid(self):
        self.client.json_data = {"foo": "bar"}
        self.assertTrue(self.client.validate_json())

    def test_validate_json_invalid(self):
        self.client.json_data = None
        self.assertFalse(self.client.validate_json())

    def test_post_no_file(self):
        self.assertIsNone(self.client.post(None, "cath-funsites", overwrite=False))

    def test_post_not_valid(self):
        with open("tmp.json", "w") as tmp:
            tmp.write('{"asd": "asd"}')
        self.client.json_data = None
        self.assertIsNone(self.client.post("tmp.json", "cath-funsites", overwrite=False))
        os.system("rm tmp.json")

    @mock.patch('requests.post', side_effect=mocked_requests_post)
    def test_post_without_overwrite(self, mock):
        response = self.client.post_without_overwrite("cath-funsites")
        self.assertEqual(response.status_code, 201)

    @mock.patch('requests.post', side_effect=mocked_requests_post)
    def test_post_with_overwrite(self, mock):
        self.client.json_data = {"pdb_id": "1abc"}
        response = self.client.post_with_overwrite("cath-funsites")
        self.assertEqual(response.status_code, 201)
        self.client.json_data["pdb_id"] = "0000"
        response = self.client.post_with_overwrite("cath-funsites")
        self.assertIsNone(response)

    def mock_validate_json_true(self):
        return True

    def test_pre_post_checks(self):
        self.client.validate_json = lambda: True
        self.assertTrue(self.client.pre_post_checks(None, "cath-funsites", True))

    def test_post(self):
        self.client.pre_post_checks = lambda x,y,z: True
        self.client.post_with_overwrite = lambda x: 1
        self.assertEqual(self.client.post(path=None, resource="cath-funsites", overwrite=True), 1)
        self.client.post_without_overwrite = lambda x: 2
        self.assertEqual(self.client.post(path=None, resource="cath-funsites", overwrite=False), 2)
        self.client.pre_post_checks = lambda x,y,z: False
        self.assertIsNone(self.client.post(path=None, resource="cath-funsites", overwrite=True))


    # @mock.patch('requests.post', side_effect=mocked_requests_post)
    # def test_post(self, mock):
    #     with open("tmp.json", "w") as tmp:
    #         tmp.write('{"foo": "bar"}')
    #     call = self.client.post("tmp.json", "cath-funsites", overwrite=False)
    #     os.system("rm tmp.json")
    #     self.assertEqual(call.status_code, 201)

    def test_post_bad_resource(self):
        self.assertIsNone(self.client.post("path", "invalid", overwrite=False))

    # @mock.patch('requests.post', side_effect=mocked_requests_post)
    # def test_put(self, mock):
    #     with open("tmp.json", "w") as tmp:
    #         tmp.write('{"foo": "bar"}')
    #     call = self.client.put("tmp.json", "1abc", "cath-funsites")
    #     os.system("rm tmp.json")
    #     self.assertEqual(call.status_code, 201)
    #
    # def test_put_no_file(self):
    #     self.assertIsNone(self.client.put(None, "1abc", "cath-funsites"))
    #
    # def test_put_not_valid(self):
    #     with open("tmp.json", "w") as tmp:
    #         tmp.write('{"asd": "asd"}')
    #     self.client.json_data = None
    #     self.assertIsNone(self.client.put("tmp.json", "1abc", "cath-funsites"))
    #     os.system("rm tmp.json")
    #
    # def test_put_bad_resource(self):
    #     self.assertIsNone(self.client.put("foo", "1abc", "bar"))
    #
    # def test_put_bad_pdb_id(self):
    #     self.assertIsNone(self.client.put("invalid", "something", "cath-funsites"))

    @mock.patch('requests.delete', side_effect=mocked_requests_delete)
    def test_delete(self, mock):
        self.assertEqual(self.client.delete_one("1abc", "nod").status_code, 301)

    @mock.patch('requests.delete', side_effect=mocked_requests_delete)
    def test_delete_not_there(self, mock):
        self.assertEqual(self.client.delete_one("2abc", "nod").status_code, 404)

    def test_delete_bad_pdb_id(self):
        self.assertIsNone(self.client.delete_one("invalid", "cath-funsites"))

    def test_delete_bad_resource(self):
        self.assertIsNone(self.client.delete_one("1abc", "invalid"))

    def test_check_pdb_id_valid(self):
        self.assertTrue(self.client.check_pdb_id("1abc"))

    def test_check_pdb_id_invalid(self):
        self.assertFalse(self.client.check_pdb_id("whatever"))

    def test_check_resource_valid(self):
        self.assertTrue(self.client.check_resource("cath-funsites"))

    def test_check_resource_invalid(self):
        self.assertFalse(self.client.check_resource("whatever"))

    def test_set_api_url(self):
        self.client.set_api_url("foo")
        self.assertEqual(self.client.api_url, "foo")

    def test_construct_get_url(self):
        self.assertIsNone(self.client.construct_get_url(None, None))