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
from funpdbe_client.user import User


class TestUser(TestCase):

    def setUp(self):
        self.user = User('foo', 'bar')

    def test_set_user(self):
        self.user.set_user()
        self.assertEqual('foo', self.user.user_name)

    def test_set_pwd(self):
        self.user.set_pwd()
        self.assertEqual('bar', self.user.user_pwd)

