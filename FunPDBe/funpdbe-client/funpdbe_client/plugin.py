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

from funpdbe_client.client import Client
from funpdbe_client.user import User
from funpdbe_client.schema import Schema


class Plugin(object):
    """
    This Plugin module can be used to import the POST functionality
    of the FunPDBe client into any Python pipeline

    Please note that you need a valid (authenticated) user name and
    password in order to use the plugin.

    You also have to use valid PDB id and resource name

    This POST functionality only pushes one JSON (data) per call
    """

    def __init__(self, user, pwd, data):
        self.user = user
        self.pwd = pwd
        self.data = data
        self.resource = self.set_value("data_resource")
        self.pdb_id = self.set_value("pdb_id")
        self.client = self.initialize_client()

    def set_value(self, key):
        """
        Helper function to check if a key
        exists, and then return the value
        from the data
        :param key: String, key in JSON data
        :return: String, PDB id or data resource
        """
        if self.check_if_key_exists(key):
            return self.data[key]
        return None

    def check_if_key_exists(self, key):
        """
        Helper function to perform simple check
        :param key: String, key in JSON data
        :return: True if key exists, False if not
        """
        if key in self.data.keys():
            return True
        return False

    def post(self, update=False):
        """
        POST data to FunPDBe deposition API
        :param update: True for overwriting existing entry, False for not
        :return: True if successful, False if not
        """
        if not self.check_if_variables_are_set():
            print("POSTing failed due to missing information")
            return False
        self.client.json_data = self.data
        return self.client.post(path=None, resource=self.resource, plugin=True, overwrite=update)

    def check_if_variables_are_set(self):
        """
        Checks if mandatory keys have values
        :return: True if all required values exist, False if not
        """
        for variable in [self.user, self.pwd, self.data, self.resource, self.pdb_id]:
            if not variable:
                print("Fatal: please check if user name, password and data are all provided!")
                return False
        return True

    def initialize_client(self):
        """
        Create Schema() and User() instances and
        create a Client() instance using them
        :return: new Client()
        """
        schema = Schema()
        user = self.initialize_user()
        return Client(schema, user)

    def initialize_user(self):
        """
        Create User() instance based on
        user name and password
        :return: new User()
        """
        return User(user=self.user, pwd=self.pwd)
