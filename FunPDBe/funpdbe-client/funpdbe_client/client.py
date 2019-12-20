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

import requests
import json
import re
from funpdbe_client.constants import PDB_ID_PATTERN, PROD_API_URL, RESOURCES, CLIENT_ERRORS
from funpdbe_client.logger_config import FunPDBeClientLogger, generic_error
from funpdbe_client.utils import check_exists, check_status


class Client(object):
    """
    The FunPDBe deposition client allows users to deposit, delete and view
    data in the FunPDBe deposition system
    """

    def __init__(self, schema, user):
        self.user = user
        self.schema = schema
        self.api_url = PROD_API_URL
        self.json_data = None
        self.logger = FunPDBeClientLogger("client")

    def __str__(self):
        return """
The FunPDBe deposition client allows users to deposit, delete and view  
data in the FunPDBe deposition system                                   

Usage parameters:
-h, --help:       Help (this is what you see now)
-u, --user:       FunPDBe user name
-p, --pwd:        FunPDBe password
-m, --mode:       Running mode (get, post, delete, validate)
-o, --overwrite:  Overwrite existing entries
-i, --pdbid:      PDB id of an entry
-r, --resource:   Name of a resource
-f, --path:       Path to JSON file (.json ending), or files (folder name)
-d, --debug:      Enable debug logging
-a, --api:        API to be used (prod (default), dev, local)
        """

    def set_api_url(self, new_url):
        """
        Setting API to other than default
        :param new_url: String, API URL
        :return: None
        """
        self.api_url = new_url

    def get_all(self, resource):
        """
        Get all FunPDBe entries based on resource name
        :param resource: String, resource name
        :return: Response
        """
        message = "GET all entries for %s" % resource
        if not self.check_resource(resource):
            return None
        self.logger.log().info(message)

        self.user_info()
        url = self.construct_get_url(resource)
        r = requests.get(url, auth=(self.user.user_name, self.user.user_pwd))

        if r.status_code == 200:
            self.logger.log().info("[%i] success" % r.status_code)
        else:
            self.log_api_error(r.status_code, r.text)

        print(r.text)
        return r

    def get_one(self, pdb_id, resource):
        """
        Get one FunPDBe entry based on PDB id and
        resource name
        :param pdb_id: String, PDB id
        :param resource: String, resource name
        :return: Response
        """
        message = "GET entry for %s" % pdb_id
        if resource and self.check_pdb_id(pdb_id):
            message += " from %s" % resource
            self.logger.log().info(message)
        else:
            return None

        self.user_info()
        url = self.construct_get_url(resource, pdb_id)
        r = requests.get(url, auth=(self.user.user_name, self.user.user_pwd))

        if r.status_code == 200:
            self.logger.log().info("[%i] success" % r.status_code)
        else:
            self.log_api_error(r.status_code, r.text)

        print(r.text)
        return r

    def post(self, path, resource, overwrite, plugin=False):
        """
        POST JSON to deposition API
        :param path: String, path to JSON file
        :param resource: String, resource name
        :param plugin: Boolean, plugin mode
        :return: None
        """
        message = "POST %s to %s (overwrite=%s)" % (path, resource, overwrite)
        self.logger.log().info(message)
        self.user_info()
        if not self.pre_post_checks(path, resource, plugin):
            return None
        if overwrite:
            return self.post_with_overwrite(resource)
        else:
            return self.post_without_overwrite(resource)

    def pre_post_checks(self, path, resource, plugin):
        if not self.check_resource(resource):
            return False
        if not self.parse_json(path) and not plugin:
            return False
        if not self.validate_json():
            return False
        return True

    def post_with_overwrite(self, resource):
        pdb_id = self.json_data["pdb_id"]
        if not self.check_pdb_id(pdb_id):
            return None
        url = self.api_url
        url += "resource/%s/%s/" % (resource, pdb_id)
        r = requests.post(url, json=self.json_data, auth=(self.user.user_name, self.user.user_pwd))
        check_status(r, 201, self.logger)
        return r

    def post_without_overwrite(self, resource):
        url = self.api_url
        url += "resource/%s/" % resource
        r = requests.post(url, json=self.json_data, auth=(self.user.user_name, self.user.user_pwd))
        check_status(r, 201, self.logger)
        return r

    def delete_one(self, pdb_id, resource):
        """
        DELETE entry based on PDB id
        :param pdb_id: String, PDB id
        :param resource: String, resource name
        :return: none
        """
        message = "DELETE entry %s from %s" % (pdb_id, resource)
        print(message)
        self.logger.log().info(message)

        if not self.check_resource(resource):
            return None
        if not self.check_pdb_id(pdb_id):
            return None
        self.user_info()
        url = self.api_url
        url += "resource/%s/%s/" % (resource, pdb_id)
        r = requests.delete(url, auth=(self.user.user_name, self.user.user_pwd))
        check_status(r, 301, self.logger)
        return r

    def construct_get_url(self, resource, pdb_id=None):
        """
        Create the GET URL based on resource
        and pdb_id
        :param resource: String
        :param pdb_id: String
        :return: String
        """
        url = self.api_url
        if resource and self.check_resource(resource):
            url += "resource/%s/" % resource
            if pdb_id:
                url += "%s/" % pdb_id
            return url
        return None

    def check_pdb_id(self, pdb_id):
        """
        Check if PDB id exists and if it matches
        the regular expression pattern of a valid
        PDB identifier
        :param pdb_id: String
        :return: Boolean
        """
        if not check_exists(pdb_id, "no_pdb", self.logger):
            return False
        if re.match(PDB_ID_PATTERN, pdb_id):
            return True
        generic_error()
        self.logger.log().error(CLIENT_ERRORS["bad_pdb"])
        return False

    def check_resource(self, resource):
        """
        Check if resource name exists and
        if it is a known (registered) resource
        :param resource: String
        :return: Boolean
        """
        if not check_exists(resource, "no_resource", self.logger):
            return False
        if resource in RESOURCES:
            return True
        self.logger.log().error(CLIENT_ERRORS["unknown_resource"])
        generic_error()
        return False

    def parse_json(self, path):
        """
        Parse user JSON file
        :param path: String, path to JSON
        :return: Boolean
        """
        if not check_exists(path, "no_path", self.logger):
            return None
        try:
            with open(path) as json_file:
                try:
                    self.json_data = json.load(json_file)
                    self.logger.log().info("JSON parsed")
                    return True
                except ValueError as valerr:
                    self.logger.log().error(valerr)
                    generic_error()
                    return False
        except IOError as ioerr:
            self.logger.log().error(ioerr)
            generic_error()
            return False

    def validate_json(self):
        """
        Validate JSON against schema
        :return: Boolean
        """
        if not self.schema.json_schema:
            self.schema.get_schema()
        if self.schema.validate_json(self.json_data):
            self.logger.log().info("JSON complies with FunPDBe schema")
            self.json_data = self.schema.clean_json(self.json_data)
            return True
        self.logger.log().error(CLIENT_ERRORS["bad_json"])
        generic_error()
        return False

    def log_api_error(self, status_code, text):
        """
        Log error based on status code and
        response text
        :param status_code: Int
        :param text: String
        :return: None
        """
        self.logger.log().error("[%s] - %s" % (status_code, text))

    def user_info(self):
        """
        Prompting user to provide info if missing
        :return: None
        """
        self.user.set_user()
        self.user.set_pwd()
