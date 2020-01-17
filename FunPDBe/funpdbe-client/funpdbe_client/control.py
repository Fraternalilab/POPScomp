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

import glob
from funpdbe_client.logger_config import FunPDBeClientLogger, generic_error
from funpdbe_client.constants import CONTROL_ERRORS
from funpdbe_client.utils import map_url


class Control(object):

    def __init__(self, opts, client):
        self.opts = opts
        self.client = client
        self.path = self.loop_options("-f", "--path")
        self.pdb_id = self.loop_options("-i", "--pdb_id")
        self.mode = self.loop_options("-m", "--mode")
        self.resource = self.loop_options("-r", "--resource")
        self.user_name = self.loop_options("-u", "--user")
        self.pwd = self.loop_options("-p", "--pwd")
        self.help = False
        self.overwrite = False
        self.api_url = self.loop_options("-a", "--api")
        self.logger = FunPDBeClientLogger("control")

    def set_help_and_overwrite_flags(self):
        for option, value in self.opts:
            if option in ["-h", "--help"]:
                self.help = True
            elif option in ["-d", "--debug"]:
                self.logger.logger.setLevel("DEBUG")
            elif option in ["-o", "--overwrite"]:
                self.overwrite = True

    def run(self):
        """
        Main entry point
        :return: Response.text or None
        """
        self.configure()
        self.set_help_and_overwrite_flags()
        if self.help:
            print(self.client)
        elif not self.mode:
            generic_error()
            self.logger.log().error(CONTROL_ERRORS["no_mode"])
        else:
            return self.action()

        return None

    def action(self):
        """
        API calls that can be performed by
        the client
        :return: Response.text or None
        """
        actions = {
            "get": self.get,
            "post": self.post,
            "delete": self.delete,
            "validate": self.validate
        }
        if self.mode in actions.keys():
            return actions[self.mode]()
        return None

    def configure(self):
        """
        Set the user name and password
        :return: None
        """
        self.client.user.user_name = self.user_name
        self.client.user.user_pwd = self.pwd
        if self.api_url:
            url_to_use = map_url(self.api_url)
            if url_to_use:
                self.client.set_api_url(url_to_use)

    def validate(self):
        """
        Validate one or more JSON files against
        FunPDBe Schema
        :return: True or False
        """
        if not self.check_path():
            return False

        if self.path.endswith(".json"):
            return self.single_validate(self.path)
        else:
            return self.multi_validate()

    def multi_validate(self):
        """
        Validate multiple JSON files against
        FunPDBe Schema
        :return: True or False
        """
        all_valid = True
        for json_path in glob.glob("%s/*.json" % self.path):
            if not self.single_validate(json_path):
                all_valid = False
        return all_valid

    def single_validate(self, path):
        """
        Validate a single JSON
        :param path: String
        :return: Boolean
        """
        print("Parsing and validating %s" % path)
        if self.client.parse_json(path):
            return self.client.validate_json()
        return False

    def get(self):
        """
        Make GET call
        :return: Response.text or None
        """
        if self.resource:
            if self.pdb_id:
                return self.client.get_one(self.pdb_id, self.resource)
            else:
                return self.client.get_all(self.resource)
        self.logger.log().error(CONTROL_ERRORS["no_pdb_or_resource"])
        print("Please provide --resource and optionally --pdb_id")
        return None

    def post(self):
        """
        Make POST call
        :return: Response.text or None
        """
        response = None
        if not self.check_path():
            return response
        if self.path.endswith(".json"):
            response = self.client.post(self.path, self.resource, self.overwrite)
        else:
            response = self.batch_post()
        return response

    def batch_post(self):
        """
        Make batch POST call
        :return: Response
        """
        response = None
        attempted = 0
        succeeded = 0
        for json_path in glob.glob("%s/*.json" % self.path):
            self.log_delimiter("start")
            response = self.client.post(json_path, self.resource, self.overwrite)
            self.client.json_data = None
            if response and response.status_code == 201:
                succeeded += 1
            attempted += 1
            self.log_delimiter("end")
        message = "Batch POSTing: %i out of %i POSTed successfully" % (succeeded, attempted)
        self.logger.log().info(message)
        print(message)
        return response

    def check_entry_exists(self):
        response = self.client.get_one(self.pdb_id, self.resource)
        if response.status_code == 200:
            return True
        return False

    def delete(self):
        """
        Make DELETE call
        :return: Response.text or None
        """
        self.log_delimiter("start")
        response = self.client.delete_one(self.pdb_id, self.resource)
        self.log_delimiter("end")
        return response

    def loop_options(self, opt1, opt2):
        """
        Loop through options
        :param opt1: String
        :param opt2: String
        :return: String or None
        """
        for option, value in self.opts:
            if option == opt1 or option == opt2:
                return value
        return None

    def log_delimiter(self, mode):
        if mode == "start":
            self.logger.log().info("BEGIN")
        elif mode == "end":
            self.logger.log().info("END\n")

    def check_path(self):
        """
        Check if path exists
        :return: Boolean
        """
        if not self.path:
            self.logger.log().error(CONTROL_ERRORS["no_path"])
            return False
        return True
