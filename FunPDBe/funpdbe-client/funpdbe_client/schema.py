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

import logging
import requests
import json
import jsonschema


class Schema(object):
    """
    Schema object to retrieve the JSON schema from GitHub
    and to validate user JSON against this schema
    """

    def __init__(self):
        self.url_base = "https://gitlab.ebi.ac.uk/pdbe-kb/funpdbe/funpdbe-schema/raw/master"
        self.json_url = "%s/funpdbe_schema.json" % self.url_base
        self.json_schema = None

    def get_schema(self):
        """
        Getting JSON schema
        :return: JSON, schema or None
        """
        logging.debug("Getting JSON schema")
        response = requests.get(self.json_url)
        if response.status_code == 404:
            logging.error(response.text)
            return None
        try:
            self.json_schema = json.loads(response.text)
        except ValueError as valerr:
            logging.warning(valerr)

    def validate_json(self, json_data):
        """
        Validating JSON data against schema
        :param json_data: JSON, user data
        :return: True if JSON is valid, False is invalid or other problems
        """
        logging.debug("Validating JSON")
        if not self.json_schema or not json_data:
            return False
        try:
            jsonschema.validate(json_data, self.json_schema)
            return True
        except jsonschema.exceptions.ValidationError as err:
            logging.warning("JSON does not comply with schema")
            logging.warning(err)
            return False

    @staticmethod
    def clean_json(json_data):
        """
        Convert source-database and pdb_id Strings to
        lower case uniformly
        :param json_data: JSON
        :return: JSON
        """
        json_copy = json_data
        for i in range(len(json_data["sites"])):
            if "source_database" in json_data["sites"][i].keys():
                json_copy["sites"][i]["source_database"] = json_data["sites"][i]["source_database"].lower()
        json_copy["pdb_id"] = json_data["pdb_id"].lower()
        return json_copy
