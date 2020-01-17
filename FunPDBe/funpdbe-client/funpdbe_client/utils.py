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

from funpdbe_client.constants import CLIENT_ERRORS, LOCAL_API_URL, DEV_API_URL, PROD_API_URL
from funpdbe_client.logger_config import generic_error


def check_exists(value, error, component_logger):
    """
    Check if a value is not none, and
    log error message if it is
    :param value: Any
    :param error: String
    :param component_logger: Logger instance
    :return: Boolean
    """
    if value:
        return True
    component_logger.log().error(CLIENT_ERRORS[error])
    return False


def check_status(response, expected, component_logger):
    """
    Check if status code is what is expected
    and log message accordingly
    :param response: Response
    :param expected: Int
    :param component_logger: Logger instance
    :return: None
    """
    if response.status_code == expected:
        component_logger.log().info("[%i] SUCCESS" % response.status_code)
    else:
        generic_error()
        component_logger.log().error("[%i] FAIL - %s" % (response.status_code, response.text))


def map_url(label):
    """
    Return the correct API URL based on user input
    parameters
    :param label: String, prod, dev or local
    :return: String, API URL
    """
    urls = {
        "prod": PROD_API_URL,
        "dev": DEV_API_URL,
        "local": LOCAL_API_URL
    }
    if label not in urls.keys():
        return None
    return urls[label]


def set_attribute(attribute, text):
    """
    Set attribute to user input
    :param attribute:
    :param text:
    :return: String, user input
    """
    value = attribute
    if not value:
        while not value:
            value = input(text)
    return value
