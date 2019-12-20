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
from funpdbe_client.constants import LOG_FILENAME


class FunPDBeClientLogger(object):
    """
    The FunPDBe client uses logging to log
    information and error to an output
    file. The file path is defined in
    constants.LOG_FILENAME
    """

    def __init__(self, name="general", write_mode="a"):
        self.write_mode = write_mode
        self.logger = logging.getLogger(name)
        self.configure()

    def configure(self):
        """
        Configure a handler for the logger
        :return: None
        """
        config = logging.FileHandler(LOG_FILENAME, mode=self.write_mode)
        config.setLevel(logging.INFO)
        formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        config.setFormatter(formatter)
        self.logger.addHandler(config)
        self.logger.setLevel("INFO")

    def log(self):
        """
        Get the logger
        :return: Logger
        """
        return self.logger


def generic_error():
    """
    Print a generic error message that points
    to the log file
    :return: None
    """
    print("FAILED - details log saved in %s" % LOG_FILENAME)
