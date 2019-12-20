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

# These are the base URLs of the API of the FunPDBe deposition system
PROD_API_URL = "https://www.ebi.ac.uk/pdbe/funpdbe/deposition/entries/"
DEV_API_URL = "https://wwwdev.ebi.ac.uk/pdbe/funpdbe/deposition/entries/"
LOCAL_API_URL = "http://127.0.0.1:8000/entries/"

CLIENT_ERRORS = {
    "no_pdb": "No PDB identifier specified",
    "bad_pdb": "Invalid PDB identifier pattern",
    "no_resource": "No resource name specified",
    "unknown_resource": "Unknown resource name",
    "no_path": "No file path to JSON(s) specified",
    "bad_json": "JSON does not comply with FunPDBe schema"
}

CONTROL_ERRORS = {
    "no_mode": "No running mode was specified (--mode=)",
    "no_path": "No path to JSON file(s) provided",
    "no_pdb_or_resource": "No PDB identifier or resource name"
}

# This is the file name of the log file that the client generates
LOG_FILENAME = 'client.log'

# This is the regular expression pattern for a PDB identifier
# as defined in the mmCIF format
PDB_ID_PATTERN = "^[1-9][a-zA-Z0-9]{3}$"

# These are the valid resource name in FunPDBe
RESOURCES = (
    "cath-funsites",
    "3dligandsite",
    "nod",
    "popscomp",
    "14-3-3-pred",
    "dynamine",
    "cansar",
    "credo",
    "depth",
    "akid"
)
