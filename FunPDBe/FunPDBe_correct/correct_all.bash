#! /bin/bash

ls -1 ~/database/JSON/??/*.json | xargs -i bash -c './correctJSONfunpdbe.sh {}'

