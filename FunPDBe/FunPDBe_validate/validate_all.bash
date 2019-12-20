#! /bin/bash
ls -1 ~/database/JSON/??/*.json | xargs -i bash -c './validateJSONfunpdbe.sh {} 2>error.log && grep complies'
