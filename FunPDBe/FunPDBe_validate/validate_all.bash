#! /bin/bash
ls -1 ~/database/JSON/c9/*.json | xargs -i bash -c './validateJSONfunpdbe.sh {} 2>error.log && grep complies'
