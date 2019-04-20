/*==============================================================================
json.c : JSON routines, using the cJSON library
Copyright (C) 2018 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef JSON_H
#define JSON_H

#include <stdio.h>

#include "arg.h"
#include "argpdb.h"
#include "cJSON.h"
#include "error.h"
#include "getpdb.h"
#include "sasa.h"
#include "sasa_const.h"

void make_resSasaJson(Arg *arg, Str *pdb, ResSasa *resSasa, cJSON *json);
void print_json(Arg *arg, cJSON *json);
void make_resbSasaJson(Arg *arg, Str *pdb, ResSasa *resSasa, cJSON *jsonb);
void print_jsonb(Arg *arg, cJSON *jsonb);

#endif


