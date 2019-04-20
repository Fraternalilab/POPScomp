/*==============================================================================
arg.h : parse command line arguments
Copyright (C) 2007 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef ARG_H
#define ARG_H

#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "argpdb.h"
#include "error.h"

/*____________________________________________________________________________*/
/* structures */

/* variables for commmand line arguments */
typedef struct  
{
    FILE *pdbInFile;
	char *pdbInFileName;
	char *pdbmlInFileName;
	int pdbml;
	int zipped;
    FILE *trajInFile;
	char *trajInFileName;
    float rProbe;
	int silent;
	char *outDirName;
    FILE *sasaOutFile;
    char *sasaOutFileName;
    FILE *sasatrajOutFile;
    char *sasatrajOutFileName;
    FILE *bsasaOutFile;
    char *bsasaOutFileName;
    FILE *bsasatrajOutFile;
    char *bsasatrajOutFileName;
    FILE *sigmaOutFile;
    char *sigmaOutFileName;
    FILE *sigmatrajOutFile;
    char *sigmatrajOutFileName;
    FILE *interfaceOutFile;
    char *interfaceOutFileName;
	int interfaceOut;
	int compositionOut;
	int topologyOut;
	int typeOut;
	int atomOut;
	int residueOut;
	int chainOut;
    FILE *neighbourOutFile;
    char *neighbourOutFileName;
	int neighbourOut;
    FILE *parameterOutFile;
    char *parameterOutFileName;
	int parameterOut;
	int noTotalOut;
	int noHeaderOut;
	int padding;
	int rout;
	int jsonOut;
	FILE *jsonOutFile;
	char *jsonOutFileName;
	FILE *jsonbOutFile;
	char *jsonbOutFileName;
} Arg;

/*____________________________________________________________________________*/
/* prototypes */
int parse_args(int argc, char **argv, Arg *arg, Argpdb *argpdb);

#endif
