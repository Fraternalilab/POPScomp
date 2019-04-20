/*==============================================================================
error.h : error msg routines
Copyright (C) 2008 Jens Kleinjung and Alessandro pandini
Read the COPYING file for license information.
==============================================================================*/

#ifndef ERRORMSG_H
#define ERRORMSG_H

#include <stdio.h>
#include <stdlib.h>

/*____________________________________________________________________________*/
/* prototypes */
void Warning(char *message);
void WarningSpec(char *message, char *spec);
void Error(char *message);
void ErrorSpec(char *message, char *spec); 
void ErrorSpecNoexit(char *message, char *spec); 
void ErrorLoc(char *message, char *file, int line);

#endif
