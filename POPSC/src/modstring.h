/*==============================================================================
modstring.h : modify strings
Copyright (C) 2008 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#if !defined MODSTRING_H
#define MODSTRING_H

#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*___________________________________________________________________________*/
void strip_char(const char *s, char *t); 
int strpos0(char *string, char *substring);
int strpos1(char *string, char *substring);

#endif
