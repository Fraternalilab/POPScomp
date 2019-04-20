/*==============================================================================
fs.h : file system routines
Copyright (C) 2008 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#ifndef FS_H
#define FS_H

#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "safe.h"

/*___________________________________________________________________________*/
typedef struct
{
   char (*file)[128];
   int nFile;
} FileList;

/*___________________________________________________________________________*/
int listfiles(char *dir, FileList *filelist);
int file_exists(char *fileName);

#endif


