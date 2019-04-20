/*==============================================================================
fs.c : file system routines
Copyright (C) 2008 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "fs.h"

#ifdef MPI
#include <mpi.h>
#endif
extern int nodes;
extern int my_rank;

/*___________________________________________________________________________*/
/* list directory content */
/* modified from http://www.gnu.org/software/libtool/manual/libc/Simple-Directory-Lister.html */
int listfiles(char *dir, FileList *filelist)
{
	DIR *dp;
	struct dirent *ep;
	unsigned int n = 0;
	unsigned int allocated = 64;

	filelist->file = safe_malloc(allocated * sizeof(char [128]));

	dp = opendir (dir);
	if (dp != 0) {
		while ((ep = readdir(dp))) {
			/*puts(ep->d_name);*/
			if (ep->d_name[0] == '.')
				continue;
			strcpy(filelist->file[n], ep->d_name);
			filelist->nFile = ++ n;
			if (n == allocated) {
				allocated += 64;
				filelist->file = safe_realloc(filelist->file, allocated * sizeof(char [128]));
			}
		}
		return closedir (dp);
    } else {
		perror ("Failed opening directory");
		return 1;
	}
}

/*___________________________________________________________________________*/
int file_exists(char *fileName)
{
   struct stat buf;
   int i = stat(fileName, &buf);
     if (i == 0)
       return 1; /* file found */

     return 0; /* file not found */
}

