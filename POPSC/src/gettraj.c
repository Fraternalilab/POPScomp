/*==============================================================================
getTRAJ.c : Routines for reading GROMOS96 trajectory file
Copyright (C) 2008-2011 Jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "gettraj.h"

/*____________________________________________________________________________*/
/** compile pattern trajectory */
static void compile_pattern_trajectory(regex_t *regex)
{
    char matchPattern[64];

	/* pattern for GROMOS trajectory coordinate line */
    strcpy (matchPattern, ".*[[:print:]]{10,15}.*[[:print:]]{10,15}.*[[:print:]]{10,15}");

    assert(regcomp(regex, matchPattern, REG_EXTENDED) == 0);
}

/*____________________________________________________________________________*/
/** compile pattern POSITIONRED */
static void compile_pattern_positionred(regex_t *regex)
{
    char matchPattern[16];

	/* pattern for POSITIONRED */
    strcpy (matchPattern, "POSITIONRED");

    assert(regcomp(regex, matchPattern, REG_EXTENDED) == 0);
}

/*____________________________________________________________________________*/
/** read GROMOS96 trajectory file */
/* Definition of GROMOS trajectory format:
 * see GROMOS96 manual, page III-42, chapter 3.4.2 Atomic coordinates, 
 * Reduced trajectory information in formatted form
 * FORMAT (3F15.9)
 * with F: float */
int read_gromos_traj(Traj *traj, Arg *arg, int protEnd)
{
	unsigned int allocated_inc = 64;
	unsigned int allocated_frame = allocated_inc;
	unsigned int allocated_atom = protEnd;
	char line[80];
    regex_t trajectory; /* regular expression of gromos coordinate line */
    regex_t positionred; /* regular expression of POSITIONRED */
	int *pnAtom; /* pointer to atom number */
	unsigned int nAtom_mem = 0;

    if (! arg->silent) fprintf(stdout, "\tGRO96 file: %s\n", arg->trajInFileName);
    arg->trajInFile = safe_open(arg->trajInFileName, "r");

	/** initialise/allocate memory for trajectory frames and atoms */
    traj->frame = safe_malloc(allocated_frame * sizeof(Frame));
    traj->nFrame = 0;
	traj->frame[traj->nFrame].trajatom = \
		safe_malloc(allocated_atom * sizeof(Trajatom));
	traj->frame[traj->nFrame].nAtom = 0;
	pnAtom = &(traj->frame[traj->nFrame].nAtom);

	/* compile trajectory coordinates pattern */
	compile_pattern_trajectory(&trajectory);
	compile_pattern_positionred(&positionred);

	/* read coordinate file */
	while(fgets(line, 80, arg->trajInFile) != 0) { /* read line */
		/* if this line (search pattern) matches the coordinate line format */
		/*fprintf(stderr, "=> %d %d\n", traj->nFrame, *pnAtom);*/
		if (match_pattern(&trajectory, line) == 0) {
			/* scan this line in and check whether the matching works */
			/*____________________________________________________________________________*/
			if (sscanf(&(line[0]), "%f%f%f\n", 
				&(traj->frame[traj->nFrame].trajatom[*pnAtom].pos.x), 
				&(traj->frame[traj->nFrame].trajatom[*pnAtom].pos.y), 
				&(traj->frame[traj->nFrame].trajatom[*pnAtom].pos.z)) 
				== 3) {

				/* convert from nm to A */
				traj->frame[traj->nFrame].trajatom[*pnAtom].pos.x *= 10;
				traj->frame[traj->nFrame].trajatom[*pnAtom].pos.y *= 10;
				traj->frame[traj->nFrame].trajatom[*pnAtom].pos.z *= 10;

				/* count atoms in this frame */
				++ *pnAtom;

				/* allocate more memory for atoms if needed */
				if (traj->frame[traj->nFrame].nAtom == allocated_atom) {
					allocated_atom += allocated_inc;
					traj->frame[traj->nFrame].trajatom = \
						safe_realloc(traj->frame[traj->nFrame].trajatom, allocated_atom * sizeof(Atom));
				}

				/*____________________________________________________________________________*/
				/* if end of protein is reached */
				if (*pnAtom >= protEnd) {
					/*while((fgets(line, 80, trajFile) != 0) && (match_pattern(&trajectory, line) != 0))
						fprintf(stderr, "==>skipped %s", line);*/
					while((fgets(line, 80, arg->trajInFile) != 0) && (match_pattern(&positionred, line) != 0))
						;

					/* assert that all frames have the same atom number */
					if (traj->nFrame == 0)
						nAtom_mem = *pnAtom;
					else
						assert(*pnAtom == nAtom_mem);

					/* count frames in this trajectory */
					++ traj->nFrame;
					
					/* allocate more memory for frames if needed */
					if (traj->nFrame == allocated_frame) {
						allocated_frame += allocated_inc;
						traj->frame = \
								safe_realloc(traj->frame, allocated_frame * sizeof(Frame));
					}

					/* initialise new frame */
					allocated_atom = protEnd;
					traj->frame[traj->nFrame].trajatom = \
						safe_malloc(allocated_atom * sizeof(Trajatom));
					pnAtom = &(traj->frame[traj->nFrame].nAtom);
					*pnAtom = 0;
				}
			}
		}
	}

	regfree(&trajectory); /* free expression structure */
	regfree(&positionred);

	fclose(arg->trajInFile);

	if (! arg->silent)
		fprintf(stdout, "\tGromos trajectory file content (water and ions excluded):\n"
						"\tnAtom = %d (per frame, taken from reference molecule file)\n\tnFrame = %d\n",
			nAtom_mem, traj->nFrame);

	return 0;
}

/*____________________________________________________________________________*/
/** copy coordinates from trajectory to template molecule */
void copy_coordinates(Str *pdb, Traj *traj, int frame)
{
	unsigned int i;

    for (i = 0; i < pdb->nAtom; ++ i)
		v_copy(&(pdb->atom[i].pos), &(traj->frame[frame].trajatom[pdb->atomMap[i]].pos));
}

