/*=============================================================================
seq.h : sequence
Copyright (C) 2008 Jens Kleinjung
Read the COPYING file for license information.
===============================================================================*/
  
#ifndef SEQ_H
#define SEQ_H

/*____________________________________________________________________________*/
/* sequence */
typedef struct  
{
    char *name; /* sequence name */
    char *res; /* array of residues = sequence */
    int length; /* length of sequence */
} Seq;

#endif
