/*==============================================================================
modstring.c : modify strings
Copyright (C) 2008 jens Kleinjung
Read the COPYING file for license information.
==============================================================================*/

#include "modstring.h"

/*___________________________________________________________________________*/
/** strip string s of white space and return resulting string to t */
void strip_char(const char *s, char *t) 
{
	unsigned int i = 0;
	unsigned int j = 0;

	while ((t[j] = s[i++]) != '\0')
		if (isspace((unsigned char) t[j]) == 0)
			j++;
}

/*___________________________________________________________________________*/
/** position (min = 0) of substring in string */
int strpos0(char *string, char *substring)
{
    char *a, *b;
    unsigned int pos = 0;

    /* first find single character match, then match rest of substring */
    b = substring;
    if (*b == 0)
        return 0;

    for ( ; *string != 0; string += 1) {
        ++ pos;
        if (*string != *b)
            continue;

		a = string;
        while (1) {
            if (*b == 0)
                return (pos - 1);

            if (*a ++ != *b ++)
                break;
        }
        b = substring;
    }
    return (pos - 1);
}

/*___________________________________________________________________________*/
/** position (min = 1) of substring in string */
int strpos1(char *string, char *substring)
{
    char *a, *b;
    unsigned int pos = 0;

    /* first find single character match, then match rest of substring */
    b = substring;
    if (*b == 0)
        return 0;

    for ( ; *string != 0; string += 1) {
        ++ pos;
        if (*string != *b)
            continue;

		a = string;
        while (1) {
            if (*b == 0)
                return pos;

            if (*a ++ != *b ++)
                break;
        }
        b = substring;
    }
    return pos;
}

