/**
 * errors.c
 * Copyright (C) Jonathan Merritt 2004.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

 /***************************************************************************\
 * Methods for handling errors.                                              *
 *                                                                           *
 * The approach taken for setting up the Tsai camera matching code as a      *
 * Python extension module is as follows.  The library originally used a     *
 * combination of sprintf() and exit() to indicate errors and exit:          *
 *      sprintf(stderr, "some error\n");                                     *
 *      exit(-1);                                                            *
 * This has been replaced by setting of an error flag, and the functions     *
 * that can fail now return 0 on failure and 1 on success rather than being  *
 * of void return type as before:                                            *
 *      pytsai_raise("some error");                                          *
 *      return 0;                                                            *
 *                                                                           *
 * To trap errors for the Python wrapper, the following is the template:     *
 *      pytsai_clear();                                                      *
 *      --- perform C-library calls here ---                                 *
 *      if (pytsai_haserror())                                               *
 *              --- raise Python exception ---                               *
 \***************************************************************************/

#include <string.h>
#include "errors.h"

/* true / false error flag */
int pytsai_error;
/* character array to store an error string */
char pytsai_string[ERROR_BUFFER_SIZE];

/**
 * Clears the error flag.
 */
void pytsai_clear()
{
        pytsai_error = 0;
}

/**
 * Sets the error flag and stores an error message.
 */
void pytsai_raise(char *message)
{
        pytsai_error = 1;
        strncpy(pytsai_string, message, ERROR_BUFFER_SIZE-1);
}

/**
 * Returns non-zero if the error flag has been set.
 */
int pytsai_haserror()
{
        return pytsai_error;
}

