/**
 * errors.h
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

/* See the file errors.c for a description. */

#ifndef ERRORS_H
#define ERRORS_H

/* Size of the buffer containing any error string. */
#define ERROR_BUFFER_SIZE 1024

/* Storage of error quantities. */
extern int pytsai_error;        /* true / false error flag  */
extern char pytsai_string[];    /* error description string */

/* Error methods. */
void pytsai_clear();
void pytsai_raise(char *message);
int pytsai_haserror();

#endif /* ERRORS_H */

