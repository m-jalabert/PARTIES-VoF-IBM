#ifndef INPUT_H
#define INPUT_H

#include "DataTypes.h"

void Input_set_parameters(Parameters *params, char *infile);
int Input_handler(void *user, const char *section, const char *name, const char *value);

void *Input_handle_array(const char *values, int data_type, const char *name);

static char* Input_lskip(const char* s);
static char* Input_find_char_next(const char* s);

// Longest array length to read in
#define MAX_ARRAY_LENGTH 50

#endif /* INPUT_H */
