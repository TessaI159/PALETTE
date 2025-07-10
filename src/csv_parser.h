#ifndef CSV_PARSER_H
#define CSV_PARSER_H

#include <stddef.h>

#define CSV_MAX_FIELDS 64

// Parse a single CSV line into fields.
// Modifies `line` in-place by inserting null terminators.
// Returns number of fields parsed, or -1 on error.
int csv_parse_line(char *line, char **fields, int max_fields);

#endif // CSV_PARSER_H
