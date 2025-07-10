#include "csv_parser.h"
#include <string.h>

int csv_parse_line(char *line, char **fields, int max_fields) {
	int field_count = 0;
	char *ptr = line;

	while (*ptr && field_count < max_fields) {
		if (*ptr == '"') {
			// Quoted field
			ptr++;
			fields[field_count++] = ptr;

			while (*ptr) {
				if (*ptr == '"' && ptr[1] == '"') {
					// Escaped quote
					memmove(ptr, ptr + 1, strlen(ptr));
					ptr++;
				} else if (*ptr == '"') {
					*ptr++ = '\0';
					break;
				}
				ptr++;
			}

			// Skip comma after closing quote
			if (*ptr == ',') ptr++;
		} else {
			// Unquoted field
			fields[field_count++] = ptr;
			while (*ptr && *ptr != ',') ptr++;
			if (*ptr == ',') *ptr++ = '\0';
		}
	}

	if (*ptr != '\0') return -1; // Too many fields
	while (field_count < max_fields) fields[field_count++] = NULL;
	return field_count;
}
