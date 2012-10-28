/*
 * error.c
 *
 *  Created on: Oct 16, 2012
 *      Author: andrei
 */

#include "error.h"

//void throw_error(bool is_exit, const char *error_template, ...)
//{
//	va_list args;
//    va_start (args, error_template);
//	vprintf(error_template, args);
//	va_end (args);
//	if (is_exit)
//		exit(EXIT_FAILURE);
//}

inline void throw_error(const char *error_template, ...)
{
	va_list args;
	    va_start (args, error_template);
		vprintf(error_template, args);
		va_end (args);
			exit(EXIT_FAILURE);
}

