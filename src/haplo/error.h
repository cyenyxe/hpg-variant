/*
 * error.h
 *
 *  Created on: Oct 16, 2012
 *      Author: andrei
 */

#ifndef ERROR_H_
#define ERROR_H_

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#define ERR_10 ""

//void throw_error(bool is_exit, const char *error_template, ...);
void throw_error(const char *error_template, ...);


#endif /* ERROR_H_ */
