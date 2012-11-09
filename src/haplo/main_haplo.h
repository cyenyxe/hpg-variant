/* 
 * File:   main_haplo.h
 * Author: andrei
 *
 * Created on October 30, 2012, 9:29 PM
 */

#ifndef MAIN_HAPLO_H
#define	MAIN_HAPLO_H

#include <argtable2.h>

#include "globals.h"
#include "ld.h"
#include "file_handling.h"

#ifdef	__cplusplus
extern "C" {
#endif

int haplo_main(int argc, char *argv[], const char *configuration_file);
haplo_options_t *new_haplo_cli_options(void) ;
haplo_options_data_t *new_haplo_options_data(haplo_options_t *options);
int read_haplo_configuration(const char *filename, haplo_options_t *haplo_options, shared_options_t *shared_options);
void **parse_haplo_options(int argc, char *argv[], haplo_options_t *haplo_options, shared_options_t *shared_options);
void **merge_haplo_options(haplo_options_t *haplo_opts, shared_options_t *shared_options, struct arg_end *arg_end);


#ifdef	__cplusplus
}
#endif

#endif	/* MAIN_HAPLO_H */

