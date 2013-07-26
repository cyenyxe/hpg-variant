/*
 * Copyright (c) 2012-2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 * Copyright (c) 2012 Andrei Alic
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * hpg-variant is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-variant. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MAIN_HAPLO_H
#define	MAIN_HAPLO_H

#include <commons/argtable/argtable2.h>

#include "haplo.h"
#include "shared_options.h"
#include "ld.h"
#include "file_handling.h"
#include "hpg_variant_utils.h"

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

