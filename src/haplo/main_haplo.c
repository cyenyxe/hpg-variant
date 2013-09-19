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

#include "main_haplo.h"
#include "haplo_runner.h"

int main(int argc, char *argv[]) {

    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_options = new_shared_cli_options(0);
    haplo_options_t *haplo_options = new_haplo_cli_options();
    void **argtable;

    // If no arguments or only -h / --help are provided, show usage
    if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
        argtable = merge_haplo_options(haplo_options, shared_options, arg_end(haplo_options->num_options + shared_options->num_options));
        show_usage(argv[0], argtable);
        arg_freetable(argtable, 32);
        return 0;
    } else if (!strcmp(argv[1], "--version")) {
        show_version("Haplotypes");
        return 0;
    }

    /* ******************************
     * 	    Execution steps	    *
     * ******************************/

    init_log_custom(LOG_DEFAULT_LEVEL, 1, "hpg-var-haplo.log", "w");
    
    array_list_t *config_search_paths = get_configuration_search_paths(argc, argv);
    const char *configuration_file = retrieve_config_file("hpg-variant.conf", config_search_paths);
    
    // Step 1: read options from configuration file
    int config_errors = read_shared_configuration(configuration_file, shared_options);
    config_errors &= read_haplo_configuration(configuration_file, haplo_options, shared_options);
    
    if (config_errors) {
        LOG_FATAL("Configuration file read with errors\n");
        return CANT_READ_CONFIG_FILE;
    }

    // Step 2: parse command-line options
    argtable = parse_haplo_options(argc, argv, haplo_options, shared_options);
    
    // Step 3: check that all options are set with valid values
    // Mandatory options that couldn't be read from the config file must be set via command-line
    // If not, return error code!
    int check_haplo_opts = verify_haplo_options(haplo_options, shared_options);
    if (check_haplo_opts > 0) {
        return check_haplo_opts;
    }
    
    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_options_data = new_shared_options_data(shared_options);
    haplo_options_data_t *haplo_options_data = new_haplo_options_data(haplo_options);
 
    init_log_custom(shared_options_data->log_level, 1, "hpg-var-effect.log", "w");
    
    // Step 5: Perform the requested task
    int ret_code = run_haplotyes_calculation(shared_options_data, haplo_options_data);

    //free_marker_array(all_markers, samples_names->size);
    free(haplo_options);
    //free_shared_options_data(shared_opts);
    arg_freetable(argtable, 32);
    
    return ret_code;
}


haplo_options_t *new_haplo_cli_options(void) {
    haplo_options_t *options = (haplo_options_t*) malloc (sizeof(haplo_options_t));
    options->alleleref = arg_int0(NULL, "alleleref", NULL, "Determine which algorithm will be used to count the alleles");
    options->hwcut = arg_dbl0(NULL, "hwcut", NULL, "Hardy-Weinberg cut");
    options->fgcut = arg_int0(NULL, "fgcut", NULL, "Failed genotype cut");
    options->mendcut = arg_int0(NULL, "mendcut", NULL, "Maximum number of mendelian errors in a variant");
    options->mafcut = arg_dbl0(NULL, "mafcut", NULL, "Limit for the minor allele frequency in a variant");
    options->num_options = NUM_HAPLO_OPTS;
    return options;
}

haplo_options_data_t *new_haplo_options_data(haplo_options_t *options) {
    haplo_options_data_t *options_data = (haplo_options_data_t*) calloc (1, sizeof(haplo_options_data_t));
    options_data->alleleref = *(options->alleleref->ival);
    options_data->hwcut = *(options->hwcut->dval);
    options_data->fgcut = *(options->fgcut->ival);
    options_data->mendcut = *(options->mendcut->ival);
    options_data->mafcut = *(options->mafcut->dval);
    return options_data;
}

void free_haplo_options(haplo_options_t *options_data) {
    free(options_data);
}