/*
 * opt_parser.h
 *
 *  Created on: Mar 23, 2018
 *  Author: shrinath
 */

#ifndef OPT_PARSER_H_
#define OPT_PARSER_H_

#include "types.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <locale>

// This should probably be a class with input_params and valid options declared here.

/*!
 *  \brief Helper function to check if option is a valid option
 */
bool is_valid_opt(const std::string& Option);

/*!
 *  \brief Function to read options file and
 */
input_params_t read_opt_file(const std::string& Optfile_name);

#endif /* OPT_PARSER_H_ */
