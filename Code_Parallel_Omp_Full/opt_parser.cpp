/*
 * opt_parser.cpp
 *
 *  Created on: Mar 23, 2018
 */

#include "opt_parser.h"


// Add more options to these if needed
// Kind of unnecessary to have a input_params struct when we already have this
// but can't think of a better way to do this right now :(
OptionsSet valid_opts({
		"N_bath",
		"N_tslice",
		"N_cut",
		"N_sample",
		"total_time",
		"beta",
		"delta",
		"w_max",
		"eta",
		"prop_tstep",
		"Regression_test"});

bool is_valid_opt(const std::string& Option){
	if (valid_opts.find(Option) != valid_opts.end()) {
		return true; // the option exists in our set
	}
	else {
		return false;
	}

}

/*!
 *	\brief Helper function to trim whitespace from strings in place.
 */
static void trim(std::string &str) {
	// trim whitespace from the start of string.
	str.erase(str.begin(),std::find_if(str.begin(),str.end(), [](int ch) {
		return !std::isspace(ch);
	}));

	// trim whitespace from end.
	str.erase(std::find_if(str.rbegin(), str.rend(), [](int ch) {
	        return !std::isspace(ch);
	    }).base(), str.end());
}

// Someone please implement this better...
input_params_t read_opt_file(const std::string& Optfile_name) {

	input_params_t input;
	std::ifstream ifs(Optfile_name, std::ios::in);
	std::string line;
	OptionsMap opts;

	while (std::getline(ifs,line)) {

		std::string key;
		std::istringstream sub_line(line);
		if(std::getline(sub_line, key, '='))
		{
			std::string value;
			if (key[0] == '#') {  //change this to whatever you
							      //want comments to start with
				continue;
			}
			trim(key);
			if (!is_valid_opt(key)) {
				throw ("Invalid option" + key); //not really bothered to implement an exception so we'll just catch const char*
			}

			if (std::getline(sub_line,value)) {
				trim(value);
				opts[key] = value;
			}
		}
	}

	// Set all values of input. It would be nice to automatically iterate over members but i don't think this is possible without reflection.
	// need a better approach.
	// Whatif some variables aren't set? won't fail any checks currently will just crash.
	input.N_bath      = std::stoi(opts["N_bath"],nullptr);
	input.N_sample    = std::stoi(opts["N_sample"],nullptr);
	input.N_tslice    = std::stoi(opts["N_tslice"],nullptr);
	input.N_cut       = std::stoi(opts["N_cut"],nullptr);
	input.total_time  = std::stod(opts["total_time"],nullptr);
	input.beta        = std::stod(opts["beta"],nullptr);
	input.delta       = std::stod(opts["delta"],nullptr);
	input.w_max       = std::stod(opts["w_max"],nullptr);
	input.eta         = std::stod(opts["eta"],nullptr);
	input.prop_tstep  = std::stod(opts["prop_tstep"],nullptr);
	input.Regression_test = std::stoi(opts["Regression_test"], nullptr);

	return input;
}



