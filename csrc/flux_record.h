#pragma once

#include <vector>
#include <string>

enum comp_type { CT_UNCOMPRESSED = 0, CT_ZLIB = 1 };

class flux_record {
	public:
		int track, side;
		std::vector<int> fluxes;

		flux_record(int track_in, int side_in,
			comp_type compression_type, std::string fluxengine_data);
};

// Read the flux records from the given FluxEngine filename. Returns a
// vector of flux records. If verbose is set to true, output stats as
// the records are read.

std::vector<flux_record> get_flux_record(std::string flux_filename,
	bool verbose);