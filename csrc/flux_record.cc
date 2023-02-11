#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

#include <zlib.h>
#include <sqlite3.h>

#include "flux_record.h"

// -- HELPER FUNCTIONS --

// Taken from a zlib example; I can't claim to fully understand zlib.
// If I really want to optimize stuff, there exist faster implementations
// of inflate than zlib's.
std::string decompress_zlib(const std::string & input) {

	// Initialize zstream
	z_stream zstream;

	zstream.opaque = Z_NULL;
	zstream.zalloc = Z_NULL;
	zstream.zfree = Z_NULL;

	zstream.avail_in = 0;
	zstream.next_in = Z_NULL;

	zstream.avail_out = 0; // pump priming

	int ret = inflateInit(&zstream);
	if (ret != Z_OK) {
		throw std::runtime_error("Could not initialize zlib stream! Error: " +
			std::string(zError(ret)));
	}

	// Set the input parameters.
	zstream.next_in = (Bytef *)input.c_str();
	zstream.avail_in = input.size();

	std::string decompressed;

	// Decompress a buffer's worth at a time.
	const int BUFLEN = 32768;

	char decompressed_buf[BUFLEN];
	while (zstream.avail_out == 0) {
		zstream.avail_out = BUFLEN;
		zstream.next_out = (Bytef *) decompressed_buf;

		ret = inflate(&zstream, Z_NO_FLUSH);

		if (ret != Z_OK && ret != Z_STREAM_END) {
			inflateEnd(&zstream);
			throw std::runtime_error("Could not decompress data! Error: " +
				std::string(zError(ret)));
		}

		decompressed += std::string(decompressed_buf,
			decompressed_buf + BUFLEN - zstream.avail_out);
	}

	// All done, do some cleanup, and then return.
	inflateEnd(&zstream);
	return decompressed;
}

// -- FLUX_RECORD --

flux_record::flux_record(int track_in, int side_in,
	comp_type compression_type, std::string fluxengine_data) {

	track = track_in;
	side = side_in;

	std::string uncompressed_data;

	// Decompress if necessary.
	switch(compression_type) {
		case CT_UNCOMPRESSED:
			uncompressed_data = fluxengine_data;
			break;
		case CT_ZLIB:
			uncompressed_data = decompress_zlib(fluxengine_data);
			break;
		default:
			// TODO: itos
			throw std::runtime_error(
				"Unknown FluxEngine compression type!");
	}

	// As far as I understand it, the FluxEngine native format
	// is simply: take each byte & 0x3F. If it is 0x3F, then that
	// just indicates that the flux delay is larger than 63, and
	// so we should accumulate them. Otherwise the flux delay
	// is just the byte & 0x3F.

	// Always skip the first such byte that we would otherwise
	// output because we don't know where the head was set down;
	// it could have been in the middle of a region between two
	// flux transitions. The last byte should be okay because
	// presumably the floppy drive doesn't report the fraction of
	// a clock before it was turned off.

	fluxes.reserve(uncompressed_data.size());

	int next_delay = 0;
	bool skip_first = true;
	for (unsigned char data_byte: uncompressed_data) {
		next_delay += data_byte & 0x3F;
		if ((data_byte & 0x3F) != 0x3F) {
			if (!skip_first) {
				fluxes.push_back(next_delay);
			}
			skip_first = false;
			next_delay = 0;
		}
	}

	if (next_delay != 0 && !skip_first) {
		fluxes.push_back(next_delay);
	}
}

// Read the flux data from the given FluxEngine filename. Returns a
// vector of flux data. If verbose is set to true, output stats as
// the data is read.

std::vector<flux_record> get_flux_record(std::string flux_filename,
	bool verbose) {

	sqlite3 * database;

	// Open FluxEngine sqlite database. Don't open if the file doesn't exist.
	if (sqlite3_open_v2(flux_filename.c_str(), &database,
		SQLITE_OPEN_READONLY | SQLITE_OPEN_FULLMUTEX, NULL) != SQLITE_OK) {
		sqlite3_close(database);
		throw std::runtime_error("Could not open flux file " + flux_filename);
	}

	// Now get every track and side. (NOTE: This differs from the MFM blocks,
	// which are categorized by sector, IIRC.)
	sqlite3_stmt * statement;
	std::string flux_query = "SELECT track, side, data, compression FROM zdata";
	if (sqlite3_prepare_v2(database, flux_query.c_str(),
		-1, &statement, NULL) != SQLITE_OK) {

		std::string err = "Could not prepare statement! sqlite error: " +
			std::string(sqlite3_errmsg(database));
		sqlite3_close(database);

		throw std::runtime_error(err);
	}

	std::vector<flux_record> flux_records;

	bool done = false;
	while (!done) {
		switch(sqlite3_step(statement)) {
			case SQLITE_DONE:
				done = true;
				break;
			case SQLITE_ROW: {
				int track = sqlite3_column_int(statement, 0);
				int side = sqlite3_column_int(statement, 1);
				comp_type compression = (comp_type)
					sqlite3_column_int(statement, 3);

				int data_size = sqlite3_column_bytes(statement, 2);
				const char * data = (const char *)sqlite3_column_blob(
					statement, 2);
				std::string data_str(data, data + data_size);

				flux_record this_row(track, side, compression, data_str);

				if (verbose) {
					std::cout << "T: " << track << ", S: "
						<< side << ", compression: " << compression
						<< ", num bytes: " << data_size << ", total bytes: " <<
						this_row.fluxes.size() << std::endl;
				}

				// Dump it into the flux_records vector. (Perhaps indexing it by head and
				// side would make more sense?)
				flux_records.push_back(this_row);
				break;
			}
			default: {
				std::cout << sqlite3_step(statement) << std::endl;
				std::string err = "Could not retrieve row! sqlite error: " +
					std::string(sqlite3_errmsg(database));
				sqlite3_close(database);

				throw std::runtime_error(err);
			}
		}
	}

	// All done, get outta here.
	sqlite3_finalize(statement);
	sqlite3_close(database);

	return flux_records;
}