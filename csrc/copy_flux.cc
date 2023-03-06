// This program takes the format [source flux file]
// [destination flux file] [track] [head].
// It will then copy the given track from the source FluxEngine file to
// the destination file, creating it (and the database) if it doesn't
// exist. This is used for creating single-sector flux files as test
// cases for certain bugs.

// Would it make sense to push this into flux_record as a write data
// function? Then there would be a pointless decompression and
// recompression of the flux data, but it would also make it possible
// to e.g. write repaired data.

// TODO? Error check if the track,head table already exists?

#include <iostream>
#include <sqlite3.h>

#include <stdio.h>
#include <sqlite3.h>

int main(int argc, char ** argv) {

	if (argc < 5) {
		std::cout << "Usage: " << argv[0] << " [source flux file] "
			"[destination flux file] [track] [head]\n";
		std::cout << "Copies the flux data for the given track and head "
			"from the source to\nthe destination file.\n";
		std::cout << "If the destination doesn't exist, it'll be created.\n";
		return -1;
	}

	sqlite3 * source_db, * dest_db;
	size_t track = atoi(argv[3]), head = atoi(argv[4]);

	std::string source_flux_name = argv[1],
		dest_flux_name = argv[2];

	std::cout << "Copying track " << track << " head " << head <<  " from "
		<< source_flux_name << " to " << dest_flux_name << std::endl;

	if (sqlite3_open_v2(source_flux_name.c_str(), &source_db,
		SQLITE_OPEN_READONLY | SQLITE_OPEN_FULLMUTEX, NULL) != SQLITE_OK) {
		sqlite3_close(source_db);
		sqlite3_close(dest_db);
		throw std::runtime_error("Could not open flux file " +
			source_flux_name);
	}

	if (sqlite3_open_v2(dest_flux_name.c_str(), &dest_db,
		SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_FULLMUTEX,
		NULL) != SQLITE_OK) {
		sqlite3_close(source_db);
		sqlite3_close(dest_db);
		throw std::runtime_error("Could not open flux file " +
			dest_flux_name);
	}

	sqlite3_exec(dest_db, "BEGIN TRANSACTION", NULL, NULL, NULL);

	// Create the tables in the destination database if they don't already
	// exist. No error handling, fix later if required

	for (std::string statement_str: {
		"CREATE TABLE IF NOT EXISTS properties (key TEXT UNIQUE NOT NULL PRIMARY KEY,  value TEXT);",
		"CREATE TABLE IF NOT EXISTS zdata (track INTEGER,  side INTEGER,  data BLOB,  compression INTEGER,  PRIMARY KEY(track, side));",
		"INSERT OR IGNORE INTO properties (key, value) VALUES (\"version\", 3);"
	}) {
		if (sqlite3_exec(dest_db, statement_str.c_str(), 0, 0,
			NULL) != SQLITE_OK) {

			std::string err = "Could not create tables: " +
				std::string(sqlite3_errmsg(dest_db));
			sqlite3_close(source_db);
			sqlite3_close(dest_db);

			throw std::runtime_error(err);
		}
	}

	sqlite3_exec(dest_db, "COMMIT", NULL, NULL, NULL);
	sqlite3_exec(dest_db, "BEGIN TRANSACTION", NULL, NULL, NULL);

	sqlite3_stmt * src_statement;
	// Transfer the relevant track.
	std::string flux_query = "SELECT * FROM zdata WHERE track=? AND side=?";
	if (sqlite3_prepare_v2(source_db, flux_query.c_str(),
		-1, &src_statement, NULL) != SQLITE_OK) {

		std::string err = "Could not prepare statement! sqlite error: " +
			std::string(sqlite3_errmsg(source_db));
		sqlite3_close(source_db);

		throw std::runtime_error(err);
	}

	sqlite3_stmt * dest_statement;
	// Transfer the relevant track.
	std::string dest_query = "INSERT INTO zdata VALUES (?, ?, ?, ?)";
	if (sqlite3_prepare_v2(dest_db, dest_query.c_str(),
		-1, &dest_statement, NULL) != SQLITE_OK) {

		std::string err = "Could not prepare statement! sqlite error: " +
			std::string(sqlite3_errmsg(dest_db));
		sqlite3_close(dest_db);

		throw std::runtime_error(err);
	}

	sqlite3_bind_int (src_statement, 1, track); // one-indexed
	sqlite3_bind_int (src_statement, 2, head);

	bool done = false, found_something = false;
	while (!done) {
		switch (sqlite3_step(src_statement)) {
			case SQLITE_DONE:
				done = true;
				break;
			case SQLITE_ROW: {
				found_something = true;
				// zero-indexed because why not? :-P
				// (There's an explanation - https://stackoverflow.com/questions/15962798 -
				//  but it's certainly an ... interesting surprise.)
				int data_size = sqlite3_column_bytes(src_statement, 2);
				const char * data = (const char *)sqlite3_column_blob(src_statement, 2);
				int compression_type = sqlite3_column_int(src_statement, 3);

				sqlite3_bind_int(dest_statement, 1, track);
				sqlite3_bind_int(dest_statement, 2, head);
				sqlite3_bind_blob(dest_statement, 3, data, data_size, SQLITE_TRANSIENT);
				sqlite3_bind_int(dest_statement, 4, compression_type);

				std::cout << "Data size " << data_size << std::endl;

				if (sqlite3_step(dest_statement) != SQLITE_DONE) {
					std::cout << "Something happened while writing..." << std::endl;
				}
				sqlite3_finalize(dest_statement);
				break;
			}
			default:
				std::string err = "Could not retrieve row! sqlite error: " +
					std::string(sqlite3_errmsg(source_db));
				throw std::runtime_error(err);
		}
	}

	sqlite3_finalize(src_statement);
	sqlite3_close(source_db);
	sqlite3_exec(dest_db, "COMMIT", NULL, NULL, NULL);
	sqlite3_close(dest_db);

	if (found_something) {
		std::cout << "Copied data for track " << track
			<< ", head " << head << "\n";
	} else {
		std::cout << "Could not find any data for track " << track
			<< ", head " << head << "\n";
	}
}