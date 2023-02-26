#pragma once

#include <cstddef>
#include <vector>
#include <array>

// Header for the different IBM floppy headers or address marks.
// The start of a chunk is marked by a synchronization header. There are
// four types:

// Abbrev. Full name                 Byte seq.  What it is
// IAM	   Index Address Mark        C2C2C2 FC, start of a track
// IDAM	   ID Address Mark           A1A1A1 FE, start of a sector
// DAM	   Data Address Mark         A1A1A1 FB, start of sector data
//		   Deleted Data Address Mark A1A1A1 F8,	start of deleted sector data

// Each address mark is followed by a different header structure,
// though the IAM has no header at all. Note that the byte sequence
// contains a deliberate error to signal out-of-band that it's an address
// mark, not ordinary data.

enum address_mark_t { A_IAM, A_IDAM, A_DAM, A_DDAM, A_UNKNOWN };

// http://nerdlypleasures.blogspot.com/2015/11/ibm-pc-floppy-disks-deeper-look-at-disk.html
const int MAX_DATALEN_IDX = 6;
const std::array<int, MAX_DATALEN_IDX+1> IDAM_datalen = {
	128, 256, 512, 1024, 2048, 4096, 8192};

// ------------- ID Address Mark ------------

class IDAM {
	public:
		// ID Address Mark info. Note these are raw data
		// so e.g. datalen would be an index into a vector. TODO:
		// Abstract that away, either here or in address_mark.

		int track, head, sector, datalen;
		unsigned short CRC;

		// The address mark doesn't include this, but since we're only
		// reading, we just check the validity of the CRC at read time.
		bool CRC_OK = false;
		bool truncated;

		// By convention, the raw bytes start at the IDAM preamble (A1A1A1FE).
		void set(std::vector<unsigned char> & raw_bytes);
		void print_info() const;

		bool operator<(const IDAM & other) const {
			if (track != other.track) { return track < other.track; }
			if (head != other.head) { return head < other.head; }
			if (sector != other.sector) { return sector < other.sector; }

			return false;
		}
};


// ------------ Data Address Mark (also used for Deleted Data) -----------

class DAM {
	public:
		std::vector<char> data;
		unsigned short CRC;

		// Auxiliary info
		bool CRC_OK;
		bool deleted; // Is it a DAM or a DDAM?

		void set(std::vector<unsigned char> & raw_bytes, int datalen);
		void print_info() const;
};

// ------------- Index Address Mark ------------

class IAM {
	public:
		// The IAM doesn't actually contain any data or metadata.

		void print_info() const;
};


// ------------- Container class ------------

// It's easier to put all the address mark types into one class
// than deal with the problems of reconciling strong typing with a
// list of different types of address marks. This is doubly true because
// a Data Address Mark's metadata is listed in the IDAM that precedes it;
// so we need to know the relative location of DAMs and IDAMs.

class address_mark {
	public:
		address_mark_t mark_type;
		// TODO: When we start doing caching, each of these should
		// probably have a one-time cache indicator so that we can
		// detect staleness if e.g. we add something to or remove
		// something from the byte stream. For now, this is just
		// an index into where the byte stream starts.
		size_t byte_stream_index;

		IDAM idam;
		IAM iam;
		DAM dam, ddam;

		// TODO? Some kind of reference of where this was located all
		// the way back to the pulse train, so that we can cross off
		// successfully decoded chunks and work by a process of
		// elimination.

		void set_address_mark_type(
			const std::vector<unsigned char> & bytes);
		void print_info() const;

		address_mark() {
			dam.deleted = false;
			ddam.deleted = true;
		}
};