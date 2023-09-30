#pragma once

#include <vector>

// Implements the CRC16 that's used on IBM MFM floppies: initial value
// 0xFFFF, polynomial 0x1021. This is sometimes called CCITT-CRC16.
// See also https://stackoverflow.com/a/44947877/6183577.

std::vector<unsigned short> create_crc16_table();
extern std::vector<unsigned short> crc16_table;

// Calculate a CCITT-CRC16 on the given region, slowly.
unsigned short crc16_bitwise(
	std::vector<unsigned char>::const_iterator start,
	std::vector<unsigned char>::const_iterator end);

// Calculate a CCITT-CRC16 on the given region with table lookup.
unsigned short crc16(
	std::vector<unsigned char>::const_iterator start,
	std::vector<unsigned char>::const_iterator end);

void test_crc16();