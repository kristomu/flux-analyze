// First stab at factoring out some full-length (flux to sectors) decoding logic.

#include <memory>

#include "decoder.h"
#include "flux_record.h"

#include "pulse_train/flux_decoder.h"
#include "rabin_karp/rabin_karp.h"

class full_decoder {
	protected:
		virtual decoded_tracks decode_track(
			const std::vector<int> & fluxes) const = 0;

		IBM_preamble preamble_info;
		rabin_karp preamble_search;

	public:
		decoded_tracks decode_track(
			const flux_record & f) const;

		decoded_tracks decode_floppy(
			const std::vector<flux_record> & flux_records) const;

		decoded_tracks decode_floppy(
			std::string flux_filename) const;

		full_decoder(): preamble_search(preamble_info.A1_sequence,
				PREAMBLE_ID_A1) {
			preamble_search.add(preamble_info.C2_sequence,
				PREAMBLE_ID_C2);
		}

		virtual std::string name() const = 0;
};

// Ordinary (complex) ordinal decoder, without any
// brute-forcing.

class ordinal_full_decoder : public full_decoder {
	private:
		decoded_tracks decode_track(
			const std::vector<int> & fluxes) const;

	public:
		std::string name() const {
			return "Ordinal full decoder";
		}
};

// Applies the given pulse decoder to the whole flux,
// then segments and returns the decoded data.

// TODO: Make pulse decoders return names too and chain
// the given decoder's name into ours.

class once_through_decoder : public full_decoder {
	private:
		decoded_tracks decode_track(
			const std::vector<int> & fluxes) const;

		std::shared_ptr<pulse_decoder> flux_to_MFM_decoder;

	public:
		std::string name() const {
			return "Once-through (full flux) decoder";
		}

		once_through_decoder(
			std::shared_ptr<pulse_decoder> p_decoder_in) : full_decoder() {
			flux_to_MFM_decoder = p_decoder_in;
		}
};