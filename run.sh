param=t150

cargo run --bin fingerprint --release -- \
	  output/$param.smiles --min-pts 15 --parameter $param --epsilon 0.72
