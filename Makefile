target := target/release/rsearch

$(target): src/main.rs src/lib.rs build.rs
	cargo build --release

RDKIT_SYS_PATH = /home/brent/Projects/rdkit-sys/include
RDKIT_PATH = /home/brent/omsf/clone/rdkit

export LD_LIBRARY_PATH=$(RDKIT_SYS_PATH):$(RDKIT_PATH)/build/lib

chembl_33.smiles: chembl_33.sdf src/bin/parse.rs
	cargo run --bin parse --release $< $@

# parsing smiles file output from main chembl search to cluster by morgan
# fingerprint and report results
fingerprint:
	cargo run --bin fingerprint --release -- t18a.smiles --min-pts 10 -p t18a

# generate a set of inchi keys for our existing training and benchmark datasets
inchis.dat: src/bin/inchis.rs
	cargo run --bin inchis --release

test:
	cargo test -- $(FLAGS)

clippy:
	cargo clippy --workspace --tests

run: src/main.rs want.params chembl_33.smiles
	cargo run --bin rsearch --release -- $(ARGS)

# Usage:
# $(call profile, bin-name, args...)
profile = RUSTFLAGS=-g cargo build --release --bin $1 \
			&& perf record --call-graph dwarf target/release/$(strip $1) $2

prof.fingerprint:
	$(call profile, fingerprint, t18a.smiles -m 10 -p t18a)

prof.parse:
	$(call profile, parse, chembl_33.sdf chembl_33.smiles)
