target := target/release/rsearch

$(target): src/main.rs src/lib.rs build.rs
	cargo build --release

RDKIT_SYS_PATH = /home/brent/Projects/rdkit-sys/include
RDKIT_PATH = /home/brent/omsf/clone/rdkit

export LD_LIBRARY_PATH=$(RDKIT_SYS_PATH):$(RDKIT_PATH)/build/lib

# parsing smiles file output from main chembl search to cluster by morgan
# fingerprint and report results
fingerprint:
	cargo build --bin fingerprint --release

# generate a set of inchi keys for our existing training and benchmark datasets
inchis.dat: src/bin/inchis.rs
	cargo run --bin inchis --release

test:
	cargo test -- $(FLAGS)

clippy:
	cargo clippy --workspace --tests $(ARGS)

run:
	cargo run --bin rsearch --release -- \
				--molecule-file data/chembl_33.sdf \
				--output-dir output \
				--search-params input/want.params \
				--forcefield input/tm.v2.offxml

all.smiles:
	find . -name '*.html' | xargs ./aggregate.awk > $@

# Usage:
# $(call profile, bin-name, args...)
profile = RUSTFLAGS=-g cargo build --release --bin $1 \
			&& perf record --call-graph dwarf target/release/$(strip $1) $2

prof.fingerprint:
	$(call profile, fingerprint, output/t48a.smiles -m 10 -p t18a)

