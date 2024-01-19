target := target/release/rsearch

$(target): src/main.rs src/lib.rs build.rs
	cargo build --release

RDKIT_SYS_PATH = /home/brent/Projects/rdkit-sys/include
RDKIT_PATH = /home/brent/omsf/clone/rdkit

export LD_LIBRARY_PATH=$(RDKIT_SYS_PATH):$(RDKIT_PATH)/build/lib

# parsing smiles file output from main chembl search to cluster by morgan
# fingerprint and report results
fingerprint:
	cargo run --bin fingerprint --release -- t18a.smiles --min-pts 10

# generate a set of inchi keys for our existing training and benchmark datasets
inchis.dat: src/bin/inchis.rs
	cargo run --bin inchis --release

test:
	cargo test -- $(FLAGS)

clippy:
	cargo clippy --workspace --tests

run:
	cargo run --release -- $(ARGS)
