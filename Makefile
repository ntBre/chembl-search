target := target/release/rsearch

$(target): src/main.rs src/lib.rs build.rs
	cargo build --release

RDKIT_SYS_PATH = /home/brent/Projects/rdkit-sys/include
RDKIT_PATH = /home/brent/omsf/clone/rdkit

export LD_LIBRARY_PATH=$(RDKIT_SYS_PATH):$(RDKIT_PATH)/build/lib

test: $(target) /tmp/rsearch.out
	cargo test
ifdef FULL
	sed -n '/t18a/,/^$$/p' /tmp/rsearch.out | sed 1d | sort \
		| awk '!/^$$/ {print $$1}' > /tmp/rsearch.got
	sort t18a.smiles | awk '{print $$2}' > /tmp/rsearch.want
	diff /tmp/rsearch.got /tmp/rsearch.want || exit 1
endif

/tmp/rsearch.out: $(target)
	$(target) | tee /tmp/rsearch.out
