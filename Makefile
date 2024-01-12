# tested with ib-dev-new environment

molecules.out: chembl_33.sdf.gz search.py
	python search.py --output $@

boost := $(addprefix -lboost_,system serialization iostreams)
rdkit := $(addprefix -lRDKit,GraphMol SmilesParse FileParsers SubstructMatch)

CPPFLAGS += -I/home/brent/omsf/clone/rdkit/Code -Wall -Wextra -pedantic -O2
LDFLAGS += -L/home/brent/omsf/clone/rdkit/build/lib $(boost) $(rdkit)

run: search
	LD_LIBRARY_PATH=/home/brent/omsf/clone/rdkit/build/lib ./search
