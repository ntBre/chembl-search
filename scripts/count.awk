#!/usr/bin/awk -f

# aggregate the clustering html results into a single data file for further
# processing by python

# Usage: find output -name '*.html' | xargs aggregate.awk

/SMILES/ {
	if (!match(FILENAME, /^.*\/([^.]+).html/, matches)) {
		printf "malformed input at line %d of %s\n", FNR, FILENAME > /dev/stderr
		exit 1
	}
	tot[matches[1]]++
}

END {
	for (p in tot)
		print p, tot[p] | "sort -nk2"
}
