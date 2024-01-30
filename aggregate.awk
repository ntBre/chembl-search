#!/usr/bin/awk -f

# aggregate the clustering html results into a single data file for further
# processing by python

# Usage: find . -name '*.html' | xargs aggregate.awk

/SMILES/ {
	gsub(/<\/?p>/, "")
	if (!match(FILENAME, /^.*\/([^.]+).html/, matches)) {
		printf "malformed input at line %d of %s\n", FNR, FILENAME > /dev/stderr
		exit 1
	}
	print matches[1], $2
}
