#!/usr/bin/awk -f

# find the slow recap calls in the logging output from fingerprint

/finished leaves/ { finished[$NF] }
/started leaves/ { all[$(NF - 3)]; started[$(NF - 3)]; smiles[$(NF-3)] = $NF}

END {
	for (a in all) {
		if (a in started && !(a in finished))
			print a, smiles[a]
	}
}
