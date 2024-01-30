#!/usr/bin/awk -f

BEGIN {
	COMPLETE = 6
	head = "cargo run --bin fingerprint --release -- output/"
	tail = " 2> /dev/null > /tmp/smiles"
}

!/^ *#/ && NF < COMPLETE {
	smiles_file = "output/" $1 ".smiles"
	param = $1
	min_pts = NF > 3 ? $4 : 1
	eps = NF > 4 ? $5 : 0.7
	system(make_cmd(param, min_pts, eps))
	cmd = "wc -l /tmp/smiles"
	cmd | getline words
	close(cmd)
	split(words, arr)
	printf "%6s %14s %8s %7s %6.2f %6s\n", $1, $2, $3, min_pts, eps, arr[1]

	next
}

# { print }

function make_cmd(param, min_pts, eps) {
	return sprintf("%s%s.smiles --parameter %s --min-pts %d --epsilon %f %s",
				   head, param, param, min_pts, eps, tail)
}
