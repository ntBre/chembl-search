#!/usr/bin/awk -f

/^t[0-9]+/ {
	print | "sort -nk2 -k1"
}
