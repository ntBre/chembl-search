#!/usr/bin/awk -f

# filter new.dat to see how many matches each parameter got

/^#/ { next }
/^t/ { param = $1; count[param] = 0; next } # explicit zero to cause output
{ count[param]++ }
END {
	for (param in count)
		printf "%5s %5d\n", param, count[param]
}
