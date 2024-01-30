#!/bin/bash

set -e

params=( t141c t141a )

default_pts=1
default_eps=7

for param in ${params[@]}; do
	echo -e "\n$param"
	firefox output/$param.html

	brk=0
	while [[ $brk != "y" ]]; do
		read -p "Min pts [$default_pts]: " pts
		pts=${pts:-$default_pts}
		read -p "Eps [$default_eps]: " eps
		eps=${eps:-$default_eps}

		default_pts=$pts
		default_eps=$eps

		time target/release/fingerprint \
			output/$param.smiles \
			--min-pts $pts \
			--parameter $param \
			--epsilon 0.$eps | wc -l

		read -p "Accept? " brk
	done
done

