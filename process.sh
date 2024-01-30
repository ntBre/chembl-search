#!/bin/bash

set -e

params=( t112 t142j t116i t164 t8 t113 t7 t155 t126 t73h t89 t74h t142k t54 t63
		 t55 t101 t133j t91 t100 t48a t141 t65 t30 t153 t142h t6 t12 t52 t56 t53
		 t19a t88 t57 t18a t37 t128 t125 t124 t136 t116j t72h t115h t71h t158
		 t122i t82a t134 t49 t123a t82h t135 t165 t99 t68 t141b t150 t81 t83h
		 t141c t141a )

default_pts=1
default_eps=7

for param in ${params[@]}; do
	echo -e "\n$param"

	brk=0
	while [[ $brk != "y" ]]; do
		read -p "Min pts [$default_pts]: " pts
		pts=${pts:-$default_pts}
		read -p "Eps [$default_eps]: " eps
		eps=${eps:-$default_eps}

		default_pts=$pts
		default_eps=$eps

		target/release/fingerprint \
			output/$param.smiles \
			--min-pts $pts \
			--parameter $param \
			--epsilon 0.$eps | wc -l

		firefox output/$param.html

		read -p "Accept? " brk
	done
done

