export RUST_LOG=debug
bin=target/release/fingerprint

eps=0.7
pts=10
fragment=-x

params=( t142j t142k t101 t141 t88 t18a t116j t72h t115h
# t165
# t141b
# t83h
# t141c
)

set -x

for param in ${params[@]}; do
	$bin output/$param.smiles -p $param -e $eps -m $pts $fragment -t 8

	if [[ $? -eq 0 ]]; then
		firefox output/$param.html
	else
		echo error in $param
	fi
done
