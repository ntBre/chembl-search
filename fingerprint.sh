export RUST_LOG=debug
bin=target/release/fingerprint

#param=t138a
eps=0.7
pts=2
fragment=-x

# params=( t118h t142j t116i t89 t142k t101 t141 t88 t18a t116j t72h t115h t165
# t141b t83h t141c )

params=( t141 )

set -x

for param in ${params[@]}; do
	$bin output/$param.smiles -p $param -e $eps -m $pts $fragment -t 8

	if [[ $? -eq 0 ]]; then
		firefox output/$param.html
	else
		echo error in $param
	fi
done
