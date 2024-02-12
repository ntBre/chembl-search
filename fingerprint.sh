export RUST_LOG=debug
bin=target/release/fingerprint

param=t61g
eps=0.7
pts=1
fragment=

set -x

$bin output/$param.smiles -p $param -e $eps -m $pts $fragment -t 1

if [[ $? ]]; then
	firefox output/$param.html
fi
