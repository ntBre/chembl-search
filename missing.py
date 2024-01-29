# compare want.params to files in output to see which parameters are not
# covered by chembl

from pathlib import Path


def load_coverage(filename):
    cover = {}
    with open(filename) as inp:
        for line in inp:
            if line.startswith("t"):
                sp = line.split()
                cover[sp[0]] = int(sp[2])
    return cover


with open("want.params") as inp:
    want = [line.strip() for line in inp]

train_cover = load_coverage("../projects/valence-fitting/coverage/v2-tm.dat")
test_cover = load_coverage(
    "../projects/valence-fitting/coverage/v2-tm.industry.dat"
)

for w in want:
    p = Path(f"output/{w}.smiles")
    if not p.exists():
        print(f"{w:<5} {train_cover[w]:5} {test_cover[w]:5}")
