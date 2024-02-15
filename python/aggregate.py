import re
from pathlib import Path

# fmt: off
params = [
    "t133g", "t18b", "t87a", "t130h", "t120h", "t119h", "t62g", "t133h",
    "t61g", "t132h", "t138a", "t118h", "t142j", "t116i", "t89",
    "t142k", "t101", "t141", "t88", "t18a", "t116j", "t72h", "t115h",
    "t165", "t141b", "t83h", "t141c",
]
# fmt: on


output = Path("output")
pat = re.compile("</?p>")


def get_smiles(path):
    ret = []
    with open(path) as inp:
        for line in inp:
            if "SMILES" in line:
                line = pat.sub("", line)
                ret.append(line.split()[1])
    return ret


with open("td.dat", "w") as td, open("opt.dat", "w") as opt:
    for param in params:
        p = output.joinpath(param).with_suffix(".html")
        smiles = get_smiles(p)
        for s in smiles[:2]:
            print(param, s, file=td)
        for s in smiles[2:5]:
            print(param, s, file=opt)
