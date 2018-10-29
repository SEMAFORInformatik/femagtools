'''
  jcf2msh Convert JMag JCF Files to Gmsh/msh

  Usage:

    jcf2msh <infile> <outfile>

  Author: Jonas Knöpfel

'''
import argparse
import re
from zipfile import ZipFile
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("outfile")

args = parser.parse_args()

infile = args.infile
outfile = args.outfile

with ZipFile(infile, 'r') as jcf:
    with jcf.open("mesh.dat") as f:
        stuff = [l.decode("utf-8") for l in f.readlines()]
    
mesh = {"Node": [], "Element": []}

def get_line():
    return re.split(" +|\n+|\*", stuff[i])
    
i = 0
line = []
while i < len(stuff):
    line = get_line()
    if line[-2] in ("Node", "Element"):
        reading = line[-2]
        i+=1
        for el in range(int(line[1])):
            line = get_line()
            i+=1
            mesh[reading].append(line)
    else:
        i+=1

ele_type = {3: "2", 4: "3"}

msh = []
msh.append("$MeshFormat\n2.2 0 8\n$EndMeshFormat")
msh.append("$Nodes")
msh.append(str(len(mesh["Node"])))

for n in mesh["Node"]:
    msh.append(" ".join([n[2], n[3], n[4], "0.0"]))
msh.append("$EndNodes")

msh.append("$Elements")
msh.append(str(len(mesh["Element"])))

for e in mesh["Element"]:
    e_id = e[3]
    num_v = int(e[2])
    n_ids = []
    for v in reversed(range(num_v)):
        n_ids.append(e[-(v+2)])
    msh.append(" ".join([e[1], ele_type[num_v], "2", "0", str(e_id)] + n_ids))
msh.append("$EndElements")
               
with open(outfile, "w") as outfile:
    outfile.write("\n".join(msh))
