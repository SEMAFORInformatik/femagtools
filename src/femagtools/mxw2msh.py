'''
  mxw2msh Convert Maxwell/msh Files to Gmsh/msh

  Usage:

    mxw2msh <infile> <outfile>

  Author: Jonas Knöpfel

'''
import argparse
import re
from collections import namedtuple

parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("outfile")

args = parser.parse_args()

infile = args.infile
outfile = args.outfile

with open(infile) as f:
    f = iter(f.readlines())
        
def next_line(f, n=1):
    for i in range(n):
        line = next(f)
    line = line.strip()
    line = re.split(" |,", line)
    return line

Node = namedtuple("Node", ["x", "y"])
nodes = []

Element = namedtuple("Element", ["n1", "n2", "n3", "entity"])
elements = []

while True:
    try:
        line = next_line(f)
        if line[0] == "$begin":
            
            if line[1] == "'points'":
                line = next_line(f, 4)
                
                while line[0] != "$end":
                    nodes.append(Node(line[3], line[5]))
                    line = next_line(f)
                    
            if line[1] == "'elements'":
                line = next_line(f)
                
                while line[0] != "$end":
                    elements.append(Element(line[3], line[5], line[7], line[9]))
                    line = next_line(f)
                    
    except StopIteration:
        break

msh = ["$MeshFormat\n2.2 0 8\n$EndMeshFormat"]

msh.append("$Nodes")
msh.append(str(len(nodes)))
for i, n in enumerate(nodes):
    msh.append("{} {} {} 0.0".format(i + 1, n.x, n.y))
msh.append("$EndNodes")

msh.append("$Elements")
msh.append(str(len(elements)))
for i, e in enumerate(elements):
    msh.append("{} 2 2 0 {} {} {} {}".format(i + 1, e.entity, e.n1, e.n2, e.n3))
msh.append("$EndElements")
               
with open(outfile, "w") as outfile:
    outfile.write("\n".join(msh))
