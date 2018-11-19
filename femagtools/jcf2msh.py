'''
  jcf2msh Convert JMag JCF Files to Gmsh/msh

  Usage:

    jcf2msh <infile> <outfile>

  Author: Jonas Knöpfel

'''
import argparse
import re
from zipfile import ZipFile
from collections import namedtuple

parser = argparse.ArgumentParser()
parser.add_argument("infile")
parser.add_argument("outfile")

args = parser.parse_args()

infile = args.infile
outfile = args.outfile

with ZipFile(infile, 'r') as jcf:
    with jcf.open("mesh.dat") as f:
        f = [l.decode("utf-8") for l in f.readlines()]
        f = iter(f)
        
def next_line(f, n=1):
    for i in range(n):
        line = next(f)
    line = line.strip()
    line = re.split(" +|,", line)
    return line


Node = namedtuple("Node", ["id", "x", "y"])
nodes = []

Element = namedtuple("Element", ["nodes", "entity"])
elements = []

while True:
    try:
        line = next_line(f)
        
        if line[-1] == "*Node":
            for n in range(int(line[0])):
                line = next_line(f)
                nodes.append(Node(line[1], line[2], line[3]))
                
            replace = {}
            for i, n1 in enumerate(nodes):
                for n2 in nodes[i+1:]:
                    if n1.x == n2.x and n1.y == n2.y and n1.id != n2.id:
                        if n1.id not in replace:
                            replace[n2.id] = n1.id
                            
            nodes = [n for n in nodes if n not in replace]
                
        if line[-1] == "*Element":
            for n in range(int(line[0])):
                line = next_line(f)
                if line[1] == "3":
                    elements.append(
                        Element([line[-3], line[-2], line[-1], None], line[2]))
                elif line[1] == "4":
                    elements.append(
                        Element([line[-4], line[-3], line[-2], line[-1]], line[2]))
                    
            for e in elements:
                for i, n in enumerate(e.nodes):
                    if n in replace:
                        e.nodes[i] = replace[n]
                    
    except StopIteration:
        break

msh = ["$MeshFormat\n2.2 0 8\n$EndMeshFormat"]

msh.append("$Nodes")
msh.append(str(len(nodes)))
for i, n in enumerate(nodes):
    msh.append("{} {} {} 0.0".format(n.id, n.x, n.y))
msh.append("$EndNodes")

msh.append("$Elements")
msh.append(str(len(elements)))
for i, e in enumerate(elements):
    if e.nodes[3]:
        msh.append("{} 3 2 0 {} {} {} {} {}".format(
            i + 1, e.entity, e.nodes[0], e.nodes[1], e.nodes[2], e.nodes[3]))
    else:
        msh.append("{} 2 2 0 {} {} {} {}".format(
            i + 1, e.entity, e.nodes[0], e.nodes[1], e.nodes[2]))
msh.append("$EndElements")
               
with open(outfile, "w") as outfile:
    outfile.write("\n".join(msh))
