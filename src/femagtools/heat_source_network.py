# -*- coding: utf-8 -*-
"""
  heat source network analysis

@author: werner b. vetter
"""

import xml.etree.ElementTree as ET
import numpy as np
import json


class HeatSourceNetwork:
    #Read the json-file with the heat source network definition

    #If you modify the netlist, you must call the command "process()" again

    #Parameter:
    #   filename (str): Filename of heat source network with extension (.hsn)

    def __init__(self, netlist):
        self.netlist = netlist
        self.process()

    def get_node_names(self):
        return [n['Name'] for n in self.netlist['Nodes']]

    def process(self):
        #Set up independent admittance matrix and the source vector

        self.nodes = []
        for branch in self.netlist['Branches']:
               branch_nodes = branch['Nodes']
               for nd in range(len(branch_nodes)):
                      node = branch_nodes[nd]
                      if not(node in self.nodes):
                             self.nodes.append(node)
        self.nodes.sort()

        self.G = np.zeros((len(self.nodes)-1,len(self.nodes)-1))
        self.P = np.zeros(len(self.nodes)-1)

        for branch in self.netlist['Branches']:
            ind1 = self.nodes.index(branch['Nodes'][0])-1
            ind2 = self.nodes.index(branch['Nodes'][1])-1
            if branch['Type'] == 'R_th':
                if ind1 != -1:
                    self.G[ind1,ind1] = self.G[ind1,ind1] + 1/branch['val']
                if ind2 != -1:
                    self.G[ind2,ind2] = self.G[ind2,ind2] + 1/branch['val']
                if ind1 != -1 and ind2 != -1:
                    self.G[ind1,ind2] = self.G[ind1,ind2] - 1/branch['val']
                    self.G[ind2,ind1] = self.G[ind2,ind1] - 1/branch['val']

            if branch['Type'] == 'Power_Source':
                if ind1 == -1:
                    self.P[ind2] = self.P[ind2] + branch['val']
                if ind2 == -1:
                    self.P[ind1] = self.P[ind1] + branch['val']

            if branch['Type'] == 'Temperature_Source':
                if ind1 == -1:
                    self.P[ind2] = self.P[ind2] + branch['T_val']/branch['R_val']
                    self.G[ind2,ind2] = self.G[ind2,ind2] + 1/branch['R_val']
                if ind2 == -1:
                    self.P[ind1] = self.P[ind1] + branch['T_val']/branch['R_val']
                    self.G[ind1,ind1] = self.G[ind1,ind1] + 1/branch['R_val']


    def solve(self):
        #Solve the system of equations

        self.T = np.linalg.solve(self.G, self.P)
        return(self.T)

    def draw(self,filename):
        #Creates an xml file of the network that can be displayed with draw.io

        #Parameter:
        #    filename (str): Filename of diagram with extension (.xml)

        HeatSourceDiagram.draw(self.netlist,filename)
        return()


def read(filename):
    with open(filename) as fp:
        netlist = json.load(fp)
    return HeatSourceNetwork(netlist)


class HeatSourceDiagram:
    #Init a diagram to create a diagram of a heat source netlist
    def __init__(self):
        self.mxfile = self.createFile()
        self.diagram = self.addDiagram(self.mxfile)
        self.graph = self.addGraphModel(self.diagram)
        self.root = self.addRoot(self.graph)
        mxCell_dict = {"id":"0"}
        self.addCell(self.root,mxCell_dict)
        mxCell_dict = {"id":"1", "parent":"0"}
        self.addCell(self.root,mxCell_dict)


    def writeFile(self,filename):
        tree = ET.ElementTree(self.mxfile)
        tree.write(filename, encoding="utf-8")


    def createFile(self):
        mxfile = ET.Element('mxfile')
        mxfile.attrib = {"host":"Electron",
                        "modified":"2023-03-08T08:47:19.445Z",
                        "agent":"5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) draw.io/20.8.16 Chrome/106.0.5249.199 Electron/21.4.0 Safari/537.36",
                        "etag":"5VmoL6vL6aGyHSSOdqAj",
                        "version":"20.8.16",
                        "type":"device"
                        }
        return mxfile

    def addDiagram(self,parent):
        diagram = ET.SubElement(parent, 'diagram')
        diagram.attrib = {"name":"Seite-1",
                         "id":"5XRmyqCaqz2rHzUy15r-"
                         }
        return diagram

    def addGraphModel(self,parent):
        mxGraphModel = ET.SubElement(parent, 'mxGraphModel')
        mxGraphModel.attrib = {"dx":"770",
                              "dy":"569",
                              "grid":"1",
                              "gridSize":"10",
                              "guides":"1",
                              "tooltips":"1",
                              "connect":"1",
                              "arrows":"1",
                              "fold":"1",
                              "page":"1",
                              "pageScale":"1",
                              "pageWidth":"827",
                              "pageHeight":"1169",
                              "math":"0",
                              "shadow":"0"
                              }
        return mxGraphModel

    def addRoot(self,parent):
        root = ET.SubElement(parent, 'root')
        root.attrib = {}
        return root

    def addCell(self,parent,cell_dict):
        mxCell = ET.SubElement(parent, 'mxCell')
        mxCell.attrib = cell_dict
        return mxCell

    def addGeometry(self,parent,geometry_dict):
        mxGeometry = ET.SubElement(parent, 'mxGeometry')
        mxGeometry.attrib = geometry_dict
        return mxGeometry

    def addPntArray(self,parent,pntArray):
        array = ET.SubElement(parent, 'Array')
        array.attrib = {"as":"points"}
        mxPoint = ET.SubElement(array, 'mxPoint')
        if isinstance(pntArray[0], list):
           for i in range(len(pntArray)):
                pnt = pntArray[i]
                mxPoint.attrib = {"x":str(pnt[0]), "y":str(pnt[1])}
        else:
           pnt = pntArray
           mxPoint.attrib = {"x":str(pnt[0]), "y":str(pnt[1])}
        return array

    def addNode(self,nodeKey,x,y):
        mxCell_dict = {"id":"node_"+str(nodeKey),
                       "value":str(nodeKey),
                       "style":"shape=waypoint;sketch=0;fillStyle=solid;size=6;pointerEvents=1;points=[];fillColor=none;resizable=0;rotatable=0;perimeter=centerPerimeter;snapToPoint=1;verticalAlign=top;spacingBottom=0;spacingTop=-5;fontFamily=Verdana;fontSize=12;",
                       "vertex":"1",
                       "parent":"1"
                       }
        cell = self.addCell(self.root,mxCell_dict)
        mxGeometry_dict = {"x":str(x-20),
                           "y":str(y-20),
                           "width":"40",
                           "height":"40",
                           "as":"geometry"
                           }
        geometry = self.addGeometry(cell,mxGeometry_dict)

        return geometry

    def addPowerSource(self,name,value,x,y,rotation=0):
        mxCell_dict = {"id":name,
                       "value":"Name = "+name+"\nLosses = "+str(value),
                       "style":"shape=stencil(vVXhboMgEH4a/i4I8QEWur7AHqCh9DZJFcyJ7fb2Q8CtteJWs42YmPsOvrvv44yEi66SLRBGjWyA8A1hbCN2okcE43bPtkfls8xvqGK2oDTG5xizMZZdC8pF8CRRy32dTnYO7RHO+uAShzYVoHZDlj8R+uj3DA8XyhrjSbQ13VXmIu/JpDb+LH2LZKn8e4oeyhi3vkIDDjC1HdEvlQkgbPvjSsW6Smy+EhceyQnmYi/V8RVtbw6z/UVPla0teiC+QzXC+DasZX043FXOw/kLng5Ajjq2llOc0cXFi0VYKZiGtdxVK4fpu0mM6caeIOdHuUg9MtTaXDCU91AEYzINLvi5SlpxV2fz4j5H4j/UQV3rtstbez2txa9N6zp3J858U37e3inJX7gbjt18cgGNf4QAfAA=);whiteSpace=wrap;html=1;labelPosition=center;verticalLabelPosition=bottom;align=center;verticalAlign=top;rotation="+str(rotation)+";",
                       "vertex":"1",
                       "parent":"1"
                       }
        cell = self.addCell(self.root,mxCell_dict)
        mxGeometry_dict = {"x":str(x-60),
                           "y":str(y-30),
                           "width":"120",
                           "height":"60",
                           "as":"geometry"
                           }
        geometry = self.addGeometry(cell,mxGeometry_dict)

        return geometry


    def addTemperatureSource(self,name,T,Rth,x,y,rotation=0):
        mxCell_dict = {"id":name,
                       "value":"Name = "+name+"\nT = "+str(T)+"\nRth = "+str(Rth),
                       "style":"shape=stencil(vVXbboMwDP2avE65jO15YusPdO9TSt0RNSQohHb7+4UEtpaSbGUVERKynRz7HJuAWN6UvAZEseIVIPaMKH2FqgbDbWvgba1bU7gwdTvKECYYB/sYbDrYvKmhsMF54EbwjexPNtboPRzF1vYYQpVghO2i7AXhJ7ene1heaKUciNCqOYucxB0YF8qdxR8BrE//2Vt3WbAdB1GBBdOXHbw/NHsHoqs/ZyLzMtHpTCx3nhhhlm94sX83ulXbyfqCpoWW2jhHePtsiLKVX2l+putVTMPpBo8HIAYdSosxjvBi+U4bmEkY+5Wuqubd9F0EhnClDxDTI0tCDwhSqBOE7BoIL0ykwISes6iRqyqbJvc9EkuwAylF3cSlPZ9WcrNpnaXu4//FJWxBcXdCyoW+KTIamvsZ2txC3lEZDzfASFNJt6jrwK9tSl+pFzen94Y/u3d8AQ==);whiteSpace=wrap;html=1;verticalAlign=top;spacing=0;labelPosition=center;verticalLabelPosition=bottom;align=center;verticalAlign=top;rotation="+str(rotation)+";",
                       "vertex":"1",
                       "parent":"1"
                       }
        cell = self.addCell(self.root,mxCell_dict)
        mxGeometry_dict = {"x":str(x-60),
                           "y":str(y-30),
                           "width":"120",
                           "height":"60",
                           "as":"geometry"
                           }
        geometry = self.addGeometry(cell,mxGeometry_dict)

        return geometry

    def addResistor(self,name,value,x,y,rotation=0):
        mxCell_dict = {"id":name,
                        "value":"Name = "+name+"\nR = "+str(value),
                        "style":"shape=stencil(vVRtbsMgDD0NfycC6wEq1h5gN6Cpt6AmEAFtt9sXMOlniNpoGoqU+Dl5fs+2QrhwjeyBMKplB4R/EMY+wSnnjQ2PAW8QfKcYHjFkNMfS9VB7BA/SKrlpATPOW7ODo9r6TKF0A1b5mOUrQpfhnXhxURutA4ky2t1krvKBTCodvqU/SJbL/+bobYFxHyp04MEiXiF68ZYBwtZPV6rmVWLjlbgISMkwFxtZ776t2evtqD7saW3aOB2K91SNML5OZ9qfjbMq9XB8wHfzLzGjspLhgi0uvoyFmX5pOtOqehmX7yExpDtzgFI72CT1wNAqfcWweIUiNaYgcKKfs6xVLykbN3feiP9wd7uod+oHIXlVqz9d1YeNTCj+JRNwAg==);whiteSpace=wrap;html=1;labelPosition=center;verticalLabelPosition=bottom;align=center;verticalAlign=top;rotation="+str(rotation)+";",
                        "vertex":"1",
                        "parent":"1"
                        }
        cell = self.addCell(self.root,mxCell_dict)
        mxGeometry_dict = {"x":str(x-60),
                           "y":str(y-10),
                           "width":"120",
                           "height":"20",
                           "as":"geometry"
                           }
        geometry = self.addGeometry(cell,mxGeometry_dict)

        return geometry

    def addCapacitor(self,name,value,x,y,rotation=0):
        mxCell_dict = {"id":name,
                       "value":"Name = "+name+"\nC = "+str(value),
                       "style":"shape=stencil(zVXRUsMgEPwaXh0C+tBHB+1/UHoapgkwBFv9+wJHtaYhajqOMplJbjfs3d6RCeFiaKUDwqiRPRD+QBgT0kmlg/XxORItorcUwwOGjJZYDg5UQHAvvZabDpAZgrc7OOhtKBLatOB1SCx/JPQ+vpMuLpQ1Jopoa4ZPzBkfxaQ2cS99RbGS/q1EN3cYu5ihhwAe8QbRD3MFIGz97UzNskxsOhMXEakZ5mIj1e7Z2xeznawPe6psl6ZD8Z6zEcbXec3782lWtR5OD3g0/5oyVlYzXLHFxZP1sNAvzWu+KifT4bsgTnRv91BrB5uVPil02pwprH4ikRtTKXCmn4usNc315t5PxL9zN2r79YP76rD/3eAWeBtL/Iq5vO3iY84o/mEycAQ=);whiteSpace=wrap;html=1;labelPosition=center;verticalLabelPosition=bottom;align=center;verticalAlign=top;rotation="+str(rotation)+";",
                       "vertex":"1",
                       "parent":"1"
                       }
        cell = self.addCell(self.root,mxCell_dict)
        mxGeometry_dict = {"x":str(x-60),
                           "y":str(y-18.75),
                           "width":"120",
                           "height":"37.5",
                           "as":"geometry"
                           }
        geometry = self.addGeometry(cell,mxGeometry_dict)

        return geometry

    def addConnection(self,name,startID,startPnt,endID,endPnt,pntArray=0):
        # Pnt = [X,Y] gibt die relative Position des Anschlusspunktes an
        # X=0 => links
        # X=1 => rechts
        # Y=0 => oben
        # Y=1 => unten
        mxCell_dict = {"id":name,
                       "value":"",
                       "style":"endArrow=none;entryX="+str(endPnt[0])+";entryY="+str(endPnt[1])+";exitX="+str(startPnt[0])+";exitY="+str(startPnt[1])+";html=1;rounded=0;entryDx=0;entryDy=0;entryPerimeter=0;exitDx=0;exitDy=0;",
                       "edge":"1",
                       "parent":"1",
                       "source":startID,
                       "target":endID
                       }
        cell = self.addCell(self.root,mxCell_dict)
        mxGeometry_dict = {"width":"50",
                           "height":"50",
                           "relative":"1",
                           "as":"geometry"
                           }
        geometry = self.addGeometry(cell,mxGeometry_dict)

        if isinstance(pntArray, list):
            self.addPntArray(geometry,pntArray)

        return

    def draw(netlist,file):
        #Creates an xml file of the netlist that can be displayed with draw.io

        #Parameter:
        #    netlist (list): List of branches
        #    file (str): Filename of diagram with extension (.xml)

        # list of nodes
        nodes = []
        for branch in netlist['Branches']:
            branch_nodes = branch['Nodes']
            for nd in range(len(branch_nodes)):
                node = branch_nodes[nd]
                if not(node in nodes):
                    nodes.append(node)
        nodes.sort()

        #list of connections
        connections = [None] * len(nodes)
        for node in nodes:
            connection = []
            for branch in netlist['Branches']:
                branch_nodes = branch['Nodes']
                if node in branch_nodes:
                    for nd in branch_nodes:
                        if nd!=node and not(nd in connection):
                            connection.append(nd)
            i = nodes.index(node)
            connections[i]=connection

        diagram = HeatSourceDiagram()

        #determine nodes position
        nd_pos = []
        xstep = 250
        ystep = 250
        xpos = 0
        ypos = ystep*len(connections)
        nd_pos.append([xpos,ypos])
        def_nd = [0]
        ypos = ypos-ystep
        for i in range(len(connections)):
            #new nodes in connection
            new_nd = []
            for nd in connections[i]:
                if nd>i and not(nd in new_nd) and not(nd in def_nd) :
                    new_nd.append(nd)
            xpos = -xstep*(1+(len(new_nd)-1)/2)
            draw = False
            for nd in new_nd:
                xpos = xpos+xstep
                nd_pos.append([xpos,ypos])
                def_nd.append(nd)
                draw = True
            #next level
            if draw:
                ypos = ypos-ystep


        #draw nodes
        for nd in def_nd:
            i = def_nd.index(nd)
            diagram.addNode(nd,nd_pos[i][0],nd_pos[i][1])

        # draw branches
        for branch in netlist['Branches']:
            #print(branch)
            xb = 0
            yb = 0
            for nd in branch['Nodes']:
                i = def_nd.index(nd)
                xb = xb+nd_pos[i][0]/len(branch['Nodes'])
                yb = yb+nd_pos[i][1]/len(branch['Nodes'])
            if len(branch['Nodes'])==2:
                i1 = def_nd.index(branch['Nodes'][0])
                i2 = def_nd.index(branch['Nodes'][1])
                dx = nd_pos[i1][0]-nd_pos[i2][0]
                dy = nd_pos[i1][1]-nd_pos[i2][1]
                angle = np.arctan2(dy,dx)*180/np.pi
            else:
                angle = 0

            if branch['Type']=='R_th':
                diagram.addResistor(branch['Name'],branch['val'],xb,yb,angle)

            if branch['Type']=='C_th':
                diagram.addCapacitor(branch['Name'],branch['val'],xb,yb,angle)

            if branch['Type']=='Power_Source':
                diagram.addPowerSource(branch['Name'],branch['val'],xb,yb,angle)

            if branch['Type']=='Temperature_Source':
                diagram.addTemperatureSource(branch['Name'],branch['T_val'],branch['R_val'],xb,yb,angle)

        # draw connections
        num_con = 0
        for branch in netlist['Branches']:
            num_con = num_con+1
            branch_id = branch['Name']
            y_offset = 0.5
            diagram.addConnection("connect"+str(num_con),\
                                  "node_"+str(branch['Nodes'][0]),[0,0],branch_id,[1,y_offset])
            num_con = num_con+1
            diagram.addConnection("connect"+str(num_con),\
                                  "node_"+str(branch['Nodes'][1]),[0,0],branch_id,[0,y_offset])


        diagram.writeFile(file)

        return
