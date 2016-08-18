#!/usr/bin/env python
import sys
import numpy as np

def findNotNone( l ):
    """return lower and upper indexes of not none values in list"""
    for i in range(len(l)):
        if l[i]:
            break
    for j in range(len(l)-1,-1,-1):
        if l[j]:
            break
    return (i,j)

def fit( f, B, losses, Bo, fo ):
    """fit coeffs of losses(f,B)=cw*(f/fo)**alfa*(B/Bo)**beta"""
    betacoeffs=[]
    # fit beta: pfe = cw_b*(B/Bo)**beta
    for i in range( len(f) ):
        if f[i]>0:
            pfe=np.array(losses).T[i]
            j,k= findNotNone(pfe)
            if j<=k:
                y=[ np.log10(p) for p in pfe[j:k+1] ]
                x=[ np.log10(b/Bo) for b in B[j:k+1]]
                A=np.vstack([x, np.ones(len(x))]).T
                beta,cw=np.linalg.lstsq(A, y)[0]
                betacoeffs.append(beta)

    # fit alfa: pfe = cw_f*(f/fo)**alfa
    alfacoeffs=[]
    for i in range( len(B) ):
        if i<len(losses):
            pfe=np.array(losses)[i]
            j,k= findNotNone(pfe)
            if f[j]<1e-2: j+=1
            if j<=k:
                y=[ np.log10(p) for p in pfe[j:k+1] ]
                x=[ np.log10(fx/fo) for fx in f[j:k+1]]
                A=np.vstack([x, np.ones(len(x))]).T
                alfa,cw=np.linalg.lstsq(A, y)[0]
                if alfa>1.2 and alfa<1.8:
                    alfacoeffs.append(alfa)

    if len(f)>1:
        alfa=np.average(alfacoeffs)
    else:
        alfa=1.3
    beta=np.average(betacoeffs)
    # fit cw: pfe = cw * (f/fo)**alfa * (B/Bo)**beta
    cw=[]
    for i in range( len(f) ):
        fx=f[i]
        for k in range( len(B) ):
            b=B[k]
            if k<len(losses) and i<len(losses[k]):
                pfe=losses[k][i]
                if pfe:
                    a=(fx/fo)**alfa*(b/Bo)**beta
                    if abs(a)>1e-3: cw.append(pfe/a)
    return (np.average(cw), alfa, beta )
