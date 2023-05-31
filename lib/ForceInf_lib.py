#!/usr/bin/env python
# coding:UTF-8

import sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg


class EDGE:
    def __init__(self):
        self.junc1 = -1
        self.junc2 = -1
        self.x1 = 0
        self.y1 = 0
        self.x2 = 0
        self.y2 = 0
        self.in_out = 'N'
        self.ncell = [-1, -1]
        self.dx = 0.0
        self.dy = 0.0
        self.dist = 0.0

    def set_distance(self,x,y):
        self.dx = self.x1 - self.x2
        self.dy = self.y1 - self.y2
        self.dist = np.sqrt( self.dx*self.dx + self.dy*self.dy )

class CELL:
    def __init__(self):
        self.jnum = 0
        self.junc = np.empty(0,dtype = np.int64)
        self.in_out = 'N'
        self.edge = []
        self.NON_EDGE = np.nan
        self.area = 0

    def initialize_edge(self, NON_EDGE = 10000000000):
        self.NON_EDGE = NON_EDGE
        self.edge = np.ones(self.jnum, dtype = int)*NON_EDGE
        self.area = 0



'''
   Functions
     Set_cellNeighbors(tedge, tcell, E_NUM)  ## set edge.ncell and cell.edge
     loaddata(filename)
     GetMatrix_ForceEstimation(x,y,edge,cell,E_NUM,CELL_NUMBER,R_NUM,INV_NUM,Rnd,ERR_MAX):
     Show_L_curve(MM,B,G,E_NUM,CELL_NUMBER,mus):
'''

def Set_cellNeighbors(tedge, tcell, non_edge):
    #  version 2: much faster than previous version
    ejs = np.vstack( ([e.junc1 for e in tedge],[e.junc2 for e in tedge]) ).T
    cis = np.empty( (4,0), int)
    for (cc, c) in enumerate(tcell):
        c.initialize_edge( NON_EDGE = non_edge )
        cjuncs = np.vstack( (c.junc, np.roll(c.junc,-1)) )
        cjuncs = np.vstack( (cjuncs, cc*np.ones(c.jnum, int)) )
        cjuncs = np.vstack( (cjuncs, np.arange(c.jnum)) )
        cis = np.hstack( (cis,cjuncs) )

    for (e, ej) in enumerate(ejs):
        ne1 = np.intersect1d( np.where( ej[0] == cis[0]), np.where( ej[1] == cis[1]) )
        ne2 = np.intersect1d( np.where( ej[1] == cis[0]), np.where( ej[0] == cis[1]) )
        if (len(ne1)!=1 and len(ne2)!=1): assert('inconsistent')
        tedge[e].ncell[0] = int(cis[2,ne1])
        tedge[e].ncell[1] = int(cis[2,ne2])
        tcell[int(cis[2,ne1])].edge[int(cis[3,ne1])] =  e
        tcell[int(cis[2,ne2])].edge[int(cis[3,ne2])] = -e



def loaddata(filename,CHECK = False):
    """
    データを読み込み、変数に格納する
    Load data and store its information to variables
    """
    print('... Load data from', filename)
    stl = 0
    R_NUM=0
    CELL_NUMBER, E_NUM, V_NUM, INV_NUM = 0, 0, 0, 0
    RndV = np.empty(0, dtype = np.int64)
    RndE = np.empty(0, dtype = np.int64)
    RndC = np.empty(0, dtype = np.int64)
    title = ''

    for line in open(filename, 'r'):

        if '#' in line:
            if line[0:2] == '# ' and title == '':
                itemList = line[:-1].split()
                title = itemList[-1]
                print('title = \"',title,'\"')

            elif 'C_NUM' in line:
                itemList = line[:-1].split()
                CELL_NUMBER = int(itemList[2])
                print('CELL_NUMBER = ',  CELL_NUMBER)
                cell = [ CELL() for i in range(CELL_NUMBER)]  # Declare cell here

            elif 'IN_CNUM' in line:
                itemList = line[:-1].split()
                IN_CNUM = int(itemList[2])
                print('IN_CNUM = ', IN_CNUM )

            elif 'EX_CNUM' in line:
                itemList = line[:-1].split()
                R_NUM = int(itemList[2])
                print('R_NUM = ', R_NUM )

            elif '# V_NUM' in line:
                itemList = line[:-1].split()
                V_NUM = int(itemList[2])
                INV_NUM=V_NUM-R_NUM
                # N = 2*V_NUM   # !! Notice: Not 2*(V_NUM+R_NUM) !!
                x = np.zeros(V_NUM)  # Declare x
                y = np.zeros(V_NUM)  # Decalre y

            elif '# E_NUM' in line:
                itemList = line[:-1].split()
                E_NUM = int(itemList[2])
                print('E_NUM = ', E_NUM)
                edge = [ EDGE() for i in range(E_NUM)]  # Declare edge here

            elif '### Scale' in line:
                itemList = line[:-1].split()
                stl = int(itemList[2])
                print('scale' .stl)

        else:
            if 'V[' in line:
                itemList = line[:-1].split()
                vid = int(itemList[0].replace('V[','').replace(']',''))
                if vid >= len(x) : assert('')
                x[vid] = np.float64(itemList[1])
                y[vid] = np.float64(itemList[2])
                if 'Ext' in line:   # IDs of verteces at boundary
                    rid = int(itemList[0].replace('V[','').replace(']','') )
                    RndV = np.append(RndV, rid )

            elif 'E[' in line:
                #print(line)
                itemList = line[:-1].split()
                eid = int(itemList[0].replace('E[','').replace(']',''))
                if eid >= len(edge):   assert('edge Number is larger than expected')
                edge[eid].junc1 = int( itemList[1] )
                edge[eid].junc2 = int( itemList[2] )
                edge[eid].x1 = x[edge[eid].junc1]
                edge[eid].y1 = y[edge[eid].junc1]
                edge[eid].x2 = x[edge[eid].junc2]
                edge[eid].y2 = y[edge[eid].junc2]
                edge[eid].set_distance(x,y)

                if 'Ext' in line:  # IDs of edges at boundary
                    eid = int(itemList[0].replace('E[','').replace(']','') )
                    RndE = np.append(RndE, eid )
                    edge[eid].in_out = 'o'
                else:
                    edge[eid].in_out = 'i'


            elif 'C[' in line:
                itemList = line[:-1].split()
                cid = int(itemList[0].replace('C[','').replace(']',''))
                if cid >= len(cell):  assert('cell Number is larger than expected')
                cell[cid].jnum = int(itemList[1])
                for i in range(cell[cid].jnum):
                    cell[cid].junc = np.append( cell[cid].junc, np.int64(itemList[i+3]) )

                if 'Ext' in line:  # IDs of cells at boundary
                    cid = int(itemList[0].replace('C[','').replace(']','') )
                    RndC = np.append(RndC, cid )
                    cell[cid].in_out = 'o'
                else:
                    cell[cid].in_out = 'i'

    ### end of for line

    #    Rnd = np.append( RndJ, RndE,  RndC)

    Rnd = (RndV, RndE,RndC)
    #print(len(RndJ),RndJ)
    #print(len(RndE),RndE)
    #print(len(RndC),RndC)

    ### cell area : For checking validity of data
    ff=0
    for c in cell:
        area = 0.0
        darea = 0.0
        xx = x[c.junc]
        sx = np.roll(xx, -1)
        yy = y[c.junc]
        sy = np.roll(yy, -1)
        area = 0.5*sum(xx * sy -yy * sx)
        c.area = area
        if area<= 0.0:
            ff=1
            print('!! error: %c %d %f\n' & (c.in_out,i,area) )

    if ff==1:
        assert('')


    ###  set neibouring cells of edge: edges of cells %
    if CHECK:
        Set_cellNeighbors(edge, cell,E_NUM)

        ### check consistency of edge-cell relation 1
        ct = np.zeros(CELL_NUMBER, int);
        for e in edge:
            if e.ncell[0] == -1 or  e.ncell[1] == -1:
                assert('edge-cell error')
            ct[e.ncell[0]] += 1
            ct[e.ncell[1]] += 1

        for i,c in enumerate(cell):
            if c.in_out == 'o' and ct[i] != c.jnum-1:
                print('cell = %d count inconsistent  %d != %d ' % (i,ct[i],c.jnum-1) )
            elif c.in_out == 'i' and ct[i] !=  c.jnum:
                print('cell = %d count inconsistent  %d != %d ' % (i,ct[i],c.jnum) )


    ## Return data
    return [x,y,edge,cell,Rnd,CELL_NUMBER,E_NUM,V_NUM,INV_NUM,R_NUM,stl,title]




def GetMatrix_ForceEstimation(x,y,edge,cell,E_NUM,CELL_NUMBER,R_NUM,INV_NUM,Rnd,ERR_MAX,SPARSE = False):
    print('...    Calculate Coefficeint Matrix        ')
    #RndJ = Rnd[0,:]
    #RndE = Rnd[1,:]
    #RndC = Rnd[2,:]

    RndJ = np.array(Rnd[0],int)
    RndE = np.array(Rnd[1],int)
    RndC = np.array(Rnd[2],int)

    C_NUM = 2*(INV_NUM);       # 条件数
    X_NUM = E_NUM+CELL_NUMBER; # 未知変数

    MX = np.zeros( (INV_NUM+R_NUM,X_NUM) );  #  V_NUM = INV_NUM+R_NUM
    MY = np.zeros( (INV_NUM+R_NUM,X_NUM) );

    for (i, te) in enumerate(edge):
        MX[te.junc1,i] = -te.dx/te.dist;
        MY[te.junc1,i] = -te.dy/te.dist;
        MX[te.junc2,i] =  te.dx/te.dist;
        MY[te.junc2,i] =  te.dy/te.dist;

    for (i,c) in enumerate(cell):
        pn = c.jnum-1
        ij = c.junc

        hx = 0.0
        hy = 0.0
        if len(ij) == pn:
            print('\'not consistent polygon class. Exit.\'');
            sys.exit()

        if  not (ij[0] in RndJ):
            MX[ij[0],E_NUM+i] = 0.5*( y[ij[1]]-y[ij[pn]] )
            MY[ij[0],E_NUM+i] = 0.5*( x[ij[pn]]-x[ij[1]] )
            hx = hx + MX[ij[0],E_NUM+i]
            hy = hy + MY[ij[0],E_NUM+i]

        for j in range(1,pn):
            MX[ij[j],E_NUM+i] = 0.5*( y[ij[j+1]]-y[ij[j-1]] )
            MY[ij[j],E_NUM+i] = 0.5*( x[ij[j-1]]-x[ij[j+1]] )
            hx = hx + MX[ij[j],E_NUM+i]
            hy = hy + MY[ij[j],E_NUM+i]

        if not (ij[pn] in RndJ):
            MX[ij[pn],E_NUM+i] = 0.5*( y[ij[0]]-y[ij[pn-1]] )
            MY[ij[pn],E_NUM+i] = 0.5*( x[ij[pn-1]]-x[ij[0]] )
            hx = hx + MX[ij[pn],E_NUM+i]
            hy = hy + MY[ij[pn],E_NUM+i]

    ##

    MX = np.delete( MX, RndJ, 0 )
    MY = np.delete( MY, RndJ, 0 )
    MM = np.concatenate( (MX, MY), 0 )

    if MM.shape != (C_NUM, X_NUM):
        print('MM.shape', MM.shape)
        print('\'Not Valid Martix: incorrect size\'');
        sys.exit()


    ###    Check the validity of Matrix MM
    ze = np.arange( E_NUM )
    ze = np.delete( ze, RndE, 0 )
    zc = np.arange( CELL_NUMBER )
    zc = np.delete( zc, RndC, 0 )

    CM = np.sum( MX, 0 )
    cc = np.hstack( (CM[ze],CM[zc+E_NUM]) )
    err_x = np.abs( np.sum(cc) )

    CM = np.sum(MY,0)
    cc =  np.hstack( (CM[ze],CM[zc+E_NUM]) )
    err_y = np.abs( np.sum(cc) )

    # not boundary
    ix = np.delete( x, RndJ, 0)
    iy = np.delete( y, RndJ, 0)
    cx = np.hstack( (iy,-ix) )
    ccx = cx.dot(MM)
    ccx = np.delete( ccx, E_NUM+RndC, 0 )
    ccx = np.delete( ccx, RndE, 0 )
    err_ang = np.abs( np.sum(ccx) )

    #     print( '%e %e %e' % (err_y, err_x,err_ang) )
    # ;;;

    ce = np.hstack( ( np.zeros((E_NUM)),np.ones((CELL_NUMBER))) )
    err_iso = sum(abs(MM.dot(ce)))

    print('## Either err_ang or err_iso should be zero ( >ERR_MAX= %.2e ) ' % (ERR_MAX) );
    print('## err_x= %e   err_y= %e   err_ang=  %e  err_iso= %e\n' % (err_x,err_y,err_ang,err_iso) );
    if max([err_x, err_y, err_ang, err_iso]) > ERR_MAX:
        print('larger than %e : errors %e %e %e %e\n' % (ERR_MAX,err_x,err_y,err_ang,err_iso) )
        sys.exit()

    if SPARSE:
        sMM = sp.csr_matrix(MM)
        return [sMM, C_NUM, X_NUM]
    else:
        return [MM,C_NUM,X_NUM]


