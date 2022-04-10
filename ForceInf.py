"""
- Bayesian Force inference in epithelial tissue
- Perform MAP estimation of cell membrane tensions and cell pressures based on force balance relationship at vertices

  As presented in
    1. Shuji Ishihara and Kaoru Sugimura
    "Bayesian inference of force dynamics during morphogenesis"
    Journal of Theoretical Biology 313, p.201-211 (2012)

    2. S Ishihara, K Sugimura, SJ Cox, Isabelle Bonnet, Y Bellaïche, F Graner
    Comparative study of non-invasive force and stress inference methods in tissue
    The European Physical Journal E 36, p.1-13 (2013)

"""

import lib.ForceInf_lib
import lib.Out_lib
import lib.EBayesSP
import sys
import time
import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt


if __name__ == '__main__':

    print('# =======  Load file  ===================')
    filename = './Sample/sample/Vertex/VDat_sample.txt'

    #outhead = filename.split('/')[-1].split('.')[0])
    outhead = 'VDat_sample'

    Outfile = outhead+'_TP.txt'
    TensionFigure = outhead+'_Tension.eps'
    PressureFigure = outhead+'_Pressure.eps'

    print(' inputfile = ', filename, '\n outputfile =', Outfile,'\n')

    print('# =======  Read file  ===================')
    [x,y,edge,cell,Rnd,CELL_NUMBER,E_NUM,V_NUM,INV_NUM,R_NUM,stl,title] = lib.ForceInf_lib.loaddata( filename )
    # lib.Out_lib.DrawCells(x,y,edge,cell)  # show cells, for check

    print('\n# =======  Generate Matrices  =============')
    ERR_MAX = 1.0e-10 # this is used for checking the matrix satisfied momentum and angular momentum conservations
    [MM,C_NUM,X_NUM] = lib.ForceInf_lib.GetMatrix_ForceEstimation(x,y,edge,cell,E_NUM,CELL_NUMBER,R_NUM,INV_NUM,Rnd,ERR_MAX,SPARSE=True)

    print('## constraint      :  C_NUM= %d   [ 2x(%d) ] ' % (C_NUM, INV_NUM) )
    print('## unknown factors :  X_NUM= %d   [ E_NUM+CELL_NUMBER= %d + %d ] ' % (X_NUM,E_NUM,CELL_NUMBER) )
    print('## rounding cells  :  R_NUM= %d   ' % (R_NUM) )

    HParameter_Number = 2
    B0 = sp.spdiags( [1.0]*E_NUM ,0,X_NUM,X_NUM, format="coo")
    G = np.hstack( ( np.ones(E_NUM,dtype=float), np.zeros(CELL_NUMBER,dtype=float)) )
    V = np.zeros( C_NUM, dtype = np.float64 )

    ## Empirical Bayesian force inference
    print('\n# =======  MAP estimation  =============')
    [mu,T,P] = lib.EBayesSP.getTP_by_EBayses(MM,V,B0,G,HParameter_Number,X_NUM,C_NUM,E_NUM,CELL_NUMBER)

    ## Outputs
    print('\n# =======  Show results  =============')
    lib.Out_lib.OutputresultsTP(filename,Outfile,x,y,edge,cell,T,P,mu)

    print('\n  Draw tensions and pressures   ')
    lib.Out_lib.Draw_Tension(x,y,T,edge,T_LINE_WIDTH = 2.0, tmin = 0.6, tmax = 1.5,savefile = TensionFigure )
    lib.Out_lib.Draw_Pressure(x,y,P,edge,cell, -0.02, 0.025, savefile = PressureFigure )

