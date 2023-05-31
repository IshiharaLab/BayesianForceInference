import lib.ForceInf_lib
import lib.Out_lib
import sys
import time
import numpy as np
from scipy.optimize import fmin
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar
import scipy.sparse as sp
import sparseqr # https://github.com/yig/PySPQR   , check help(sparseqr.qr)



# ---------------------------------------------------
#     Functions
# ---------------------------------------------------

def Residue(MM,ep):
    v = MM.dot(ep)
    return v.dot(v)


def get_ABIC( x, MM, V, B, G, Parameter_Number ):
#def get_ABIC( x, args ):
    Mu = x

    C_NUM = MM.shape[0] #  condition number
    X_NUM = MM.shape[1] # unknown vaiable number

    smu = np.sqrt(Mu)
    K = 1              # number of zero-eigen values of A
    UN = np.linalg.matrix_rank( np.transpose(B).dot(B) )  # rank of B'B
    M0 = X_NUM-UN;    # Number of zero-eigen values of B
    NKM = C_NUM+K-M0;

    sSa = sp.hstack((MM,V),format="coo")  # this matrix can be defined outside
    sSb = sp.hstack( (sp.coo_matrix(smu*B),sp.coo_matrix(smu*G)),dtype=float)
    sS = sp.vstack((sSa,sSb),format="coo")  # coo is better for SPQR used in get_ABIC, see help(sparseqr.qr)

    #Q, R, E, rank = sparseqr.qr( sS )
    Q, R, E, rank = sparseqr.qr( sS, economy = True ) # QR decomposition

    dhs = np.diag(R.toarray())
    odh = np.zeros(len(dhs))
    odh[E] = dhs[:]  # need to reoderding following E

    # h = R[0:,-1]
    F = odh[-1]*odh[-1]
    odhh = odh[:-1]
    dh =np.abs(odhh[odhh!=0]) #
    detlA = 2*np.sum( np.log(dh) )
    detlB = UN*np.log(Mu)

    ABIC = NKM+NKM*np.log(2.0*np.pi*F/(NKM))+detlA-detlB+2*Parameter_Number;

    print(' fmin: %e %e' % (x, ABIC) )
    return ABIC


def get_ABIC_for_fmin( x, *args ):
#def get_ABIC( x, args ):
    Mu = x
    Parameter_Number, UN, NKM, sSa, sSb0 = args

    if Mu <=0:
        return 1.0e10

    smu = np.sqrt(Mu)
    sSb = smu*sSb0   # coo
    sS  = sp.vstack((sSa,sSb),format="coo")  # coo is better for SPQR used in get_ABIC, see help(sparseqr.qr)

    #Q, R, E, rank = sparseqr.qr( sS, economy = True ) # QR decomposition
    Q, R, E, rank = sparseqr.qr( sS ) # QR decomposition

    dhs = np.diag(R.toarray())
    odh = np.zeros(len(dhs))
    odh[E] = dhs[:]  # need to reoderding following E

    # h = R[0:,-1]
    F = odh[-1]*odh[-1]
    odhh = odh[:-1]
    dh =np.abs(odhh[odhh!=0]) #
    detlA = 2*np.sum( np.log(dh) )
    detlB = UN*np.log(Mu)

    ABIC = NKM+NKM*np.log(2.0*np.pi*F/(NKM))+detlA-detlB+2*Parameter_Number;

    print(' fmin: %e %e %e' % (x, F, ABIC) )
    return ABIC



# ---------------------------------------------------
#      Main
# ---------------------------------------------------

def getTP_by_EBayses(MM,V,B,G,HParameter_Number,X_NUM,C_NUM,E_NUM,CELL_NUMBER):
    # B is nparray, not sparray, this should be improved.. 2019/04/23
    print('# =======  EBayes: Minimization of ABIC  =============')
    mu = 1.0  # first guess

    C_NUM, X_NUM = MM.shape[0], MM.shape[1] #  condition and unknown vaiable numbers
    UN = E_NUM        # rank of B*B
    K = 1             # number of zero-eigen values of A
    M0 = X_NUM-UN;    # Number of zero-eigen values of B
    NKM = C_NUM+K-M0;

    sSa  = sp.hstack( (MM, V[:,None]), format="coo" )
    sSb0 = sp.hstack( (B,  G[:,None]), format="coo" )
    pargs = ( HParameter_Number, UN, NKM, sSa, sSb0 )

    print('# === starf fmin. ===')
    print('#     scipy.optimize.minize_salr method = \'bounded\' ')
    #print('#     scipy.optimize.minize_salr method = \'brent\' ')
    start = time.time()
    
    res = minimize_scalar( get_ABIC_for_fmin,  args=pargs, bounds= (0,10), method='bounded', tol=1e-2)
    #res = minimize_scalar( get_ABIC_for_fmin,  args=pargs, bounds= (0,10), method='brent', tol=1e-2)

    elapsed_time = time.time() - start
    print ("elapsed_time for optimization:{0}".format(elapsed_time) + "[sec]")

    print('# =======  EBayes: MAP estimation  =============')
    mu = res.x
    smu = np.sqrt(mu)

    sS = sp.vstack( (MM, smu*B), format="coo" )
    sb = sp.vstack( (V[:,None], smu*G[:,None]), format="coo" )

    ep = sparseqr.solve( sS, sb, tolerance=None).toarray()
    ep = ep.reshape(len(ep))

    print(ep)

    ep[E_NUM:] = ep[E_NUM:] - sum(ep[E_NUM:])/CELL_NUMBER
    norm = np.mean( ep[0:E_NUM] )
    ep = ep/norm
    print( "Residue = %f " % (Residue(MM,ep),) )

    T, P = ep[0:E_NUM], ep[E_NUM:X_NUM] # infered tension and pressure
    return [mu,T,P]




