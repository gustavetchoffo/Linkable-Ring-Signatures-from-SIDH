from two_isogenies.Theta_SageMath.utilities.supersingular import *
from two_isogenies.Theta_SageMath.montgomery_isogenies.kummer_isogeny import *
from two_isogenies.Theta_SageMath.utilities.discrete_log import *

#...............................................................................|
#................  FUNCTIONS RELATED TO ISOGENIES ..............................|
#...............................................................................|
def my_torsion_basis(E,d,T=None):
        if T==None:
            return torsion_basis(E,d)
        R=compute_linearly_independent_point(E, T, d, x_start=0)
        return R,T

def kumer_isogeny(E,xP,d):
    '''Compute the isogeny of domain E and kernel <P> wher P=(xP,yP) using Kummer line'''
    L = KummerLine(E)
    xP=L(xP)
    phi = KummerLineIsogeny(L, xP, d)
    return phi



def CGL(E,c,d,k,return_kernel=False):
    M=c.digits(d,padto=k)
    #print('c in CGL=',c)
    P,Q=torsion_basis(E,d) 
    #print('P=',P,'\n Q=',Q)
    if d%2==0:
        f=factor(d)[0][1]
        P,Q=fix_torsion_basis_renes(P, Q, f)
    K=P+M[0]*Q
    L = KummerLine(E)
    xK=L(K[0])
    phi_i=KummerLineIsogeny(L,xK,d)
    vec=[]
    #phi_i=Ei.isogeny(K,algorithm="factored")
    if return_kernel:
        vec.append(xK)
        
    for i in range(1,k):
        # Compute the curve from the Kummer Line
        Li=phi_i.codomain()
        Ei=Li.curve()
        # Speed up SageMath by setting the order of the curve
        p = Ei.base_ring().characteristic()
        Ei.set_order((p+1)**2, num_checks=0)
        xQ=Li(Q[0])
        x=phi_i(xQ)
        phi_Q=x.curve_point()
        P,Q=my_torsion_basis(Ei,d,phi_Q)
        #print('P=',P,'\n Q=',Q)
        K=P+M[i]*Q
        xK=Li(K[0])
        phi_i=KummerLineIsogeny(Li,xK,d)
        if return_kernel:
            vec.append(xK)
        #PHI.append(phi_i) 
    E1=phi_i.codomain().curve()
    E1.set_order((p+1)**2, num_checks=0)
    if return_kernel:
        return vec,E1
    return E1

def dual_isogeny(phi):
    L0=phi.domain()
    E0=L0.curve()
    L1=phi.codomain()
    E1=L1.curve
    d=phi.degree()
    '''
    P,Q=torsion_basis(E0,d)
    xP=L0(P[0])
    xR=phi(xP)
    if not xR:
        xQ=L0(Q[0])
        xR=phi(xQ)
    '''
    Qs = generate_point_order_D(E0, d, x_start=0)
    for Q in Qs:
        xQ=L0(Q[0])
        xR=phi(xQ)
        if xR:
            if d%2==0:
                f=factor(d)[0][1]
                R=xR.curve_point()
                if ((2**f)*R)[0]==0:
                    raise NotImplementedError('stil to understand how to compute the dual isogeny of power 2 Kummer isogeny')
            break
        continue
    phi_tild=KummerLineIsogeny(L1, xR, d)
    E=phi_tild.codomain().curve()
    #print(E0.is_isomorphic(E))
    #eta=E.isomorphism_to(E0)
    #print('E0:',E0.based_field(), '\n E:', E.based_fild())
    #phi_dual=eta*psi
    #T=isogeny_kernel(phi_dual)
    return phi_tild

def kernel_dual(xK):
    ''' 
    imput: a Kummer point genarating ker(phi)
    output: a Kummer point generating ker(phi_dual)
    ''' 
    L0=xK.parent()
    E0=L0.curve()
    K=xK.curve_point()
    d=K.order()
    phi=KummerLineIsogeny(L0,xK,d)
    L1=phi.codomain()
    Q=compute_linearly_independent_point(E0, K, d, x_start=0)
    xQ=L0(Q[0])
    xR=L1(phi(xQ).curve_point()[0])
    return xR

def SIDH_diagram(xK_phi,d1,xK_psi,d2):
    '''
    Construction of SIDH diagram (with rational isogenies)
    Input: a generator of Ker(phi) and a generator of Ker(psi_prim)
    Output: a generator of Ker(psi) and a generator of Ker(phi_prim)
    '''
    L0=xK_phi.parent()
    E0=L0.curve()
    if L0!=xK_psi.parent():
        raise ValueError('The codomaine of phi chould be the domaine of psi_prim')
    psi=KummerLineIsogeny(L0,xK_psi,d2)
    phi=KummerLineIsogeny(L0,xK_phi,d1)
    L1=phi.codomain()
    #compute the right isogeny
    xK_psi_prim=phi(xK_psi)
    psi_prim=KummerLineIsogeny(L1,xK_psi_prim,d2)
    #botom isogeny
    xK_phi_prim=psi(xK_phi)
    L2=psi.codomain()
    phi_prim=KummerLineIsogeny(L2,xK_phi_prim,d1)
    assert phi_prim.codomain()==psi_prim.codomain()
    return xK_phi_prim,xK_psi_prim
    
    
def SIDH_lider_2(vect_xK_phi,d1,vect_xK_psi,d2):
    '''
    For 2x2 SIDH lider
    Imput: vectors of points generating Ker(phi) and Ker(psi_prim) respectively
    Output: vectors of points generating Ker(psi) and Ker(phi_prim)
    '''
    xK={}
    xT={}
    #Ker(phi)=:K
    #Ker(psi)=:T
    xK[0,0]=vect_xK_phi[0]
    xK[0,1]=vect_xK_phi[1]
    xT[0,0]=vect_xK_psi[0]
    xT[1,0]=vect_xK_psi[1]
    xK[1,0],xT[0,1]=SIDH_diagram(xK[0,0],d1,xT[0,0],d2)
    xK[1,1],xT[0,2]=SIDH_diagram(xK[0,1],d1,xT[0,1],d2)
    xK[2,0],xT[1,1]=SIDH_diagram(xK[1,0],d1,xT[1,0],d2)
    xK[2,1],xT[1,2]=SIDH_diagram(xK[1,1],d1,xT[1,1],d2)
    
    return [xK[2,0],xK[2,1]],[xT[0,2],xT[1,2]] #phi_pri, psi_prim
