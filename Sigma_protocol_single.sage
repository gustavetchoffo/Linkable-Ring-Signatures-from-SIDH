load('/home/gustave/sage/Ring_Signature_SIDH/setup.sage')  

class KeyGen():
    def __init__(self,pp,rd):
        self.rd=rd
        self.pp=pp
        d1=pp.d1()
        k1=pp.k1
        E0=pp.initial_curve()
        A=pp.deg_phi()
        Z_A=IntegerModRing(A)
        set_random_seed(self.rd)
        self.sk=Integer(Z_A.random_element())
        self.pk=CGL(E0,self.sk,d1,k1,return_kernel=False)

    def __repr__(self):
        return f'pk={self.pk} and sk={self.sk}'

class Commitment():
    def __init__(self,pp,sk,t,Ring,seed=None,rd0=None):
        
        B=pp.deg_psi()
        ZB=IntegerModRing(B)
        self.seed=set_seed(seed)
        
        self.a=0
        self.rd=[]
        self.path=[]
        if rd0==None:
            self.rd0=ZZ.random_element(2^(pp.lamda))
        else:
            self.rd0=rd0
        self.pp=pp
        self.sk=sk
        self.t=t
        self.Ring=Ring
        self.vect_K_psi=[]
        self.vect_K_phi_prim=[] # actually phi_dual!!
    '''    
    def com_R(self):
        n=len(self.Ring)
        d2=self.pp.d2()
        k2=self.pp.k2
        Fq=pp.Fq
        self.a,self.rd=PRNG(self.seed,n,B,pp.lamda)
        vect_com_R=[]
        for i in range(n):
            E1=self.Ring[i]
            E3=CGL(E1,self.a,d2,k2,return_kernel=False)
            j3=E3.j_invariant()
            
            vect_com_R.append(self.pp.C(j3.to_bytes(),self.rd[i]))
            print('j3=',j3)
        root,self.path=self.pp.hidden_merkle_tree(vect_com_R,vect_com_R[self.t])
        
        return root
        
    def com_L(self):
        pp=self.pp
        d2=self.pp.d2()
        k2=self.pp.k2
        d1=self.pp.d1()
        k1=self.pp.k1
        
        E0=self.pp.initial_curve()
        vect_K_phi,E1=CGL(E0,sk,d1,k1,return_kernel=True)
        assert E1==Ring[self.t]
        vect_K_phi_dual=[kernel_dual(xK) for xK in vect_K_phi]
        vect_K_phi_dual.reverse()
        print('test1',vect_K_phi_dual[0].parent().curve()==E1)
        
        vect_K_psi_prim,E3=CGL(E1,self.a,d2,k2,return_kernel=True)
        
        self.vect_K_phi_prim,self.vect_K_psi=SIDH_lider_2(vect_K_phi_dual,d1,vect_K_psi_prim,d2)
        E2=self.vect_K_phi_prim[1].parent()
        j2=E2.j_invariant()
        hx=(self.rd0).str(base=16)
        bytes_rd=bytes.fromhex(hx)
        C0=pp.C(j2.to_bytes(),bytes_rd)
        #E2=vect_K_phi_prim[0].curve()
        return C0
    '''   
    def com(self):
        pp=self.pp
        n=len(self.Ring)
        d1=pp.d1()
        d2=pp.d2()
        k1=pp.k1
        k2=pp.k2
        t=self.t
        Fq=pp.Fq
        self.a,self.rd=PRNG(self.seed,n,B,pp.lamda)
        E0=self.pp.initial_curve()
        
        vect_com_R=[]
        for i in range(n):
            if i==t:
                vect_K_phi,E1=CGL(E0,sk,d1,k1,return_kernel=True)
                assert E1==Ring[self.t]
                vect_K_phi_dual=[kernel_dual(xK) for xK in vect_K_phi]
                vect_K_phi_dual.reverse()
                
                vect_K_psi_prim,E3=CGL(E1,self.a,d2,k2,return_kernel=True)
                
                j3=E3.j_invariant()
                C1=self.pp.C(j3.to_bytes(),self.rd[i].to_bytes(pp.lamda/8))
                vect_com_R.append(C1)
                self.vect_K_phi_prim,self.vect_K_psi=SIDH_lider_2(vect_K_phi_dual,d1,vect_K_psi_prim,d2)
                #print('test degree',self.vect_K_phi_prim[0].curve_point().order()==self.vect_K_phi_prim[0].curve_point().order()==d1)
                xK=self.vect_K_phi_prim[1]
                L1=xK.parent()
                xT=self.vect_K_psi[1]
                L2=xT.parent()
                phi1=KummerLineIsogeny(L1,xK,d1)
                psi1=KummerLineIsogeny(L2,xT,d2)
                assert phi1.codomain()==psi1.codomain()
                
                E2=phi1.codomain().curve()
                self.j2=E2.j_invariant()
                #hx=(self.rd0).str(base=16)
                #bytes_rd=bytes.fromhex(hx)
                C0=pp.C(self.j2.to_bytes(),self.rd0.to_bytes(pp.lamda/8))  
            else:
                E1=self.Ring[i]
                E3=CGL(E1,self.a,d2,k2,return_kernel=False)
                j3=E3.j_invariant()
                
                C1=pp.C(j3.to_bytes(),self.rd[i].to_bytes(pp.lamda/8))
                vect_com_R.append(C1)
        root,self.path=self.pp.hidden_merkle_tree(vect_com_R,vect_com_R[self.t])
        root1=self.pp.hidden_merkle_tree(vect_com_R)
        assert root==root1
        
        return [C0,root]


def response(pp,sk,R,t,com,ch,with_rd=False):
    assert ch in [2,0,1]
    if ch==2:
        if with_rd:
            resp=[com.vect_K_psi,com.rd0]
        else:
            resp=[com.vect_K_psi]
    else:
        if ch==1:
            resp=com.seed
        else:
            if with_rd:
                resp=[[com.j2,com.rd0],[com.vect_K_phi_prim,com.rd[t]],com.path]
            else:
                resp=[[com.j2],[com.vect_K_phi_prim,com.rd[t]],com.path]
    return resp
'''
TO DO: represent the isogenies with integers c=(c1,c2) so that it can be recovred by CGL
'''

def verivication(pp,R,com,ch,resp):
    [C0,C1]=com
    d1=pp.d1()
    d2=pp.d2()
    if ch==2:
        [vect_K_psi,rd0]=resp
        xT0=vect_K_psi[0]
        xT1=vect_K_psi[1]
        L0=xT0.parent()
        L1=xT1.parent()
        if L0.curve()!=pp.initial_curve():
            print('incorrect domain for psi')
            return False
        psi0=KummerLineIsogeny(L0,xT0,d2)
        if psi0.codomain()!=L1:
            print('incorrect intermediate domain for psi')
            return False
        psi1=KummerLineIsogeny(L1,xT1,d2)
        j2=(psi1.codomain().curve()).j_invariant()
        if pp.C(j2.to_bytes(),rd0.to_bytes(pp.lamda/8))!=C0:
            print('incorrect codomain for psi')
            return False
        return True
    if ch==1:
        seed=resp
        n=len(R)
        a,rd=PRNG(seed,n,B,pp.lamda)
        vect_com=[]
        d=pp.d2()
        for i in range(n):
            E3=CGL(R[i],a,d,2,return_kernel=False)
            j3=E3.j_invariant()
            C2=pp.C(j3.to_bytes(),rd[i].to_bytes(pp.lamda/8))
            vect_com.append(C2)
        root=pp.hidden_merkle_tree(vect_com)
        if root!=C1:
            print('incorrect rigth seed')
            return False
        return True
    if ch==0:
        [[j2,rd0],[vect_K_phi_prim,rd],path]=resp
        # Authantication of the left commitment
        if pp.C(j2.to_bytes(),rd0.to_bytes(pp.lamda/8))!=C0:
            print('Error authanticating the  domain of phi_prim')
            return False
        xK0=vect_K_phi_prim[0]
        xK1=vect_K_phi_prim[1]
        L0=xK0.parent()
        E3=L0.curve()
        j3=E3.j_invariant()
        Ct=pp.C(j3.to_bytes(),rd.to_bytes(pp.lamda/8))
        if pp.reconstructRoot(Ct,path)!=C1:
            print('Error authanticating the  codomain of phi_prim')
            return False
        L1=xK1.parent()
        phi0=KummerLineIsogeny(L0,xK0,d1)
        if phi0.codomain()!=L1:
            print('incorrect intermediate codomain for phi_prim')
            return False
        phi1=KummerLineIsogeny(L1,xK1,d1)
        L2=phi1.codomain()
        E2=L2.curve()
        if E2.j_invariant()!=j2:
            print('incorrect domain for phi_prim')
            return False
        return True

'''
#.....................................................
#................ test .............................
#.....................................................

#setup test
pp=SetUp(32)
p=pp.prime()
d1=pp.d1()
d2=pp.d2()
print('p=',factor(p+1),'-1')
E0=pp.initial_curve()
A=pp.deg_phi()
B=pp.deg_psi()
keys=[KeyGen(pp,rd) for rd in range(8)]
#print('keys=',keys)
Ring=[K.pk for K in keys]
#print('Ring=',Ring)
t=2
sk=keys[t].sk
#sk=116852
k1=2
ch=2
vect_K_phi,E1=CGL(E0,sk,d1,k1,return_kernel=True)

t1=time.time()
Com=Commitment(pp,sk,t,Ring,seed=None,rd0=None)
com=Com.com()
resp=response(pp,sk,Ring,t,Com,ch,with_rd=True)
print(verivication(pp,Ring,com,ch,resp))
t2=time.time()
xR=vect_K_phi[0]
E=xR.parent().curve()
R=xR.curve_point()
P,Q=torsion_basis(E,d1)
c=point_to_int(xR,P,Q)
#print(R==c*P+Q)
print(R)
#print(P+c*Q)



Z_A=IntegerModRing(A)
Z_B=IntegerModRing(B)
a=Integer(Z_A.random_element())
b=Integer(Z_B.random_element())

vect_xK_phi,E1=CGL(E0,a,d1,2,return_kernel=True)
L1=KummerLine(E1)
xK1,xK2=vect_xK_phi
Lp1=xK1.parent()
Lp2=xK2.parent()
phi1=KummerLineIsogeny(Lp1,xK1,d1)
phi2=KummerLineIsogeny(Lp2,xK2,d1)

#
print('test',phi2.codomain().curve()==E1)
t1=time.time()
vect_xK_phi_dual=[kernel_dual(xK) for xK in vect_xK_phi]
t2=time.time()
vect_xK_phi_dual.reverse()
[xK1,xK2]=vect_xK_phi_dual


vect_xK_psi_prim,E3=CGL(E1,b,d2,2,return_kernel=True)
xT1,xT2=vect_xK_psi_prim
print('test1:',xT1.parent()==L1)
print('test1:',xK1.parent()==L1)
#print(SIDH_diagram(xK1,d1,xT1,d2))

t5=time.time()

print(SIDH_lider_2(vect_xK_phi_dual,d1,vect_xK_psi_prim,d2))
t6=time.time()
print('t3=',t2-t1,'s')
'''

