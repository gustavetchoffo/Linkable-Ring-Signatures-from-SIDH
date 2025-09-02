import inspect
import os
import sys
import sage.all
import time
import hashlib
load('/home/gustave/sage/Ring_Signature_SIDH/util/generalities.sage')
#load('/home/gustave/sage/SIDH_RS/isogenies_based_functions.sage')
load('/home/gustave/sage/Ring_Signature_SIDH/util/isogenies.sage')



#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#:::::::::::::::::::::::::::::::::: DESCRIPTION OF THE PROTOOL  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

class SetUp():
    def __init__(self,lamda,linkeable=False):
        assert lamda>0
        self.l1=3
        self.l2=5
        self.lamda=lamda
        self.k1=2
        self.k2=2
        self.N=1
        self.e1=ceil(self.lamda/log(self.l1,2))
        self.e2=ceil(self.lamda/log(self.l2,2))
        if linkeable:
            self.N=7^(ceil(lamda/log(7,2)))
        
    '''the functions d1, d2 return d1 and d2 such that p=d1.d2.f-1'''
    def d1(self):
        return self.l1^(self.e1)
    
    def d2(self):
        return self.l2^(self.e2)
    
    def prime(self):
        #p_128=2^128 * 3^81-1, 2^2 * 3^81 * 5^57 -1
        return good_prime(self.d1(),self.d2(),self.N)
        
    def torsion_basis(self,E,d,T=None):
        if T==None:
            return torsion_basis(E,d)
        R=compute_linearly_independent_point(E, T, d, x_start=0)
        return T,R
       
    def deg_phi(self): #deg(phi)=A
        return (self.d1())^(self.k1)
    def deg_psi(self): #deg(psi)=B
        return (self.d2())^(self.k2)
    def initial_curve(self):
        p=self.prime()
        Fp=GF(p)
        R=Fp["x"]
        self.Fq=GF(p^2,'i',modulus=x**2+1)
        return EllipticCurve(self.Fq,[0,6,0,1,0])
        
    def H(self,data):  #The hash function
        h= hashlib.shake_256()
        
        h.update(data)
        
        hashed=h.digest(self.lamda/4)
        
        
        return hashed
    
    def C(self,data,rd): #Commitment function
        d=data+rd
        h=self.H(d)
        return h
    
    def hidden_merkle_tree(self,leaves,sibling_leaf=None):
        '''
        Return root and path where root is the root for merkel tree of leaves
        When sibling_leaf is not None, it also return the path to sibling_leaf
        '''
        tree=leaves
        leaf=sibling_leaf
        if sibling_leaf!= None:
            path=[]
            j=tree.index(leaf)
            if j%2==0:
                path.append([tree[j+1],0])
            else:
                path.append([tree[j-1],1])
        # Build the tree
        while len(tree) > 1:
            if len(tree) % 2 == 1:
                tree.append(tree[-1])  # Duplicate the last element if odd
            new_tree = []
            for i in range(0, len(tree), 2):
                h=self.H(bytes(tree[i]) + bytes(tree[i + 1]))
                new_tree.append(h)
                if leaf!= None:
                    if j in [i,i+1] and len(tree)>2:
                        leaf=h
            tree = new_tree
            if sibling_leaf!= None:
                if len(tree)>=2:
                    j=tree.index(leaf)
                    if j%2==0:
                        path.append([tree[j+1],0])
                    else:
                        path.append([tree[j-1],1])

        if sibling_leaf!= None:
            return tree[0],path  # Return the root of the Merkle tree
        else:
            return tree[0]
    

    def reconstructRoot(self,leaf,path):
        '''Reconstruct a root from a leaf and a path'''
        root=leaf
        for i in range(len(path)):
            node=path[i]
            if node[1]==0:
                root=self.H(bytes(root)+bytes(node[0]))
            else:
                root=self.H(bytes(node[0])+bytes(root))
        return root

    def hash_in_chaset(self,data,w,m):
            """hash m in ChaSet=ch in {-1,0,1}^n having w non-one bits """
            #import hashlib
            #counter=-1
            #h= hashlib.shake_256()
            #h.update(bytes(data))
            #ch= int.from_bytes(h.digest(ceil(m*log(3,2)/8)))
            #while self.count_in_base3(ch,0)>k0:
             #   counter+=1
              #  h.update(bytes(data)+bytes(counter))
            digest=bytes(data)#h.digest(self.lamda//8)
            #return Integer(ch).digits(3)
            return parse_hashs_t_w(digest, 3, m, w)

    def __repr__(self):
        return f"SetUp for (Linkeable) Ring Signature with security parameter {self.lamda}:\n p={factor(self.p+1)}-1;\n E0={self.initial_curve()}"


print('hello')
'''
       
#.....................................................
#................ test .............................
#.....................................................

#setup test
pp=SetUp(128)
p=pp.prime()
d1=pp.d1()
d2=pp.d2()
print('p=',factor(p+1),'-1')
E0=pp.initial_curve()
#print(type(E0))
P,Q=pp.torsion_basis(E0,d1)

# isogenies file tes
t0=time.time()
L0=KummerLine(E0)

A=pp.deg_phi()
B=pp.deg_psi()
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
