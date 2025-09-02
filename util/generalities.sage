#from sage.libs.gmp.all import mpz
from hashlib import shake_256
import random

def good_prime(d1,d2,N):
    for f in range(100):
        p=d1*d2*N*f-1
        if is_prime(p) and p%4==3:
            break
    return p
        
def set_seed(seed=None):
    if seed==None:
        set_random_seed()
        rseed=initial_seed()
    else: 
        set_random_seed(seed)
        rseed=seed
    return rseed

def PRNG(seed,n,d,lamda):
    Zd=IntegerModRing(d)
    set_random_seed(seed)
    a=Integer(Zd.random_element())
    rd=sample(range(2^lamda),n)
    return a,rd

def PRF(root,lamda):  
    '''Input: A lambda bits string
    output: a vector of 2 lambda bits strings'''
    import hashlib
    h= hashlib.shake_256()
    h.update(bytes(root))
    seed=h.digest(lamda//4)
    split = lamda//8
    return [seed[:split], seed[split:]]

def power_of_two(a):
    vect=Integer(a).digits(2)    
    vect2=[]
    vect2=[2^i for i, val in enumerate(vect) if val!=0]
    vect2.reverse()
    return vect2

#.............................................................
#......... CONSTRUCTION OF Seed Tree ........................
#............................................................
class Node_seed():
    def __init__(self,value,parent=None,left_child=None,right_child=None,h=0,i=0,nb_leaves=1):
        self.value=value
        self.parent=parent
        self.left_child=left_child
        self.right_child=right_child
        self.height=h
        self.position=i
        self.nb_leaves=nb_leaves
        self.cov=False      #relevant for checking whether a node can be considered as "internal node" (releaseSeed)
    
    def set_cover(self,b=bool):
        self.cov=b
        return None
        
    def set_nb_leaves(self,nb):
        self.nb_leaves=nb
        return None
        
    def update(self,left_child,right_child,h=None,i=None,val=None):
        if val!=None:
            self.value=val
        self.left_child=left_child
        self.right_child=right_child
        if h!=None:
            self.height=h
        if i!=None:
            self.position=i
    
    def other_child(self,node):
        if self.left_child==node:
            return self.right_child
        return self.left_child
    
    def is_leaf(self):
        return self.left_child==self.right_child==None
    
    def __repr__(self):
        return f"Node of binary tree of (height, position)=({self.height},{self.position}) with value {self.value}"


class SeedTree():
    def __init__(self,root,m,lamda):
        '''construction of seed tree with m leaves and rooting at root'''
        self.root=root
        self.m=m
        self.lamda=lamda
        self.nodes={}
        self.tree(root,m)

    def leaves(self): 
        return [self.nodes[nd] for nd in self.nodes if self.nodes[nd].is_leaf()]
        
    def tree_2(self,root,t):
        ''' construct a seed tree of t leaves, rooting at self.root where t is a power of 2'''
        lamda=self.lamda
        self.nodes[root.height,root.position]=root
        if t==1:
            self.nodes[root.height,root.position]=root
            return [root] 
        else:
            nb=t/2
            [lc,rc]=PRF(root.value,lamda)
            height=root.height+1
            pos_l=2*root.position
            pos_r=2*root.position+1
            #print(f'({height},{pos_l})---({height},{pos_r})\n')
            left_c=Node_seed(lc,root,left_child=None,right_child=None,h=height,i=pos_l,nb_leaves=nb)
            right_c=Node_seed(rc,root,left_child=None,right_child=None,h=height,i=pos_r,nb_leaves=nb)
            self.nodes[height,pos_l]=left_c
            self.nodes[height,pos_r]=right_c
            root.update(left_c,right_c)
            return self.tree_2(left_c,nb)+self.tree_2(right_c,nb)
        
    
    def tree_T(self,root,T):
        '''Construction of seed tree of t leaves where t is the sum of t_i in T, t_i being a power of 2'''
        lamda=self.lamda
        nb=sum(T)
        root.set_nb_leaves(nb)
        self.nodes[root.height,root.position]=root
        if len(T)==1:
            return self.tree_2(root,T[0])
        else:
            [lc,rc]=PRF(root.value,lamda)
            height=root.height+1
            pos_l=2*root.position
            pos_r=pos_l+1
            left_c=Node_seed(lc,root,left_child=None,right_child=None,h=height,i=pos_l)
            right_c=Node_seed(rc,root,left_child=None,right_child=None,h=height,i=pos_r)
            root.update(left_c,right_c)
            root.left_child.set_nb_leaves(T[0])
            L=self.tree_2(root.left_child,T[0])
            T.remove(T[0])
            return L+self.tree_T(root.right_child,T)
    
    def tree(self,root,t):
        '''construct a tree of leaves rooted at root, for any t'''
        T=power_of_two(t)
        return self.tree_T(root,T)
    
    def releaseSeed(self,ch,j):
        assert len(ch)==self.m
        nodes=self.nodes
        leaves=self.leaves()
        target_leaves=[leaves[i] for i, bit in enumerate(ch) if bit == j]
        if len(target_leaves)==0:
            return []
        if len(target_leaves)==1:
            return [[target_leaves[0].value,target_leaves[0].nb_leaves]]
        test=all([leaf.cov for leaf in target_leaves])
        leafs=target_leaves
        while not test:
            parents=[]
            for i in range(len(leafs)):
                if leafs[i].cov:
                    parents+=[leafs[i]]
                    continue
                par=leafs[i].parent
                if par in parents:
                    continue
                if par.other_child(leafs[i]) in leafs:
                    parents+=[par]
                else:
                    leafs[i].set_cover(True)
                    parents+=[leafs[i]]
            leafs=parents  
            test=all([leaf.cov for leaf in leafs])
        int_nodes=leafs 
        #print([node.value for node in target_leaves])
        int_seeds=[[node.value,node.nb_leaves] for node in int_nodes]
        return int_seeds
             
def recover_leaves(int_seeds,ch,j,lamda):
    seeds=[]
    for s in int_seeds:
        if len(s)==2:
            [val,nb]=s
            root=Node_seed(value=val,parent=None,left_child=None,right_child=None,h=0,i=0,nb_leaves=nb)
            int_tree=SeedTree(root,nb,lamda)
            leaves=int_tree.leaves()
            seeds+=[leaf.value for leaf in leaves]
        else:
            print('error')
    return seeds
    
#...................................................
#........ HACHING A DATA TO CHALLENGE SET ..........
#...................................................
def parse_hashs_t_w(d, s, t, w):
    
    '''
    IMPUT: d, a byte string digest
            s in {2,...,256}, t, w positive integers, w << t
    OUTPUT: h = [h_0, ..., h_{t-1}] in {0,..., s-1}^t with exactly w non-one elements
    We use a modified variant of Alg 9 in https://www.meds-pqc.org/spec/MEDS-2023-07-26.pdf
    '''
    
    h = [1]*t
    fb = 0
    bitlen_t = t.nbits()    # bit length of t
    bitlen_s = s.nbits()    # bit length of s
    bytelen_t = (bitlen_t+7)//8
    xof = shake_256()
    xof.update(d)

    def get_next_byte():
        nonlocal fb
        b = xof.digest(fb+1)[fb]
        fb += 1
        return b

    # Loop to set exactly w nonzero elements
    for _ in range(w):
        while True:
            fh = 0
            # compute fh from next bytelen_t bytes (little endian)
            for j in range(bytelen_t):
                b_byte = get_next_byte()
                fh += b_byte << (8*j)
            fh %= 2^bitlen_t
            if fh >= t or h[fh] != 1:
                continue
            # get h_fh for this location
            hfh = get_next_byte()
            hfh = hfh % (2^bitlen_s)
            if hfh == 1 or hfh >= s:
                continue
            h[fh] = hfh
            if 0 in h and 1 in h:
                break
    return h
