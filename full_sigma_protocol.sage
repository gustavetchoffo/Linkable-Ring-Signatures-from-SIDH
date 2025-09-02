import os
load('/home/gustave/sage/Ring_Signature_SIDH/Sigma_protocol_single.sage') 
class full_Commitment():
    def __init__(self,pp,sk,t,Ring):
        self.pp=pp
        self.sk=sk
        self.t=t
        self.Ring=Ring
        self.vect_com=[] #vector of object commitments
        self.com=[]
        self.seed_tree=None
        self.rd0_tree=None
        
        
    def commit(self):
        ''' Algorithm tild{P_1}'''
        lamda=self.pp.lamda
        size=pari(Integer(2**lamda))
        [val_seed_root,val_rd0_root]=sample(range(size),2)
        #seed_root=ZZ.random_element(2^lamda)
        #rd0_root=ZZ.random_element(2^lamda)
        m=ceil(lamda*1.71)
        
        seed_root=Node_seed(val_seed_root,parent=None,left_child=None,right_child=None,h=0,i=0,nb_leaves=m)
        self.seed_tree=SeedTree(seed_root,m,lamda)
        seed_leaves=self.seed_tree.leaves()
        
        seed_bytes=[leaf.value for leaf in seed_leaves]
        seed=[int.from_bytes(val) for val in seed_bytes]
        
        rd0_root=Node_seed(val_rd0_root,parent=None,left_child=None,right_child=None,h=0,i=0,nb_leaves=m)
        self.rd0_tree=SeedTree(rd0_root,m,lamda)
        rd0_leaves=self.rd0_tree.leaves()
        rd0_bytes=[leaf.value for leaf in rd0_leaves]
        rd0=[int.from_bytes(val) for val in rd0_bytes]
        
        for i in range(m):
            self.vect_com.append(Commitment(self.pp,self.sk,self.t,self.Ring,seed=seed[i],rd0=rd0[i]))
        self.com=[Co.com() for Co in self.vect_com]
        #print(f'...nb_com={len(self.vect_com)}')
        return self.com
        
def Challenge(pp,com,msg=None):
    lamda=pp.lamda
    m=ceil(pp.lamda*1.71)
    w=ceil(pp.lamda/7)
    data=b''
    [data:=data+c[0]+c[1] for c in com]
    if msg!=None:
        msg=msg.encode()
        data+=msg
        #d=randint(1,2**lamda-1)
        #data=d.to_bytes(lamda/8)
    ch=parse_hashs_t_w(data, 3, m, w)
    return ch

def Response(Commitment,ch):
    '''Commitent is a full_Commitment object'''
    sk=Commitment.sk
    pp=Commitment.pp
    t=Commitment.t
    Ring=Commitment.Ring
    com=Commitment.com
    seed_tree=Commitment.seed_tree
    rd0_tree=Commitment.rd0_tree
    vect_com=Commitment.vect_com
    m=len(com)
    resp={}

    #resp for ch=-1=2
    for i in range(m):
        if ch[i]==2:
            resp[i]=response(pp,sk,Ring,t,vect_com[i],2)
            
    #[resp[i]=response(pp,sk,Ring,t,vect_com[i],2) for i in range(m) if ch[i]==2]
    rd_internal_2=rd0_tree.releaseSeed(ch,2)
    
    #resp for ch=0
    for i in range(m):
        if ch[i]==0:
            resp[i]=response(pp,sk,Ring,t,vect_com[i],0)
            
    #[resp[i]=response(pp,sk,Ring,t,vect_com[i],0) for i in range(m) if ch[i]==0]
    rd_internal_0=rd0_tree.releaseSeed(ch,0)
    
    # resp for ch=1
    seed_internal=seed_tree.releaseSeed(ch,1)
    
    return [resp,rd_internal_2,rd_internal_0,seed_internal]
    
def verify(pp,Ring,com,ch,rsp):
    
    lamda=pp.lamda
    [resp,rd_internal_2,rd_internal_0,seed_internal]=rsp
    
    leaves_b=recover_leaves(seed_internal,ch,1,lamda)
    leaves=[int.from_bytes(leaf) for leaf in leaves_b]
    
    rd0_2_b=recover_leaves(rd_internal_2,ch,2,lamda)
    rd0_2=[int.from_bytes(leaf) for leaf in rd0_2_b]
    rd0_0_b=recover_leaves(rd_internal_0,ch,0,lamda)
    rd0_0=[int.from_bytes(leaf) for leaf in rd0_0_b]
    for i in range(len(ch)):
        if ch[i]==1:
            resp[i]=leaves[0]
            leaves.remove(leaves[0])
        else:
            if ch[i]==2:
                resp[i].append(rd0_2[0])
                rd0_2.remove(rd0_2[0])
            else:
                resp[i][0].append(rd0_0[0])
                rd0_0.remove(rd0_0[0])
    return all([verivication(pp,Ring,com[i],ch[i],resp[i]) for i in range(len(ch))])   
    
    
    
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
t1=time.time()
fcom=full_Commitment(pp,sk,t,Ring)
com=fcom.commit()
t2=time.time()
print('commitment:',t2-t1,'s')
ch=Challenge(pp,com)
t3=time.time()
print('challenge:',t3-t2,'s')
resp=Response(fcom,ch)
t4=time.time()
print('response:',t4-t3,'s')
verif=verify(pp,Ring,com,ch,resp)
t5=time.time()
print('verify:',t5-t4,'s',verif)

ch=2
vect_K_phi,E1=CGL(E0,sk,d1,k1,return_kernel=True)

t1=time.time()
Com=Commitment(pp,sk,t,Ring,seed=None,rd0=None)
com=Com.com()
resp=response(pp,sk,Ring,t,Com,ch,with_rd=True)
print(verivication(pp,Ring,com,ch,resp))
t2=time.time()
#xR=vect_K_phi[0]
#E=xR.parent().curve()
#R=xR.curve_point()
#P,Q=torsion_basis(E,d1)
#c=point_to_int(xR,P,Q)
#print(R==c*P+Q)
#print(R)
#print(P+c*Q)
'''
