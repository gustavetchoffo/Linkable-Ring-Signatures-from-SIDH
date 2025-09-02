load('/home/gustave/sage/Ring_Signature_SIDH/Sigma_protocol_single.sage') 



seed=86903692069173005082753832226369567526
n=len(Ring)
a,rd=PRNG(seed,n,B,pp.lamda)
vect_com=[]
d=pp.d2()
for i in range(n):
    E3=CGL(Ring[i],a,d,2,return_kernel=False)
    j3=E3.j_invariant()
    vect_com.append(pp.C(j3.to_bytes(),rd[i].to_bytes(pp.lamda/8)))
root1,path=pp.hidden_merkle_tree(vect_com,vect_com[t])

seed=86903692069173005082753832226369567526
n=len(Ring)
a,rd=PRNG(seed,n,B,pp.lamda)
vect_com=[]
d=pp.d2()
for i in range(n):
    E3=CGL(Ring[i],a,d,2,return_kernel=False)
    j3=E3.j_invariant()
    vect_com.append(pp.C(j3.to_bytes(),rd[i].to_bytes(pp.lamda/8)))
root2=pp.hidden_merkle_tree(vect_com)
assert root1==root2

