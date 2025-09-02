load('/home/gustave/sage/Ring_Signature_SIDH/full_sigma_protocol.sage') 

def RS_sign(pp,sk,t,Ring,msg):
    fcom=full_Commitment(pp,sk,t,Ring)
    com=fcom.commit()
    ch=Challenge(pp,com,msg)
    resp=Response(fcom,ch)
    return com,resp
    
def RS_verify(pp,Ring,msg,sigma):
    com,resp=sigma
    ch=Challenge(pp,com,msg)
    verif=verify(pp,Ring,com,ch,resp)
    if not verif:
        raise ValueError('Non authentic message')
    return verif
    


#.....................................................
#................ test .............................
#.....................................................

#setup test
pp=SetUp(32)
keys=[KeyGen(pp,rd) for rd in range(8)]
#print('keys=',keys)
Ring=[K.pk for K in keys]
#print('Ring=',Ring)
t=2
sk=keys[t].sk
msg='isogeny-based ring signature'
print('.....................................................')
for _ in range(10):
    t1=time.time()
    sigma=RS_sign(pp,sk,t,Ring,msg)
    t2=time.time()
    print('signature: time=',t2-t1,'s')
    t3=time.time()
    verif=RS_verify(pp,Ring,msg,sigma)
    t4=time.time()
    print('verification: time=',t4-t3,'s ',verif)
