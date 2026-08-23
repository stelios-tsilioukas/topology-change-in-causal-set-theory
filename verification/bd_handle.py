import numpy as np
from math import comb

# ---------- smeared BD weight (4D) ----------
def f4(n,eps):
    r=eps/(1.0-eps)
    return (1-eps)**n*(1-9*r*n+8*r**2*n*(n-1)-(4.0/3.0)*r**3*n*(n-1)*(n-2))

def sprinkle_box(N,T,Lx,seed):
    rng=np.random.default_rng(seed)
    t=rng.random(N)*T
    x=rng.random((N,3))*Lx-Lx/2
    return t,x

def causal_flat(t,x):
    dt=t[None,:]-t[:,None]
    dx2=((x[None,:,:]-x[:,None,:])**2).sum(-1)
    R=(dt>0)&(dt**2>dx2); np.fill_diagonal(R,False); return R

def add_wormhole(t,x,R,p,q,a):
    """Glue a ball of radius a at p to a ball of radius a at q:
       x~p and y~q become causally related as if adjacent (throat length 0)."""
    inA=np.linalg.norm(x-p,axis=1)<a
    inB=np.linalg.norm(x-q,axis=1)<a
    R2=R.copy()
    dt=t[None,:]-t[:,None]
    # travel through throat: proper distance = |x-p| + |y-q|
    dp=np.linalg.norm(x-p,axis=1); dq=np.linalg.norm(x-q,axis=1)
    dwh=dp[:,None]+dq[None,:]
    new=(dt>0)&(dt**2>dwh**2)&(inA[:,None]&inB[None,:])
    R2|=new
    new2=(dt>0)&(dt**2>(dq[:,None]+dp[None,:])**2)&(inB[:,None]&inA[None,:])
    R2|=new2
    n=R.shape[0]
    for k in range(n): R2|=np.outer(R2[:,k],R2[k,:])
    return np.triu(R2,1), int(inA.sum()), int(inB.sum())

def bd_core_smeared(R,eps):
    """ -N + eps * sum_{related pairs} f4(n,eps) ; n = # elements strictly between """
    n=R.shape[0]
    Ri=R.astype(np.int16); inter=Ri@Ri
    ii,jj=np.nonzero(R)
    cards=inter[ii,jj]
    tab={}
    tot=0.0
    for c in cards:
        c=int(c)
        if c not in tab: tab[c]=f4(c,eps)
        tot+=tab[c]
    return -float(n)+eps*tot

# ---------------------------------------------------------------
N=700; T=2.0; Lx=2.0; EPS=0.35
print("="*76); print(" BENINCASA-DOWKER ACTION OF A HANDLE CONFIGURATION"); print("="*76)
print(f" N={N}, smearing eps={EPS} (= retention f, by the thinning theorem)\n")

p=np.array([-0.55,0,0]); q=np.array([0.55,0,0])
print(f"{'a/l':>8} {'<n_A>':>7} {'dS_core':>12} {'+/-':>9} {'dS_EH = dS/(2sqrt6 pi)':>24}")
rows=[]
rng=np.random.default_rng(0)
ells=(N/(T*Lx**3))**(-0.25)     # discreteness length from density
for a in (0.16,0.22,0.28,0.34,0.40):
    ds=[]; nA=[]
    for s in range(9):
        t,x=sprinkle_box(N,T,Lx,s)
        R0=causal_flat(t,x)
        s0=bd_core_smeared(R0,EPS)
        R1,na,nb=add_wormhole(t,x,R0,p,q,a)
        s1=bd_core_smeared(R1,EPS)
        ds.append(s1-s0); nA.append(na)
    ds=np.array(ds)
    rows.append((a/ells,np.mean(ds),np.std(ds)/3,np.mean(nA)))
    print(f"{a/ells:>8.2f} {np.mean(nA):>7.1f} {np.mean(ds):>12.2f} {np.std(ds)/3:>9.2f}"
          f" {np.mean(ds)/(2*np.sqrt(6)*np.pi):>24.3f}")

print(f"\n discreteness length l = {ells:.4f} (box units)")
A=np.array([r[0] for r in rows]); D=np.array([r[1] for r in rows])
m=D>0
if m.sum()>=3:
    sl,ic=np.polyfit(np.log(A[m]),np.log(D[m]),1)
    print(f"\n FIT  dS_core ~ (a/l)^p  with  p = {sl:.2f}")
    print(f"      continuum expectation for a gravitational handle: p = 2")
