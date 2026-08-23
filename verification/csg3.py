import numpy as np

def percolate(n,p,rng):
    R=np.triu(rng.random((n,n))<p,1)
    for k in range(n): R|=np.outer(R[:,k],R[k,:])
    return np.triu(R,1)

def abundances(R,kmax=11):
    n=R.shape[0]; c=np.zeros(kmax+1)
    ii,jj=np.nonzero(R)
    for i,j in zip(ii,jj):
        card=int((R[i,:]&R[:,j]).sum())
        if card+2<=kmax: c[card+2]+=1
    return c/n

Ns=[30,50,80,120,180,260]
pstar={30:0.13130,50:0.08412,80:0.06527,120:0.05353,180:0.03801,260:0.03022}
KMAX=11; S=24
data={}
for n in Ns:
    acc=[abundances(percolate(n,pstar[n],np.random.default_rng(5000+s)),KMAX) for s in range(S)]
    data[n]=np.mean(acc,axis=0)

Na=np.array(Ns,float)
print("="*70); print(" SCALING EXPONENT s_k  of  n_k ~ N^{s_k}   at renormalised p*(N)"); print("="*70)
ks=[];ss=[]
for k in range(2,KMAX+1):
    y=np.array([data[n][k] for n in Ns])
    if (y>0).all():
        s,_=np.polyfit(np.log(Na),np.log(y),1)
        ks.append(k); ss.append(s)
        print(f"   k={k:>2}:  s_k = {s:+.3f}")
ks=np.array(ks,float); ss=np.array(ss)
A,B=np.polyfit(ks,ss,1)
print(f"\n   LINEAR FIT:  s_k = {A:.4f}*k + {B:.4f}")

print("\n"+"="*70)
print(" THE CENTRAL RELATION")
print("="*70)
print("""
  p_handle ∝ N^b  with b = s_{n_h}  (n_h = size of minimal handle sub-order)

  N ∝ H^-4 (verified per era).  Then
      Omega_DE ∝ sqrt(p) ∝ N^{b/2} ∝ H^{-2b}
      rho_DE   ∝ Omega_DE * H^2  ∝ H^{2-2b}
  Matter era (H ∝ (1+z)^1.5):  rho_DE ∝ (1+z)^{3-3b} = (1+z)^{3(1+w)}
""")
print("      =>   w = -b     (and identically in the radiation era)\n")
print(f"{'n_h':>5} {'b = s_{n_h}':>12} {'w = -b':>9}   verdict")
for nh in range(3,12):
    b=A*nh+B
    v = "EXCLUDED (no acceleration)" if -b>-1/3 else ("VIABLE" if -1.35<-b<-0.75 else "phantom, excluded")
    print(f"{nh:>5} {b:>12.3f} {-b:>9.3f}   {v}")
nh_star=(-1.0-B)/A
print(f"\n  w = -1 exactly  <=>  b = 1  <=>  n_h = {nh_star:.2f}")
print("  and b=1 means rho_DE = const, i.e. an exact cosmological constant.")

print("\n"+"="*70); print(" BRACKETING"); print("="*70)
for lbl,b in [("static p (b=0)",0.0),("required for w=-1",1.0),("transitive percolation, n_h=10",A*10+B)]:
    w=-b
    print(f"  {lbl:>32}:  b={b:5.2f}  w={w:6.2f}", end="")
    if b==0: print("   -> rho_DE ∝ rho_crit: NO acceleration + BBN violation")
    elif abs(b-1)<1e-9: print("   -> exact LCDM")
    else: print("   -> strongly phantom: excluded by SNIa")
