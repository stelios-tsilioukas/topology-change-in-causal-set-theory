import sympy as sp

# ---- Product metric on S^2(a) x S^2(b), Euclidean ----
th1,ph1,th2,ph2 = sp.symbols('theta1 phi1 theta2 phi2')
a,b,Lam,G = sp.symbols('a b Lambda G', positive=True)
x=[th1,ph1,th2,ph2]
g = sp.diag(a**2, a**2*sp.sin(th1)**2, b**2, b**2*sp.sin(th2)**2)
ginv = g.inv()
n=4
Gam=[[[sum(ginv[l,s]*(sp.diff(g[s,m],x[k])+sp.diff(g[s,k],x[m])-sp.diff(g[m,k],x[s])) for s in range(n))/2
       for k in range(n)] for m in range(n)] for l in range(n)]
Gam=[[[sp.simplify(Gam[l][m][k]) for k in range(n)] for m in range(n)] for l in range(n)]
def Ric(m,k):
    e=0
    for l in range(n):
        e+=sp.diff(Gam[l][m][k],x[l])-sp.diff(Gam[l][m][l],x[k])
        for s in range(n):
            e+=Gam[l][l][s]*Gam[s][m][k]-Gam[l][k][s]*Gam[s][m][l]
    return sp.simplify(e)
R_mn=sp.Matrix(n,n,lambda i,j:Ric(i,j))
Rs=sp.simplify(sum(ginv[i,j]*R_mn[i,j] for i in range(n) for j in range(n)))
print("Ricci diag =", [sp.simplify(R_mn[i,i]) for i in range(n)])
print("Ricci/g    =", [sp.simplify(R_mn[i,i]/g[i,i]) for i in range(n)])
print("R          =", Rs)

# ---- 1. Pure gravity + Lambda:  R_mn = Lam g_mn ----
print("\n--- Vacuum Einstein with Lambda: R_mn = Lam g_mn ---")
sol=sp.solve([sp.Eq(1/a**2,Lam), sp.Eq(1/b**2,Lam)],[a,b],dict=True)
print("  forces a = b = 1/sqrt(Lambda):", sol)
Vol = (4*sp.pi*a**2)*(4*sp.pi*b**2)
SE = -sp.Rational(1,16)/sp.pi/G*(Rs-2*Lam)*Vol
print("  S_E(a=b=1/sqrt(Lam)) =", sp.simplify(SE.subs({a:1/sp.sqrt(Lam),b:1/sp.sqrt(Lam)})))
# de Sitter S^4 cross-check
r=sp.sqrt(3/Lam); VolS4=sp.Rational(8,3)*sp.pi**2*r**4
print("  [check] S^4:", sp.simplify(-(4*Lam-2*Lam)*VolS4/(16*sp.pi*G)), " (expect -3pi/(G Lam))")

# ---- 2. Einstein-Maxwell with 2-form flux on the two S^2 cycles ----
print("\n--- Einstein-Maxwell on S^2 x S^2, magnetic flux q on S^2_1, p on S^2_2 ---")
q,p=sp.symbols('q p', real=True)
F1sq = 2*q**2/a**4      # F_{mu nu}F^{mu nu} from block 1
F2sq = 2*p**2/b**4
Fsq  = F1sq+F2sq
# T_mn = F_ma F_n^a - 1/4 g_mn F^2 ; block-proportional
T1 = sp.Rational(1,4)*(F1sq-F2sq)     # coefficient of g_mn on block 1
T2 = sp.Rational(1,4)*(F2sq-F1sq)     # block 2
# Einstein: R_mn - 1/2 R g_mn + Lam g_mn = 8 pi G T_mn   (coefficients of g_mn)
E1 = sp.simplify(1/a**2 - Rs/2 + Lam - 8*sp.pi*G*T1)
E2 = sp.simplify(1/b**2 - Rs/2 + Lam - 8*sp.pi*G*T2)
print("  block1:", sp.Eq(sp.simplify(E1),0))
print("  block2:", sp.Eq(sp.simplify(E2),0))
print("  SUM of the two equations:", sp.Eq(sp.simplify(E1+E2),0))
print("  => Lambda =", sp.solve(sp.Eq(sp.simplify(E1+E2),0),Lam))
print("  Flux terms cancel in the sum: Lambda is forced POSITIVE, independent of q,p.")

# ---- 3. Same test for the axion (3-form) on S^2 x S^2 ----
print("\n--- Axion 3-form flux on S^2 x S^2 ---")
print("  H_3 flux needs a 3-cycle. Betti numbers of S^2 x S^2: b = (1,0,2,0,1) -> b_3 = 0.")
print("  => no 3-cycle, no quantised axion charge. The axion cannot support S^2 x S^2.")
