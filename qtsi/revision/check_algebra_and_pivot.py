import numpy as np
I=np.eye(2,dtype=complex)
s1=np.array([[0,1],[1,0]],dtype=complex);s2=np.array([[0,-1j],[1j,0]]);s3=np.array([[1,0],[0,-1]],dtype=complex)
def Th(M): return s1@M.conj()@s1                      # Theta M Theta^-1, Theta=s1 K   (5.7)
n=lambda A: np.abs(A).max()

print("=== algebra checks ===")
E,Dl=2.0,1.0; H=E*s3+1j*Dl*s1
print(f" (1.6) H^2=(E^2-D^2)I      resid {n(H@H-(E**2-Dl**2)*I):.1e}")
print(f" (1.4) s3 H s3 = H^dag     resid {n(s3@H@s3-H.conj().T):.1e}")
K3,K1,K2=s3/2,1j*s1/2,1j*s2/2
print(f" (1.5) [K1,K2]=-iK3        resid {n(K1@K2-K2@K1+1j*K3):.1e}")
print(f" (1.5) [K2,K3]=+iK1        resid {n(K2@K3-K3@K2-1j*K1):.1e}")
print(f" (1.5) [K3,K1]=+iK2        resid {n(K3@K1-K1@K3-1j*K2):.1e}")

print("\n=== Sec.5.3: the two vertices ===")
chi=0.4; Vt=np.cosh(chi)*I+np.sinh(chi)*s1          # w real   (5.10)
ph =0.4; Vs=np.cos(ph)*I+1j*np.sin(ph)*s1           # w imag   (5.10)
th=0.3; Dm=np.diag([np.exp(-1j*th),np.exp(1j*th)])
for nm,V in (("V_time (w real)",Vt),("V_space (w imag)",Vs)):
    print(f" {nm:18} |V^d s3 V - s3|={n(V.conj().T@s3@V-s3):.1e}   |V^d V - I|={n(V.conj().T@V-I):.1e}   |Theta-cov|={n(Th(V)-V):.1e}")
print(f" {'D (propagation)':18} {'':34}   |Theta-cov|={n(Th(Dm)-Dm):.1e}    <- (5.8)")

print("\n=== THE PIVOT: can a vertex be BOTH unitary AND Theta-covariant? ===")
print(" V = I + w*s1 :  unitarity  => w + w* = 0  (5.6)  => w imaginary")
print("                 Theta-cov  => w = w*     (5.9)  => w real")
print(" both at once => w = 0 : no folds, no coupling, H = E*s3 alone.")
for w in (0.3, 0.3j, 0.0):
    V=I+w*s1
    print(f"   w={str(w):>6}  |V^dV-I|={n(V.conj().T@V-I):.3f}   |Theta(V)-V|={n(Th(V)-V):.3f}")
print(" => the two conditions are mutually exclusive. Sec.5.3 does not ADD a")
print("    postulate to unitarity, it SWAPS unitarity out for Theta-covariance.")

print("\n=== Sec.7: do conditions (7.4) really force Gamma2 anti-Hermitian? ===")
for b in (1.0, 2.0, 0.5):
    G1=s3; G2=np.array([[0,b],[-1/b,0]],dtype=complex)
    ok=(n(G1@G1-I)<1e-12) and (n(G2@G2+I)<1e-12) and (n(G1@G2+G2@G1)<1e-12)
    anti=n(G2+G2.conj().T)<1e-12
    H2=E*G1+Dl*G2
    print(f"   b={b:<4} satisfies (7.4): {ok}   anti-Hermitian: {anti}   H^2=(E^2-D^2)I: {n(H2@H2-(E**2-Dl**2)*I)<1e-12}")
print(" => (7.4) admits Gamma2 that is NOT anti-Hermitian (any |b| != 1).")
print("    So eta = Gamma1 in (7.6) needs an extra Hermiticity choice; it is not")
print("    forced by the signature conditions alone.")

print("\n=== Sec.1.7: the positive metric eta+ at E=2D ===")
ev,R=np.linalg.eig(H); L=np.linalg.inv(R)
etap=L.conj().T@L                                    # biorthogonal eta+
etap/= np.abs(np.linalg.det(etap))**0.5
print(f"   H eigenvalues {np.sort(ev.real)}  (paper: +-sqrt(3)={np.sqrt(3):.4f})")
print(f"   intertwining |eta+ H - H^dag eta+| = {n(etap@H-H.conj().T@etap):.1e}")
print(f"   eta+ eigenvalues {np.sort(np.linalg.eigvalsh(etap))}   (paper: 2/3 and 2 -> product {2/3*2:.4f}, ratio 3)")
print(f"   ratio = {max(np.linalg.eigvalsh(etap))/min(np.linalg.eigvalsh(etap)):.4f}")
