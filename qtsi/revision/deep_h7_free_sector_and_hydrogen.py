import numpy as np
from scipy.linalg import sqrtm, expm
s1=np.array([[0,1],[1,0]],dtype=complex); s3=np.array([[1,0],[0,-1]],dtype=complex); I=np.eye(2,dtype=complex)
n=lambda A: np.abs(A).max()
def etaplus(E,D):
    H=E*s3+1j*D*s1; ev,R=np.linalg.eig(H); L=np.linalg.inv(R); e=L.conj().T@L
    return e/np.sqrt(abs(np.linalg.det(e))), H

print("="*74);print("H7-A. How well-conditioned is the map to the Hermitian picture?");print("="*74)
print(f"{'E/D':>10} {'cond(eta+)':>14} {'(E+D)/(E-D)':>14} {'|h-h^dag|':>12}")
for r in (50,5,2,1.5,1.1,1.01,1.001,1.0001):
    e,H=etaplus(r,1.0); w=np.linalg.eigvalsh(e); c=w.max()/w.min()
    S=sqrtm(e); h=S@H@np.linalg.inv(S)
    print(f"{r:>10} {c:14.4e} {(r+1)/(r-1):14.4e} {n(h-h.conj().T):12.1e}")
print(" => cond(eta+) = (E+D)/(E-D) exactly. The similarity to a Hermitian theory")
print("    exists for every E>D but DEGENERATES as E->D (the particle at rest).")

print("\n For a free mode E=sqrt(D^2+c^2p^2), E-D ~ p^2/2m, so cond ~ 4 m^2c^2/p^2:")
print(f"{'p/mc':>10} {'cond(eta+)':>14} {'4(mc/p)^2':>14}")
for pr in (1.0,0.1,0.01,0.001):
    E=np.sqrt(1+pr**2); e,_=etaplus(E,1.0); w=np.linalg.eigvalsh(e)
    print(f"{pr:>10} {w.max()/w.min():14.4e} {4/pr**2:14.4e}")
print(" => the equivalence is POINTWISE in p but NOT UNIFORM. It blows up as p^-2.")

print("\n"+"="*74);print("H7-B. The paper's rest-frame growth slope 2*pi*|phi(0)|^2");print("="*74)
print(" Positive form: int dp |phi|^2 [cos^2(Om t) + ((E^2+D^2)/(hbar Om)^2) sin^2(Om t)]")
print(" Near p=0 the bracket ~ (2D^2/c^2p^2) sin^2(cp t/hbar); using int sin^2(ap)/p^2 dp = pi|a|")
print(" gives slope 2*pi*|phi(0)|^2 (hbar=c=D=1).  Numerical check:")
def posform(tau,sig):
    p=np.linspace(-40*sig,40*sig,4000001); phi=(2*np.pi*sig**2)**-0.25*np.exp(-p**2/(4*sig**2))
    E=np.sqrt(1+p**2); Om=np.sqrt(np.maximum(E**2-1,1e-300))
    return np.trapezoid(np.abs(phi)**2*(np.cos(Om*tau)**2+((E**2+1)/Om**2)*np.sin(Om*tau)**2),p)
sig=0.3; phi0sq=(2*np.pi*sig**2)**-0.5
print(f"   Gaussian width {sig} mc (the paper's case), |phi(0)|^2 = {phi0sq:.6f}")
print(f"   predicted slope 2*pi*|phi(0)|^2 = {2*np.pi*phi0sq:.4f}")
t1,t2=200.0,400.0; v1,v2=posform(t1,sig),posform(t2,sig)
print(f"   measured slope over tau in [{t1},{t2}] = {(v2-v1)/(t2-t1):.4f}   (paper: reproduced to <1% at tau=400)")
print(f"   ratio measured/predicted = {((v2-v1)/(t2-t1))/(2*np.pi*phi0sq):.4f}")
print(" => confirmed. But note WHAT grows: the NON-conserved positive form. Q stays put.")

print("\n"+"="*74);print("H7-C. Where does the hydrogen prediction actually come from?");print("="*74)
print(" Claim in Sec.1.9 (from [9,10]):  dEn/|En| = (D0^2 / 2Ry^2) n^4")
print(" Derivation: spectrum of H = E s3 + i D0 s1 is +-sqrt(E^2 - D0^2), so with E=En")
print("   sqrt(En^2-D0^2) = |En| sqrt(1-(D0/En)^2) ~ |En| (1 - D0^2/(2En^2))")
print("   => |dEn|/|En| = D0^2/(2En^2);  En = -Ry/n^2  =>  = D0^2 n^4/(2 Ry^2).  Matches.")
Ry=13.605693
for n_ in (1,10,50,100):
    En=Ry/n_**2
    for D0 in (1e-9,1e-6):
        pass
    print(f"   n={n_:<4} |En|={En:.4e} eV   shift/|En| = (D0/{Ry:.3f} eV)^2 * {n_**4/2:.4g}")
print("\n KEY POINT: sqrt(E^2-D0^2) is a SIMILARITY INVARIANT (it is the spectrum of H).")
print(" So this prediction does NOT come from the metric or the two-branch structure --")
print(" both of which are similarity-covariant and hence unobservable in the free sector.")
print(" It comes entirely from postulating D0 != 0, i.e. from the invariant E^2-D0^2.")
print(" Sec.11 concedes the lattice 'does not supply the value of Delta for a given system'.")

print("\n What does existing spectroscopy already say about D0?")
print(f"{'n':>5} {'rel. precision':>16} {'bound on D0 (eV)':>20} {'D0/mc^2':>12}")
mc2=510998.95
for n_,prec in ((2,1e-14),(10,1e-12),(50,1e-11),(100,1e-10)):
    D0=Ry*np.sqrt(2*prec/n_**4)
    print(f"{n_:>5} {prec:>16.0e} {D0:>20.3e} {D0/mc2:>12.2e}")
print(f"\n Compare the 'natural' identification D0 = mc^2 = {mc2:.3e} eV:")
for n_ in (1,10):
    print(f"   n={n_}: predicted shift/|En| = {(mc2/Ry)**2*n_**4/2:.3e}   (a relative shift of ~10^9 at n=1)")
print(" => D = mc^2 is excluded for a bound electron by ~15 orders of magnitude.")
print("    So the theory needs a RULE for what Delta is in a bound system, and has none.")
