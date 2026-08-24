import numpy as np
from scipy.linalg import sqrtm
I=np.eye(2,dtype=complex);s1=np.array([[0,1],[1,0]],dtype=complex);s3=np.array([[1,0],[0,-1]],dtype=complex)
n=lambda A: np.abs(A).max()

print("=== Sec.1.7 + 1.8 together: is the FREE sector just ordinary QM in disguise? ===")
print(" Sec.1.8: free modes have E >= D always.  Sec.1.7: for E > D a positive eta+ exists.")
print(" So for every free mode, h = eta+^(1/2) H eta+^(-1/2) should be HERMITIAN:")
for E in (1.001,1.5,2.0,5.0,50.0):
    H=E*s3+1j*1.0*s1
    ev,R=np.linalg.eig(H); L=np.linalg.inv(R); etap=L.conj().T@L
    S=sqrtm(etap); h=S@H@np.linalg.inv(S)
    print(f"   E={E:>6}D  |h - h^dag| = {n(h-h.conj().T):.2e}   spectrum {np.sort(np.linalg.eigvals(h).real)}")
print(" => yes. Every free mode is a similarity transform of a Hermitian 2-level system.")
print("    No free-particle observable can separate QTSI from ordinary QM.")

print("\n=== physical numbers quoted in Sec.8 ===")
hbar=1.054571817e-34;c=2.99792458e8;G=6.67430e-11;me=9.1093837015e-31
lp=np.sqrt(hbar*G/c**3);tp=lp/c;mp=np.sqrt(hbar*c/G);Ep=mp*c**2
tauC=hbar/(me*c**2)
print(f"   l_p = {lp:.4e} m      (paper 1.6e-35)")
print(f"   t_p = {tp:.4e} s      (paper 5.391e-44)")
print(f"   m_p = {mp:.4e} kg     (paper 2.2e-8)")
print(f"   E_p = {Ep/1.602176634e-10:.4e} GeV (paper 1.2e19)")
print(f"   tau_C(e-) = {tauC:.4e} s  (paper 1.288e-21)")
print(f"   w(e-) = t_p/tau_C = {tp/tauC:.4e}   = m_e/m_p = {me/mp:.4e}   (paper 4.185e-23)")
print(f"   (D/Ep)^2/3 for electron = {(me/mp)**2/3:.4e}   (paper 'about 1e-45')")
lC=hbar/(me*c); rs=2*G*me/c**2
print(f"   lambda_C * r_s = {lC*rs:.4e}   2 l_p^2 = {2*lp**2:.4e}   ratio {lC*rs/(2*lp**2):.6f}")

print("\n=== where is the UPPER band edge, i.e. the Sec.8 'collapse' point? ===")
for y in (me/mp, 1e-3, 0.15):
    Eup=np.arccos(-np.sqrt(1-y**2))
    print(f"   D/Ep={y:.3e}  ->  E_upper = {Eup:.4f} Ep = {Eup*Ep/1.602176634e-10:.3e} GeV")
print(" (LHC reaches ~1.4e4 GeV.)")
