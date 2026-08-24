import numpy as np
s1=np.array([[0,1],[1,0]],dtype=complex); s3=np.array([[1,0],[0,-1]],dtype=complex)
hbar=1.0; D=1.0; eps=1e-3; N=10000            # eps*D/hbar=1e-3, tau*D/hbar=10
w=D*eps/hbar; chi=np.arctanh(w)
V=np.cosh(chi)*np.eye(2)+np.sinh(chi)*s1

print("=== Sec.10 trace expansion:  |tr T| - 2  ~=  -(E^2-D^2) eps^2/hbar^2 ===")
for E in (2.0,0.5):
    th=E*eps/hbar; T=V@np.diag([np.exp(-1j*th),np.exp(1j*th)])
    print(f"  E={E}D : |trT|-2 = {abs(np.trace(T).real)-2: .4e}   predicted {-(E**2-D**2)*eps**2: .4e}   detT-1={abs(np.linalg.det(T)-1):.1e}")

print("\n=== Sec.10 evolution, Psi0=(1,0), tau*D/hbar=10 ===")
print(f"{'E':>6} {'regime':>10} {'norm growth':>14} {'paper':>10} {'max|dQ|':>11} {'paper':>9}")
papers={2.0:('elliptic','1.66','2.6e-12'),1.0:('parabolic','197','3.2e-12'),0.5:('hyperbolic','1.87e7','1.6e-8')}
lyap={}
for E in (2.0,1.0,0.5):
    th=E*eps/hbar; T=V@np.diag([np.exp(-1j*th),np.exp(1j*th)])
    psi=np.array([1,0],dtype=complex); Q0=(psi.conj()@s3@psi).real
    nmin=nmax=1.0; dQ=0.0
    for k in range(N):
        psi=T@psi
        n=(psi.conj()@psi).real; nmin=min(nmin,n); nmax=max(nmax,n)
        dQ=max(dQ,abs((psi.conj()@s3@psi).real-Q0))
    g=nmax/nmin if E==2.0 else n
    ev=np.linalg.eigvals(T); lyap[E]=np.log(max(abs(ev)))/eps
    print(f"{E:>5}D {papers[E][0]:>10} {g:14.4g} {papers[E][1]:>10} {dQ:11.2e} {papers[E][2]:>9}")

print(f"\n=== Lyapunov at E=0.5D ===\n  measured {lyap[0.5]:.6f}   predicted sqrt(D^2-E^2)/hbar = {np.sqrt(1-0.25):.6f}   (paper: 0.866026 vs 0.866025)")

print("\n=== Band edges of (8.3): cos(Om*eps)=cosh(chi)cos(E/Ep) ===")
print("  upper edge needs cos(E/Ep) = -sech(chi); E_upper in units of Ep:")
for y in (0.001,0.15,0.5,0.9):
    ch=1/np.sqrt(1-y**2)
    print(f"   D/Ep={y:<6} E_lower/Ep={np.arccos(np.sqrt(1-y**2)):.4f} (=D/Ep? {y:.4f})   E_upper/Ep={np.arccos(-np.sqrt(1-y**2)):.4f}")
print("\n  max(Om*eps) over the band:")
for y in (0.001,0.15,0.5,0.9):
    ch=1/np.sqrt(1-y**2); th=np.linspace(0,np.pi,2000001)
    r=ch*np.cos(th); ok=np.abs(r)<=1
    print(f"   D/Ep={y:<6} max(Om*eps)={np.arccos(r[ok]).max():.6f}   (paper: 3.141593 = pi)")
