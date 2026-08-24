import numpy as np
s1=np.array([[0,1],[1,0]],dtype=complex);s3=np.array([[1,0],[0,-1]],dtype=complex)
eps=1e-3;N=10000;w=eps;chi=np.arctanh(w)
V=np.cosh(chi)*np.eye(2)+np.sinh(chi)*s1
def run(E,order):
    D=np.diag([np.exp(-1j*E*eps),np.exp(1j*E*eps)])
    T={'VD':V@D,'DV':D@V,'sym':np.diag([np.exp(-1j*E*eps/2),np.exp(1j*E*eps/2)])@V@np.diag([np.exp(-1j*E*eps/2),np.exp(1j*E*eps/2)])}[order]
    p=np.array([1,0],dtype=complex);mn=mx=1.0
    for _ in range(N):
        p=T@p;nn=(p.conj()@p).real;mn=min(mn,nn);mx=max(mx,nn)
    return mx/mn if E==2.0 else (p.conj()@p).real
print(f"{'order':>6} {'E=2D (osc range)':>18} {'E=D (poly)':>12} {'E=0.5D (exp)':>14}")
for o in ('VD','DV','sym'):
    print(f"{o:>6} {run(2.0,o):18.4g} {run(1.0,o):12.4g} {run(0.5,o):14.4g}")
print(f"{'paper':>6} {'1.66':>18} {'197':>12} {'1.87e7':>14}")
print()
print("continuum closed forms at tau*D/hbar=10, Psi0=(1,0):")
for E,lab in ((2.0,'E=2D'),(1.0,'E=D'),(0.5,'E=0.5D')):
    H=E*s3+1j*s1; ev=np.linalg.eigvals(H)
    from scipy.linalg import expm
    p=expm(-1j*H*10)@np.array([1,0],dtype=complex)
    print(f"   {lab:>7}: |psi|^2 at tau=10 -> {(p.conj()@p).real:.4g}")
