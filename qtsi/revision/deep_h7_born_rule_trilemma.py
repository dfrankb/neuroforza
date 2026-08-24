import numpy as np
from scipy.linalg import sqrtm, expm
s1=np.array([[0,1],[1,0]],dtype=complex); s3=np.array([[1,0],[0,-1]],dtype=complex); I=np.eye(2,dtype=complex)
def etaplus(E,D):
    H=E*s3+1j*D*s1; ev,R=np.linalg.eig(H); L=np.linalg.inv(R); e=L.conj().T@L
    return e/np.sqrt(abs(np.linalg.det(e))), H

print("="*74);print("H7-D. What is the Born rule in QTSI? Three candidates, three failures.");print("="*74)
E,D=2.0,1.0; e,H=etaplus(E,D)
states={"(1,0)":np.array([1,0],dtype=complex),
        "(1,1)/sqrt2":np.array([1,1],dtype=complex)/np.sqrt(2),
        "(1,1.3)/norm":np.array([1,1.3],dtype=complex)/np.linalg.norm([1,1.3])}
print("\n (i) <psi|psi>  -- positive, but NOT conserved:")
for lbl,p in states.items():
    v=[float((expm(-1j*H*t)@p).conj()@(expm(-1j*H*t)@p)) for t in (0,1,2,3)]
    print(f"   {lbl:14} <psi|psi> at tau=0,1,2,3: {['%.4f'%x for x in v]}")
print("\n (ii) Q = <psi|s3|psi>  -- conserved, but INDEFINITE (zero and negative allowed):")
for lbl,p in states.items():
    q0=float(p.conj()@s3@p); q3=float((expm(-1j*H*3)@p).conj()@s3@(expm(-1j*H*3)@p))
    print(f"   {lbl:14} Q(0)={q0:+.4f}  Q(3)={q3:+.4f}  conserved={abs(q0-q3)<1e-10}   usable as probability={q0>0}")
print("   -> a perfectly normalisable state can have Q = 0 or Q < 0. Q is a charge, as the paper says.")
print("\n (iii) <psi|eta+|psi>  -- conserved AND positive, but:")
for lbl,p in states.items():
    n0=float(p.conj()@e@p); n3=float((expm(-1j*H*3)@p).conj()@e@(expm(-1j*H*3)@p))
    print(f"   {lbl:14} <eta+>(0)={n0:.4f}  <eta+>(3)={n3:.4f}  conserved={abs(n0-n3)<1e-10}")
print("   -> works, but exists only for E>D and degenerates as (E+D)/(E-D) at rest.")
print("   -> AND it makes the free sector exactly ordinary QM. Check transition probs:")
S=sqrtm(e); h=S@H@np.linalg.inv(S)
psi=np.array([1,0],dtype=complex)
for t in (0.5,1.0,2.0):
    U=expm(-1j*H*t); pq=abs((psi.conj()@e@(U@psi)))**2/ (float(psi.conj()@e@psi)**2)
    Uh=expm(-1j*h*t); chi=S@psi; chi/=np.linalg.norm(chi)
    ph=abs(chi.conj()@(Uh@chi))**2
    print(f"      tau={t}: |<psi|eta+|U psi>|^2 = {pq:.10f}   Hermitian picture = {ph:.10f}   equal={abs(pq-ph)<1e-9}")
print("\n => TRILEMMA: no bilinear on the two-branch space is simultaneously")
print("    (a) positive, (b) conserved by the free dynamics, and (c) capable of")
print("    distinguishing QTSI from ordinary QM. The paper never states a Born rule.")
