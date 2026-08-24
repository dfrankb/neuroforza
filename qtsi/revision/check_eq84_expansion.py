from mpmath import mp, mpf, cos, acos, sqrt, cosh, atanh, tanh
mp.dps = 50

def u2_exact(x, y):
    """u = Omega*eps from cos(u) = cosh(chi) cos(x), tanh(chi)=y."""
    chi = atanh(y)
    return acos(cosh(chi)*cos(x))**2

def u2_paper(x, y):
    # (8.4): (hbar Om)^2 = E^2 - D^2 + D^2 E^2/(3Ep^2) - 2 D^4/(3Ep^2)
    # divided by Ep^2 with x=E/Ep, y=D/Ep
    return x**2 - y**2 + x**2*y**2/3 - 2*y**4/3

def u2_leading(x, y):
    return x**2 - y**2

print("Scaling test: set x=2s, y=s and shrink s. Residual should be O(s^6) if (8.4) is right.")
print(f"{'s':>10} {'exact-paper':>16} {'ratio':>10} {'exact-leading':>16} {'ratio':>10}")
prev_p = prev_l = None
for k in range(1, 7):
    s = mpf(10)**(-k)
    x, y = 2*s, s
    rp = u2_exact(x,y) - u2_paper(x,y)
    rl = u2_exact(x,y) - u2_leading(x,y)
    rat_p = float(prev_p/rp) if prev_p else float('nan')
    rat_l = float(prev_l/rl) if prev_l else float('nan')
    print(f"{float(s):10.0e} {float(rp):16.4e} {rat_p:10.3e} {float(rl):16.4e} {rat_l:10.3e}")
    prev_p, prev_l = rp, rl
print()
print("ratio ~1e6 => residual is O(s^6) => the 4th-order terms of (8.4) are CORRECT")
print("ratio ~1e4 => residual is O(s^4) => leading term only, as expected")
