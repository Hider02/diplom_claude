k = 1.4
R = 287.0
cp = k * R / (k - 1)
def gdf_tau(l): return 1.0 - (k-1.0)/(k+1.0) * l**2
def gdf_pi(l): return (1.0 - (k-1.0)/(k+1.0) * l**2)**(k/(k-1.0))
def gdf_q(l):
    C = ((k+1.0)/2.0)**(1.0/(k-1.0))
    return C * (1.0 - (k-1.0)/(k+1.0) * l**2)**(1.0/(k-1.0)) * l

def solve_lambda(q_val, tol=1e-7):
    low, high = 0.001, 0.999
    for _ in range(60):
        mid = (low + high) / 2.0
        if gdf_q(mid) < q_val: low = mid
        else: high = mid
        if high - low < tol: break
    return (low + high) / 2.0