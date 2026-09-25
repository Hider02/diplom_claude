"""Газодинамические функции для воздуха (k = 1,4)"""
import math

k = 1.4
R = 287.0
cp = k * R / (k - 1)
# Константа расхода s_в = sqrt(k·(2/(k+1))^((k+1)/(k-1))/R) = 0,0404 (кг·К/Дж)^0,5
s_v = math.sqrt(k * (2 / (k + 1)) ** ((k + 1) / (k - 1)) / R)


def a_cr(T_star):
    """Критическая скорость звука, м/с"""
    return math.sqrt(2 * k / (k + 1) * R * T_star)


def gdf_tau(lam):
    return 1.0 - (k - 1.0) / (k + 1.0) * lam ** 2


def gdf_pi(lam):
    return gdf_tau(lam) ** (k / (k - 1.0))


def gdf_eps(lam):
    return gdf_tau(lam) ** (1.0 / (k - 1.0))


def gdf_q(lam):
    return ((k + 1.0) / 2.0) ** (1.0 / (k - 1.0)) * lam * gdf_eps(lam)


def lambda_from_q(q_val, tol=1e-9):
    """Дозвуковое решение q(λ) = q_val (0 < λ < 1)"""
    if not 0 < q_val <= 1:
        raise ValueError(f"q(λ) = {q_val:.3f} вне диапазона (0; 1] — сечение не пропускает расход")
    low, high = 0.0, 1.0
    while high - low > tol:
        mid = (low + high) / 2.0
        if gdf_q(mid) < q_val:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0
