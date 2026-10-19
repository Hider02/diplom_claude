"""
Расчёт выходного устройства (ВУ) за свободной турбиной: кольцевой осевой диффузор → поворот в радиальный
канал → выхлопной (улиточный) патрубок с одним выходом вверх, как у аналога ВК-800С.
Вал СТ и его опоры проходят внутри кожуха — внутренней стенки диффузора (втулочный диаметр выхода РК СТ).

Исходные данные — из проекта: расход и параметры за РК СТ (free_turbine/gasdynamiks.py), давление за СТ,
принятое в термодинамике (termo/termodynamiks.py: P5 = P0/σ_вых·P8/P0). Назначаются (отмечены «принято»):
скорость на выходе патрубка, эквивалентный угол и длина кольцевого участка, коэффициенты потерь.

Потери (Идельчик; Холщевников): кольцевой диффузор — ζ = φ_расш·(1 − 1/n)² + ζ_тр, φ_расш = 3,2·tg^1,25(θ/2);
поворот на 90° с радиальным каналом и улиткой — ζ_пов ≈ 0,3 от скоростного напора за кольцевым участком;
закрутка за РК (c_2u) в патрубке не восстанавливается. Выходная скорость теряется в атмосферу.

Запуск из корня проекта: python konstrukciya/vyhodnoe_ustr.py
"""
import contextlib
import io
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# ─── Принятые величины ───
C_OUT = 80.0  # скорость на выходе патрубка, м/с (принято: 60…100 м/с)
THETA_EQ = 10.0  # эквивалентный угол раскрытия кольцевого участка, град (принято: ≤ 8…10°)
L_AX = 150.0  # длина кольцевого участка, мм (принято по компоновке: опоры СТ внутри кожуха)
N_RAD = 2.0  # степень расширения на входе в радиальный канал (до улитки) F_рад/F_вх (принято)
R_TURN = 40.0  # радиус поворота наружной стенки в радиальный канал, мм (принято)
A_FLANGE = 180.0  # осевой размер выходного фланца патрубка, мм (принято)
FLANGE_TILT = 10.0  # наклон плоскости выходного фланца к оси двигателя, град (как у аналога)
ZETA_FR = 0.02  # трение на кольцевом участке (принято)
ZETA_TURN = 0.3  # поворот + радиальный канал + улитка, от напора за кольцевым участком (принято: 0,25…0,4)


def _models():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    show = plt.show
    plt.show = lambda *a, **k: None
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            import termo.termodynamiks as t
            import free_turbine.gasdynamiks as ct
    finally:
        plt.show = show
        plt.close("all")
    return t, ct


def calc():
    t, ct = _models()
    R, k, cp = t.R_g, ct.k_g, ct.cp_prime
    G = t.g_gas
    T_s, p_s = ct.T_2_full, ct.P_2_full  # полные параметры за РК СТ
    c_a, c_u, rho_in = ct.c_2a, ct.c_2u, ct.rho_2
    d_h, D_in = ct.d_rk * 1000, ct.D_rk * 1000  # мм
    F_in = math.pi / 4 * (D_in ** 2 - d_h ** 2) * 1e-6

    # Кольцевой участок: втулка постоянна (кожух опор), наружная стенка раскрывается
    r_eq_in = math.sqrt(F_in / math.pi)
    r_eq_1 = r_eq_in + L_AX / 1000 * math.tan(math.radians(THETA_EQ / 2))
    F_1 = math.pi * r_eq_1 ** 2
    n_1 = F_1 / F_in
    D_1 = math.sqrt(4 * F_1 / math.pi * 1e6 + d_h ** 2)
    gamma_out = math.degrees(math.atan((D_1 - D_in) / 2 / L_AX))

    # Радиальный канал за поворотом
    r_rad = D_1 / 2 + R_TURN
    F_rad = N_RAD * F_in
    b_rad = F_rad / (2 * math.pi * r_rad / 1000) * 1000

    # Выход патрубка (статическое давление — атмосферное)
    T_out = T_s - C_OUT ** 2 / (2 * cp)
    rho_out = t.P0 / (R * T_out)
    F_out = G / (rho_out * C_OUT)
    c_flange = F_out * 1e6 / A_FLANGE
    D_out_eq = math.sqrt(4 * F_out / math.pi) * 1000

    # Потери
    q_a = rho_in * c_a ** 2 / 2
    q_u = rho_in * c_u ** 2 / 2
    phi = 3.2 * math.tan(math.radians(THETA_EQ / 2)) ** 1.25
    zeta_d = phi * (1 - 1 / n_1) ** 2 + ZETA_FR
    c_1 = c_a / n_1
    q_1 = rho_in * c_1 ** 2 / 2
    q_out = rho_out * C_OUT ** 2 / 2
    dp_d, dp_turn, dp_sw = zeta_d * q_a, ZETA_TURN * q_1, q_u
    dp_sum = dp_d + dp_turn + dp_sw
    p_out_full = t.P0 + q_out
    P5_req = p_out_full + dp_sum  # полное давление за СТ, нужное для выхода в атмосферу
    sigma_vu = p_out_full / P5_req  # σ ВУ без потери выходной скорости
    sigma_tot = t.P0 / P5_req

    # Влияние на мощность СТ (формула termo/termodynamiks.py с P5 = P5_req)
    def N_ct(P5):
        pit = t.P45 / P5
        L = t.cp_ct * t.T4 * (1 - pit ** (-(t.k_ct - 1) / t.k_ct)) * t.eta_ct * t.eta_mex_ct
        return L * t.g_gas, pit

    N0, pit0 = N_ct(t.P5)
    N1, pit1 = N_ct(P5_req)

    return dict(G=G, T_s=T_s, p_s=p_s, c_a=c_a, c_u=c_u, c_2=ct.c_2, lam_in=ct.c_2 / math.sqrt(2 * k / (k + 1) * R * T_s),
                d_h=d_h, D_in=D_in, F_in=F_in, L_ax=L_AX, theta=THETA_EQ, n_1=n_1, D_1=D_1, F_1=F_1, gamma_out=gamma_out,
                R_turn=R_TURN, r_rad=r_rad, F_rad=F_rad, b_rad=b_rad, n_rad=N_RAD,
                C_out=C_OUT, T_out=T_out, rho_out=rho_out, F_out=F_out, a_flange=A_FLANGE, c_flange=c_flange,
                D_out_eq=D_out_eq, tilt=FLANGE_TILT, n_tot=F_out / F_in,
                zeta_d=zeta_d, dp_d=dp_d, dp_turn=dp_turn, dp_sw=dp_sw, dp_sum=dp_sum, q_out=q_out,
                P5_cycle=t.P5, P5_req=P5_req, sigma_vu=sigma_vu, sigma_tot=sigma_tot, sigma_cycle=t.sigma_exh,
                pit_cycle=pit0, pit_real=pit1, N_cycle=N0, N_real=N1, P0=t.P0)


def report(r):
    hp = 735.499
    rows = [
        ("Расход газа G", f"{r['G']:.3f} кг/с"),
        ("Полные T*, p* за РК СТ", f"{r['T_s']:.1f} К; {r['p_s'] / 1e3:.2f} кПа"),
        ("Скорость за РК СТ c2 (c2a; c2u)", f"{r['c_2']:.1f} м/с ({r['c_a']:.1f}; {r['c_u']:.1f}), λ = {r['lam_in']:.3f}"),
        ("Вход ВУ: ⌀вт / ⌀нар, F_вх", f"{r['d_h']:.0f} / {r['D_in']:.0f} мм; {r['F_in'] * 1e4:.1f} см²"),
        ("Кольцевой участок: L, θ_экв (принято)", f"{r['L_ax']:.0f} мм; {r['theta']:.1f}°"),
        ("  выход: ⌀нар, n1 = F1/F_вх", f"{r['D_1']:.1f} мм; {r['n_1']:.3f}"),
        ("  угол наружной стенки", f"{r['gamma_out']:.1f}°"),
        ("Радиальный канал: r, ширина b, n", f"{r['r_rad']:.0f} мм; {r['b_rad']:.1f} мм; {r['n_rad']:.2f}"),
        ("Выход патрубка: c_вых (принято)", f"{r['C_out']:.0f} м/с, T = {r['T_out']:.1f} К, ρ = {r['rho_out']:.3f} кг/м³"),
        ("  F_вых, фланец a × b", f"{r['F_out'] * 1e4:.0f} см²; {r['a_flange']:.0f} × {r['c_flange']:.0f} мм "
                                   f"(экв. ⌀{r['D_out_eq']:.0f} мм), наклон {r['tilt']:.0f}°"),
        ("Степень расширения n = F_вых/F_вх", f"{r['n_tot']:.2f}"),
        ("Потери: диффузор (ζ), поворот+улитка, закрутка", f"{r['dp_d']:.0f} (ζ={r['zeta_d']:.3f}); {r['dp_turn']:.0f}; "
                                                          f"{r['dp_sw']:.0f} Па"),
        ("Сумма потерь Δp*, выходной напор", f"{r['dp_sum']:.0f} Па; {r['q_out']:.0f} Па"),
        ("Нужное p* за СТ / принятое в цикле", f"{r['P5_req'] / 1e3:.2f} / {r['P5_cycle'] / 1e3:.2f} кПа"),
        ("Располагаемое p* за РК СТ", f"{r['p_s'] / 1e3:.2f} кПа — "
                                       f"{'запас есть' if r['p_s'] >= r['P5_req'] else 'НЕ хватает'}"),
        ("σ ВУ (без вых. скорости) / с ней / в цикле", f"{r['sigma_vu']:.4f} / {r['sigma_tot']:.4f} / {r['sigma_cycle']:.3f}"),
        ("π_СТ: в цикле → с расчётным ВУ", f"{r['pit_cycle']:.4f} → {r['pit_real']:.4f}"),
        ("N_СТ: в цикле → с расчётным ВУ", f"{r['N_cycle'] / 1e3:.1f} → {r['N_real'] / 1e3:.1f} кВт "
                                            f"({r['N_cycle'] / hp:.0f} → {r['N_real'] / hp:.0f} л.с., "
                                            f"{(r['N_real'] / r['N_cycle'] - 1) * 100:+.2f} %)"),
    ]
    w = max(len(a) for a, _ in rows)
    print("\n═══ ВЫХОДНОЕ УСТРОЙСТВО (за СТ) ═══")
    for a, b in rows:
        print(f"{a:<{w}} | {b}")


if __name__ == "__main__":
    report(calc())
