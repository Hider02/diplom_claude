"""
Расчёт радиального входного устройства (ВхУ) перед ЦБК, как у аналога ВК-800С: воздух входит радиально снаружи
внутрь через кольцевое окно с защитной сеткой, канал плавно поворачивает поток на 90° к входу в колесо (конфузор),
передняя опора — в ступице ВхУ, канал пересекают полые радиальные стойки.

Исходные данные — из проекта (termo/termodynamiks.py, kompressor/cbk.py через t.komp): расход воздуха, диаметры
входа в колесо D1вт/D1, осевая скорость c1a, статические параметры на входе в колесо. Атмосфера — стандартная.
Назначено (отмечено «принято»): радиус окна, скорость в окне, радиус поворота покрывной стенки, коэффициенты
потерь сетки и поворота, число и толщина стоек.

Потери: сетка ζ_с ≈ 0,4 от напора в окне (Идельчик, сетки с живым сечением ~0,7); поворот на 90° с конфузорностью
и стойками ζ_пов ≈ 0,06 от напора на входе в колесо. Влияние σ_вх на мощность — пересчётом цикла
termo.termodynamiks.cycle() с другим σ_вх (в памяти, файлы проекта не меняются).

Запуск из корня проекта: python konstrukciya/vhodnoe_ustr.py
"""
import contextlib
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

K_RW = 2.2  # r_окна / r1 (периферия входа в колесо), принято по аналогу (2,1…2,4)
C_W = 50.0  # скорость в окне, м/с (принято: 30…50)
R_S = 20.0  # радиус поворота покрывной (задней) стенки, мм (принято)
X_C = -1.0  # осевая координата центра поворота (вход в колесо — x = 0), мм
Z_STRUT, T_STRUT = 6, 10.0  # число и толщина полых стоек, мм (принято)
ZETA_GRID = 0.4  # сетка, от напора в окне (принято)
ZETA_TURN = 0.06  # поворот 90° + стойки, от напора на входе в колесо (принято)


def _models():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    show = plt.show
    plt.show = lambda *a, **k: None
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            import termo.termodynamiks as t
    finally:
        plt.show = show
        plt.close("all")
    return t


def calc():
    t = _models()
    k = t.komp
    G, P0, T0, R = t.g_v, t.P0, t.T0, t.R_v
    rho0 = P0 / (R * T0)
    r1, r_vt = k["D1t"] / 2 * 1000, k["D1h"] / 2 * 1000
    F1 = math.pi * (r1 ** 2 - r_vt ** 2) * 1e-6
    c1a, rho1 = k["c1a"], k["rho1"]

    # Окно
    r_w = K_RW * r1
    F_w = G / (rho0 * C_W)
    b_w = F_w / (2 * math.pi * r_w / 1000) * 1000
    blk = 1 - Z_STRUT * T_STRUT / (2 * math.pi * r_w)  # загромождение стойками
    n_konf = F_w / F1

    # Канал: концентрический поворот (центр X_C, r_c), радиальная часть — от r_c до окна, ширина от h1 до b_w
    h1 = r1 - r_vt
    r_c = r1 + R_S
    R_H = R_S + h1
    x_r = X_C - R_S  # задняя (покрывная) стенка радиальной части
    x_f_c = X_C - R_H  # передняя стенка в начале поворота
    x_f_w = x_r - b_w  # передняя стенка у окна

    # Площадь по средней линии: радиальная часть (2π r w), поворот (по нормали между дугами), осевая (F1)
    s, F = [], []
    for r in np.linspace(r_w, r_c, 30):
        w = h1 + (b_w - h1) * (r - r_c) / (r_w - r_c)
        F.append(2 * math.pi * r * w * 1e-6)
        s.append(r_w - r)
    L_rad = r_w - r_c
    R_m = (R_S + R_H) / 2
    for a in np.linspace(0, math.pi / 2, 20)[1:]:  # от радиального направления к осевому
        r_m = r_c - R_m * math.sin(a)  # средний радиус в сечении поворота (сечение — отрезок по радиусу дуги)
        F.append(2 * math.pi * r_m * h1 * 1e-6)
        s.append(L_rad + R_m * a)
    F.append(F1)
    s.append(s[-1] + 1.0)
    monot = all(F[i + 1] <= F[i] * 1.001 for i in range(len(F) - 1))

    # Потери и σ_вх
    q_w = rho0 * C_W ** 2 / 2
    q_1 = rho1 * c1a ** 2 / 2
    dp = ZETA_GRID * q_w + ZETA_TURN * q_1
    sigma = 1 - dp / P0

    # Влияние на цикл (пересчёт в памяти)
    s0 = t.sigma_vh
    with contextlib.redirect_stdout(io.StringIO()):
        c0 = t.cycle()
        t.sigma_vh = sigma
        try:
            c1 = t.cycle()
        finally:
            t.sigma_vh = s0

    return dict(G=G, rho0=rho0, r1=r1, r_vt=r_vt, F1=F1, c1a=c1a, r_w=r_w, C_W=C_W, F_w=F_w, b_w=b_w, blk=blk,
                c_w_eff=C_W / blk, n_konf=n_konf, h1=h1, r_c=r_c, R_S=R_S, R_H=R_H, X_C=X_C, x_r=x_r, x_f_c=x_f_c,
                x_f_w=x_f_w, s=s, F=F, monot=monot, q_w=q_w, q_1=q_1, dp=dp, sigma=sigma, sigma_cycle=s0,
                N0=c0["N_ct"], N1=c1["N_ct"], G0=c0["g_v"], G1=c1["g_v"], z=Z_STRUT, t_strut=T_STRUT)


def report(r):
    hp = 735.499
    rows = [
        ("Расход воздуха G_в", f"{r['G']:.3f} кг/с (атмосфера: ρ = {r['rho0']:.3f} кг/м³)"),
        ("Вход в колесо: ⌀вт / ⌀пер, F1, c1a", f"{2 * r['r_vt']:.0f} / {2 * r['r1']:.0f} мм; {r['F1'] * 1e4:.1f} см²; "
                                                 f"{r['c1a']:.0f} м/с"),
        ("Окно: r_окна = 2,2·r1 (принято), ⌀", f"{r['r_w']:.0f} мм; ⌀{2 * r['r_w']:.0f} мм"),
        ("Скорость в окне (принято) / со стойками", f"{r['C_W']:.0f} / {r['c_w_eff']:.0f} м/с "
                                                    f"({r['z']} стоек по {r['t_strut']:.0f} мм)"),
        ("Площадь окна F_ок, ширина b", f"{r['F_w'] * 1e4:.0f} см²; {r['b_w']:.1f} мм"),
        ("Степень конфузорности F_ок/F1", f"{r['n_konf']:.2f}"),
        ("Поворот: R покрывной / R ступицы", f"{r['R_S']:.0f} / {r['R_H']:.1f} мм (концентрический, h = {r['h1']:.1f} мм)"),
        ("Площадь по средней линии убывает", "да ✓" if r["monot"] else "НЕТ ✗"),
        ("Потери: сетка + поворот/стойки", f"{r['dp']:.0f} Па (q_ок = {r['q_w']:.0f} Па, q1 = {r['q_1']:.0f} Па)"),
        ("σ_вх расчётное / принятое в цикле", f"{r['sigma']:.4f} / {r['sigma_cycle']:.3f}"),
        ("Расход G_в при расчётном σ_вх", f"{r['G0']:.3f} → {r['G1']:.3f} кг/с (G_пр = const)"),
        ("Мощность СТ по циклу", f"{r['N0'] / 1e3:.1f} → {r['N1'] / 1e3:.1f} кВт ({r['N0'] / hp:.0f} → {r['N1'] / hp:.0f} л.с., "
                                 f"{(r['N1'] / r['N0'] - 1) * 100:+.1f} %)"),
    ]
    w = max(len(a) for a, _ in rows)
    print("\n═══ ВХОДНОЕ УСТРОЙСТВО (радиальное) ═══")
    for a, b in rows:
        print(f"{a:<{w}} | {b}")


if __name__ == "__main__":
    report(calc())
