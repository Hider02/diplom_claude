"""
Дополнительные расчёты конструкции (по схеме konstrukciya/konstr.py, аналог ВК-800С):
  1. валы газогенератора и СТ на кручение;
  2. подшипники по параметру d·n;
  3. болтовые соединения: колесо ЦБК — диск ТВД, диск СТ — валы (фрикционные, по крутящему моменту),
     фланцы корпуса КС (по внутреннему давлению), стяжка газогенератора;
  4. первые критические частоты роторов (жёсткие опоры, сосредоточенные массы дисков, масса вала не учтена,
     ступицы колеса/дисков, стянутые с валом, работают на изгиб заодно с ним — принято);
  5. радиальное входное устройство (konstrukciya/vhodnoe_ustr.py), выходное устройство
     (konstrukciya/vyhodnoe_ustr.py) и итоговые удельные параметры двигателя.
Принятые величины (материалы, коэффициенты) — в константах ниже.

Запуск из корня проекта: python konstrukciya/raschet.py
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from konstrukciya import konstr  # noqa: E402
from konstrukciya.vhodnoe_ustr import report as report_vh  # noqa: E402
from konstrukciya.vyhodnoe_ustr import _models, report as report_vu  # noqa: E402

HP = 735.499
E_SHAFT = 2.1e11  # модуль упругости стали вала, Па
TAU_T = 0.58 * 835e6  # предел текучести на сдвиг, 40ХН2МА (σ_т = 835 МПа), Па
N_MIN_SHAFT = 1.5  # требуемый запас по текучести вала
DN_LIM = {"ball": 2.2e6, "roller": 2.0e6}  # допустимое d·n, мм·об/мин (авиационные, струйная смазка)
SIGMA_T_BOLT = 900e6  # предел текучести болтов (30ХГСА, термообработка), Па
F_FRICTION, K_SLIP, K_TIGHT = 0.15, 1.5, 1.3  # трение в стыке, запас от сдвига, учёт кручения при затяжке
A_S = {6: 20.1e-6, 8: 36.6e-6, 10: 58.0e-6}  # площадь по резьбе, м²
TIE_PRELOAD = 25e3  # усилие затяжки стяжки газогенератора, Н (принято)


def shaft_torsion(name, N, n_rpm, D, d):
    w = n_rpm * math.pi / 30
    M = N / w
    Wp = math.pi * (D / 1000) ** 3 / 16 * (1 - (d / D) ** 4)
    tau = M / Wp
    return dict(name=name, N=N, n=n_rpm, M=M, D=D, d=d, tau=tau, n_t=TAU_T / tau)


def friction_joint(name, M, z, r_mm, d_bolt):
    F0 = K_SLIP * M / (F_FRICTION * z * r_mm / 1000)
    sigma = K_TIGHT * F0 / A_S[int(d_bolt)]
    return dict(name=name, M=M, z=z, r=r_mm, d=d_bolt, F0=F0, sigma=sigma, n=SIGMA_T_BOLT / sigma)


def flange_joint(name, dp, D_mm, z, d_bolt, k_tight=1.5):
    F = dp * math.pi / 4 * (D_mm / 1000) ** 2
    Fb = k_tight * F / z
    sigma = K_TIGHT * Fb / A_S[int(d_bolt)]
    return dict(name=name, F=F, z=z, d=d_bolt, Fb=Fb, sigma=sigma, n=SIGMA_T_BOLT / sigma,
                pitch=math.pi * (D_mm + 22) / z)


def EI_steps(steps, D, d):
    """Изгибная жёсткость по длине, Н·м²: базовый вал ⌀D/⌀d (мм) и участки [(x0, x1, D_экв, d_экв)] ступиц
    и цапф, стянутых с валом (считаются работающими на изгиб заодно с валом — принято)."""
    I = lambda Do, di: math.pi * ((Do / 1000) ** 4 - (di / 1000) ** 4) / 64

    def f(x):
        for x0, x1, Do, di in steps:
            if x0 <= x <= x1:
                return E_SHAFT * I(Do, di)
        return E_SHAFT * I(D, d)
    return f


def crit_beam(xs, ms, xa, xb, EI):
    """Первая критическая частота (об/мин) вала на двух шарнирных опорах xa, xb (мм) с сосредоточенными массами
    ms (кг) в точках xs (в т. ч. на консоли). Коэффициенты влияния — интеграл Мора по сетке, EI(x) переменная."""
    L = xb - xa
    lo, hi = min(xs + [xa]), max(xs + [xb])
    x = np.linspace(lo, hi, 4001)
    dx = (x[1] - x[0]) / 1000
    ei = np.array([EI(v) for v in x])

    def moment(p):  # эпюра от единичной силы в точке p (м·Н на 1 Н), реакции опор
        forces = [(xa, (xb - p) / L), (xb, (p - xa) / L), (p, -1.0)]
        return sum(F * np.clip(x - s, 0, None) / 1000 for s, F in forces)

    Ms = [moment(p) for p in xs]
    A = np.array([[np.sum(Mi * Mj / ei) * dx for Mj in Ms] for Mi in Ms])
    lam = np.linalg.eigvals(A @ np.diag(ms)).real.max()
    return 30 / math.pi / math.sqrt(lam)


def main():
    t, ct = _models()
    P = konstr.build()
    vu = P["vu"]
    GG, CT = konstr.GG, konstr.CT
    n_gg, n_ct = t.KOMP["n_rpm"], ct.n_ct
    N_tvd = t.L_tvd_actual * t.g_gas
    N_ct_max = max(ct.N_st, 800 * HP)

    report_vh(P["vh"])
    report_vu(vu)

    # 1. Кручение
    sh = [shaft_torsion("Вал газогенератора (ТВД → ЦБК)", N_tvd, n_gg, 2 * GG["r_sh"], 2 * GG["r_bore"]),
          shaft_torsion("Вал СТ (на выходе к редуктору)", N_ct_max, n_ct, 2 * CT["r_sh"], 2 * CT["r_bore"])]
    print("\n═══ ВАЛЫ НА КРУЧЕНИЕ (40ХН2МА, τ_т = {:.0f} МПа) ═══".format(TAU_T / 1e6))
    print(f"{'Вал':<34} | {'N, кВт':>7} | {'n, об/мин':>9} | {'M, Н·м':>7} | {'⌀D/⌀d, мм':>10} | {'τ, МПа':>7} | {'n_т':>5}")
    for s in sh:
        ok = "✓" if s["n_t"] >= N_MIN_SHAFT else "✗"
        print(f"{s['name']:<34} | {s['N'] / 1e3:>7.1f} | {s['n']:>9.0f} | {s['M']:>7.1f} | "
              f"{s['D']:>4.0f}/{s['d']:<5.0f} | {s['tau'] / 1e6:>7.1f} | {s['n_t']:>4.1f} {ok}")

    # 2. Подшипники
    print("\n═══ ПОДШИПНИКИ (d·n) ═══")
    print(f"{'Опора':<10} | {'тип':<8} | {'d×D×B, мм':<12} | {'n, об/мин':>9} | {'d·n, 10⁶':>8} | {'допуск':>6}")
    for b in P["bearings"]:
        n = n_gg if b["name"] in ("Опора 1", "Опора 2") else n_ct
        dn = b["d"] * n
        kind = "шарик." if b["kind"] == "ball" else "ролик."
        print(f"{b['name']:<10} | {kind:<8} | {b['d']:.0f}×{b['D']:.0f}×{b['B']:.0f}{'':<3} | {n:>9.0f} | "
              f"{dn / 1e6:>8.2f} | {DN_LIM[b['kind']] / 1e6:>4.1f} {'✓' if dn <= DN_LIM[b['kind']] else '✗'}")

    # 3. Болты
    M_k = t.N_k / (n_gg * math.pi / 30)
    M_ct = N_ct_max / (n_ct * math.pi / 30)
    bj = [friction_joint("Колесо ЦБК — диск ТВД", M_k, GG["bolt_n"], GG["bolt_r"], GG["bolt_d"]),
          friction_joint("Диск СТ — валы СТ", M_ct, CT["bolt_n"], CT["bolt_r"], CT["bolt_d"])]
    fj = [flange_joint("Корпус ЦБК — корпус КС", t.P2 - t.P0, 2 * P["geo"]["r_out_tip"], 24, 6),
          flange_joint("Корпус СТ — выходное устройство", vu["p_s"] - t.P0, vu["D_in"], 16, 6)]
    print(f"\n═══ БОЛТОВЫЕ СОЕДИНЕНИЯ (σ_т болтов = {SIGMA_T_BOLT / 1e6:.0f} МПа) ═══")
    for j in bj:
        print(f"{j['name']:<34}: M = {j['M']:.0f} Н·м, {j['z']}×М{j['d']:.0f} на r = {j['r']:.0f} мм; затяжка "
              f"F₀ = {j['F0'] / 1e3:.2f} кН (f = {F_FRICTION}, запас {K_SLIP}); σ = {j['sigma'] / 1e6:.0f} МПа, "
              f"n = {j['n']:.1f}")
    for j in fj:
        print(f"{j['name']:<34}: F = {j['F'] / 1e3:.1f} кН, {j['z']}×М{j['d']:.0f}, на болт {j['Fb'] / 1e3:.2f} кН; "
              f"σ = {j['sigma'] / 1e6:.0f} МПа, n = {j['n']:.1f}; шаг {j['pitch']:.0f} мм")
    A_tie = math.pi / 4 * ((2 * GG["r_tie"]) ** 2 - (2 * GG["r_tie_in"]) ** 2) * 1e-6
    print(f"{'Стяжка газогенератора':<34}: труба ⌀{2 * GG['r_tie']:.0f}/⌀{2 * GG['r_tie_in']:.0f}, затяжка "
          f"{TIE_PRELOAD / 1e3:.0f} кН (принято), σ = {K_TIGHT * TIE_PRELOAD / A_tie / 1e6:.0f} МПа, "
          f"n = {SIGMA_T_BOLT / (K_TIGHT * TIE_PRELOAD / A_tie):.1f}")

    # 4. Критические частоты
    m_w, x_w = P["masses"]["ЦБК"]
    m_t, x_t = P["masses"]["ТВД"]
    xa, xb = P["masses"]["опоры ГГ"]
    g = P["geo"]
    # жёсткие участки ГГ: ступица и задняя цапфа колеса ЦБК (⌀92) и ступица диска ТВД (⌀88), внутри — вал ⌀36
    ei_gg = EI_steps([(0.0, 118.0, 92.0, 36.0), (g["x_tvd"] - 13, g["x_tvd"] + 13, 88.0, 36.0)],
                     2 * GG["r_sh"], 2 * GG["r_bore"])
    n_cr_gg = crit_beam([x_w, x_t], [m_w, m_t], xa, xb, ei_gg)
    m_c, x_c = P["masses"]["СТ"]
    x3, x4 = P["masses"]["опоры СТ"]
    # СТ: сплошная ступица диска (⌀92) и фланец вала (⌀88/28)
    ei_ct = EI_steps([(x_c - 15, x_c + 15, 92.0, 0.0), (x_c + 15, x_c + 21, 88.0, 2 * CT["r_bore"])],
                     2 * CT["r_sh"], 2 * CT["r_bore"])
    n_cr_ct = crit_beam([x_c], [m_c], x3, x4, ei_ct)
    print("\n═══ КРИТИЧЕСКИЕ ЧАСТОТЫ (жёсткие опоры, массы дисков сосредоточены) ═══")
    for name, ncr, n, extra in (
            ("Ротор ГГ (ЦБК + ТВД, вал на 2 опорах)", n_cr_gg, n_gg,
             f"m_ЦБК = {m_w:.2f} кг, m_ТВД = {m_t:.2f} кг, пролёт {xb - xa:.0f} мм"),
            ("Ротор СТ (консольный диск)", n_cr_ct, n_ct,
             f"m_СТ = {m_c:.2f} кг, вылет {x3 - x_c:.0f} мм, пролёт {x4 - x3:.0f} мм")):
        margin = ncr / n - 1
        verdict = "докритический, запас ✓" if margin >= 0.2 else (
            "надкритический (рабочая частота выше критической) — нужны упругодемпферные опоры" if margin < -0.2
            else "✗ критическая частота близко к рабочей")
        print(f"{name:<38}: n_кр = {ncr:.0f} об/мин, n_раб = {n:.0f}, запас {margin * 100:+.0f} % — {verdict}")
        print(f"{'':<38}  ({extra})")

    # 5. Итоговые удельные параметры
    N = vu["N_real"]
    Ce = t.W_fuel * 3600 / (N / 1e3)
    print("\n═══ ИТОГОВЫЕ ПАРАМЕТРЫ ДВИГАТЕЛЯ (с расчётным выходным устройством) ═══")
    rows = [("Мощность на валу СТ N_e", f"{N / 1e3:.1f} кВт = {N / HP:.0f} л.с. (задано 800 л.с.)"),
            ("Расход воздуха G_в", f"{t.g_v:.3f} кг/с"),
            ("Удельная мощность N_уд = N_e/G_в", f"{N / 1e3 / t.g_v:.1f} кВт·с/кг = {N / HP / t.g_v:.0f} л.с.·с/кг"),
            ("Расход топлива G_т", f"{t.W_fuel_hour:.1f} кг/ч"),
            ("Удельный расход топлива C_e", f"{Ce:.3f} кг/(кВт·ч) = {Ce * HP / 1e3:.3f} кг/(л.с.·ч)"),
            ("Эффективный КПД η_e", f"{N / (t.W_fuel * t.H_u):.3f}"),
            ("π_к*, T_г*", f"{t.pik:.1f}; {t.T_g:.0f} К"),
            ("Температура выхлопа", f"{vu['T_out']:.0f} К ({vu['T_out'] - 273.15:.0f} °C)")]
    w = max(len(a) for a, _ in rows)
    for a, b in rows:
        print(f"{a:<{w}} | {b}")
    return dict(shafts=sh, joints=bj + fj, n_cr_gg=n_cr_gg, n_cr_ct=n_cr_ct)


if __name__ == "__main__":
    main()
