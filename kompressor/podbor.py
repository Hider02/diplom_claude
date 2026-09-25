"""
Подбор π*_к, расхода воздуха и конструктивных параметров ЦБК под заданную мощность двигателя.

1. Для каждой π*_к конструктивные параметры ЦБК подбираются в пределах рекомендаций книги так, чтобы получить
   максимальный η_к при выполнении ограничений (M_w1 ≤ 0,9; M3 ≤ 0,95; u2 ≤ 550 м/с; b2 ≥ 5 мм; D1/D2 = 0,55…0,70;
   c̄1a = 0,25…0,35; z_ЛД = 9…25; c_к < c4; η_к < η_РК). Коэффициенты потерь (ξ1, ξ2, β, n_БЛД, ξ_ЛД, ξ_ВУ)
   не подбираются — это опытные данные книги, берутся средние значения.
2. Расход воздуха подбирается так, чтобы мощность свободной турбины была N_target_hp·(1 + ZAPAS).
3. Выводится таблица вариантов; параметры лучшего варианта (минимальный удельный расход топлива)
   нужно перенести в termo/termodynamiks.py (pik, g_v_pr, KOMP).

Запуск из корня проекта: python kompressor/podbor.py
"""
import random
import termo.termodynamiks as t
from kompressor.cbk import compressor

PIK_VARIANTS = [3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0]
ZAPAS = 0.05  # Запас по мощности к N_target_hp
B2_MIN = 5e-3  # Минимальная ширина лопатки на выходе из РК, м

# Подбираемые параметры и их пределы по книге
BOUNDS = dict(
    G_bar=(0.5, 0.7),  # с. 151, п. 2
    d_vt=(0.5, 0.7),  # с. 152, п. 7
    alpha1=(60.0, 90.0),  # с. 152, п. 8 (50…90°)
    beta2l=(50.0, 90.0),  # с. 136
    D3_D2=(1.05, 1.15),  # с. 159
    D4_D3=(1.25, 1.35),  # с. 161
    d_alpha4=(12.0, 18.0),  # с. 161
    z_vna=(19, 25),  # с. 156
)


def violation(r):
    """Суммарное нарушение ограничений (0 — все выполнены)"""
    v = 0.0
    v += max(0.0, r["M_w1"] - 0.9) * 10
    v += max(0.0, r["M3"] - 0.95) * 10
    v += max(0.0, r["u2"] - 550) / 50
    v += max(0.0, B2_MIN - r["b2"]) * 1000
    v += max(0.0, 0.55 - r["D1_D2"]) * 10 + max(0.0, r["D1_D2"] - 0.70) * 10
    v += max(0.0, 0.25 - r["c1a_bar"]) * 10 + max(0.0, r["c1a_bar"] - 0.35) * 10
    v += max(0.0, 9 - r["z_ld_raw"]) / 5 + max(0.0, r["z_ld_raw"] - 25) / 5
    v += max(0.0, r["c_k"] - r["c4"]) / 20
    v += max(0.0, r["eta_k"] - r["eta_rk"]) * 10
    return v


def evaluate(pik, G, x):
    try:
        r = compressor(pik, G, t.P1, t.T1, **x)
    except (ValueError, ZeroDivisionError, OverflowError, TypeError):
        return None, -1e9
    return r, r["eta_k"] - violation(r)


def random_point(rnd):
    x = {}
    for name, (lo, hi) in BOUNDS.items():
        x[name] = rnd.randint(lo, hi) if isinstance(lo, int) else rnd.uniform(lo, hi)
    return x


def optimize_design(pik, G, n_random=400, n_refine=6, seed=1):
    """Случайный поиск + покоординатное уточнение"""
    rnd = random.Random(seed)
    best_x, best_f, best_r = None, -1e18, None
    for _ in range(n_random):
        x = random_point(rnd)
        r, f = evaluate(pik, G, x)
        if f > best_f:
            best_x, best_f, best_r = x, f, r
    for _ in range(n_refine):
        improved = False
        for name, (lo, hi) in BOUNDS.items():
            step = 1 if isinstance(lo, int) else (hi - lo) / 20
            for sign in (-1, 1):
                x = dict(best_x)
                x[name] = min(max(x[name] + sign * step, lo), hi)
                r, f = evaluate(pik, G, x)
                if f > best_f + 1e-7:
                    best_x, best_f, best_r, improved = x, f, r, True
        if not improved:
            break
    return best_x, best_r


def podbor_variant(pik):
    """Подбор конструкции ЦБК и расхода под заданную мощность при данной π*_к"""
    N_target = t.N_target_hp * (1 + ZAPAS) * t.HP
    g_v_pr = t.g_v_pr
    x, r = None, None
    for _ in range(4):
        G = g_v_pr * t.sigma_vh
        x, r = optimize_design(pik, G)
        c = t.cycle(pik=pik, g_v_pr=g_v_pr, komp_params=x)
        g_v_pr_new = g_v_pr * N_target / c["N_ct"]  # мощность пропорциональна расходу
        if abs(g_v_pr_new - g_v_pr) / g_v_pr < 0.002:
            break
        g_v_pr = g_v_pr_new
    c = t.cycle(pik=pik, g_v_pr=g_v_pr, komp_params=x)
    return dict(pik=pik, g_v_pr=g_v_pr, x=x, komp=c["komp"], cycle=c)


if __name__ == "__main__":
    print(f"Подбор под N = {t.N_target_hp} л.с. + {ZAPAS * 100:.0f}% запаса, T_г* = {t.T_g} К\n")
    results = []
    hdr = (f"{'π*к':<5} | {'G_пр, кг/с':<10} | {'η_к':<6} | {'η_РК':<6} | {'N, л.с.':<8} | {'C_e, кг/(л.с.·ч)':<16} | "
           f"{'n, об/мин':<9} | {'D2, мм':<7} | {'u2, м/с':<7} | {'b2, мм':<6} | {'M_w1':<5} | {'M3':<5} | ограничения")
    print(hdr)
    print("-" * len(hdr))
    for pik in PIK_VARIANTS:
        v = podbor_variant(pik)
        r, c = v["komp"], v["cycle"]
        ok = violation(r) < 1e-9
        v["ok"] = ok
        results.append(v)
        print(f"{pik:<5} | {v['g_v_pr']:<10.3f} | {r['eta_k']:<6.3f} | {r['eta_rk']:<6.3f} | {c['N_ct'] / t.HP:<8.1f} | "
              f"{c['SFC'] * t.HP / 1000:<16.4f} | {r['n_rpm']:<9.0f} | {r['D2'] * 1000:<7.1f} | {r['u2']:<7.1f} | "
              f"{r['b2'] * 1000:<6.2f} | {r['M_w1']:<5.2f} | {r['M3']:<5.2f} | {'✅' if ok else '❌ нарушены'}")

    feasible = [v for v in results if v["ok"]] or results
    best = min(feasible, key=lambda v: v["cycle"]["SFC"])
    print(f"\nЛучший вариант по удельному расходу топлива: π*к = {best['pik']}, G_пр = {best['g_v_pr']:.3f} кг/с")
    print("Перенесите в termo/termodynamiks.py:")
    print(f"g_v_pr = {best['g_v_pr']:.3f}")
    print(f"pik = {best['pik']}")
    x = {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in best["x"].items()}
    print(f"KOMP = dict({', '.join(f'{kk}={vv}' for kk, vv in x.items())})")
