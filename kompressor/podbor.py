"""
Подбор π*_к и конструктивных параметров ЦБК при заданном приведённом расходе воздуха.

1. Для каждой π*_к конструктивные параметры ЦБК подбираются на максимум η_к при выполнении всех проверок
   kompressor/cbk.py (M_w1к ≤ 1,3; u2 ≤ 620 м/с; w2/w1к ≥ 0,6; b2 ≥ 5 мм; D1к/D2 ≤ 0,7; α2 = 60…75° от радиуса;
   M на входе в осевой ЛД ≤ 0,85; D_eq ≤ 2 и D ≤ 0,6 в венцах ЛД; c_к = 100…140 м/с).
   Коэффициенты моделей потерь не подбираются.
2. Для каждого варианта считается цикл (termo.termodynamiks.cycle) и мощность свободной турбины.
3. Лучший вариант — минимальный удельный расход топлива среди вариантов с N ≥ N_target_hp;
   его параметры нужно перенести в termo/termodynamiks.py (pik, KOMP).

Запуск из корня проекта: python kompressor/podbor.py
"""
import random
import termo.termodynamiks as t
from kompressor.cbk import compressor

PIK_VARIANTS = [6.0, 6.5, 7.0, 7.5, 8.0]

BOUNDS = dict(
    n_rpm=(35000, 50000),
    nu_hub=(0.30, 0.50),
    k_D1=(1.0, 1.10),
    beta2l=(45.0, 65.0),
    z_main=(11, 17),
    cr_ratio=(0.70, 1.00),
    b3_b2=(0.80, 1.00),
    Rb_b3=(1.5, 3.0),
    AR_bend=(1.0, 1.2),
    n_rows=(3, 4),
    sigma_ld=(1.4, 2.0),
    AR_ld=(1.0, 1.2),
)


def violation(r):
    """Мера нарушения ограничений: число замечаний расчёта ЦБК"""
    return len(r["warnings"])


def evaluate(pik, G, x):
    # D3/D2 не подбирается: высота радиального диффузора задана стыковкой с противоточной КС
    # (kamera_sgoraniya/podbor_diffuzor.py)
    x = dict(x, z_split=x["z_main"], D3_D2=t.KOMP["D3_D2"])
    try:
        r = compressor(pik, G, t.P1, t.T1, **x)
    except (ValueError, ZeroDivisionError, OverflowError, TypeError):
        return None, -1e9
    return r, r["eta_k"] - 0.1 * violation(r)


def random_point(rnd):
    return {name: (rnd.randint(lo, hi) if isinstance(lo, int) else rnd.uniform(lo, hi))
            for name, (lo, hi) in BOUNDS.items()}


def optimize_design(pik, G, n_random=300, n_refine=8, seed=1):
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
    return dict(best_x, z_split=best_x["z_main"], D3_D2=t.KOMP["D3_D2"]), best_r


if __name__ == "__main__":
    G = t.g_v_pr * t.sigma_vh
    print(f"Подбор при G_пр = {t.g_v_pr} кг/с (G = {G:.3f} кг/с), T_г* = {t.T_g} К, "
          f"требуется N ≥ {t.N_target_hp} л.с.\n")
    hdr = (f"{'π*к':<5} | {'η_к':<6} | {'η_РК':<6} | {'N, л.с.':<8} | {'C_e, кг/(л.с.·ч)':<16} | {'n, об/мин':<9} | "
           f"{'D2, мм':<7} | {'u2, м/с':<7} | {'b2, мм':<6} | {'M_w1к':<5} | {'венцов ЛД':<9} | замечания")
    print(hdr)
    print("-" * len(hdr))
    results = []
    for pik in PIK_VARIANTS:
        x, r = optimize_design(pik, G)
        c = t.cycle(pik=pik, komp_params=x)
        results.append(dict(pik=pik, x=x, komp=r, cycle=c))
        print(f"{pik:<5} | {r['eta_k']:<6.3f} | {r['eta_rk']:<6.3f} | {c['N_ct'] / t.HP:<8.1f} | "
              f"{c['SFC'] * t.HP / 1000:<16.4f} | {r['n_rpm']:<9.0f} | {r['D2'] * 1000:<7.1f} | {r['u2']:<7.1f} | "
              f"{r['b2'] * 1000:<6.2f} | {r['M_w1t']:<5.2f} | {x['n_rows']:<9} | {len(r['warnings']) or '✅'}")
    ok = [v for v in results if not v["komp"]["warnings"] and v["cycle"]["N_ct"] / t.HP >= t.N_target_hp]
    if not ok:
        print(f"\n⚠️ Ни один вариант не даёт {t.N_target_hp} л.с. без замечаний — повысьте T_г или расход")
        ok = results
    best = min(ok, key=lambda v: v["cycle"]["SFC"])
    print(f"\nЛучший вариант по удельному расходу топлива: π*к = {best['pik']}")
    print("Перенесите в termo/termodynamiks.py:")
    print(f"pik = {best['pik']}")
    x = {kk: (round(vv, 3) if isinstance(vv, float) else vv) for kk, vv in best["x"].items()}
    print(f"KOMP = dict({', '.join(f'{kk}={vv}' for kk, vv in x.items())})")
