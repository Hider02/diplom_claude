"""
Расчёт противоточной камеры сгорания, увязанный с уже выполненными расчётами ЦБК и ТВД:
  - вход в КС — сечение за осевым ЛД компрессора (kompressor/cbk.py, через termo/termodynamiks.py):
    втулочный D4m − h5 и наружный D4m + h5 диаметры; радиальный диффузор ЦБК поднят так, что ЛД выходит
    прямо в наружный канал КС (D3/D2 — kamera_sgoraniya/podbor_diffuzor.py);
  - выход из КС — сечение на входе в СА ТВД (turbine/gasdynamiks.py): втулочный d_sa и наружный D_sa.
Методика — Митрофанов, books/Mitrofanov.pdf, с.48…57, рис.18 (см. kamera_sgoraniya/mitrofanov.py).

Запуск из корня проекта: python kamera_sgoraniya/main.py
"""
import contextlib
import io
import sys

import matplotlib.pyplot as plt

from kamera_sgoraniya import mitrofanov as ks

MM = 1000.0


def compute(t, tvd, komp=None, **params):
    """Расчёт КС по результатам термодинамического расчёта t (с ЦБК t.komp) и ТВД tvd.
    komp — другой вариант расчёта ЦБК (для подбора диффузора, см. podbor_diffuzor.py)."""
    komp = komp or t.komp
    diff = dict(r3=komp["D3"] / 2, b3=komp["b3"], Rb=komp["Rb"], h4=komp["h4"], L_ld=komp["L_ld"])
    return ks.combustor(
        G_v=t.g_v * t.v_ks,  # расход воздуха в КС (за вычетом отбора на нужды ЛА)
        P_k=t.P2, T_k=t.T2, T_z=t.T_g, alpha=t.alpha, L_o=t.Km0,
        D_k_vn=komp["D4m"] - komp["h5"], D_k_n=komp["D4m"] + komp["h5"],
        D_vyh_vn=tvd.d_sa, D_vyh_n=tvd.D_sa, diff=diff, **params)


def turbine_mm(tvd):
    return dict(d_sa=tvd.d_sa * MM, D_sa=tvd.D_sa * MM, d_rk=tvd.d_rk * MM, D_rk=tvd.D_rk * MM,
                S_sa=tvd.S_sa * MM, S_rk=tvd.S_rk * MM)


def main():
    from kamera_sgoraniya import plot as ksplot

    # Импорт расчётов ЦБК и ТВД тихо (подавляем их собственные print/plt.show)
    _show = plt.show
    plt.show = lambda *a, **k: None
    with contextlib.redirect_stdout(io.StringIO()):
        import termo.termodynamiks as t
        import turbine.gasdynamiks as tvd
    plt.show = _show
    plt.close("all")

    r = compute(t, tvd)
    ks.print_report(r)

    g, tr, komp = r["geo"], turbine_mm(tvd), t.komp
    L_tvd = tr["S_sa"] + 0.25 * tr["S_rk"] + tr["S_rk"]
    print("\n─── Согласование с ЦБК и ТВД ───")
    print(f"  ЦБК: радиальный диффузор до D3 = {komp['D3'] * MM:.1f} мм (D3/D2 = {t.KOMP['D3_D2']}), колено, "
          f"осевой ЛД ⌀{r['D_k_vn'] * MM:.1f}/⌀{r['D_k_n'] * MM:.1f} мм → малый диффузор → наружный канал КС")
    print(f"  Корпус КС D_ксн = {r['D_ksn'] * MM:.1f} мм = наружная стенка малого диффузора за ЛД ЦБК")
    print(f"  ТВД, вход СА: ⌀{tr['d_sa']:.1f}/⌀{tr['D_sa']:.1f} мм → выход газосборника (стык по обоим диаметрам)")
    print(f"  Осевая длина ступени ТВД (СА + зазор + РК) {L_tvd:.1f} мм, место под жаровой трубой "
          f"{g['L_under'] * MM:.1f} мм → " + ("ступень целиком под КС" if L_tvd <= g["L_under"] * MM
                                               else "ступень выходит за задний торец КС"))

    ksplot.plot_combustor(r, komp, tr)


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
