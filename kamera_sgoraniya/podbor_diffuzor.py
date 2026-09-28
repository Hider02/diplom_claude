"""
Подбор высоты радиального диффузора ЦБК (D3/D2) под петлевую КС по Мингазову (kamera_sgoraniya/mingazov.py).

В петлевой КС наружный диаметр КС соответствует наружному диаметру ЦБК (Мингазов, п. 8.10): осевой ЛД ЦБК
выходит в диффузор КС (п. 8.3), и через внезапное расширение воздух уходит в наружный канал над ЖТ.
ЖТ по п. 8.10 имеет высоту H_ж = 0,7·H_к, H_к = (D_н − D_т)/2, и стоит по среднему диаметру КС. Чем выше
диффузор ЦБК, тем больше D_н и H_к, тем ниже наружная стенка ЖТ относительно зоны внезапного расширения.
Скрипт находит наименьшее D3/D2 (самую компактную КС и наименьшие потери в радиальном диффузоре), при
котором:
  - наружная стенка ЖТ не выше внутренней стенки зоны внезапного расширения (D_ж.н ≤ D_2вн) — поток из
    диффузора КС уходит в наружный канал без поджатия;
  - газосборник и ЖТ не задевают стенку колена и ЛД ЦБК (зазор ≥ CLEAR_MIN).
Остальные параметры ЦБК — из termo/termodynamiks.py (KOMP).

Изменение D3/D2 меняет КПД ЦБК, а через цикл — и размеры ТВД, поэтому после переноса значения в
termo/termodynamiks.py скрипт нужно запустить ещё раз и убедиться, что подобранное значение не изменилось.

Запуск из корня проекта: python kamera_sgoraniya/podbor_diffuzor.py
"""
import contextlib
import io
import math
import sys

import matplotlib.pyplot as plt

TOL = 0.002  # допуск согласования текущего D3/D2 с подобранным


def main():
    _show = plt.show
    plt.show = lambda *a, **k: None
    with contextlib.redirect_stdout(io.StringIO()):
        import termo.termodynamiks as t
        import turbine.gasdynamiks as tvd
    plt.show = _show
    plt.close("all")
    from kompressor.cbk import compressor
    from kamera_sgoraniya.main import compute
    from kamera_sgoraniya.mingazov import CLEAR_MIN

    def variant(D3_D2):
        komp = compressor(t.pik, t.g_v, t.P1, t.T1, **dict(t.KOMP, D3_D2=D3_D2))
        return komp, compute(t, tvd, komp=komp)

    def excess(D3_D2):
        _, r = variant(D3_D2)
        return min((r["D2_vn"] - r["D_zh_n"]) / 2, r["geo"]["clear"] - CLEAR_MIN)

    lo, hi = 1.05, 2.5
    for _ in range(50):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if excess(mid) < 0 else (lo, mid)
    D3_D2 = math.ceil(hi * 1000) / 1000
    komp, r = variant(D3_D2)
    print(f"Текущее D3/D2 = {t.KOMP['D3_D2']}, требуется D3/D2 = {D3_D2}")
    print(f"  D3 = {komp['D3'] * 1000:.1f} мм, выход ЛД ⌀{r['D_k_vn'] * 1000:.1f}/⌀{r['D_k_n'] * 1000:.1f} мм, "
          f"D_н = {r['D_n'] * 1000:.1f} мм")
    print(f"  D_ж.н = {r['D_zh_n'] * 1000:.1f} мм ≤ D_2вн = {r['D2_vn'] * 1000:.1f} мм; зазор до стенки ЦБК "
          f"{r['geo']['clear'] * 1000:.1f} мм")
    print(f"  η_к = {komp['eta_k']:.4f} (сейчас {t.komp['eta_k']:.4f}), c_к = {komp['c5']:.1f} м/с")
    print(f"  Замечания ЦБК: {komp['warnings'] or 'нет'}")
    if abs(t.KOMP["D3_D2"] - D3_D2) > TOL:
        print(f"\nПеренесите в termo/termodynamiks.py: KOMP[...] D3_D2={D3_D2} и запустите скрипт ещё раз")
    else:
        print(f"\n✅ Текущее D3/D2 = {t.KOMP['D3_D2']} согласовано с КС")


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
