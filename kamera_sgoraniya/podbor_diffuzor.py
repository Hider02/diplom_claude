"""
Подбор высоты радиального диффузора ЦБК (D3/D2) под противоточную КС.

В противоточной КС (рис.18 Митрофанова, чертёж books/гг.png) осевой ЛД ЦБК выходит прямо в наружный
кольцевой канал КС, поэтому корпус КС D_ксн = D_к + h_д1 задаётся компрессором, а жаровая труба
(п.4: Джтрвн = √(D_ксн² − 4/π·(F_кк + F_ж))) должна пройти над СА ТВД: Джтрвн ≥ D_са + 70 мм (с.57)
и l_гс/h_ср ≥ X3 (с.50). Скрипт находит наименьшее D3/D2, при котором это выполняется, при остальных
параметрах ЦБК из termo/termodynamiks.py (KOMP).

Изменение D3/D2 меняет КПД ЦБК, а через цикл — и размеры ТВД, поэтому после переноса значения в
termo/termodynamiks.py скрипт нужно запустить ещё раз и убедиться, что подобранное значение не изменилось.

Запуск из корня проекта: python kamera_sgoraniya/podbor_diffuzor.py
"""
import contextlib
import io
import math
import sys

import matplotlib.pyplot as plt

MARGIN = 0.001  # запас по Джтрвн, м


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

    def variant(D3_D2):
        komp = compressor(t.pik, t.g_v, t.P1, t.T1, **dict(t.KOMP, D3_D2=D3_D2))
        return komp, compute(t, tvd, komp=komp)

    def excess(D3_D2):
        try:
            _, r = variant(D3_D2)
        except ValueError:  # корпус КС слишком мал даже для жаровой трубы
            return -1.0
        return r["D_zhtr_vn"] - r["D_zhtr_vn_req"] - MARGIN

    lo, hi = 1.05, 2.5
    if excess(lo) >= 0:
        hi = lo
    else:
        for _ in range(50):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if excess(mid) < 0 else (lo, mid)
    D3_D2 = math.ceil(hi * 1000) / 1000
    komp, r = variant(D3_D2)
    print(f"Текущее D3/D2 = {t.KOMP['D3_D2']}, требуется D3/D2 = {D3_D2}")
    print(f"  D3 = {komp['D3'] * 1000:.1f} мм, выход ЛД ⌀{r['D_k_vn'] * 1000:.1f}/⌀{r['D_k_n'] * 1000:.1f} мм, "
          f"D_ксн = {r['D_ksn'] * 1000:.1f} мм")
    print(f"  Джтрвн = {r['D_zhtr_vn'] * 1000:.1f} мм (требуется ≥ {r['D_zhtr_vn_req'] * 1000:.1f}: {r['req_why']})")
    print(f"  η_к = {komp['eta_k']:.4f} (сейчас {t.komp['eta_k']:.4f}), c_к = {komp['c5']:.1f} м/с")
    print(f"  Замечания ЦБК: {komp['warnings'] or 'нет'}")
    print(f"  Замечания КС:  {r['warnings'] or 'нет'}")
    r_cur = compute(t, tvd)
    if r_cur["warnings"] or t.KOMP["D3_D2"] < D3_D2:
        print(f"\nПеренесите в termo/termodynamiks.py: KOMP[...] D3_D2={D3_D2} и запустите скрипт ещё раз")
    else:
        print(f"\n✅ Текущее D3/D2 = {t.KOMP['D3_D2']} согласовано с КС (Джтрвн = "
              f"{r_cur['D_zhtr_vn'] * 1000:.1f} мм, l_гс/h_ср = {r_cur['l_gs_rel']:.2f}, замечаний нет)")


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
