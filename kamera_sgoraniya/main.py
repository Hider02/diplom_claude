"""
Расчёт противоточной камеры сгорания, увязанный с уже выполненными расчётами ЦБК и ТВД:
вход в КС берётся из сечения за осевым ЛД компрессора (kompressor/cbk.py, через
termo/termodynamiks.py), выход из КС — из сечения на входе в СА ТВД (turbine/gasdynamiks.py).
Методика — Митрофанов, books/Mitrofanov.pdf, с.48…57 (см. kamera_sgoraniya/mitrofanov.py).

Запуск из корня проекта: python kamera_sgoraniya/main.py
"""
import math
import contextlib
import io

import matplotlib.pyplot as plt

# Импорт расчётов ЦБК и ТВД тихо (подавляем их собственные print/plt.show)
_show = plt.show
plt.show = lambda *a, **k: None
_buf = io.StringIO()
with contextlib.redirect_stdout(_buf):
    import termo.termodynamiks as t
    import turbine.gasdynamiks as tvd
plt.show = _show
plt.close("all")

from kamera_sgoraniya import mitrofanov as ks
from kamera_sgoraniya import plot as ksplot

komp = t.komp

# ── Вход в КС: сечение за осевым ЛД компрессора ──
D_k = komp["D4m"]
h_k = komp["h4"]
F_k = math.pi * D_k * h_k

# ── Выход из КС: сечение на входе в СА ТВД ──
D_vyh = (tvd.d_sa + tvd.D_sa) / 2
h_vyh = tvd.l_c
F_vyh = math.pi / 4 * (tvd.D_sa ** 2 - tvd.d_sa ** 2)

# ── Термодинамика: расход, давление и температура за компрессором, состав смеси ──
G_v = t.g_v * t.v_ks  # расход воздуха, поступающего в КС (за вычетом отбора на нужды ЛА)
P_k = t.P2
T_k = t.T2
T_z = t.T_g
alpha = t.alpha
L_o = t.Km0

r = ks.combustor(G_v, P_k, T_k, T_z, alpha, L_o, D_k, F_k, h_k, D_vyh, F_vyh, h_vyh)
ks.print_report(r)

D_komp_naruzh = komp["D4m"] + komp["h5"]
print("\n─── Согласование с габаритом компрессора ───")
print(f"  Наружный диаметр ЦБК (выход осевого ЛД): {D_komp_naruzh * 1000:.1f} мм")
print(f"  Наружный диаметр корпуса КС D_ксн:        {r['D_ksn'] * 1000:.1f} мм "
      f"({(r['D_ksn'] / D_komp_naruzh - 1) * 100:+.0f}% к диаметру компрессора)")
if r["D_ksn"] > D_komp_naruzh:
    print("  ⚠️ КС шире компрессора — на компоновке потребуется плавный переход по наружному "
          "диаметру между компрессором и КС (см. визуализацию)")

ksplot.plot_combustor(r)
