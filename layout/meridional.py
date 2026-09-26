"""
Совместная меридиональная схема тракта двигателя: ЦБК (центробежный компрессор) + КС (схематично) +
ТВД (турбина высокого давления, на одном валу с ЦБК) + переходной канал + СТ (свободная турбина).

Геометрия узлов берётся из уже существующих расчётов без их изменения:
  - ЦБК: kompressor/cbk.py (через termo/termodynamiks.py, где результат уже посчитан как t.komp);
  - ТВД: turbine/gasdynamiks.py;
  - СТ:  free_turbine/gasdynamiks.py.
Все три модуля при импорте сами печатают протокол расчёта и строят свои графики — здесь их собственные
print и plt.show() временно подавляются, чтобы не открывать лишние окна и не засорять вывод.

Расстояния между узлами (длина КС и переходного канала) заданы условно — по задаче важны только
радиальные размеры (они взяты из расчёта без изменений), а не реальная компоновка по длине.

Запуск из корня проекта: python layout/meridional.py
Результат: layout/meridional_layout.png (и окно, если доступен интерактивный backend).
"""
import contextlib
import io
import math

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

MM = 1000.0  # перевод м -> мм для всех расчётных величин ниже

# Условные (не из расчёта) осевые зазоры между узлами, мм — только чтобы разнести блоки на рисунке.
# КС принята противоточной: поток разворачивается внутри жаровой трубы, поэтому по факту КС не растягивает
# тракт компрессор-турбина так, как показал бы прямоточный вариант — зазор здесь взят небольшим.
GAP_KS = 60.0      # камера сгорания (противоточная, схематично)
GAP_MID = 80.0     # переходной канал между ТВД и СТ


def _get_geometry():
    """Импортирует существующие расчёты, подавляя их собственные print/plt.show(), и возвращает геометрию в мм."""
    _show = plt.show
    plt.show = lambda *a, **k: None
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            import termo.termodynamiks as t
            import turbine.gasdynamiks as tvd
            import free_turbine.gasdynamiks as ct
    finally:
        plt.show = _show
        plt.close("all")

    komp = t.komp
    scale_keys = {"S", "b2", "b3", "h4", "D1h", "D1t", "D2", "D3", "D4m", "Rb", "h5"}
    r = {k: (v * MM if k in scale_keys else v) for k, v in komp.items()}

    turb = dict(D_sa=tvd.D_sa * MM, d_sa=tvd.d_sa * MM, D_rk=tvd.D_rk * MM, d_rk=tvd.d_rk * MM,
                S_sa=tvd.S_sa * MM, S_rk=tvd.S_rk * MM, D_sr=tvd.D_sr * MM)
    st = dict(D_sa=ct.D_sa * MM, d_sa=ct.d_sa * MM, D_rk=ct.D_rk * MM, d_rk=ct.d_rk * MM,
              S_sa=ct.S_sa * MM, S_rk=ct.S_rk * MM, D_sr=ct.D_sr * MM)
    return r, turb, st


def draw_compressor(ax, r, x0):
    """Меридиональное сечение ЦБК (копия геометрии kompressor/plot.py:plot_meridional, сдвинутая на x0).
    Возвращает (x_in, hub_in, tip_in, x_out, hub_out, tip_out)."""
    S, b2, b3, h4 = r["S"], r["b2"], r["b3"], r["h4"]
    r_vt, r1, r2, r3, r4 = r["D1h"] / 2, r["D1t"] / 2, r["D2"] / 2, r["D3"] / 2, r["D4m"] / 2
    Rb = r["Rb"]
    t = np.linspace(0, math.pi / 2, 60)
    hub_x, hub_r = x0 + S * np.sin(t), r_vt + (r2 - r_vt) * (1 - np.cos(t))
    shr_x, shr_r = x0 + (S - b2) * np.sin(t), r1 + (r2 - r1) * (1 - np.cos(t))
    x_c = x0 + S - b2 / 2
    phi = np.linspace(0, math.pi / 2, 60)
    b_phi = b3 + (h4 - b3) * phi / (math.pi / 2)
    cx, cr = x_c + Rb * (1 - np.cos(phi)), r3 + Rb * np.sin(phi)
    inner = (cx + b_phi / 2 * np.cos(phi), cr - b_phi / 2 * np.sin(phi))
    outer = (cx - b_phi / 2 * np.cos(phi), cr + b_phi / 2 * np.sin(phi))
    x4 = x_c + Rb
    n_rows = len(r["rows"])
    row_len, gap = 1.3 * h4, 0.4 * h4
    L_ld = n_rows * row_len + (n_rows + 1) * gap
    r5_in, r5_out = r4 - r["h5"] / 2, r4 + r["h5"] / 2
    x_out = x4 + L_ld

    ax.fill(np.concatenate([hub_x, shr_x[::-1]]), np.concatenate([hub_r, shr_r[::-1]]),
            color="lightcoral", alpha=0.6, label="Рабочее колесо ЦБК")
    ax.fill([x0 + S - b2, x0 + S, x_c + b3 / 2, x_c - b3 / 2], [r2, r2, r3, r3],
            color="lightyellow", label="Радиальный БЛД")
    ax.fill(np.concatenate([inner[0], outer[0][::-1]]), np.concatenate([inner[1], outer[1][::-1]]),
            color="honeydew", label="Переходное колено")
    ax.fill([x4, x_out, x_out, x4], [r4 - h4 / 2, r5_in, r5_out, r4 + h4 / 2], color="lightblue",
            alpha=0.8, label="Осевой ЛД")
    for j in range(n_rows):
        xr0 = x4 + gap + j * (row_len + gap)
        ax.add_patch(patches.Rectangle((xr0, r4 - h4 / 2), row_len, h4, color="steelblue", alpha=0.8))
    for xs, rs in ((hub_x, hub_r), (shr_x, shr_r), inner, outer):
        ax.plot(xs, rs, "k-", lw=1.4)
    ax.plot([x0 + S - b2, x_c - b3 / 2], [r2, r3], "k-", lw=1.4)
    ax.plot([x0 + S, x_c + b3 / 2], [r2, r3], "k-", lw=1.4)
    ax.plot([x4, x_out], [r4 - h4 / 2, r5_in], "k-", lw=1.4)
    ax.plot([x4, x_out], [r4 + h4 / 2, r5_out], "k-", lw=1.4)
    ax.text((x0 + x_out) / 2, r4 + h4 * 3.2, "ЦБК", ha="center", fontsize=13, fontweight="bold")

    # Размерные линии: диаметры на входе (D1) и на выходе из осевого ЛД, осевая длина колеса и ЛД
    dx1, dx2 = 0.16 * r1, 0.34 * r1
    off = 0.12 * r1
    _dim_diam(ax, x0, x0 - dx1, r_vt, f"⌀{2 * r_vt:.0f}", -1)
    _dim_diam(ax, x0, x0 - dx2, r1, f"⌀{2 * r1:.0f}", -1)
    _dim_diam(ax, x_out, x_out + dx1, r5_in, f"⌀{2 * r5_in:.0f}", 1)
    _dim_diam(ax, x_out, x_out + dx2, r5_out, f"⌀{2 * r5_out:.0f}", 1)
    _dim_width(ax, x0, x0 + S, r1, r2, r2 + off, f"S={S:.0f}")
    _dim_width(ax, x4, x_out, r4 + h4 / 2, r4 + h4 / 2, r4 + h4 / 2 + off, f"L={L_ld:.0f}")

    return x0, r_vt, r1, x_out, r5_in, r5_out


def _dim_diam(ax, x_block, x_dim, r, label, side):
    """Размерная линия диаметра: от оси (y=0) до поверхности (y=r), с линией-выноской от контура.
    side = -1 (линия слева от узла) или +1 (справа) — определяет, куда сдвинута подпись."""
    ax.plot([x_block, x_dim], [r, r], ":", color="dimgray", lw=0.7)
    ax.annotate("", xy=(x_dim, r), xytext=(x_dim, 0),
                arrowprops=dict(arrowstyle="<->", color="dimgray", lw=0.9))
    ax.text(x_dim + side * 3.5, r / 2, label, rotation=90, ha="center", va="center", fontsize=7.5,
            color="dimgray")


def _dim_width(ax, x0, x1, y0, y1, y_dim, label):
    """Размерная линия осевой ширины узла (СА или РК), под втулочной линией."""
    ax.plot([x0, x0], [y0, y_dim], ":", color="dimgray", lw=0.7)
    ax.plot([x1, x1], [y1, y_dim], ":", color="dimgray", lw=0.7)
    ax.annotate("", xy=(x1, y_dim), xytext=(x0, y_dim),
                arrowprops=dict(arrowstyle="<->", color="dimgray", lw=0.9))
    ax.text((x0 + x1) / 2, y_dim, label, ha="center", va="top", fontsize=7.5, color="dimgray")


def draw_axial_turbine(ax, x0, geo, axial_gap, label, sa_color="lightblue", rk_color="lightcoral"):
    """Меридиональное сечение осевой ступени СА+РК (копия geometрии turbine/plot_turb.py:plot_meridional_section,
    сдвинутая на x0). Возвращает (x_in, hub_in, tip_in, x_out, hub_out, tip_out)."""
    S_sa, S_rk = geo["S_sa"], geo["S_rk"]
    D_ca_hub, D_ca_tip = geo["d_sa"], geo["D_sa"]
    D_rk_hub, D_rk_tip = geo["d_rk"], geo["D_rk"]

    z_ca_in = x0
    z_ca_out = z_ca_in + S_sa
    z_gap_out = z_ca_out + axial_gap
    z_rk_in = z_gap_out
    z_rk_out = z_rk_in + S_rk

    r_ca_hub, r_ca_tip = D_ca_hub / 2, D_ca_tip / 2
    r_rk_hub_out, r_rk_tip_out = D_rk_hub / 2, D_rk_tip / 2

    ax.add_patch(patches.Polygon([(z_ca_in, r_ca_hub), (z_ca_in, r_ca_tip), (z_ca_out, r_ca_tip),
                                  (z_ca_out, r_ca_hub)], closed=True, lw=1.4, edgecolor="black",
                                 facecolor=sa_color, alpha=0.85))
    ax.add_patch(patches.Polygon([(z_ca_out, r_ca_hub), (z_ca_out, r_ca_tip), (z_gap_out, r_ca_tip),
                                  (z_gap_out, r_ca_hub)], closed=True, lw=1.0, edgecolor="gray",
                                 facecolor="lightyellow", alpha=0.6, linestyle="--"))
    ax.add_patch(patches.Polygon([(z_rk_in, r_ca_hub), (z_rk_in, r_ca_tip), (z_rk_out, r_rk_tip_out),
                                  (z_rk_out, r_rk_hub_out)], closed=True, lw=1.4, edgecolor="black",
                                 facecolor=rk_color, alpha=0.85))
    ax.plot([z_ca_in, z_ca_out, z_gap_out, z_rk_out], [r_ca_hub, r_ca_hub, r_ca_hub, r_rk_hub_out], "k-", lw=1.6)
    ax.plot([z_ca_in, z_ca_out, z_gap_out, z_rk_out], [r_ca_tip, r_ca_tip, r_ca_tip, r_rk_tip_out], "k-", lw=1.6)
    ax.text((z_ca_in + z_ca_out) / 2, r_ca_tip + 6, "СА", ha="center", fontsize=10, fontweight="bold",
            color="darkblue")
    ax.text((z_rk_in + z_rk_out) / 2, r_rk_tip_out + 6, "РК", ha="center", fontsize=10, fontweight="bold",
            color="darkred")
    ax.text((z_ca_in + z_rk_out) / 2, max(r_ca_tip, r_rk_tip_out) + 22, label, ha="center", fontsize=13,
            fontweight="bold")

    # Размерные линии: диаметры втулки/корпуса на входе в СА и на выходе из РК, осевая ширина СА и РК
    dx1, dx2 = 0.16 * r_ca_tip, 0.34 * r_ca_tip
    off = 0.14 * r_ca_tip
    _dim_diam(ax, z_ca_in, z_ca_in - dx1, r_ca_hub, f"⌀{D_ca_hub:.0f}", -1)
    _dim_diam(ax, z_ca_in, z_ca_in - dx2, r_ca_tip, f"⌀{D_ca_tip:.0f}", -1)
    _dim_diam(ax, z_rk_out, z_rk_out + dx1, r_rk_hub_out, f"⌀{D_rk_hub:.0f}", 1)
    _dim_diam(ax, z_rk_out, z_rk_out + dx2, r_rk_tip_out, f"⌀{D_rk_tip:.0f}", 1)
    _dim_width(ax, z_ca_in, z_ca_out, r_ca_hub, r_ca_hub, r_ca_hub - off, f"S={S_sa:.0f}")
    _dim_width(ax, z_rk_in, z_rk_out, r_ca_hub, r_rk_hub_out, r_rk_hub_out - off, f"S={S_rk:.0f}")
    _dim_width(ax, z_ca_out, z_gap_out, r_ca_hub, r_ca_hub, r_ca_hub - 2.1 * off, f"зазор={axial_gap:.0f}")

    return z_ca_in, r_ca_hub, r_ca_tip, z_rk_out, r_rk_hub_out, r_rk_tip_out


def label_gap(ax, x0, x1, y, text):
    """Текстовая подпись условного зазора между узлами (без отрисовки переходного канала).
    va="top" — точка y задаёт верх текстового блока, чтобы многострочная подпись росла вниз."""
    ax.text((x0 + x1) / 2, y, text, ha="center", va="top", fontsize=11, fontweight="bold", style="italic")


def main():
    r, tvd, ct = _get_geometry()

    fig, ax = plt.subplots(figsize=(16, 7))

    x_in_c, hub_in_c, tip_in_c, x_out_c, hub_out_c, tip_out_c = draw_compressor(ax, r, 0.0)

    # Верхняя граница по радиусу известна заранее — подписи зазоров размещаются у неё, выше всех
    # размерных линий диаметров (те не поднимаются выше макс. радиуса узлов) во избежание наложения
    r_max = max(tip_out_c, tvd["D_sa"] / 2, ct["D_sa"] / 2, ct["D_rk"] / 2) * 1.28
    y_label = r_max * 0.85

    x_ks_end = x_out_c + GAP_KS
    label_gap(ax, x_out_c, x_ks_end, y_label, "КС\n(противоточная,\nсхематично)")

    axial_gap_tvd = 0.25 * tvd["S_rk"]
    (x_in_t, hub_in_t, tip_in_t, x_out_t, hub_out_t, tip_out_t) = draw_axial_turbine(
        ax, x_ks_end, tvd, axial_gap_tvd, "ТВД")

    x_mid_end = x_out_t + GAP_MID
    label_gap(ax, x_out_t, x_mid_end, y_label, "переходной\nканал")

    axial_gap_ct = 0.25 * ct["S_rk"]
    (x_in_s, hub_in_s, tip_in_s, x_out_s, hub_out_s, tip_out_s) = draw_axial_turbine(
        ax, x_mid_end, ct, axial_gap_ct, "СТ (свободная турбина)")

    x_total = x_out_s
    # Запас под размерные линии диаметров на входе в ЦБК (слева) и на выходе из СТ (справа) —
    # те же формулы dx2, что и в draw_compressor/draw_axial_turbine
    x_left = -(0.34 * (r["D1t"] / 2) + 12)
    x_right = x_total + 0.34 * (ct["D_sa"] / 2) + 12
    ax.plot([x_left, x_right], [0, 0], "k-.", lw=0.9)

    ax.set_aspect("equal")
    ax.set_xlim(x_left, x_right)
    ax.set_ylim(-15, r_max)
    ax.set_xlabel("Осевая координата, мм (контуры узлов — в масштабе; расстояния между узлами — условные)")
    ax.set_ylabel("Радиус, мм")
    ax.set_title("Меридиональная схема тракта: ЦБК → КС → ТВД → СТ (общий масштаб по радиусу)",
                fontsize=13, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)
    plt.tight_layout()

    out_path = __file__.replace("meridional.py", "meridional_layout.png")
    fig.savefig(out_path, dpi=150)
    print(f"Сохранено: {out_path}")

    print("\n─── Радиальные размеры узлов (мм) ───")
    print(f"{'Узел':<28} | {'D_вт':>8} | {'D_нар':>8} | {'D_ср':>8}")
    print("-" * 62)
    print(f"{'ЦБК: вход (D1)':<28} | {r['D1h']:>8.1f} | {r['D1t']:>8.1f} | {'-':>8}")
    print(f"{'ЦБК: выход РК (D2)':<28} | {'-':>8} | {'-':>8} | {r['D2']:>8.1f}")
    print(f"{'ЦБК: выход БЛД (D3)':<28} | {'-':>8} | {'-':>8} | {r['D3']:>8.1f}")
    print(f"{'ЦБК: выход осевого ЛД':<28} | {2 * hub_out_c:>8.1f} | {2 * tip_out_c:>8.1f} | {r['D4m']:>8.1f}")
    print(f"{'ТВД: вход СА':<28} | {tvd['d_sa']:>8.1f} | {tvd['D_sa']:>8.1f} | {tvd['D_sr']:>8.1f}")
    print(f"{'ТВД: выход РК':<28} | {tvd['d_rk']:>8.1f} | {tvd['D_rk']:>8.1f} | {'-':>8}")
    print(f"{'СТ: вход СА':<28} | {ct['d_sa']:>8.1f} | {ct['D_sa']:>8.1f} | {ct['D_sr']:>8.1f}")
    print(f"{'СТ: выход РК':<28} | {ct['d_rk']:>8.1f} | {ct['D_rk']:>8.1f} | {'-':>8}")

    plt.show()


if __name__ == "__main__":
    main()
