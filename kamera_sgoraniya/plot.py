"""
Визуализация меридионального сечения противоточной КС по рис.18 методики Митрофанова
(books/Mitrofanov.pdf, с.49) и чертежу реального двигателя books/гг.png; контур —
kamera_sgoraniya/mitrofanov.py:layout().

draw_combustor() рисует КС на уже созданной оси со сдвигом x0 = выход осевого ЛД ЦБК (используется и
здесь, и в общей схеме тракта layout/meridional.py); plot_combustor() — отдельный рисунок
ЦБК + КС + ТВД (ТВД стоит под жаровой трубой).
"""
import math

import matplotlib.pyplot as plt
import numpy as np

MM = 1000.0


def _arc(xc, rc, R, phi0, phi1, n=60):
    phi = np.linspace(math.radians(phi0), math.radians(phi1), n)
    return xc + R * np.cos(phi), rc + R * np.sin(phi)


def draw_combustor(ax, r, x0=0.0, dims=True):
    """Меридиональное сечение противоточной КС на оси ax, x0 — выход осевого ЛД ЦБК, мм.
    Возвращает словарь координат (мм) для стыковки с турбиной."""
    g = {k: v * MM for k, v in r["geo"].items() if isinstance(v, float)}
    X = lambda v: x0 + v
    xc, xd, x_end, x_d1, x_back = X(g["x_c"]), X(g["x_dome"]), X(g["x_end"]), X(g["x_d1"]), X(g["x_back"])
    hub_x = np.array([X(p[0] * MM) for p in r["geo"]["hub"]])
    hub_r = np.array([p[1] * MM for p in r["geo"]["hub"]])
    h_zh, r_zh_sr = r["h_zh"] * MM, r["D_zhtr_sr"] / 2 * MM

    out_x, out_r = _arc(xc, g["rc_o"], g["Ro"], 90, 270)  # наружная стенка газосборника: Джтрн → втулка СА
    in_x, in_r = _arc(xc, g["rc_i"], g["Ri"], 90, 270)  # внутренняя: Джтрвн → корпус СА
    dome_x, dome_r = _arc(xd, r_zh_sr, h_zh / 2, -90, 90)

    # Малый «безотрывный» диффузор (пп.14…16)
    ax.fill([x0, x_d1, x_d1, x0], [g["r_k_vn"], g["r_d1_vn"], g["r_d1_n"], g["r_k_n"]], facecolor="lightyellow",
            edgecolor="black", lw=1.2, label="Малый диффузор КС (пп.14…16)")

    # Воздушная полость: наружный канал F_кк над жаровой трубой, за куполом вниз, внутренняя полость над ТВД,
    # спереди — полость между газосборником и задней стенкой диффузора ЦБК
    air_x = np.concatenate([[x_d1, x_d1, x_end, x_end, xc], in_x[::-1], [xd], dome_x, [xc], out_x,
                            [x_back, x_back], hub_x, [x0]])
    air_r = np.concatenate([[g["r_d1_vn"], g["r_ksn"], g["r_ksn"], g["r_vyh_n"], g["r_vyh_n"]], in_r[::-1],
                            [g["r_zh_vn"]], dome_r, [g["r_zh_n"]], out_r,
                            [g["r_vyh_vn"], g["r3"]], hub_r, [g["r_k_vn"]]])
    ax.fill(air_x, air_r, color="lightblue", alpha=0.45, edgecolor="none", label="Воздух КС (F_кк и полости)")

    # Газ: жаровая труба + газосборник
    gas_x = np.concatenate([out_x[::-1], [xd], dome_x[::-1], [xc], in_x, [xc]])
    gas_r = np.concatenate([out_r[::-1], [g["r_zh_n"]], dome_r[::-1], [g["r_zh_vn"]], in_r, [g["r_vyh_vn"]]])
    ax.fill(gas_x, gas_r, color="mistyrose", edgecolor="none", label="Жаровая труба + газосборник")

    k = dict(color="black", lw=1.5)
    ax.plot(out_x, out_r, **k)
    ax.plot(in_x, in_r, **k)
    ax.plot([xc, xd], [g["r_zh_n"]] * 2, **k)
    ax.plot([xc, xd], [g["r_zh_vn"]] * 2, **k)
    ax.plot(dome_x, dome_r, **k)
    # корпус: наружный, задний торец, внутренний (= корпус СА ТВД)
    ax.plot([x_d1, x_end, x_end, xc], [g["r_ksn"], g["r_ksn"], g["r_vyh_n"], g["r_vyh_n"]], color="black", lw=1.8)
    # задняя стенка радиального диффузора ЦБК (передняя стенка КС) и опора СА по втулке
    ax.plot([x_back, x_back, xc], [g["r3"], g["r_vyh_vn"], g["r_vyh_vn"]], color="black", lw=1.8)
    ax.plot([x0, x_d1], [g["r_k_vn"], g["r_d1_vn"]], **k)
    ax.plot([x0, x_d1], [g["r_k_n"], g["r_d1_n"]], **k)
    ax.plot([xc, xc], [g["r_vyh_vn"], g["r_vyh_n"]], color="gray", lw=0.8, ls="--")

    # Форсунка через задний торец корпуса к куполу
    ax.plot([xd + h_zh / 2, x_end + 8], [r_zh_sr] * 2, color="dimgray", lw=3, solid_capstyle="butt")

    # Направления потоков
    arr = dict(arrowstyle="-|>", lw=1.3, mutation_scale=11)
    r_kk = (g["r_zh_n"] + g["r_ksn"]) / 2
    ax.annotate("", xy=(xd, r_kk), xytext=(x_d1 + 2, r_kk), arrowprops=dict(color="royalblue", **arr))
    ax.annotate("", xy=(xd + h_zh / 2 + 0.6 * r["Delta"] * MM, r_zh_sr), xytext=(xd + 3, g["r_ksn"] - 1),
                arrowprops=dict(color="royalblue", connectionstyle="arc3,rad=-0.4", **arr))
    ax.annotate("", xy=(xc + 3, r_zh_sr), xytext=(xd - 3, r_zh_sr), arrowprops=dict(color="firebrick", **arr))
    ax.annotate("", xy=(xc + 2, (g["r_vyh_vn"] + g["r_vyh_n"]) / 2),
                xytext=(xc - g["Ro"] * 0.9, g["rc_o"] - 0.3 * g["Ro"]),
                arrowprops=dict(color="firebrick", connectionstyle="arc3,rad=0.4", **arr))

    c = dict(xc=xc, xd=xd, x_end=x_end, x_back=x_back, x0=x0, g=g)
    if dims:
        _dimensions(ax, c)
    return c


def _dim_diam(ax, x_from, x_dim, rr, label, side=1):
    ax.plot([x_from, x_dim], [rr, rr], ":", color="dimgray", lw=0.7)
    ax.annotate("", xy=(x_dim, rr), xytext=(x_dim, 0), arrowprops=dict(arrowstyle="<->", color="dimgray", lw=0.8))
    ax.text(x_dim + side * 3.5, rr * 0.55, label, rotation=90, ha="center", va="center", fontsize=7.5,
            color="dimgray")


def _dim_len(ax, x0, x1, y, label):
    ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="<->", color="dimgray", lw=0.8))
    ax.text((x0 + x1) / 2, y + 2, label, ha="center", va="bottom", fontsize=7.5, color="dimgray")


def _dimensions(ax, c):
    g, xc, xd, x_end = c["g"], c["xc"], c["xd"], c["x_end"]
    step = 8.0
    for i, (x_from, rr, name) in enumerate([
            (c["x0"], g["r_k_vn"], "D_к.вн"), (c["x0"], g["r_k_n"], "D_к.н"),
            (xc, g["r_vyh_vn"], "D_вых.вн"), (x_end, g["r_vyh_n"], "D_вых.н"),
            (xd, g["r_zh_vn"], "Джтрвн"), (xd, g["r_zh_n"], "Джтрн"), (x_end, g["r_ksn"], "D_ксн")]):
        _dim_diam(ax, x_from, x_end + step * (i + 1), rr, f"{name} ⌀{2 * rr:.0f}", 1)
    y_top = g["r_ksn"] + 12
    _dim_len(ax, c["x_back"], x_end, y_top + 16, f"L_кс = {x_end - c['x_back']:.0f}")
    _dim_len(ax, xc, xd, y_top, f"l_зг = {xd - xc:.0f}")


def plot_combustor(r, komp, tvd):
    """Отдельный рисунок ЦБК → противоточная КС → ТВД. komp — результат расчёта ЦБК (м),
    tvd — словарь размеров ТВД в мм (d_sa, D_sa, d_rk, D_rk, S_sa, S_rk)."""
    from layout.meridional import draw_compressor, draw_axial_turbine, scale_komp

    fig, ax = plt.subplots(figsize=(13, 8))
    _, _, _, x_out, _, _ = draw_compressor(ax, scale_komp(komp), 0.0, out_dims=False)
    c = draw_combustor(ax, r, x_out)
    draw_axial_turbine(ax, c["xc"], tvd, 0.25 * tvd["S_rk"], "ТВД", in_dims=False)
    x_right = c["x_end"] + 70
    ax.plot([-40, x_right], [0, 0], "k-.", lw=0.9)

    ax.set_aspect("equal")
    ax.set_xlim(-40, x_right)
    ax.set_ylim(-5, c["g"]["r_ksn"] + 45)
    ax.set_xlabel("Осевая координата от входа в ЦБК, мм")
    ax.set_ylabel("Радиус, мм")
    ax.set_title(f"ЦБК + противоточная КС + ТВД: D3 = {komp['D3'] * MM:.0f} мм, D_ксн = {r['D_ksn'] * MM:.0f} мм, "
                 f"Джтрн/вн = {r['D_zhtr_n'] * MM:.0f}/{r['D_zhtr_vn'] * MM:.0f} мм", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper left", fontsize=8)
    plt.tight_layout()

    out_path = __file__.replace("plot.py", "combustor_layout.png")
    fig.savefig(out_path, dpi=150)
    print(f"\nСохранено: {out_path}")
    plt.show()
