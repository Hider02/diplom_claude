"""
Визуализация меридионального сечения петлевой (противоточной) КС по Мингазову (рис. 8.3,
books/Мингазов.pdf) и чертежу реального двигателя books/гг.png; контур — kamera_sgoraniya/mingazov.py:layout().

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


def draw_combustor(ax, r, x0=0.0, dims=True, duct=None):
    """Меридиональное сечение петлевой (противоточной) КС на оси ax, x0 — выход осевого ЛД ЦБК, мм.
    duct — переходный канал ТВД → СТ (free_turbine/perehod.py, м, координаты КС): задняя стенка
    КС садится на его наружную стенку; без него — на корпус СА ТВД. Возвращает координаты (мм)."""
    g = {k: v * MM for k, v in r["geo"].items() if isinstance(v, float)}
    X = lambda v: x0 + v
    xc, xd, x_end, x_d1, x_back = X(g["x_c"]), X(g["x_dome"]), X(g["x_end"]), X(g["x_d1"]), X(g["x_back"])
    hub_x = np.array([X(p[0] * MM) for p in r["geo"]["hub"]])
    hub_r = np.array([p[1] * MM for p in r["geo"]["hub"]])
    h_zh, r_zh_sr = r["H_zh"] * MM, r["D_sr"] / 2 * MM
    # Скошенная задняя стенка корпуса и «пол» КС: корпус СА/РК ТВД и наружная стенка канала к СТ
    x_top = X(g["x_top"])
    if duct is not None:
        floor_x = [X(duct["x_b"] * MM), X(duct["x_rk"] * MM), xc]
        floor_r = [duct["r_b"] * MM, duct["r_rk_tip"] * MM, g["r_vyh_n"]]
    else:
        floor_x = [X(g["x_wall_t"]), xc]
        floor_r = [g["r_vyh_n"], g["r_vyh_n"]]

    out_x, out_r = _arc(xc, g["rc_o"], g["Ro"], 90, 270)  # наружная стенка газосборника: Джтрн → втулка СА
    in_x, in_r = _arc(xc, g["rc_i"], g["Ri"], 90, 270)  # внутренняя: Джтрвн → корпус СА
    dome_x, dome_r = _arc(xd, r_zh_sr, h_zh / 2, -90, 90)

    # Безотрывный участок диффузора КС (Мингазов, п. 8.3), за ним — внезапное расширение до корпуса D_н
    ax.fill([x0, x_d1, x_d1, x0], [g["r_k_vn"], g["r_d1_vn"], g["r_d1_n"], g["r_k_n"]], facecolor="lightyellow",
            edgecolor="black", lw=1.2, label="Диффузор КС (п. 8.3)")

    # Воздушная полость: наружный канал F_кк над жаровой трубой, за куполом вниз, внутренняя полость над ТВД,
    # спереди — полость между газосборником и задней стенкой диффузора ЦБК
    air_x = np.concatenate([[x_d1, x_d1, x_top], floor_x, in_x[::-1], [xd], dome_x, [xc], out_x,
                            [x_back, x_back], hub_x, [x0]])
    air_r = np.concatenate([[g["r_d1_vn"], g["r_ksn"], g["r_ksn"]], floor_r, in_r[::-1],
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
    # корпус: наружный, задняя стенка, внутренний (= корпус ТВД и наружная стенка канала к СТ)
    ax.plot([x_d1, x_d1, x_top] + floor_x, [g["r_d1_n"], g["r_ksn"], g["r_ksn"]] + floor_r, color="black", lw=1.8)
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
    ax.annotate("", xy=(xd + h_zh / 2 + 0.6 * r["dl"] * MM, r_zh_sr), xytext=(xd + 3, g["r_ksn"] - 1),
                arrowprops=dict(color="royalblue", connectionstyle="arc3,rad=-0.4", **arr))
    ax.annotate("", xy=(xc + 3, r_zh_sr), xytext=(xd - 3, r_zh_sr), arrowprops=dict(color="firebrick", **arr))
    ax.annotate("", xy=(xc + 2, (g["r_vyh_vn"] + g["r_vyh_n"]) / 2),
                xytext=(xc - g["Ro"] * 0.9, g["rc_o"] - 0.3 * g["Ro"]),
                arrowprops=dict(color="firebrick", connectionstyle="arc3,rad=0.4", **arr))

    c = dict(xc=xc, xd=xd, x_end=floor_x[0], x_top=x_top, x_back=x_back, x0=x0, g=g)
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
    g, xc, xd, x_end, x_top = c["g"], c["xc"], c["xd"], c["x_end"], c["x_top"]
    step = 8.0
    for i, (x_from, rr, name) in enumerate([
            (c["x0"], g["r_k_vn"], "D_к.вн"), (c["x0"], g["r_k_n"], "D_к.н"),
            (xc, g["r_vyh_vn"], "D_т.вн"), (x_end, g["r_vyh_n"], "D_т.н"),
            (xd, g["r_zh_vn"], "D_ж.вн"), (xd, g["r_zh_n"], "D_ж.н"), (x_end, g["r_ksn"], "D_н")]):
        _dim_diam(ax, min(x_from, x_top), x_top + step * (i + 1), rr, f"{name} ⌀{2 * rr:.0f}", 1)
    y_top = g["r_ksn"] + 12
    _dim_len(ax, c["x_back"], x_end, y_top + 16, f"L_кс = {x_end - c['x_back']:.0f}")
    _dim_len(ax, xc, xd, y_top, f"L_п.з+L_з.д = {xd - xc:.0f}")


def plot_combustor(r, komp, tvd, duct=None):
    """Отдельный рисунок ЦБК → противоточная КС → ТВД. komp — результат расчёта ЦБК (м),
    tvd — словарь размеров ТВД в мм (d_sa, D_sa, d_rk, D_rk, S_sa, S_rk)."""
    from layout.meridional import draw_compressor, draw_axial_turbine, scale_komp, draw_duct

    fig, ax = plt.subplots(figsize=(13, 8))
    _, _, _, x_out, _, _ = draw_compressor(ax, scale_komp(komp), 0.0, out_dims=False)
    c = draw_combustor(ax, r, x_out, duct=duct)
    draw_axial_turbine(ax, c["xc"], tvd, 0.25 * tvd["S_rk"], "ТВД", in_dims=False, label_y=0.2 * tvd["d_sa"] / 2)
    x_right = c["x_top"] + 8 * 8 + 15  # за размерными линиями КС
    if duct is not None:
        x_ct = draw_duct(ax, duct, x_out)
        ax.text(x_ct + 2, duct["r_in_tip"] * MM, "→ СТ", fontsize=10, fontweight="bold", va="center")
        x_right = max(x_right, x_ct + 20)
    ax.plot([-40, x_right], [0, 0], "k-.", lw=0.9)

    ax.set_aspect("equal")
    ax.set_xlim(-40, x_right)
    ax.set_ylim(-5, c["g"]["r_ksn"] + 45)
    ax.set_xlabel("Осевая координата от входа в ЦБК, мм")
    ax.set_ylabel("Радиус, мм")
    ax.set_title(f"ЦБК + петлевая КС (Мингазов) + ТВД: D3 = {komp['D3'] * MM:.0f} мм, D_н = {r['D_n'] * MM:.0f} мм, "
                 f"D_ж.н/вн = {r['D_zh_n'] * MM:.0f}/{r['D_zh_vn'] * MM:.0f} мм", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper left", fontsize=8)
    plt.tight_layout()

    out_path = __file__.replace("plot.py", "combustor_layout.png")
    fig.savefig(out_path, dpi=150)
    print(f"\nСохранено: {out_path}")
    plt.show()
