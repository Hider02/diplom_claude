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
    # «Пол» КС: корпус СА/РК ТВД и наружная стенка канала к СТ
    x_top = X(g["x_top"])
    # Угол корпус D_н → задняя стенка срезан прямой (как на рис. 8.3 Мингазова)
    arc_x = [X(g["ch_x_top"]), X(g["ch_x_w"])]
    arc_r = [g["r_ksn"], g["ch_r_w"]]
    # Передняя стенка корпуса у ЦБК: вертикально до носка газосборника, ниже — дугой вокруг него (рис. 8.3)
    fr_x, fr_r = _arc(xc, g["rc_o"], g["R_front"], 180, 270)
    fr_x, fr_r = list(fr_x), list(fr_r)
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
    air_x = np.concatenate([[x_d1, x_d1], arc_x, floor_x, in_x[::-1], [xd], dome_x, [xc], out_x,
                            fr_x[::-1], [x_back], hub_x, [x0]])
    air_r = np.concatenate([[g["r_d1_vn"], g["r_ksn"]], arc_r, floor_r, in_r[::-1],
                            [g["r_zh_vn"]], dome_r, [g["r_zh_n"]], out_r,
                            fr_r[::-1], [g["r3"]], hub_r, [g["r_k_vn"]]])
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
    ax.plot([x_d1, x_d1] + arc_x + floor_x, [g["r_d1_n"], g["r_ksn"]] + arc_r + floor_r, color="black", lw=1.8)
    # передняя стенка КС: задняя стенка радиального диффузора ЦБК, ниже носка газосборника — дуга вокруг него
    ax.plot([x_back] + fr_x + [xc], [g["r3"]] + fr_r + [g["r_vyh_vn"]], color="black", lw=1.8)
    ax.plot([x0, x_d1], [g["r_k_vn"], g["r_d1_vn"]], **k)
    ax.plot([x0, x_d1], [g["r_k_n"], g["r_d1_n"]], **k)
    ax.plot([xc, xc], [g["r_vyh_vn"], g["r_vyh_n"]], color="gray", lw=0.8, ls="--")

    # Форсунка — через срез корпуса, перпендикулярно ему, к куполу (как на рис. 8.3)
    nx, nr = r["geo"]["ch_nx"], r["geo"]["ch_nr"]
    r_f0, r_f1 = h_zh / 2, h_zh / 2 + r["dl"] * MM + 15
    ax.plot([xd + r_f0 * nx, xd + r_f1 * nx], [r_zh_sr + r_f0 * nr, r_zh_sr + r_f1 * nr], color="dimgray",
            lw=3, solid_capstyle="butt")

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

    c = dict(xc=xc, xd=xd, x_end=floor_x[0], x_top=x_top, x_wall=x_end, x_back=x_back, x0=x0, g=g)
    if dims:
        _dimensions(ax, c)
    return c


def _spread(segments, n):
    """n точек, равномерно разнесённых по отрезкам segments [(a, b), ...] (как по одному склеенному отрезку)."""
    segments = [(a, b) for a, b in segments if b >= a]
    total = sum(b - a for a, b in segments)
    pts = []
    for k in range(n):
        t = total * k / (n - 1) if n > 1 else 0.0
        for a, b in segments:
            if t <= b - a + 1e-9:
                pts.append(a + t)
                break
            t -= b - a
    return pts


def _belt(ax, P, t, nv, lip, step, gap):
    """Условный пояс охлаждения в точке стенки P: t — направление течения газа вдоль стенки, nv — нормаль внутрь
    ЖТ (единичные векторы). Окно gap в стенке, уступ step внутрь, козырёк lip по потоку, плёнка воздуха."""
    P, t, nv = np.asarray(P), np.asarray(t), np.asarray(nv)
    pts = lambda *q: np.array(q).T
    ax.plot(*pts(P, P + gap * t), color="white", lw=2.5, solid_capstyle="butt", zorder=3)
    ax.plot(*pts(P, P + step * nv, P + step * nv + lip * t), color="black", lw=1.3, zorder=4)
    deg = lambda v: math.degrees(math.atan2(v[1], v[0]))
    ax.annotate("", xy=P + (lip + 2.5) * t + 0.5 * step * nv, xytext=P + 0.5 * gap * t - 4.0 * nv, zorder=5,
                arrowprops=dict(arrowstyle="-|>", color="royalblue", lw=0.9, mutation_scale=7,
                                connectionstyle=f"angle,angleA={deg(nv):.1f},angleB={deg(t):.1f},rad=2"))


def draw_liner_holes(ax, r, c, n_turn=2, lip=4.0, step=1.8, gap=1.8, margin=1.0, legend_kw=None):
    """Отверстия и пояса охлаждения ЖТ (п. 8.6); газ на прямом участке идёт в −x, расстояния u отсчитываются от
    стыка купола с прямым участком:
      - отверстия зоны горения (ОЗГ, d_озг) — ряд на границе первичной зоны, u = L_п.з;
      - отверстия смесителя (ОЗС, d_озс) — ряд в конце зоны догорания, вплотную к газосборнику (смесителю);
      - пояса охлаждения (п. 8.6.7, n_ohl): n_ohl − n_turn — на обеих стенках прямого участка, равномерно
        в промежутках между куполом и рядами отверстий; n_turn — на наружной стенке газосборника (поворот
        на 180°, самая длинная и горячая стенка смесителя), равномерно по углу. Внутренняя стенка газосборника —
        дуга R ≈ 6 мм — сразу переходит в корпус СА, поясов на ней нет.
    В сечении отверстие — разрыв стенки шириной d со струёй глубиной m·H_ж.т (п. 8.6.1). Щель пояса h_s ≪ 1 мм
    в масштабе не видна, поэтому пояс условный (мм): секция стенки выше по потоку уступом step уходит внутрь ЖТ
    и козырьком lip перекрывает начало следующей; воздух входит в окно gap и выходит плёнкой вдоль стенки.
    Осевые положения рядов и поясов в расчёте не задаются — приняты здесь только для рисунка.
    legend_kw — положение легенды отверстий (по умолчанию — под основной легендой combustor_layout.png)."""
    from matplotlib.lines import Line2D

    g, xc, xd = c["g"], c["xc"], c["xd"]
    walls = ((g["r_zh_n"], -1), (g["r_zh_vn"], 1))  # s — направление внутрь ЖТ
    d_g, d_s = r["d_ozg"] * MM, r["d_ozs"] * MM
    u_g = r["L_pz"] * MM
    u_s = (xd - xc) - d_s / 2 - margin
    jet = r["params"]["m"] * r["H_zh"] * MM

    rows = (("navy", d_g, u_g), ("teal", d_s, u_s))
    for color, d, u in rows:
        xh = xd - u
        for rw, s in walls:
            ax.plot([xh - d / 2, xh + d / 2], [rw, rw], color="white", lw=2.5, solid_capstyle="butt", zorder=3)
            ax.annotate("", xy=(xh, rw + s * jet), xytext=(xh, rw - s * 2.5), zorder=5,
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=1.4, mutation_scale=10))

    n = r["n_ohl"]
    n_turn = min(n_turn, n)
    free = [(0.0, u_g - d_g / 2 - margin - lip), (u_g + d_g / 2 + margin, u_s - d_s / 2 - margin - lip)]
    for u in _spread(free, n - n_turn):
        for rw, s in walls:
            _belt(ax, (xd - u, rw), (-1.0, 0.0), (0.0, s), lip, step, gap)
    # наружная стенка газосборника: дуга от D_ж.н (φ = 90°) через носок (180°) к втулке СА (270°)
    for k in range(n_turn):
        phi = math.radians(90 + 180 * (k + 1) / (n_turn + 1))
        cs, sn = math.cos(phi), math.sin(phi)
        _belt(ax, (xc + g["Ro"] * cs, g["rc_o"] + g["Ro"] * sn), (-sn, cs), (-cs, -sn), lip, step, gap)

    handles = [
        Line2D([], [], color="navy", lw=1.4, marker="^", ms=5,
               label=f"ОЗГ ⌀{d_g:.1f} мм: {r['n_ozg_n']} + {r['n_ozg_vn']} = {r['n_ozg']} шт., t = {r['t_o'] * MM:.0f} мм"),
        Line2D([], [], color="teal", lw=1.4, marker="^", ms=5,
               label=f"ОЗС ⌀{d_s:.1f} мм: {r['n_ozs_n']} + {r['n_ozs_vn']} = {r['n_ozs']} шт., t = {r['t_ozs'] * MM:.0f} мм"),
        Line2D([], [], color="royalblue", lw=0.9, marker="<", ms=3,
               label=f"Пояса охлаждения: {n} шт., щель {r['h_s'] * MM:.2f} мм\n"
                     f"({n - n_turn} на прямом участке, {n_turn} на газосборнике)"),
    ]
    legend_kw = dict(dict(loc="upper left", bbox_to_anchor=(0.0, 0.79), fontsize=7.5, title_fontsize=8),
                     **(legend_kw or {}))  # по умолчанию — под основной легендой рисунка
    ax.add_artist(ax.legend(handles=handles, title="Стенки ЖТ (п. 8.6), шт. = нар. + вн.", **legend_kw))


def _dim_diam(ax, x_from, x_dim, rr, label, side=1):
    ax.plot([x_from, x_dim], [rr, rr], ":", color="dimgray", lw=0.7)
    ax.annotate("", xy=(x_dim, rr), xytext=(x_dim, 0), arrowprops=dict(arrowstyle="<->", color="dimgray", lw=0.8))
    ax.text(x_dim + side * 3.5, rr * 0.55, label, rotation=90, ha="center", va="center", fontsize=7.5,
            color="dimgray")


def _dim_len(ax, x0, x1, y, label):
    ax.annotate("", xy=(x1, y), xytext=(x0, y), arrowprops=dict(arrowstyle="<->", color="dimgray", lw=0.8))
    ax.text((x0 + x1) / 2, y + 2, label, ha="center", va="bottom", fontsize=7.5, color="dimgray")


def _dimensions(ax, c):
    g, xc, xd, x_end, x_top = c["g"], c["xc"], c["xd"], c["x_end"], c["x_wall"]  # размеры — за задней стенкой
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
    draw_liner_holes(ax, r, c)
    draw_axial_turbine(ax, c["xc"], tvd, 0.25 * tvd["S_rk"], "ТВД", in_dims=False, label_y=0.2 * tvd["d_sa"] / 2)
    x_right = c["x_wall"] + 8 * 8 + 15  # за размерными линиями КС
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
