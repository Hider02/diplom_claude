"""Визуализация ЦБК: меридиональное сечение, треугольники скоростей РК, решётки осевого лопаточного диффузора"""
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


def _dim_diam(ax, x_block, x_dim, r, label, side):
    """Размерная линия диаметра: от оси (y=0) до поверхности (y=r), с линией-выноской от контура.
    side = -1 (линия слева от узла) или +1 (справа) — определяет, куда сдвинута подпись."""
    ax.plot([x_block, x_dim], [r, r], ":", color="dimgray", lw=0.7)
    ax.annotate("", xy=(x_dim, r), xytext=(x_dim, 0),
                arrowprops=dict(arrowstyle="<->", color="dimgray", lw=0.9))
    ax.text(x_dim + side * 0.0035, r / 2, label, rotation=90, ha="center", va="center", fontsize=7.5,
            color="dimgray")


def _dim_width(ax, x0, x1, y0, y1, y_dim, label):
    """Размерная линия осевого размера, с выносками от контура."""
    ax.plot([x0, x0], [y0, y_dim], ":", color="dimgray", lw=0.7)
    ax.plot([x1, x1], [y1, y_dim], ":", color="dimgray", lw=0.7)
    ax.annotate("", xy=(x1, y_dim), xytext=(x0, y_dim),
                arrowprops=dict(arrowstyle="<->", color="dimgray", lw=0.9))
    ax.text((x0 + x1) / 2, y_dim, label, ha="center", va="top", fontsize=7.5, color="dimgray")


def plot_meridional(r):
    """
    Меридиональное сечение: РК, радиальный БЛД (D2…D3), переходное колено (радиальное → осевое),
    осевой ЛД (горизонтальный). Контуры колеса — дуги эллипсов (эскиз).
    """
    S, b2, b3, h4 = r["S"], r["b2"], r["b3"], r["h4"]
    r_vt, r1, r2, r3, r4 = r["D1h"] / 2, r["D1t"] / 2, r["D2"] / 2, r["D3"] / 2, r["D4m"] / 2
    Rb = r["Rb"]
    t = np.linspace(0, math.pi / 2, 60)
    # Колесо: втулка (задний диск) и покрывающая поверхность
    hub_x, hub_r = S * np.sin(t), r_vt + (r2 - r_vt) * (1 - np.cos(t))
    shr_x, shr_r = (S - b2) * np.sin(t), r1 + (r2 - r1) * (1 - np.cos(t))
    # БЛД: канал между x = S − b2 (покрывной) и x = S (задний), сужение b2 → b3
    x_c = S - b2 / 2
    # Колено: средняя линия — четверть окружности радиуса Rb с центром (x_c + Rb, r3)
    phi = np.linspace(0, math.pi / 2, 60)
    b_phi = b3 + (h4 - b3) * phi / (math.pi / 2)
    cx, cr = x_c + Rb * (1 - np.cos(phi)), r3 + Rb * np.sin(phi)
    inner = (cx + b_phi / 2 * np.cos(phi), cr - b_phi / 2 * np.sin(phi))
    outer = (cx - b_phi / 2 * np.cos(phi), cr + b_phi / 2 * np.sin(phi))
    x4 = x_c + Rb
    # Осевой ЛД
    n_rows = len(r["rows"])
    row_len, gap = 1.3 * h4, 0.4 * h4
    L_ld = n_rows * row_len + (n_rows + 1) * gap
    r5_in, r5_out = r4 - r["h5"] / 2, r4 + r["h5"] / 2

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.fill(np.concatenate([hub_x, shr_x[::-1]]), np.concatenate([hub_r, shr_r[::-1]]),
            color="lightcoral", alpha=0.6, label="Рабочее колесо")
    ax.fill([S - b2, S, x_c + b3 / 2, x_c - b3 / 2], [r2, r2, r3, r3], color="lightyellow", label="Радиальный БЛД")
    ax.fill(np.concatenate([inner[0], outer[0][::-1]]), np.concatenate([inner[1], outer[1][::-1]]),
            color="honeydew", label="Переходное колено")
    ax.fill([x4, x4 + L_ld, x4 + L_ld, x4], [r4 - h4 / 2, r5_in, r5_out, r4 + h4 / 2], color="lightblue",
            alpha=0.8, label="Осевой ЛД")
    for j in range(n_rows):
        x0 = x4 + gap + j * (row_len + gap)
        ax.add_patch(patches.Rectangle((x0, r4 - h4 / 2), row_len, h4, color="steelblue", alpha=0.8))
        ax.text(x0 + row_len / 2, r4 + h4 * 0.9, f"{j + 1}", ha="center", fontsize=9)
    for xs, rs in ((hub_x, hub_r), (shr_x, shr_r), inner, outer):
        ax.plot(xs, rs, "k-", lw=1.8)
    ax.plot([S - b2, x_c - b3 / 2], [r2, r3], "k-", lw=1.8)
    ax.plot([S, x_c + b3 / 2], [r2, r3], "k-", lw=1.8)
    ax.plot([x4, x4 + L_ld], [r4 - h4 / 2, r5_in], "k-", lw=1.8)
    ax.plot([x4, x4 + L_ld], [r4 + h4 / 2, r5_out], "k-", lw=1.8)
    ax.plot([-0.3 * S, x4 + L_ld * 1.2], [0, 0], "k-.", lw=0.8)
    for x, rr, name in ((0, (r_vt + r1) / 2, "1"), (S + b2, r2, "2"), (S + b2, r3, "3"), (x4, r4 - h4, "4"),
                        (x4 + L_ld, r4 - h4, "5")):
        ax.text(x, rr, name, fontsize=11, fontweight="bold")

    # Размерные линии: диаметры на входе (D1) и на выходе из осевого ЛД, осевая длина колеса и ЛД
    dx1, dx2 = 0.16 * r1, 0.34 * r1
    off = 0.12 * r1
    _dim_diam(ax, 0.0, -dx1, r_vt, f"⌀{2 * r_vt * 1000:.0f}", -1)
    _dim_diam(ax, 0.0, -dx2, r1, f"⌀{2 * r1 * 1000:.0f}", -1)
    _dim_diam(ax, x4 + L_ld, x4 + L_ld + dx1, r5_in, f"⌀{2 * r5_in * 1000:.0f}", 1)
    _dim_diam(ax, x4 + L_ld, x4 + L_ld + dx2, r5_out, f"⌀{2 * r5_out * 1000:.0f}", 1)
    _dim_width(ax, 0.0, S, r1, r2, r2 + off, f"S={S * 1000:.0f}")
    _dim_width(ax, x4, x4 + L_ld, r4 + h4 / 2, r4 + h4 / 2, r4 + h4 / 2 + off, f"L={L_ld * 1000:.0f}")

    ax.set_aspect("equal")
    ax.set_xlim(min(-0.3 * S, -dx2 - 0.010), max(x4 + L_ld * 1.25, x4 + L_ld + dx2 + 0.010))
    ax.set_ylim(0, max(r4 + max(h4, r["h5"]) * 3, r2 + off + 0.015, r4 + h4 / 2 + off + 0.015))
    ax.set_xlabel("Осевая координата, м")
    ax.set_ylabel("Радиус, м")
    ax.set_title(f"Меридиональное сечение ЦБК: D2 = {r['D2'] * 1000:.0f} мм, b2 = {b2 * 1000:.1f} мм, "
                 f"D_ЛД = {r['D4m'] * 1000:.0f} мм", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="lower right")
    plt.tight_layout()
    plt.show()


def _arc(ax, theta1, theta2, rad, label):
    ax.add_patch(patches.Arc((0, 0), 2 * rad, 2 * rad, theta1=theta1, theta2=theta2, color="black", lw=1.1))
    mid = math.radians((theta1 + theta2) / 2)
    ax.text(1.25 * rad * math.cos(mid), 1.25 * rad * math.sin(mid), label, ha="center", va="center", fontsize=11)


def _triangle(ax, c_u, c_m, u, alpha, beta, beta_l, title, m_name):
    """c и w из общей точки, u замыкает треугольник; углы от окружного направления; штрих — лопатка β_л"""
    C, W = (c_u, c_m), (c_u - u, c_m)
    for (x0, y0), (x1, y1), color in (((0, 0), C, "red"), ((0, 0), W, "green"), (W, C, "blue")):
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=2.2, mutation_scale=16, shrinkA=0, shrinkB=0))
    L = 1.15 * math.hypot(*W)
    ax.plot([0, -L * math.cos(math.radians(beta_l))], [0, L * math.sin(math.radians(beta_l))], "--", color="gray",
            lw=1.3, label=f"лопатка β_л = {beta_l:.1f}°")
    rad = 0.28 * c_m
    if alpha < 89.9:
        _arc(ax, 0, alpha, rad, "α")
    _arc(ax, 180 - beta, 180, rad * 1.3, "β")
    ax.plot([], [], color="red", lw=2.2, label=f"c = {math.hypot(c_u, c_m):.0f} м/с, α = {alpha:.1f}°")
    ax.plot([], [], color="green", lw=2.2, label=f"w = {math.hypot(*W):.0f} м/с, β = {beta:.1f}°")
    ax.plot([], [], color="blue", lw=2.2, label=f"u = {u:.0f} м/с")
    xs = [0, C[0], W[0]]
    pad = 0.1 * (max(xs) - min(xs) + c_m)
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(-pad, c_m * 1.9)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlabel("Окружное направление, м/с (→ вращение)")
    ax.set_ylabel(f"{m_name}, м/с")


def plot_velocity_triangles(r):
    """Треугольники скоростей на входе в РК (втулка, средний, периферия) и на выходе из РК"""
    secs = r["sections"]
    chosen = [(secs[0], "втулка"), (secs[len(secs) // 2], "средний диаметр"), (secs[-1], "периферия")]
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), layout="constrained")
    for ax, (s, name) in zip(axes.flat[:3], chosen):
        _triangle(ax, s["cu"], s["ca"], s["u"], s["alpha"], s["beta"], s["beta_l"],
                  f"Вход в РК, {name} (D = {s['D'] * 1000:.0f} мм), i = {s['i']:+.1f}°", "Осевая составляющая")
    _triangle(axes.flat[3], r["c2u"], r["c2r"], r["u2"], r["alpha2"], r["beta2"], r["params"]["beta2l"],
              f"Выход из РК (D2 = {r['D2'] * 1000:.0f} мм), σ = {r['slip']:.3f}", "Радиальная составляющая")
    fig.suptitle("Треугольники скоростей рабочего колеса (углы — от окружного направления)", fontsize=14,
                 fontweight="bold")
    plt.show()


def plot_diffuser_cascades(r):
    """
    Развёртка венцов осевого ЛД на среднем диаметре: средние линии лопаток (дуги окружности между лопаточными
    углами входа и выхода) и направления потока. Ось x — осевое направление, y — окружное. Углы — от оси.
    """
    rows = r["rows"]
    fig, axes = plt.subplots(1, len(rows), figsize=(4 * len(rows), 6), layout="constrained")
    axes = np.atleast_1d(axes)
    for ax, row, j in zip(axes, rows, range(1, len(rows) + 1)):
        chord = 1.0
        pitch = chord / row["sigma"]
        a1, a2 = math.radians(row["blade_in"]), math.radians(row["blade_out"])
        s = np.linspace(0, 1, 50)
        ang = a1 + (a2 - a1) * s  # угол касательной к средней линии меняется линейно по длине (дуга окружности)
        dx, dy = np.cos(ang), np.sin(ang)
        xs = np.concatenate([[0], np.cumsum(dx[:-1] * np.diff(s))])
        ys = np.concatenate([[0], np.cumsum(dy[:-1] * np.diff(s))])
        L = math.hypot(xs[-1], ys[-1])
        xs, ys = xs / L * chord, ys / L * chord
        for kk in range(3):
            ax.plot(xs, ys + kk * pitch, "k-", lw=3)
        # Поток: на входе — к передней кромке средней лопатки, на выходе — от задней кромки
        a_in, a_out, La = math.radians(row["a_in"]), math.radians(row["a_out"]), 0.5
        le, te = (0.0, pitch), (xs[-1], ys[-1] + pitch)
        ax.annotate("", xy=le, xytext=(le[0] - La * math.cos(a_in), le[1] - La * math.sin(a_in)),
                    arrowprops=dict(arrowstyle="-|>", color="red", lw=2))
        ax.annotate("", xy=(te[0] + La * math.cos(a_out), te[1] + La * math.sin(a_out)), xytext=te,
                    arrowprops=dict(arrowstyle="-|>", color="red", lw=2))
        ax.text(le[0] - La, le[1] - La * math.sin(a_in) - 0.15, f"α_вх = {row['a_in']:.1f}°", color="red", fontsize=9)
        ax.text(te[0], te[1] + La * math.sin(a_out) + 0.1, f"α_вых = {row['a_out']:.1f}°", color="red", fontsize=9)
        ax.set_xlim(-La - 0.3, xs[-1] + La + 0.3)
        ax.set_ylim(-0.3, 2 * pitch + ys[-1] + 0.5)
        ax.set_title(f"Венец {j}: θ = {row['camber']:.1f}°, δ = {row['dev']:.1f}°,\n"
                     f"лопатка {row['blade_in']:.1f}° → {row['blade_out']:.1f}°, D_eq = {row['D_eq']:.2f}",
                     fontsize=10)
        ax.set_aspect("equal")
        ax.set_xlabel("Осевое направление")
        ax.set_ylabel("Окружное направление")
        ax.grid(True, alpha=0.3, linestyle="--")
    fig.suptitle("Решётки осевого лопаточного диффузора на среднем диаметре (углы — от оси)", fontsize=13,
                 fontweight="bold")
    plt.show()
