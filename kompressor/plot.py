"""Визуализация ЦБК: меридиональное сечение и треугольники скоростей с лопаточными углами"""
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np


def plot_meridional(r):
    """Меридиональное сечение: ВНА/РК, БЛД (D2…D3), ЛД (D3…D4). Контуры колеса — дуги эллипсов (эскиз)."""
    S, b2 = r["S"], r["b2"]
    r_vt, r1, r2, r3, r4 = r["D_vt"] / 2, r["D1"] / 2, r["D2"] / 2, r["D3"] / 2, r["D4"] / 2
    t = np.linspace(0, math.pi / 2, 60)
    hub_x, hub_r = S * np.sin(t), r_vt + (r2 - r_vt) * (1 - np.cos(t))
    shr_x, shr_r = (S - b2) * np.sin(t), r1 + (r2 - r1) * (1 - np.cos(t))

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.fill(np.concatenate([hub_x, shr_x[::-1]]), np.concatenate([hub_r, shr_r[::-1]]),
            color="lightcoral", alpha=0.6, label="Рабочее колесо")
    ax.fill([S - b2, S, S, S - b2], [r2, r2, r3, r3], color="lightyellow", alpha=0.8, label="БЛД")
    ax.fill([S - b2, S, S, S - b2], [r3, r3, r4, r4], color="lightblue", alpha=0.8, label="ЛД")
    ax.plot(hub_x, hub_r, "k-", lw=2)
    ax.plot(shr_x, shr_r, "k-", lw=2)
    ax.plot([S, S], [r2, r4], "k-", lw=2)
    ax.plot([S - b2, S - b2], [r2, r4], "k-", lw=2)
    ax.plot([-0.3 * S, S * 1.3], [0, 0], "k-.", lw=0.8)  # ось вращения
    for x0, rr, name in ((0, (r_vt, r1), "1"), (S - b2, (r2, r2), "2"), (S - b2, (r3, r3), "3"), (S - b2, (r4, r4), "4")):
        if rr[0] == rr[1]:
            ax.plot([S - b2, S], [rr[0], rr[0]], "k--", lw=0.8)
            ax.text(S + 0.03 * S, rr[0], name, va="center", fontsize=11)
        else:
            ax.plot([x0, x0], rr, "k--", lw=0.8)
            ax.text(x0 - 0.08 * S, (rr[0] + rr[1]) / 2, name, va="center", fontsize=11)
    ax.set_aspect("equal")
    ax.set_xlim(-0.4 * S, 1.5 * S)
    ax.set_ylim(0, r4 * 1.08)
    ax.set_xlabel("Осевая координата, м")
    ax.set_ylabel("Радиус, м")
    ax.set_title(f"Меридиональное сечение ЦБК: D1 = {r['D1'] * 1000:.0f} мм, D2 = {r['D2'] * 1000:.0f} мм, "
                 f"b2 = {b2 * 1000:.1f} мм", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper left")
    plt.tight_layout()
    plt.show()


def _arc(ax, theta1, theta2, rad, label):
    ax.add_patch(patches.Arc((0, 0), 2 * rad, 2 * rad, theta1=theta1, theta2=theta2, color="black", lw=1.1))
    mid = math.radians((theta1 + theta2) / 2)
    ax.text(1.25 * rad * math.cos(mid), 1.25 * rad * math.sin(mid), label, ha="center", va="center", fontsize=11)


def _triangle(ax, c_u, c_m, u, alpha, beta, beta_l, title, m_name):
    """
    Треугольник скоростей: c и w выходят из общей точки, u замыкает треугольник (c = w + u).
    Горизонталь — окружное направление (u вправо), вертикаль — меридиональное (осевое на входе, радиальное на выходе).
    α — угол c, β — угол w, оба от окружного направления; штрих — направление лопатки β_л.
    """
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
    """Треугольники скоростей на входе в ВНА (втулка, средний, периферия) и на выходе из РК"""
    secs = r["sections"]
    chosen = [(secs[0], "втулка"), (secs[len(secs) // 2], "средний диаметр"), (secs[-1], "периферия")]
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), layout="constrained")
    for ax, (s, name) in zip(axes.flat[:3], chosen):
        _triangle(ax, s["cu"], s["ca"], s["u"], s["alpha"], s["beta"], s["beta_l"],
                  f"Вход в ВНА, {name} (D = {s['D'] * 1000:.0f} мм), i = {s['i']:+.1f}°", "Осевая составляющая")
    _triangle(axes.flat[3], r["c2u"], r["c2r"], r["u2"], r["alpha2"], r["beta2"], r["params"]["beta2l"],
              f"Выход из РК (D2 = {r['D2'] * 1000:.0f} мм), μ = {r['mu']:.3f}", "Радиальная составляющая")
    fig.suptitle("Треугольники скоростей ЦБК (углы — от окружного направления)", fontsize=14, fontweight="bold")
    plt.show()
