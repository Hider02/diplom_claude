"""
Визуализация меридионального сечения противоточной КС (схематично).

Условности, принятые при построении (там, где методика задаёт только характерные диаметры/
площади, а не полный контур):
  - купол (переход от входного диффузора к жаровой трубе) и переходный участок на выходе
    из разворота — плавные скруглённые вставки схематичной (не расчётной) длины;
  - разворот потока на 180° (газосборник) показан одной дугой средней линии радиусом R_г_ср
    с шириной канала, линейно меняющейся от h_ж (вход в разворот) до h_вых (выход) — методика
    даёт ещё R_г_н и R_г_вн (радиусы по наружной/внутренней стенке жаровой трубы), но они
    сходятся в одну точку на выходе (в источнике не detalized дальше) и для наглядной картинки
    не используются;
  - наружный канал (F_кк, между корпусом D_ксн и наружной стенкой жаровой трубы) показан вдоль
    прямого участка жаровой трубы как отдельная кольцевая полость.
"""
import math
import matplotlib.pyplot as plt
import numpy as np

MM = 1000.0


def plot_combustor(r):
    D_k, h_k = r["D_k"] * MM, r["h_k"] * MM
    D_vyh, h_vyh = r["D_vyh"] * MM, r["h_vyh"] * MM
    D_ksn = r["D_ksn"] * MM
    D_zhtr_n, D_zhtr_vn, D_zhtr_sr = r["D_zhtr_n"] * MM, r["D_zhtr_vn"] * MM, r["D_zhtr_sr"] * MM
    h_zh = r["h_zh"] * MM
    R_g_sr = r["R_g_sr"] * MM
    l_zg, l_d1, Delta = r["l_zg"] * MM, r["l_d1"] * MM, r["Delta"] * MM
    h_d1 = r["h_d1"] * MM

    r_k_hub, r_k_tip = (D_k - h_k) / 2, (D_k + h_k) / 2
    r_d1_hub, r_d1_tip = (D_k - h_d1) / 2, (D_k + h_d1) / 2
    r_zh_vn, r_zh_n = D_zhtr_vn / 2, D_zhtr_n / 2
    r_vyh_hub, r_vyh_tip = (D_vyh - h_vyh) / 2, (D_vyh + h_vyh) / 2
    r_ksn = D_ksn / 2

    fig, ax = plt.subplots(figsize=(11, 8))

    # 1. Малый «безотрывный» диффузор: 0 -> l_d1, расширение от h_к до h_д1 (центр на D_к/2)
    x0, x1 = 0.0, l_d1
    ax.fill([x0, x1, x1, x0], [r_k_hub, r_d1_hub, r_d1_tip, r_k_tip], color="lightyellow",
            edgecolor="black", lw=1.2, label="Вх. диффузор")

    # 2. Зазор Δ до форсунок: постоянная высота h_д1
    x2 = x1 + Delta
    ax.fill([x1, x2, x2, x1], [r_d1_hub, r_d1_hub, r_d1_tip, r_d1_tip], color="honeydew",
            edgecolor="gray", lw=1.0, linestyle="--", label="Зазор Δ (до форсунок)")

    # 3. Купол: плавный переход от (h_д1 около D_к/2) к стенкам жаровой трубы (Джтрн/Джтрвн)
    dome_len = 0.5 * h_zh
    t = np.linspace(0, math.pi / 2, 40)
    dome_out_x = x2 + dome_len * np.sin(t)
    dome_out_r = r_d1_tip + (r_zh_n - r_d1_tip) * (1 - np.cos(t))
    dome_in_x = x2 + dome_len * np.sin(t)
    dome_in_r = r_d1_hub + (r_zh_vn - r_d1_hub) * (1 - np.cos(t))
    ax.fill(np.concatenate([dome_in_x, dome_out_x[::-1]]), np.concatenate([dome_in_r, dome_out_r[::-1]]),
            color="wheat", edgecolor="black", lw=1.2, label="Купол")
    x3 = x2 + dome_len

    # 4. Жаровая труба (зона горения l_зг): кольцевой канал Джтрвн…Джтрн
    x4 = x3 + l_zg
    ax.fill([x3, x4, x4, x3], [r_zh_vn, r_zh_vn, r_zh_n, r_zh_n], color="mistyrose",
            edgecolor="black", lw=1.4, label="Жаровая труба (зона горения)")

    # Наружный кольцевой канал F_кк (между корпусом и наружной стенкой жаровой трубы),
    # вдоль купола и зоны горения
    ax.fill(np.concatenate([dome_out_x, [x4, x3]]), np.concatenate([dome_out_r, [r_ksn, r_ksn]]),
            color="lightblue", alpha=0.5, edgecolor="none", label="Наружный канал (F_кк)")

    # 5. Газосборник: разворот потока на 180° по средней линии радиусом R_г_ср,
    # ширина канала линейно меняется от h_ж до h_вых (см. docstring)
    phi = np.linspace(0, math.pi, 80)
    center_r = D_zhtr_sr / 2 - R_g_sr
    cx = x4 + R_g_sr * np.sin(phi)
    cr = center_r + R_g_sr * np.cos(phi)
    b_phi = h_zh + (h_vyh - h_zh) * phi / math.pi
    outer_x, outer_r = cx + b_phi / 2 * np.sin(phi), cr + b_phi / 2 * np.cos(phi)
    inner_x, inner_r = cx - b_phi / 2 * np.sin(phi), cr - b_phi / 2 * np.cos(phi)
    ax.fill(np.concatenate([inner_x, outer_x[::-1]]), np.concatenate([inner_r, outer_r[::-1]]),
            color="lightcoral", alpha=0.6, edgecolor="black", lw=1.2, label="Газосборник (разворот 180°)")

    x5 = x4  # разворот не даёт смещения по оси (см. docstring), выход — там же по x, что и вход в разворот

    # Корпус КС (наружный диаметр D_ксн) — контур вдоль наружного канала
    ax.plot([x2, x4], [r_ksn, r_ksn], "k-", lw=1.8)
    ax.plot([x2, x2], [r_d1_tip, r_ksn], "k--", lw=1.0)

    # Контуры стенок
    for xs, rs in ((dome_out_x, dome_out_r), (dome_in_x, dome_in_r)):
        ax.plot(xs, rs, "k-", lw=1.4)
    ax.plot([x3, x4], [r_zh_n, r_zh_n], "k-", lw=1.4)
    ax.plot([x3, x4], [r_zh_vn, r_zh_vn], "k-", lw=1.4)
    ax.plot(outer_x, outer_r, "k-", lw=1.4)
    ax.plot(inner_x, inner_r, "k-", lw=1.4)
    ax.plot([x0, x1], [r_k_hub, r_d1_hub], "k-", lw=1.4)
    ax.plot([x0, x1], [r_k_tip, r_d1_tip], "k-", lw=1.4)

    ax.plot([-0.15 * D_ksn, x5 + R_g_sr * 1.3], [0, 0], "k-.", lw=0.9)

    ax.text(x0, r_k_tip + 6, "к", fontsize=11, fontweight="bold", ha="center")
    ax.text(x5, r_vyh_hub - 10, "вых", fontsize=11, fontweight="bold", ha="center")

    ax.set_aspect("equal")
    ax.set_xlim(-0.15 * D_ksn, x5 + R_g_sr * 1.5)
    ax.set_ylim(-0.05 * r_ksn, r_ksn * 1.15)
    ax.set_xlabel("Осевая координата, мм")
    ax.set_ylabel("Радиус, мм")
    ax.set_title(f"Меридиональное сечение противоточной КС: D_ксн = {D_ksn:.0f} мм, "
                f"Джтрн/вн = {D_zhtr_n:.0f}/{D_zhtr_vn:.0f} мм", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper right", fontsize=8.5)
    plt.tight_layout()

    out_path = __file__.replace("plot.py", "combustor_layout.png")
    fig.savefig(out_path, dpi=150)
    print(f"\nСохранено: {out_path}")
    plt.show()
