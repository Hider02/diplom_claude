import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np



def _vector(ax, start, end, color):
    """Вектор скорости со стрелкой от start к end"""
    ax.annotate('', xy=end, xytext=start,
                arrowprops=dict(arrowstyle='-|>', color=color, lw=2.2,
                                mutation_scale=16, shrinkA=0, shrinkB=0))


def _angle_arc(ax, center, theta1, theta2, r, label):
    """Дуга угла между вектором и фронтом решётки (горизонталью) с подписью"""
    ax.add_patch(patches.Arc(center, 2 * r, 2 * r, theta1=theta1, theta2=theta2,
                             color='black', lw=1.2))
    mid = np.radians((theta1 + theta2) / 2)
    ax.text(center[0] + 1.3 * r * np.cos(mid), center[1] + 1.3 * r * np.sin(mid), label,
            fontsize=11, ha='center', va='center')


def _draw_triangle(ax, P, C, W, labels, arcs, legend, title):
    """
    P — общая вершина векторов c и w (сверху), C и W — их концы на фронте решётки.
    u идёт от конца w к концу c (c = w + u) и направлена влево — по вращению.
    """
    P, C, W = (np.asarray(p, dtype=float) for p in (P, C, W))
    _vector(ax, P, C, 'red')     # c
    _vector(ax, P, W, 'green')   # w
    _vector(ax, W, C, 'blue')    # u

    # Подписи векторов — снаружи треугольника
    G = (P + C + W) / 3
    for (a, b), text, color in zip(((P, C), (P, W), (W, C)), labels, ('red', 'green', 'blue')):
        M = (a + b) / 2
        d = M - G
        d = d / (np.hypot(*d) or 1.0)
        ax.annotate(text, xy=M, xytext=(16 * d[0], 16 * d[1]), textcoords='offset points',
                    fontsize=12, color=color, fontweight='bold', ha='center', va='center')

    # Дуги углов в концах векторов на фронте решётки
    c_a = P[1]
    r = 0.3 * c_a
    for center, theta1, theta2, text in arcs:
        _angle_arc(ax, center, theta1, theta2, r, text)

    # Горизонталь — фронт решётки
    xs = [P[0], C[0], W[0]]
    pad = 0.12 * (max(xs) - min(xs) + c_a)
    ax.plot([min(xs) - pad, max(xs) + pad], [0, 0], color='black', lw=0.8)

    for text, color in zip(legend, ('r', 'b', 'g')):
        ax.plot([], [], color=color, lw=2.2, label=text)
    ax.legend(loc='upper center', fontsize=10, framealpha=0.9)

    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(-pad, 1.9 * c_a + pad)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel('Окружное направление, м/с  (← направление вращения)', fontsize=10)
    ax.set_ylabel('Осевое направление, м/с  (поток ↓)', fontsize=10)


def plot_velocity_triangles(c_1, alpha_1_rad, w_1, beta_1_rad,
                            c_2, alpha_2_rad, w_2, beta_2_rad,
                            u_sr, section_name="Среднее сечение"):
    """
    Треугольники скоростей в виде, принятом в отечественной литературе
    (Шестаков, рис. 7; Батурин, рис. 2.2):
    - поток идёт сверху вниз, векторы c и w выходят из общей вершины;
    - окружная скорость u лежит на фронте решётки и направлена влево (по вращению);
    - α1, β1 отсчитываются от направления вращения,
      α2, β2 — от направления, противоположного вращению.
    """
    # Проекции. Окружные: на входе — по вращению, на выходе — против вращения
    c_1u = c_1 * np.cos(alpha_1_rad)
    c_1a = c_1 * np.sin(alpha_1_rad)
    w_1u = w_1 * np.cos(beta_1_rad)
    w_1a = w_1 * np.sin(beta_1_rad)
    c_2u = c_2 * np.cos(alpha_2_rad)
    c_2a = c_2 * np.sin(alpha_2_rad)
    w_2u = w_2 * np.cos(beta_2_rad)
    w_2a = w_2 * np.sin(beta_2_rad)

    a1, b1 = np.degrees(alpha_1_rad), np.degrees(beta_1_rad)
    a2, b2 = np.degrees(alpha_2_rad), np.degrees(beta_2_rad)

    fig, axes = plt.subplots(1, 2, figsize=(15, 7), layout='constrained')

    # Вход: вращение влево, поэтому концы c1 и w1 левее вершины на свои окружные проекции
    C1, W1 = (-c_1u, 0.0), (-w_1u, 0.0)
    _draw_triangle(
        axes[0], P=(0.0, c_1a), C=C1, W=W1,
        labels=('c₁', 'w₁', 'u'),
        arcs=[(C1, 0.0, a1, 'α₁'), (W1, 0.0, b1, 'β₁')],
        legend=(f'c₁ = {c_1:.1f} м/с, α₁ = {a1:.1f}°',
                f'u = {u_sr:.1f} м/с',
                f'w₁ = {w_1:.1f} м/с, β₁ = {b1:.1f}°'),
        title='Треугольник скоростей на входе в РК')

    # Выход: окружные проекции c2 и w2 направлены против вращения — вправо от вершины
    C2, W2 = (c_2u, 0.0), (w_2u, 0.0)
    _draw_triangle(
        axes[1], P=(0.0, c_2a), C=C2, W=W2,
        labels=('c₂', 'w₂', 'u'),
        arcs=[(C2, 180.0 - a2, 180.0, 'α₂'), (W2, 180.0 - b2, 180.0, 'β₂')],
        legend=(f'c₂ = {c_2:.1f} м/с, α₂ = {a2:.1f}°',
                f'u = {u_sr:.1f} м/с',
                f'w₂ = {w_2:.1f} м/с, β₂ = {b2:.1f}°'),
        title='Треугольник скоростей на выходе из РК')

    fig.suptitle(f'Треугольники скоростей осевой турбины — {section_name}',
                 fontsize=14, fontweight='bold')
    plt.show()

    # Вывод числовых значений
    print("\n" + "=" * 70)
    print(f"ТРЕУГОЛЬНИКИ СКОРОСТЕЙ - {section_name}")
    print("=" * 70)
    print(f"\nНА ВХОДЕ В РАБОЧЕЕ КОЛЕСО (углы от направления вращения):")
    print(f"  Абсолютная скорость:         c₁ = {c_1:.1f} м/с")
    print(f"    - окружная составляющая:   c₁ᵤ = {c_1u:.1f} м/с (по вращению)")
    print(f"    - осевая составляющая:     c₁ₐ = {c_1a:.1f} м/с")
    print(f"  Угол абсолютной скорости:    α₁ = {a1:.2f}°")
    print(f"  Относительная скорость:      w₁ = {w_1:.1f} м/с")
    print(f"    - окружная составляющая:   w₁ᵤ = {w_1u:.1f} м/с (по вращению)")
    print(f"    - осевая составляющая:     w₁ₐ = {w_1a:.1f} м/с")
    print(f"  Угол относительной скорости: β₁ = {b1:.2f}°")
    print(f"  Окружная скорость:           u  = {u_sr:.1f} м/с")

    print(f"\nНА ВЫХОДЕ ИЗ РАБОЧЕГО КОЛЕСА (углы от направления, противоположного вращению):")
    print(f"  Абсолютная скорость:         c₂ = {c_2:.1f} м/с")
    print(f"    - окружная составляющая:   c₂ᵤ = {c_2u:.1f} м/с (против вращения)")
    print(f"    - осевая составляющая:     c₂ₐ = {c_2a:.1f} м/с")
    print(f"  Угол абсолютной скорости:    α₂ = {a2:.2f}°")
    print(f"  Относительная скорость:      w₂ = {w_2:.1f} м/с")
    print(f"    - окружная составляющая:   w₂ᵤ = {w_2u:.1f} м/с (против вращения)")
    print(f"    - осевая составляющая:     w₂ₐ = {w_2a:.1f} м/с")
    print(f"  Угол относительной скорости: β₂ = {b2:.2f}°")
    print(f"  Окружная скорость:           u  = {u_sr:.1f} м/с")
    print("=" * 70)


def plot_meridional_section(S_ca, S_rk,
                            D_ca_hub, D_ca_tip,  # СА (втулка и корпус)
                            D_rk_hub, D_rk_tip,  # РК (втулка и корпус на выходе)
                            axial_gap=0.005):
    fig, ax = plt.subplots(figsize=(10, 6))

    # 1. Осевые координаты (Z)
    z_ca_in = 0.0
    z_ca_out = S_ca
    z_gap_in = z_ca_out
    z_gap_out = z_gap_in + axial_gap
    z_rk_in = z_gap_out
    z_rk_out = z_rk_in + S_rk

    # 2. Радиусы (делим диаметры на 2)
    r_ca_hub = D_ca_hub / 2
    r_ca_tip = D_ca_tip / 2

    # На входе в РК диаметры такие же, как на выходе из СА
    r_rk_hub_in = D_ca_hub / 2
    r_rk_tip_in = D_ca_tip / 2

    # На выходе из РК
    r_rk_hub_out = D_rk_hub / 2
    r_rk_tip_out = D_rk_tip / 2

    # 3. СА (прямоугольник)
    sa_verts = [
        (z_ca_in, r_ca_hub), (z_ca_in, r_ca_tip),
        (z_ca_out, r_ca_tip), (z_ca_out, r_ca_hub)
    ]
    ax.add_patch(patches.Polygon(sa_verts, closed=True, linewidth=1.5,
                                 edgecolor='black', facecolor='lightblue', alpha=0.8, label='СА'))

    # 4. Осевой зазор (прямоугольник)
    gap_verts = [
        (z_gap_in, r_ca_hub), (z_gap_in, r_ca_tip),
        (z_gap_out, r_ca_tip), (z_gap_out, r_ca_hub)
    ]
    ax.add_patch(patches.Polygon(gap_verts, closed=True, linewidth=1.0,
                                 edgecolor='gray', facecolor='lightyellow', alpha=0.6, linestyle='--', label='Зазор'))

    # 5. РК (трапеция)
    rk_verts = [
        (z_rk_in, r_rk_hub_in),  # Вход (втулка)
        (z_rk_in, r_rk_tip_in),  # Вход (корпус)
        (z_rk_out, r_rk_tip_out),  # Выход (корпус)
        (z_rk_out, r_rk_hub_out)  # Выход (втулка)
    ]
    ax.add_patch(patches.Polygon(rk_verts, closed=True, linewidth=1.5,
                                 edgecolor='black', facecolor='lightcoral', alpha=0.8, label='РК'))

    # 6. Непрерывные линии втулки и корпуса
    hub_z = [z_ca_in, z_ca_out, z_gap_out, z_rk_out]
    hub_r = [r_ca_hub, r_ca_hub, r_rk_hub_in, r_rk_hub_out]
    ax.plot(hub_z, hub_r, 'k-', linewidth=2.5, label='Втулка')

    tip_z = [z_ca_in, z_ca_out, z_gap_out, z_rk_out]
    tip_r = [r_ca_tip, r_ca_tip, r_rk_tip_in, r_rk_tip_out]
    ax.plot(tip_z, tip_r, 'k-', linewidth=2.5, label='Корпус')

    # 7. Подписи
    ax.text(z_ca_in + S_ca / 2, (r_ca_hub + r_ca_tip) / 2, 'СА', ha='center', va='center', fontsize=12,
            fontweight='bold', color='darkblue')
    ax.text(z_rk_in + S_rk / 2, (r_rk_hub_in + r_rk_tip_in) / 2, 'РК', ha='center', va='center', fontsize=12,
            fontweight='bold', color='darkred')

    # 8. Настройки графика
    ax.set_xlabel('Осевая координата Z, м', fontsize=12)
    ax.set_ylabel('Радиус R, м', fontsize=12)
    ax.set_title('Меридиональное сечение проточной части', fontsize=14, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(loc='upper right')

    all_r = [r_ca_hub, r_ca_tip, r_rk_hub_in, r_rk_tip_in, r_rk_hub_out, r_rk_tip_out]
    ax.set_xlim(z_ca_in - 0.03, z_rk_out + 0.03)
    ax.set_ylim(min(all_r) - 0.015, max(all_r) + 0.015)

    plt.tight_layout()
    plt.show()



