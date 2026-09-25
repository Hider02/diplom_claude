import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'FreeSans']
plt.rcParams['axes.unicode_minus'] = False


def plot_comp(S, r_vt, r1, r2, h1, r3=None, h3=None, r4=None, h4=None, h2=None):
    fig, ax = plt.subplots(figsize=(6, 6))

    # Ось вращения
    ax.plot([-0.15 * S, S * 1.15], [0, 0], 'k-')

    # === РАБОЧЕЕ КОЛЕСО ===
    ax.plot([0, S], [r_vt, r_vt], 'k-')  # втулка
    ax.plot([0, S - h2], [r1, r1], 'k-')  # обод (горизонтальный)
    ax.plot([S - h2, S - h2], [r1, r2], 'k-')  # обод (вертикальный переход)
    ax.plot([S - h2, S], [r2, r2], 'k-')  # обод (выходной участок)
    ax.plot([0, 0], [r_vt, r1], 'k-')  # входное сечение
    ax.plot([S, S], [r_vt, r2], 'k-')  # выходное сечение РК

    # === ДИФФУЗОРЫ (начинаются строго от x = S) ===
    if r3 is not None and r3 > r2:
        # БЛД: радиальное расширение от r2 до r3
        ax.plot([S, S], [r2, r3], 'k-')
        # Условная верхняя граница канала БЛД
        ax.plot([S * 0.98, S * 1.02], [r3, r3], 'k--', alpha=0.4)

        if r4 is not None and r4 > r3:
            # ЛД: продолжение радиально от r3 до r4
            ax.plot([S, S], [r3, r4], 'k-')
            ax.plot([S * 0.98, S * 1.02], [r4, r4], 'k--', alpha=0.4)

    # Настройки
    ax.set_aspect('equal')
    ax.set_xlim(-0.15 * S, S * 1.2)
    ax.set_ylim(0, max(r2, r3 or 0, r4 or 0) * 1.1)
    ax.set_xlabel('x, м')
    ax.set_ylabel('r, м')
    ax.grid(True, linestyle=':', alpha=0.3)

    plt.tight_layout()
    plt.show()