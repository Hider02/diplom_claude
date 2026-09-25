"""
Вспомогательные расчёты по методике Батурина О.В. «Проектирование свободной турбины в составе
многокаскадной турбины ГТД» (СГАУ, 2013). Номера пунктов (п. NN) — по пособию.
Используется в turbine/gasdynamiks.py (ТВД) и free_turbine/gasdynamiks.py (СТ).
"""
import math


def mu_gas(T):
    """Динамическая вязкость продуктов сгорания, Па·с (п. 57)"""
    return 44.3e-6 * (T / 1073) ** 0.678


def velocity_coefficient(angle_in_deg, angle_out_deg, t, b, h, r_edge, lambda_s, lambda_s_opt, Re):
    """
    Коэффициент скорости венца (φ для СА, ψ для РК) по модели потерь ЦИАМ (п. 57…64 для СА, п. 65…72 для РК).
    angle_in_deg, angle_out_deg — углы входа и выхода потока (α0, α1 для СА; β1, β2 для РК), град;
    t — шаг, b — хорда, h — высота лопатки, r_edge — радиус выходной кромки, м;
    lambda_s — приведённая изоэнтропическая скорость на выходе, lambda_s_opt — её оптимальное значение
    (0,9 для СА, 0,8 для РК); Re — число Рейнольдса по хорде.
    Возвращает (коэффициент скорости, словарь составляющих потерь).
    """
    a_in = math.radians(angle_in_deg)
    a_out = math.radians(angle_out_deg)

    # Профильные потери трения в автомодельной области (п. 58)
    x1 = min(angle_in_deg + angle_out_deg, 110.0)
    y1 = min(math.sin(a_in) / math.sin(a_out), 1.7)
    xi_tr = 0.02185 * (0.01065 * x1 ** 2 - 2.295 * x1 + 160.5) * (0.1055 * y1 ** 2 - 0.3427 * y1 + 0.295)
    # Кромочные потери (п. 59)
    xi_kr = 0.4 * r_edge / (t * math.sin(a_out))
    # Поправка на число Рейнольдса при Re <= 10^6 (п. 60)
    d_xi_re = max(2100 / Re - 0.0021, 0.0) if Re <= 1e6 else 0.0
    # Поправка на уровень скорости (п. 61). В пособии напечатано «−3,87x»: тогда при x = 1 поправка
    # равна −7,74·ξ и потери становятся отрицательными. Принят знак «+3,87x» (РТМ 1509-75):
    # поправка ≈ 0 при оптимальной скорости и растёт при отклонении от неё.
    x = lambda_s / lambda_s_opt
    d_xi_lambda = (3.61 * x ** 3 - 7.22 * x ** 2 + 3.87 * x - 0.26) * (xi_tr + xi_kr + d_xi_re)
    # Вторичные потери (п. 62)
    xi_vt = 2 * (xi_tr + xi_kr + d_xi_re + d_xi_lambda) * t * math.sin(a_out) / h

    xi_sum = xi_tr + xi_kr + xi_vt + d_xi_re + d_xi_lambda
    parts = {"тр": xi_tr, "кр": xi_kr, "вт": xi_vt, "Re": d_xi_re, "λ": d_xi_lambda, "Σ": xi_sum}
    return math.sqrt(1 - xi_sum), parts


def zweifel_blade_number(D, angle_in_deg, angle_out_deg, S, C_L):
    """Потребное число лопаток по числу Цвайфеля (п. 49, 53), нецелое"""
    a_in = math.radians(angle_in_deg)
    a_out = math.radians(angle_out_deg)
    return 2 * math.pi * D * math.sin(a_out) ** 2 * (1 / math.tan(a_out) + 1 / math.tan(a_in)) / (S * C_L)


def zweifel_number(D, angle_in_deg, angle_out_deg, S, Z):
    """Фактическое число Цвайфеля при выбранном числе лопаток Z"""
    return zweifel_blade_number(D, angle_in_deg, angle_out_deg, S, 1.0) / Z


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


def prime_blade_number(z, exclude=()):
    """
    Ближайшее к z простое число лопаток (п. 49: «стремиться к простым и нечётным числам», не кратным
    числу лопаток соседних венцов). При равном удалении берётся большее — меньше нагрузка на профиль.
    """
    z0 = max(round(z), 3)
    for d in range(0, z0):
        for cand in (z0 + d, z0 - d):
            if is_prime(cand) and cand not in exclude:
                return cand
    return z0


def larson_miller(T_blade, tau_hours):
    """Параметр длительной прочности P·10³ для рис. 2.7 (ось абсцисс): P = T·(20,3 + lg τ)·10⁻³.
    Константа 20,3 подобрана по примеру пособия (T = 1050 К, τ = 1000 ч → P ≈ 24,5)."""
    return T_blade * (20.3 + math.log10(tau_hours)) * 1e-3


def form_factor(F_ratio):
    """Коэффициент формы лопатки К_ф = (F_к/F_вт + 1)/2 (п. 29), F_к/F_вт = 0,25…0,55"""
    return (F_ratio + 1) / 2


def blade_root_stress(u_sr, D_sr_to_h, K_f, rho_mat):
    """Напряжение растяжения от центробежных сил в корне лопатки РК, Па (п. 30, 32):
    σ = 2·u²·К_ф·ρ_м / (D_ср/h)"""
    return 2 * u_sr ** 2 * K_f * rho_mat / D_sr_to_h


def D_to_h_min(u_sr, K_f, rho_mat, sigma_dl, k_m):
    """Минимально допустимое по прочности отношение (D_ср/h)_2 (п. 30)"""
    return 2 * u_sr ** 2 * K_f * rho_mat / (sigma_dl / k_m)


def d_eta_tip_clearance(delta, h2, D2_sr, beta_2_rad, t_rk, b_rk):
    """Потери КПД от перетекания в радиальном зазоре без бандажа (п. 76)"""
    return (delta / h2) * (1 + h2 / D2_sr) * (1 + 0.3 / math.sin(beta_2_rad) * t_rk / b_rk)


def d_eta_disk(G, L_ts, D1_hub, u1_hub, rho_sr):
    """Потери КПД на трение диска (п. 77)"""
    return 1701 / (G * L_ts) * D1_hub ** 2 * (u1_hub / 100) ** 3 * rho_sr


def blade_temperature(T_w, T_blade_allow, T_blade_cooled):
    """Температура рабочей лопатки (п. 27): T_л = T*_w1, если она не выше допустимой,
    иначе лопатка охлаждается и T_л задаётся из диапазона 950…1100 К"""
    return T_w if T_w <= T_blade_allow else T_blade_cooled


def print_strength_check(name, T_w, T_blade, T_blade_allow, tau, material, rho_mat, sigma_dl, K_f, u_sr,
                         D_sr_to_h, k_m_min):
    """Проверка температуры и прочности рабочей лопатки (п. 27…32). Возвращает запас k_м."""
    print(f"\n--- Прочность рабочей лопатки {name} (Батурин, п. 27…32) ---")
    print(f"T*_w1 = {T_w:.0f} К (допустимая температура лопатки без охлаждения {T_blade_allow:.0f} К)")
    if T_w > T_blade_allow:
        print(f"⚠️ T*_w1 > {T_blade_allow:.0f} К: рабочую лопатку необходимо охлаждать (п. 27); "
              f"принято T_л = {T_blade:.0f} К. При охлаждении φ и ψ снижаются (п. 11, 19)")
    sigma_p = blade_root_stress(u_sr, D_sr_to_h, K_f, rho_mat)
    k_m = sigma_dl / sigma_p
    dh_min = D_to_h_min(u_sr, K_f, rho_mat, sigma_dl, k_m_min)
    print(f"Материал {material}, T_л = {T_blade:.0f} К, τ = {tau} ч, P·10³ = {larson_miller(T_blade, tau):.1f} "
          f"→ σ_дл = {sigma_dl / 1e6:.0f} МПа (по рис. 2.7)")
    print(f"Напряжение растяжения в корне σ_р = {sigma_p / 1e6:.0f} МПа (К_ф = {K_f:.3f})")
    print(f"Запас прочности k_м = σ_дл/σ_р = {k_m:.2f} (требуется ≥ {k_m_min})")
    print(f"D_ср/h = {D_sr_to_h:.2f}, минимально допустимое (D_ср/h)_min = {dh_min:.2f}")
    if k_m < k_m_min:
        print("⚠️ Запас прочности недостаточен: уменьшите u_ср (n или D_ср), увеличьте D_ср/h "
              "или выберите более прочный материал (ЖС30, ЖС32)")
    return k_m


# ---------------------------------------------------------------------------
# Углы лопаток. [3] Кожухов Ю.В. и др. «Расчет авиационных газотурбинных двигателей», СПбПУ, 2020, гл. 4
# ---------------------------------------------------------------------------
def deviation_angle(lam, angle_out_deg):
    """Угол отставания потока в косом срезе решётки, град ([3], ф. 4.37 для СА, 4.69 для РК).
    lam — приведённая скорость на выходе (λ_c1 для СА, λ_w2 для РК), angle_out_deg — угол выхода потока."""
    a = angle_out_deg
    return (0.0213 * lam * a * (a - 54) - 0.0255 * lam ** 2 * a * (a - 46.5)
            + 14 * lam * (1 - lam) + 0.19 * (a - 1.21))


def stagger_angle(angle_in_deg, angle_out_deg):
    """Угол установки профиля в решётке, град (Батурин, п. 51, 55)"""
    d = angle_in_deg - angle_out_deg
    return 85.95 - 1.277 * d + 0.0084 * d ** 2


def print_blade_angles(sections_ca, sections_rk, alpha_0_deg, alpha_1_deg, lambda_c1_arr, beta1_arr, beta2_arr,
                       lambda_w2_arr, S_sa, S_rk, Z_ca, Z_rk, i_ca, i_rk):
    """
    Лопаточные углы и углы установки по сечениям.
    Вход: угол лопатки = угол потока + угол атаки (i = β1л − β1, отсчёт от фронта решётки);
    рекомендуемый расчётный угол атаки небольшой отрицательный, −5…0°.
    Выход: эффективный угол лопатки = угол потока − угол отставания ([3], ф. 4.38, 4.70).
    Угол установки γ и хорда b = S/sin γ — по углам потока (Батурин, п. 51…56).
    Углы входа отсчитываются от направления вращения, выхода — от противоположного (как в таблицах выше).
    """
    print(f"\n--- 3.2. Углы лопаток и углы установки по сечениям (угол атаки i_СА = {i_ca:+.1f}°, "
          f"i_РК = {i_rk:+.1f}°) ---")
    print(f"\n{'Сечение':<15} | {'D, мм':<7} | {'α0л, °':<7} | {'α1, °':<7} | {'δα, °':<6} | {'α1л, °':<7} | "
          f"{'γ, °':<6} | {'b, мм':<6} | {'t/b':<5}")
    print("-" * 95)
    for sec, lam in zip(sections_ca, lambda_c1_arr):
        d_a = deviation_angle(lam, alpha_1_deg)
        gamma = stagger_angle(alpha_0_deg, alpha_1_deg)
        b = S_sa / math.sin(math.radians(gamma))
        t = math.pi * sec["D"] / Z_ca
        print(f"{sec['name']:<15} | {sec['D'] * 1000:<7.1f} | {alpha_0_deg + i_ca:<7.2f} | {alpha_1_deg:<7.2f} | "
              f"{d_a:<6.2f} | {alpha_1_deg - d_a:<7.2f} | {gamma:<6.2f} | {b * 1000:<6.1f} | {t / b:<5.3f}")

    print(f"\n{'Сечение':<15} | {'D, мм':<7} | {'β1, °':<7} | {'β1л, °':<7} | {'β2, °':<7} | {'δβ, °':<6} | "
          f"{'β2л, °':<7} | {'θ, °':<6} | {'γ, °':<6} | {'b, мм':<6} | {'t/b':<5}")
    print("-" * 110)
    beta1_bl, beta2_bl = [], []
    for sec, b1, b2, lam in zip(sections_rk, beta1_arr, beta2_arr, lambda_w2_arr):
        d_b = deviation_angle(lam, b2)
        b1l, b2l = b1 + i_rk, b2 - d_b
        theta = 180 - b1l - b2l  # угол изгиба профиля
        gamma = stagger_angle(b1, b2)
        b = S_rk / math.sin(math.radians(gamma))
        t = math.pi * sec["D"] / Z_rk
        beta1_bl.append(b1l)
        beta2_bl.append(b2l)
        print(f"{sec['name']:<15} | {sec['D'] * 1000:<7.1f} | {b1:<7.2f} | {b1l:<7.2f} | {b2:<7.2f} | {d_b:<6.2f} | "
              f"{b2l:<7.2f} | {theta:<6.2f} | {gamma:<6.2f} | {b * 1000:<6.1f} | {t / b:<5.3f}")
    print("β1л = β1 + i; β2л = β2 − δβ (эффективный угол выхода); θ = 180° − β1л − β2л — угол изгиба профиля")
    for name, i in (("СА", i_ca), ("РК", i_rk)):
        if i > 0:
            print(f"⚠️ Угол атаки {name} i = {i:+.1f}° > 0: на нерасчётных режимах возможен отрыв на спинке "
                  f"(рекомендуется −5…0°)")
        elif i < -10:
            print(f"⚠️ Угол атаки {name} i = {i:+.1f}° слишком мал: возможен отрыв у входной кромки со стороны корыта "
                  f"(рекомендуется −5…0°)")
    return beta1_bl, beta2_bl
