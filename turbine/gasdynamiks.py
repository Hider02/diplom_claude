import math
from turbine import plot_turb as p
from turbine import baturin as bt
from termo.termodynamiks import (k_g, R_g, cp_g, T_g, alpha, eta_tvd, g_gas, P3, N_k, P4, L_tvd, T4,
                                     eta_mex, l_agr, komp)

# Методика: [1] ЦБК+ОС+ЦБТ, п. 1.4 «Газодинамический расчет осевой турбины» (с. 31…37) — расчёт на среднем
# диаметре; [2] Батурин О.В. «Проектирование свободной турбины…» — уточнение φ и ψ по модели потерь ЦИАМ
# (п. 57…72), число лопаток (п. 49…56), прочность (п. 27…32), расчёт по высоте лопатки (п. 79…112).

# ==========================================
# 1. ИСХОДНЫЕ ДАННЫЕ (из расчета компрессора)
# ==========================================
# ТВД и ЦБК на одном валу: частота вращения и диаметр колеса — из газодинамического расчёта ЦБК (kompressor/cbk.py)
n_tk = komp["n_rpm"] if komp else 45000  # Частота вращения ротора, об/мин
D_k = komp["D2"] if komp else 0.23  # Наружный диаметр колеса компрессора, м

# Параметры газа (продукты сгорания)
k_ratio = k_g / (k_g - 1)  # k/(k-1) для газа (k_g согласован с cp_g, см. termodynamiks.py)
cp_prime = cp_g(T_g, alpha)  # Изобарная теплоемкость газа, Дж/(кг·К)

# Предварительно принимаемые коэффициенты
eta_i_guess = eta_tvd  # Начальное приближение внутреннего КПД турбины
phi = 0.98  # Коэффициент скорости соплового аппарата (начальное приближение, [1]: 0,95…0,98)
psi = 0.96  # Коэффициент скорости рабочего колеса (начальное приближение, [1]: 0,93…0,97)
alpha_1_deg = 20.0  # Угол выхода газа из сопел, град
alpha_1_rad = math.radians(alpha_1_deg)
rho_st = 0.5  # Степень реактивности
alpha_0_deg = 90  # Угол входа в СА
# Угол beta_2 не задаётся, а находится из уравнения неразрывности ([1], с. 35)
# при высоте лопатки РК l_рк = l_c + (0…2 мм), чтобы канал не сужался по тракту
dl_rk = 0.001  # Превышение высоты лопатки РК над высотой лопатки СА, м

C_l_CA = 0.85  # Число цвайфеля для СА
C_l_RK = 0.85  # Число цвайфеля для РК

# Углы атаки на входе: i = β1л − β1 (угол лопатки − угол потока), по Лапшину К.Л. (СПбПУ, 2018).
# В [1] и [2] для турбин рекомендаций нет. Выбор по литературе:
#  - при i > 0 отрыв развивается на спинке, при i < 0 — у входной кромки со стороны корыта; в конфузорных
#    (реактивных) решётках положительные углы атаки опаснее равных по модулю отрицательных (Лапшин);
#  - профильные потери почти не меняются при i = −20…+5° и растут при i > +5° (Мамаев Б.И. и др.,
#    Теплоэнергетика, 2022, № 1);
#  - РК с большим углом поворота проектируют на отрицательные углы атаки (КиберЛенинка, малорасходные турбины).
# Поэтому расчётный угол атаки берётся небольшим отрицательным: на нерасчётных режимах (дросселирование,
# неравномерность по высоте) угол атаки смещается, и запас до +5° защищает от отрыва на спинке.
i_ca = -3.0  # Угол атаки СА, град
i_rk = -3.0  # Угол атаки РК, град (рекомендуемо −5…0°)

# Уточнение φ и ψ по модели потерь ЦИАМ ([2], п. 57…72)
UTOCHNENIE_PHI_PSI = True
r_kr_uncooled = 0.4e-3  # Радиус выходной кромки неохлаждаемой лопатки, м ([2], п. 59: 0,25…1,5 мм)
r_kr_cooled = 0.6e-3  # Радиус кромки охлаждаемой лопатки, м ([2], п. 59: толщина кромки свыше 1 мм)
T_blade_allow = 1100  # Допустимая температура неохлаждаемой лопатки, К ([2], п. 27)
T_blade_cooled = 1050  # Температура охлаждаемой лопатки, К ([2], п. 27: 950…1100 К)
# Лопатки, обтекаемые газом с температурой выше допустимой, охлаждаются ([2], п. 27), что снижает
# коэффициенты скорости ([2], п. 11 и 19, конвективное охлаждение)
cooled_CA = T_g > T_blade_allow
dphi_cool = 0.006  # φ = φ0 − (5…7)·10⁻³
dpsi_cool = 0.005  # ψ = ψ0 − (4…6)·10⁻³

# Прочность рабочей лопатки ([2], п. 28…32)
material = "ЖС6К"
rho_mat = 8100  # Плотность материала, кг/м³ ([2], табл. 2.1)
tau = 1000  # Эквивалентная наработка на взлётном режиме, ч ([2]: 1000 ч для гражданской авиации)
# Разрушающее напряжение: снимается с рис. 2.7 [2] по T_л, τ и материалу (ЖС6К, T_л = 1050 К, τ = 1000 ч → 450 МПа)
sigma_dl = 450e6
F_ratio = 0.4  # Отношение площадей сечения лопатки на периферии и втулке ([2]: 0,25…0,55)
k_m_min = 1.8  # Минимальный запас прочности ([2], п. 30)

delta_rz = 0.5e-3  # Радиальный зазор РК для оценки потерь по [2], п. 76 (0,3…1 мм)
# ==========================================
# 2. ИТЕРАЦИОННЫЙ РАСЧЕТ ТУРБИНЫ (1D)
# ==========================================
print("--- НАЧАЛО РАСЧЕТА ОСЕВОЙ ТУРБИНЫ ---\n")

n_iter = 10
for iteration in range(n_iter):
    L_tad = L_tvd / eta_i_guess
    print(f"=== Итерация {iteration + 1} (eta_i_guess = {eta_i_guess:.3f}, phi = {phi:.4f}, psi = {psi:.4f}) ===")
    print(f"Адиабатный теплоперепад L_tad* = {L_tad:.2f} Дж/кг")

    c_ad = math.sqrt(2 * L_tad)
    print(f"Условная скорость c_ad = {c_ad:.1f} м/с")
    print(f"Давление перед турбиной p_T* = {P3:.1f} Па")

    D_sr = 0.9 * D_k
    print(f"Средний диаметр D_cp = {D_sr * 1000:.0f} мм")

    u_sr = (math.pi * D_sr * n_tk) / 60
    print(f"Окружная скорость u_cp = {u_sr:.1f} м/с")

    x_1 = u_sr / c_ad
    H = 1 / (x_1 ** 2)
    print(f"Отношение скоростей x_1 = {x_1:.2f} ([1]: 0,45…0,55), Коэффициент напора H = {H:.2f}")

    L_sa = (1 - rho_st) * L_tad
    c_1 = phi * math.sqrt(2 * L_sa)
    c_1u = c_1 * math.cos(alpha_1_rad)
    c_1a = c_1 * math.sin(alpha_1_rad)  # Осевая составляющая
    T_3_4_stat = T_g - (c_1 ** 2) / (2 * cp_prime)  # статическая
    P_3_4_stat = P3 * (1 - L_sa / (cp_prime * T_g)) ** k_ratio
    rho_1 = P_3_4_stat / (R_g * T_3_4_stat)

    print(f"Теплоперепад в соплах L_sa = {L_sa:.1f} Дж/кг")
    print(f"Скорость на выходе из сопел c_1 = {c_1:.1f} м/с")
    print(f"Статическая Температура после сопел T_1 = {T_3_4_stat:.1f} К, Давление p_1 = {P_3_4_stat:.1f} Па")

    l_c = round(g_gas / (math.pi * D_sr * c_1 * math.sin(alpha_1_rad) * rho_1),3)
    S_sa = round(l_c / 1.3, 3)  # Осевая ширина СА: (h/S)_са = 1,3 ([2], п. 42: 1…3 для малорасходных турбин)

    # Число лопаток СА по параметру Цвайфеля — ближайшее простое число ([2], п. 49)
    Z_ca_bl = bt.prime_blade_number(bt.zweifel_blade_number(D_sr, alpha_0_deg, alpha_1_deg, S_sa, C_l_CA))
    C_l_CA_fact = bt.zweifel_number(D_sr, alpha_0_deg, alpha_1_deg, S_sa, Z_ca_bl)

    t_ca = math.pi * D_sr / Z_ca_bl
    y_sa = round(85.95 - 1.277 * (alpha_0_deg - alpha_1_deg) + 0.0084 * (alpha_0_deg - alpha_1_deg) ** 2)
    b_sa = S_sa / math.sin(math.radians(y_sa))
    print(f"t/b {t_ca / b_sa}")
    d_sa = D_sr - l_c
    D_sa = D_sr + l_c
    print(f"Высота лопаток сопел l_c = {l_c * 1000:.1f} мм, Число лопаток Z_ca = {Z_ca_bl} (C_L = {C_l_CA_fact:.3f})")

    L_rk = rho_st * L_tad
    w_1 = math.sqrt(c_1 ** 2 + u_sr ** 2 - 2 * c_1 * u_sr * math.cos(alpha_1_rad))

    T1_w = T_3_4_stat + w_1 ** 2 / (2 * cp_prime)
    lambda_c1 = c_1 / math.sqrt(2 * k_g * R_g * T_g / (k_g + 1))
    lambda_w1 = w_1 / math.sqrt(2 * k_g * R_g * T1_w / (k_g + 1))
    if lambda_w1 >= 1 or lambda_c1 >= 1:
        print("Запирание потока!")
    # atan2 даёт верный угол и при u > c1u (beta_1 > 90°), в отличие от asin
    beta_1_rad = math.atan2(c_1a, c_1u - u_sr)
    beta_1_deg = math.degrees(beta_1_rad)

    w_2 = psi * math.sqrt(w_1 ** 2 + 2 * L_rk)
    # Статические параметры за РК — из расчёта самой ступени ([1], с. 34):
    # T2 — из энергии в относительном движении, p2 — по адиабатному теплоперепаду ступени L_tad
    T_4_stat = T_3_4_stat - (w_2 ** 2 - w_1 ** 2) / (2 * cp_prime)
    P4_stat = P3 * (1 - L_tad / (cp_prime * T_g)) ** k_ratio
    rho_2 = P4_stat / (R_g * T_4_stat)

    l_rk = l_c + dl_rk
    F_2 = l_rk * math.pi * D_sr
    sin_beta_2 = g_gas / (w_2 * F_2 * rho_2)
    if sin_beta_2 >= 1:
        raise ValueError("Расход не проходит через РК: увеличьте dl_rk или D_sr")
    beta_2_rad = math.asin(sin_beta_2)
    beta_2_deg = math.degrees(beta_2_rad)
    if beta_2_deg < 10 or beta_1_deg <= beta_2_deg:
        print("⚠️ Угол beta_2 вне рекомендаций Батурина (beta_2 >= 10°, beta_1 > beta_2)")
    w_2u = w_2 * math.cos(beta_2_rad)
    w_2a = w_2 * math.sin(beta_2_rad)
    c_2u = u_sr - w_2u
    c_2a = w_2a
    c_2 = math.sqrt(c_2u ** 2 + c_2a ** 2)
    # alpha_2 отсчитывается от направления, противоположного вращению (Батурин, Шестаков):
    # alpha_2 < 90° — закрутка на выходе против вращения, > 90° — по вращению (c_2u здесь — по вращению)
    alpha_2_rad = math.atan2(c_2a, -c_2u)
    alpha_2_deg = math.degrees(alpha_2_rad)
    print(f"Статические параметры за РК: p2 = {P4_stat:.1f} Па, T2 = {T_4_stat:.1f} К")
    print(f"Высота лопатки РК l_рк = {l_rk * 1000:.1f} мм, отношение D_ср/l_рк = {D_sr / l_rk:.2f}")
    S_rk = round(l_rk/1.5, 3)

    d_rk = D_sr - l_rk
    D_rk = D_sr + l_rk

    Z_rk_bl = bt.prime_blade_number(bt.zweifel_blade_number(D_sr, beta_1_deg, beta_2_deg, S_rk, C_l_RK),
                                    exclude=(Z_ca_bl,))
    C_l_RK_fact = bt.zweifel_number(D_sr, beta_1_deg, beta_2_deg, S_rk, Z_rk_bl)
    t_rk = math.pi * D_sr / Z_rk_bl
    y_rk = 85.95 - 1.277 * (beta_1_deg - beta_2_deg) + 0.0084 * (beta_1_deg - beta_2_deg) ** 2
    b_rk = S_rk / math.sin(math.radians(y_rk))
    print(f"t/b rk {t_rk / b_rk}")
    print(f"Число лопаток РК Z_rk = {Z_rk_bl} (C_L = {C_l_RK_fact:.3f})")
    print(f"Относительная скорость на входе w_1 = {w_1:.1f} м/с, Угол beta_1 = {beta_1_deg:.1f}°")
    print(f"Относительная скорость на выходе w_2 = {w_2:.1f} м/с, Угол beta_2 = {beta_2_deg:.1f}°")
    print(f"Абсолютная скорость на выходе c_2 = {c_2:.1f} м/с, Угол alpha_2 = {alpha_2_deg:.1f}°")

    Z_c = ((1 / (phi ** 2)) - 1) * ((c_1 ** 2) / 2)
    Z_rk = ((1 / (psi ** 2)) - 1) * ((w_2 ** 2) / 2)
    L_l = L_tad - Z_c - Z_rk
    eta_ad = L_l / L_tad

    Z_v = (c_2 ** 2) / 2
    L_u = L_l - Z_v
    eta_u = L_u / L_tad

    Z_ut = 0.01 * L_tad
    betta_koeff = 1.2
    Z_tv = (betta_koeff * (D_sr ** 2) * ((u_sr / 10) ** 3) * (rho_1 + rho_2) / (2 * g_gas))

    L_i = L_u - Z_ut - Z_tv
    eta_i = L_i / L_tad

    print(f"Коэффициент нагрузки L_u/u² = {L_u / u_sr ** 2:.2f} ([2]: 1…2,5), "
          f"коэффициент расхода c_a/u = {c_1a / u_sr:.2f}")
    print(f"Внутренний КПД eta_i = {eta_i:.3f}")

    # Уточнение коэффициентов скорости по модели потерь ЦИАМ ([2], п. 57…72)
    T2_w = T_4_stat + w_2 ** 2 / (2 * cp_prime)
    lambda_w2 = w_2 / math.sqrt(2 * k_g * R_g * T2_w / (k_g + 1))
    cooled_RK = T1_w > T_blade_allow
    Re_ca = c_1 * b_sa * rho_1 / bt.mu_gas(T_3_4_stat)
    # В РК определяющая скорость — относительная w2 (в [2], п. 65, напечатано c2)
    Re_rk = w_2 * b_rk * rho_2 / bt.mu_gas(T_4_stat)
    phi_new, xi_ca = bt.velocity_coefficient(alpha_0_deg, alpha_1_deg, t_ca, b_sa, l_c,
                                             r_kr_cooled if cooled_CA else r_kr_uncooled,
                                             lambda_c1 / phi, 0.9, Re_ca)
    psi_new, xi_rk = bt.velocity_coefficient(beta_1_deg, beta_2_deg, t_rk, b_rk, l_rk,
                                             r_kr_cooled if cooled_RK else r_kr_uncooled,
                                             lambda_w2 / psi, 0.8, Re_rk)
    phi_new -= dphi_cool if cooled_CA else 0.0
    psi_new -= dpsi_cool if cooled_RK else 0.0
    print(f"Модель потерь ЦИАМ: phi = {phi_new:.4f} (Re = {Re_ca:.2e}, ξ = {xi_ca['Σ']:.4f}), "
          f"psi = {psi_new:.4f} (Re = {Re_rk:.2e}, ξ = {xi_rk['Σ']:.4f})"
          + (" — с учётом охлаждения" if cooled_CA or cooled_RK else ""))
    if not UTOCHNENIE_PHI_PSI:
        phi_new, psi_new = phi, psi

    error = abs(eta_i - eta_i_guess) / eta_i_guess
    err_phi = abs(phi_new - phi) / phi
    err_psi = abs(psi_new - psi) / psi
    if error <= 0.02 and err_phi <= 0.01 and err_psi <= 0.01:
        print(f"\n>>> РАСЧЕТ СХОДИТСЯ! Погрешность eta_i {error * 100:.2f}% <= 2%, "
              f"phi {err_phi * 100:.2f}%, psi {err_psi * 100:.2f}% <= 1%.\n")
        break
    else:
        print(f"Погрешность eta_i {error * 100:.2f}%, phi {err_phi * 100:.2f}%, psi {err_psi * 100:.2f}%. "
              f"Повторяем расчет...\n")
        eta_i_guess = eta_i
        phi, psi = phi_new, psi_new
else:
    print(f"⚠️ Расчёт НЕ сошёлся за {n_iter} итераций (погрешность {error * 100:.2f}%)\n")

N_t = L_i * g_gas * eta_mex      # мощность, передаваемая на вал компрессора (с учётом мех. КПД)
N_k_treb = N_k * (1 + l_agr)     # потребная мощность: компрессор + коробка приводов
print("--- ИТОГОВЫЕ РЕЗУЛЬТАТЫ ---")
print(f"Мощность турбины на валу N_t = {N_t:.1f} Вт (с учётом eta_mex = {eta_mex})")
print(f"Потребная мощность N_k*(1+l_agr) = {N_k_treb:.1f} Вт")
print(f"Отличие N_t от N_k*(1+l_agr): {abs(N_t - N_k_treb) / N_k_treb * 100:.2f}% ([1]: не более 2%)")

# Сопоставление с термодинамическим расчётом.
# eta_i ([1]) включает потерю с выходной скоростью, т.е. отнесён к статическому давлению p2.
# При eta_i = eta_tvd давление p2 совпадает с P4 из termodynamiks.py: термодинамика в этом случае
# считает кинетическую энергию на выходе из ТВД потерянной (консервативная оценка).
# Если выходная скорость используется в СТ, с eta_tvd надо сравнивать КПД по параметрам торможения eta_т*.
T_4_full = T_4_stat + c_2 ** 2 / (2 * cp_prime)
P_4_full = P4_stat * (T_4_full / T_4_stat) ** k_ratio
L_tad_tt = cp_prime * T_g * (1 - (P_4_full / P3) ** (1 / k_ratio))
eta_tt = L_i / L_tad_tt
print(f"\nСогласование с термодинамикой (eta_tvd = {eta_tvd}, P4 = {P4 / 1000:.1f} кПа):")
print(f"  eta_i  = {eta_i:.3f} (с потерей выходной скорости), p2  = {P4_stat / 1000:.1f} кПа")
print(f"  eta_т* = {eta_tt:.3f} (по параметрам торможения),    p2* = {P_4_full / 1000:.1f} кПа")
if P_4_full < P4:
    print(f"⚠️ p2* < P4: перед СТ давление ниже принятого в термодинамике — мощность СТ будет меньше расчётной")

# Потери в радиальном зазоре и на трение диска по [2] (п. 76, 77) — для сравнения с [1]
d_eta_rz = bt.d_eta_tip_clearance(delta_rz, l_rk, D_sr, beta_2_rad, t_rk, b_rk)
d_eta_d = bt.d_eta_disk(g_gas, L_tad, d_rk, u_sr * d_rk / D_sr, (rho_1 + rho_2) / 2)
print(f"Потери в радиальном зазоре: [1] Z_ут/L = {Z_ut / L_tad * 100:.2f}%, "
      f"[2] п. 76 при δ = {delta_rz * 1000:.1f} мм — {d_eta_rz * 100:.2f}%")
print(f"Потери на трение диска:     [1] Z_тв/L = {Z_tv / L_tad * 100:.2f}%, [2] п. 77 — {d_eta_d * 100:.2f}%")

# Числа Маха
M_c1 = c_1 / math.sqrt(k_g * R_g * T_3_4_stat)
M_c2 = c_2 / math.sqrt(k_g * R_g * T_4_stat)
print(f"M_c1 = {M_c1:.3f}, M_c2 = {M_c2:.3f}, lambda_w2 = {lambda_w2:.3f}")

if M_c1 > 0.95 or M_c2 > 0.95 or lambda_w2 > 0.95:
    print("⚠️ ВНИМАНИЕ: скорость близка к звуковой!")
    print("   Увеличьте rho_st или T_g, уменьшите pik/phi.")

# Температура и прочность рабочей лопатки ([2], п. 27…32)
T_blade = bt.blade_temperature(T1_w, T_blade_allow, T_blade_cooled)
k_m = bt.print_strength_check("ТВД", T1_w, T_blade, T_blade_allow, tau, material, rho_mat, sigma_dl,
                              bt.form_factor(F_ratio), u_sr, D_sr / l_rk, k_m_min)
if cooled_CA:
    print(f"⚠️ T*_г = {T_g} К > {T_blade_allow} К: лопатки СА также необходимо охлаждать")

# ==========================================
# 3. РАСЧЕТ ПО ВЫСОТЕ ПРОТОЧНОЙ ЧАСТИ (МЕТОДИКА БАТУРИНА О.В.)
# ==========================================
print("\n" + "=" * 80)
print(" РАСЧЕТ ПАРАМЕТРОВ ПОТОКА ПО ВЫСОТЕ ЛОПАТКИ (Закон закрутки α1 = const)")
print("=" * 80)

# Батурин п. 79: 5…10 сечений по высоте лопатки
h_rel = [0.0, 0.25, 0.5, 0.75, 1.0]
sec_names = ["Втулка", "0,25h", "Средний", "0,75h", "Периферия"]
sections_ca = [{"name": f"{n} СА", "D": d_sa + (D_sa - d_sa) * h} for n, h in zip(sec_names, h_rel)]
sections_rk = [{"name": f"{n} РК", "D": d_rk + (D_rk - d_rk) * h} for n, h in zip(sec_names, h_rel)]

print("\n--- 3.1. Профилирование СА и РК по сечениям ---")
betta1_angle_arr = []
lambda_c1_arr = []
c1ui_arr = []
ui_arr = []
w_1i_arr = []
c_1i_arr = []
t_ca_i_arr = []
print(f"\n{'Сечение':<15} | {'D, мм':<8} | {'u, м/с':<8} | {'c1, м/с':<9} | {'c1u, м/с':<9} | {'w1, м/с':<9} | {'β1, град':<9} | {'λc1':<6} | {'t, мм':<8}")
print("-" * 105)


for sec in sections_ca:
    D_i = sec["D"] / D_sr
    u_i = u_sr * D_i

    alpha_1i = alpha_1_rad
    c_1i = c_1 / (D_i ** ((phi ** 2) * (math.cos(alpha_1i) ** 2)))
    c_1ai = c_1i * math.sin(alpha_1i)
    c_1ui = c_1i * math.cos(alpha_1i)
    if u_i < c_1ui:
        w_1ui = c_1ui - u_i
        beta_1i_rad = math.atan2(c_1ai, w_1ui)
    else:
        w_1ui = u_i - c_1ui
        beta_1i_rad = math.pi - math.atan2(c_1ai, w_1ui)
    lambda_c1i = c_1i / math.sqrt(2 * k_g * R_g * T_g / (k_g + 1))
    lambda_c1_arr.append(lambda_c1i)
    betta1_angle_arr.append(math.degrees(beta_1i_rad))
    c1ui_arr.append(c_1ui)
    ui_arr.append(u_i)
    w_1i = math.sqrt(w_1ui ** 2 + c_1ai ** 2)
    w_1i_arr.append(w_1i)
    c_1i_arr.append(c_1i)
    t_ca_i = math.pi * sec["D"] / Z_ca_bl
    t_ca_i_arr.append(t_ca_i)
    print(
        f"{sec['name']:<15} | {sec['D'] * 1000:<8.1f} | {u_i:<8.1f} | {c_1i:<9.1f} | {c_1ui:<9.1f} | {w_1i:<9.1f} | {math.degrees(beta_1i_rad):<9.2f} | {lambda_c1i:<6.3f} | {t_ca_i * 1000:<8.2f}")


c2_ui_arr = []
c2_ai_arr = []
c2_i_arr = []
w2_i_arr = []
betta2_i_arr = []
t_rk_i_arr = []
rho_i_arr = []
alpha2_i_arr = []
lambda_w2_arr = []

print(f"\n{'Сечение':<15} | {'D, мм':<8} | {'u, м/с':<8} | {'c2, м/с':<9} | {'c2u, м/с':<9} | {'w2, м/с':<9} | {'β2, град':<9} | {'ρ':<8} | {'λw2':<6} | {'t, мм':<8} | {'α2, град'}")
print("-" * 120)

# Постоянство работы по высоте (Батурин п. 106): u1*c1u - u2*c2u = L_u_ср
L_u_sr = u_sr * (c_1u - c_2u)  # Удельная работа на среднем радиусе

for i, sec in enumerate(sections_rk):
    D2_i = sec["D"] / D_sr
    u2_i = u_sr * D2_i

    c2_ui = (c1ui_arr[i] * ui_arr[i] - L_u_sr) / u2_i  # Окружная составляющая c2 на выходе из РК

    w2_ui = u2_i - c2_ui  # Окружная составляющая относительной скорости

    # ДОПУЩЕНИЕ: c2a = const по высоте. У Батурина (п. 108) принято alpha2 = const (c2a = c2u*tg alpha2),
    # но при alpha2 ≈ 90° эта формула плохо обусловлена (tg alpha2 → ∞), поэтому принято c2a = const.
    c2_ai = c_2a
    c2_i = math.sqrt(c2_ui ** 2 + c2_ai ** 2)  # Теорема Пифагора
    w2_i = math.sqrt(w2_ui ** 2 + c2_ai ** 2)
    betta2_i = math.atan2(c2_ai, w2_ui)  # Корректное определение угла
    alpha2i = math.degrees(math.atan2(c2_ai, -c2_ui))  # от направления, противоположного вращению
    alpha2_i_arr.append(alpha2i)

    t_rk_i = math.pi * sec["D"] / Z_rk_bl
    # Батурин п. 112: w2s = w2/psi, учитывается разница окружных скоростей u1 и u2
    rho_i = ((w2_i / psi) ** 2 - w_1i_arr[i] ** 2 + ui_arr[i] ** 2 - u2_i ** 2) / (2 * L_tad)
    # λ_w2 по относительной температуре торможения своего сечения (п. 105: T*_w2 = T*_w1 − (u1² − u2²)/2cp)
    T2_wi = T2_w + (u2_i ** 2 - u_sr ** 2) / (2 * cp_prime)
    lambda_w2i = w2_i / math.sqrt(2 * k_g * R_g * T2_wi / (k_g + 1))

    c2_ui_arr.append(c2_ui)
    c2_ai_arr.append(c2_ai)
    c2_i_arr.append(c2_i)
    w2_i_arr.append(w2_i)
    betta2_i_arr.append(math.degrees(betta2_i))
    t_rk_i_arr.append(t_rk_i)
    rho_i_arr.append(rho_i)
    lambda_w2_arr.append(lambda_w2i)

    print(
        f"{sec['name']:<15} | {sec['D'] * 1000:<8.1f} | {u2_i:<8.1f} | {c2_i:<9.1f} | {-c2_ui:<9.1f} | {w2_i:<9.1f} | {math.degrees(betta2_i):<9.2f} | {rho_i:<8.3f} | {lambda_w2i:<6.3f} | {t_rk_i * 1000:<8.2f} | {alpha2i:<8.2f}")
print("c2u, β2, α2 на выходе — относительно направления, противоположного вращению (как у Батурина)")
if rho_i_arr[0] < 0.1:
    print(f"⚠️ Реактивность на втулке {rho_i_arr[0]:.3f} < 0,1 — увеличьте rho_st (Батурин п. 44)")
if max(lambda_w2_arr) > 0.95:
    print(f"⚠️ lambda_w2 на периферии {max(lambda_w2_arr):.3f} > 0,95 — близко к запиранию РК "
          f"(Батурин п. 10 допускает до 1,1)")

# 3.2. ЛОПАТОЧНЫЕ УГЛЫ И УГЛЫ УСТАНОВКИ ПО СЕЧЕНИЯМ
beta1_bl_arr, beta2_bl_arr = bt.print_blade_angles(
    sections_ca, sections_rk, alpha_0_deg, alpha_1_deg, lambda_c1_arr, betta1_angle_arr, betta2_i_arr,
    lambda_w2_arr, S_sa, S_rk, Z_ca_bl, Z_rk_bl, i_ca, i_rk)

# 3.3. ПРОВЕРКА УГЛОВ РАСКРЫТИЯ МЕРИДИОНАЛЬНОГО КАНАЛА (Батурин п. 43, стр. 63)
# ==========================================
print("\n--- 3.3. Проверка углов раскрытия проточной части (Батурин п. 43) ---")

# Батурин п. 43: периферия (D2к − D1к)/(2Sрк), втулка (D1вт − D2вт)/(2Sрк).
# Положительный угол — канал расширяется, отрицательный — сужается.
# СА выполнен с цилиндрическими обводами (D0 = D1), поэтому углы раскрытия СА равны 0°.
gamma_rk_k = math.degrees(math.atan((D_rk - D_sa) / (2 * S_rk)))
gamma_rk_vt = math.degrees(math.atan((d_sa - d_rk) / (2 * S_rk)))

print(f"Угол раскрытия РК на периферии: {gamma_rk_k:.2f}° (Должно быть меньше: 15...20°)")
print(f"Угол раскрытия РК на втулке:    {gamma_rk_vt:.2f}° (Должно быть меньше: 8...12°)")

if gamma_rk_vt < 0 or gamma_rk_k < 0:
    print("⚠️ ВНИМАНИЕ: канал РК сужается (отрицательный угол раскрытия) — лопатка РК ниже лопатки СА.")
    print("   Рекомендация: увеличить dl_rk (l_рк = l_c + 0…2 мм).")
elif gamma_rk_vt > 12 or gamma_rk_k > 20:
    print("⚠️ ВНИМАНИЕ: Угол раскрытия превышает норму! Возможны отрывы потока.")
    print("   Рекомендация: Увеличить осевую ширину венцов (S = l / 1.5 или меньше) в 1D расчете.")
else:
    print("✅ Геометрия меридионального канала обеспечивает безотрывное течение.")

# Осевой зазор δ = (0,2…0,3)·S_рк (Батурин, п. 42)
p.plot_meridional_section(S_sa, S_rk, d_sa, D_sa, d_rk, D_rk, axial_gap=0.25 * S_rk)
p.plot_velocity_triangles(
    c_1=c_1,
    alpha_1_rad=alpha_1_rad,
    w_1=w_1,
    beta_1_rad=beta_1_rad,
    c_2=c_2,
    alpha_2_rad=alpha_2_rad,
    w_2=w_2,
    beta_2_rad=beta_2_rad,
    u_sr=u_sr,
    section_name=f"Среднее сечение (D = {D_sr * 1000:.0f} мм)"
)
