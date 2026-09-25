import math
from turbine import plot_turb as p
from termo.termodynamiks import (R_g, cp_g, alpha, eta_ct, eta_mex_ct, g_gas, T4, P45, P5, P0, T5,
                                 L_ct, N_ct)

# ==========================================
# 1. ИСХОДНЫЕ ДАННЫЕ (из термодинамического расчёта)
# ==========================================
# Свободная турбина (СТ) на отдельном валу: частота вращения не связана с ротором турбокомпрессора,
# выбирается из условия оптимального u/c_ad (далее — через редуктор на нагрузку)
n_ct = 28000  # Частота вращения ротора СТ, об/мин
D_sr = 0.300  # Средний диаметр СТ, м (больше, чем у ТВД 207 мм: объёмный расход выше; переходный канал с подъёмом радиуса)

# Параметры газа на входе в СТ (полные) — сечение 45 термодинамического расчёта
T_in = T4  # Полная температура перед СТ, К
P_in = P45  # Полное давление перед СТ (после промежуточного тракта), Па
cp_prime = cp_g(T_in, alpha)  # Изобарная теплоемкость газа, Дж/(кг·К)
k_g = cp_prime / (cp_prime - R_g)  # k газа, согласованный с cp (как k_ct в termodynamiks.py)
k_ratio = k_g / (k_g - 1)

# Потребная работа газа на лопатках: L_ct в termodynamiks.py — работа на валу (с eta_mex_ct)
L_ct_gas = L_ct / eta_mex_ct

# Предварительно принимаемые коэффициенты
eta_i_guess = eta_ct  # Начальное приближение внутреннего КПД турбины
phi = 0.98  # Коэффициент скорости соплового аппарата
psi = 0.96  # Коэффициент скорости рабочего колеса (Шестаков 0,93…0,97; Батурин 0,96…0,97)
alpha_1_deg = 20.0  # Угол выхода газа из сопел, град
alpha_1_rad = math.radians(alpha_1_deg)
rho_st = 0.5  # Степень реактивности
alpha_0_deg = 90  # Угол входа в СА: за ТВД поток практически осевой (alpha_2 ТВД ≈ 90°)
# Угол beta_2 не задаётся, а находится из уравнения неразрывности (Шестаков, п. 1.4)
# при высоте лопатки РК l_рк = l_c + (0…2 мм), чтобы канал не сужался по тракту
dl_rk = 0.001  # Превышение высоты лопатки РК над высотой лопатки СА, м

C_l_CA = 0.85  # Число цвайфеля для СА
C_l_RK = 0.85  # Число цвайфеля для РК
# ==========================================
# 2. ИТЕРАЦИОННЫЙ РАСЧЕТ СВОБОДНОЙ ТУРБИНЫ (1D)
# ==========================================
print("--- НАЧАЛО РАСЧЕТА СВОБОДНОЙ ТУРБИНЫ ---\n")
print(f"Потребная работа газа L_ст = {L_ct_gas:.1f} Дж/кг (на валу {L_ct:.1f} Дж/кг, eta_мех_ст = {eta_mex_ct})")

n_iter = 20
for iteration in range(n_iter):
    L_tad = L_ct_gas / eta_i_guess
    print(f"=== Итерация {iteration + 1} (eta_i_guess = {eta_i_guess:.3f}) ===")
    print(f"Адиабатный теплоперепад L_tad* = {L_tad:.2f} Дж/кг")
    if L_tad >= cp_prime * T_in:
        raise ValueError("Теплоперепад больше энтальпии торможения на входе — проверьте исходные данные")

    c_ad = math.sqrt(2 * L_tad)
    print(f"Условная скорость c_ad = {c_ad:.1f} м/с")
    print(f"Давление перед СТ p_СТ* = {P_in:.1f} Па")

    u_sr = (math.pi * D_sr * n_ct) / 60
    print(f"Средний диаметр D_cp = {D_sr * 1000:.0f} мм, окружная скорость u_cp = {u_sr:.1f} м/с")
    print(f'Коэффициент нагрузки {L_tad / u_sr ** 2}')

    x_1 = u_sr / c_ad
    H = 1 / (x_1 ** 2)
    print(f"Отношение скоростей x_1 = {x_1:.2f}, Коэффициент напора H = {H:.2f}")

    L_sa = (1 - rho_st) * L_tad
    c_1 = phi * math.sqrt(2 * L_sa)
    c_1u = c_1 * math.cos(alpha_1_rad)
    c_1a = c_1 * math.sin(alpha_1_rad)  # Осевая составляющая
    T_1_stat = T_in - (c_1 ** 2) / (2 * cp_prime)  # статическая
    P_1_stat = P_in * (1 - L_sa / (cp_prime * T_in)) ** k_ratio
    rho_1 = P_1_stat / (R_g * T_1_stat)

    print(f"Теплоперепад в соплах L_sa = {L_sa:.1f} Дж/кг")
    print(f"Скорость на выходе из сопел c_1 = {c_1:.1f} м/с")
    print(f"Статическая Температура после сопел T_1 = {T_1_stat:.1f} К, Давление p_1 = {P_1_stat:.1f} Па")

    l_c = round(g_gas / (math.pi * D_sr * c_1a * rho_1), 3)
    S_sa = round(l_c / 1.3, 3)  # Осевая ширина СА: (h/S)_са = 1,3

    # Параметр Цвайфеля для СА
    Z_ca_bl = round(2 * math.pi * D_sr * (math.sin(alpha_1_rad) ** 2) * (
                1 / math.tan(alpha_1_rad) + 1 / math.tan(math.radians(alpha_0_deg))) / (S_sa * C_l_CA))

    t_ca = math.pi * D_sr / Z_ca_bl
    y_sa = round(85.95 - 1.277 * (alpha_0_deg - alpha_1_deg) + 0.0084 * (alpha_0_deg - alpha_1_deg) ** 2)
    b_sa = S_sa / math.sin(math.radians(y_sa))
    print(f"t/b {t_ca / b_sa}")
    d_sa = D_sr - l_c
    D_sa = D_sr + l_c
    print(f"Высота лопаток сопел l_c = {l_c * 1000:.1f} мм, Число лопаток Z_ca ≈ {Z_ca_bl}")

    L_rk = rho_st * L_tad
    w_1 = math.sqrt(c_1 ** 2 + u_sr ** 2 - 2 * c_1 * u_sr * math.cos(alpha_1_rad))

    T1_w = T_1_stat + w_1 ** 2 / (2 * cp_prime)
    lambda_c1 = c_1 / math.sqrt(2 * k_g * R_g * T_in / (k_g + 1))
    lambda_w1 = w_1 / math.sqrt(2 * k_g * R_g * T1_w / (k_g + 1))
    if lambda_w1 >= 1 or lambda_c1 >= 1:
        print("Запирание потока!")
    beta_1_rad = math.atan2(c_1a, c_1u - u_sr)
    beta_1_deg = math.degrees(beta_1_rad)

    w_2 = psi * math.sqrt(w_1 ** 2 + 2 * L_rk)
    T_2_stat = T_1_stat - (w_2 ** 2 - w_1 ** 2) / (2 * cp_prime)
    P2_stat = P_in * (1 - L_tad / (cp_prime * T_in)) ** k_ratio
    rho_2 = P2_stat / (R_g * T_2_stat)

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
    # alpha_2 отсчитывается от направления, противоположного вращению (Батурин, Шестаков)
    alpha_2_rad = math.atan2(c_2a, -c_2u)
    alpha_2_deg = math.degrees(alpha_2_rad)
    print(f"Статические параметры за РК: p2 = {P2_stat:.1f} Па, T2 = {T_2_stat:.1f} К")
    print(f"Высота лопатки РК l_рк = {l_rk * 1000:.1f} мм, отношение D_ср/l_рк = {D_sr / l_rk:.2f}")
    S_rk = round(l_rk / 1.5, 3)

    d_rk = D_sr - l_rk
    D_rk = D_sr + l_rk

    Z_rk_bl = round(2 * math.pi * D_sr * (math.sin(beta_2_rad) ** 2) * (
            1 / math.tan(beta_2_rad) + 1 / math.tan(beta_1_rad)) / (S_rk * C_l_RK))
    t_rk = math.pi * D_sr / Z_rk_bl
    y_rk = 85.95 - 1.277 * (beta_1_deg - beta_2_deg) + 0.0084 * (beta_1_deg - beta_2_deg) ** 2
    b_rk = S_rk / math.sin(math.radians(y_rk))
    print(f"t/b rk {t_rk / b_rk}")
    print(f"Число лопаток РК Z_rk ≈ {Z_rk_bl}")
    print(f"Относительная скорость на входе w_1 = {w_1:.1f} м/с, Угол beta_1 = {beta_1_deg:.1f}°")
    print(f"Относительная скорость на выходе w_2 = {w_2:.1f} м/с, Угол beta_2 = {beta_2_deg:.1f}°")
    print(f"Абсолютная скорость на выходе c_2 = {c_2:.1f} м/с, Угол alpha_2 = {alpha_2_deg:.1f}°")

    Z_c = ((1 / (phi ** 2)) - 1) * ((c_1 ** 2) / 2)
    Z_rk = ((1 / (psi ** 2)) - 1) * ((w_2 ** 2) / 2)
    L_l = L_tad - Z_c - Z_rk
    eta_ad = L_l / L_tad

    # Выходная скорость СТ последней ступени в выхлоп не используется — это потеря
    Z_v = (c_2 ** 2) / 2
    L_u = L_l - Z_v
    eta_u = L_u / L_tad

    Z_ut = 0.01 * L_tad
    betta_koeff = 1.2
    Z_tv = (betta_koeff * (D_sr ** 2) * ((u_sr / 10) ** 3) * (rho_1 + rho_2) / (2 * g_gas))

    L_i = L_u - Z_ut - Z_tv
    eta_i = L_i / L_tad

    print(f"Внутренний КПД eta_i = {eta_i:.3f}")

    error = abs(eta_i - eta_i_guess) / eta_i_guess
    if error <= 0.005:
        print(f"\n>>> РАСЧЕТ СХОДИТСЯ! Погрешность {error * 100:.2f}% <= 0,5%.\n")
        break
    else:
        print(f"Погрешность {error * 100:.2f}% > 0,5%. Повторяем расчет...\n")
        eta_i_guess = eta_i
else:
    print(f"⚠️ Расчёт НЕ сошёлся за {n_iter} итераций (погрешность {error * 100:.2f}%)\n")

# Высота лопатки РК, при которой выход осевой (c_2u = 0, alpha_2 = 90°): w_2u = u, плотность rho_2 от beta_2 не зависит
if w_2 > u_sr:
    l_rk_axial = g_gas / (math.pi * D_sr * math.sqrt(w_2 ** 2 - u_sr ** 2) * rho_2)
    print(f"Для осевого выхода (alpha_2 = 90°) нужна l_рк = {l_rk_axial * 1000:.1f} мм "
          f"(dl_rk = {(l_rk_axial - l_c) * 1000:.1f} мм)")

N_st = L_i * g_gas * eta_mex_ct  # мощность на выходном валу СТ
print("--- ИТОГОВЫЕ РЕЗУЛЬТАТЫ ---")
print(f"Мощность СТ на валу N_ст = {N_st / 1000:.1f} кВт = {N_st / 1000 * 1.36:.1f} л.с. (eta_мех_ст = {eta_mex_ct})")
print(f"Мощность по термодинамике N_ст = {N_ct / 1000:.1f} кВт = {N_ct / 1000 * 1.36:.1f} л.с.")
print(f"Отличие: {abs(N_st - N_ct) / N_ct * 100:.2f}%")
if N_st / 1000 * 1.36 < 800:
    print("⚠️ Мощность меньше заданных 800 л.с.!")

# Сопоставление с термодинамикой.
# eta_i включает потерю с выходной скоростью (отнесён к статическому давлению за СТ).
# В termodynamiks.py eta_ct отнесён к полному давлению P5 за СТ, поэтому давление торможения
# за РК p2* должно быть не ниже P5 — иначе располагаемого перепада не хватит на выхлопной тракт.
T_2_full = T_2_stat + c_2 ** 2 / (2 * cp_prime)
P_2_full = P2_stat * (T_2_full / T_2_stat) ** k_ratio
L_tad_tt = cp_prime * T_in * (1 - (P_2_full / P_in) ** (1 / k_ratio))
eta_tt = L_i / L_tad_tt
print(f"\nСогласование с термодинамикой (eta_ct = {eta_ct}, P5 = {P5 / 1000:.1f} кПа, T5 = {T5:.1f} К):")
print(f"  eta_i  = {eta_i:.3f} (с потерей выходной скорости), p2  = {P2_stat / 1000:.1f} кПа")
print(f"  eta_т* = {eta_tt:.3f} (по параметрам торможения),    p2* = {P_2_full / 1000:.1f} кПа, T2* = {T_2_full:.1f} К")
print(f"  pi_ст* = {P_in / P_2_full:.4f} (по термодинамике pi_ст = {P_in / P5:.4f})")
if P_2_full < P5:
    print("⚠️ p2* < P5: располагаемого перепада не хватает — уменьшите потери (c_2, psi) или поднимите pik/T_g")
if P2_stat < P0:
    print(f"   Статическое давление за СТ ниже атмосферного ({P0 / 1000:.1f} кПа) — нужен выходной диффузор")

# Числа Маха
M_c1 = c_1 / math.sqrt(k_g * R_g * T_1_stat)
M_c2 = c_2 / math.sqrt(k_g * R_g * T_2_stat)
T2_w = T_2_stat + w_2 ** 2 / (2 * cp_prime)
a_cr_w2 = math.sqrt(2 * k_g * R_g * T2_w / (k_g + 1))
lambda_w2 = w_2 / a_cr_w2
print(f"M_c1 = {M_c1:.3f}, M_c2 = {M_c2:.3f}, lambda_w2 = {lambda_w2:.3f}")
if M_c1 > 0.95 or M_c2 > 0.95 or lambda_w2 > 0.95:
    print("⚠️ ВНИМАНИЕ: скорость близка к звуковой!")

# Прочность: напряжение растяжения в корне лопатки РК (без бандажа, коэф. формы 0,55 для сужающейся лопатки)
rho_mat = 8200  # Плотность материала лопатки (жаропрочный никелевый сплав), кг/м³
F_ann = math.pi * D_sr * l_rk  # Ометаемая площадь РК, м²
sigma_p = 0.55 * rho_mat * (math.pi * n_ct / 30) ** 2 * F_ann / (2 * math.pi)
print(f"Напряжение растяжения в корне лопатки РК sigma_p ≈ {sigma_p / 1e6:.0f} МПа")

# ==========================================
# 3. РАСЧЕТ ПО ВЫСОТЕ ПРОТОЧНОЙ ЧАСТИ (МЕТОДИКА БАТУРИНА О.В.)
# ==========================================
print("\n" + "=" * 80)
print(" РАСЧЕТ ПАРАМЕТРОВ ПОТОКА ПО ВЫСОТЕ ЛОПАТКИ (Закон закрутки α1 = const)")
print("=" * 80)

h_rel = [0.0, 0.25, 0.5, 0.75, 1.0]
sec_names = ["Втулка", "0,25h", "Средний", "0,75h", "Периферия"]
sections_ca = [{"name": f"{n} СА", "D": d_sa + (D_sa - d_sa) * h} for n, h in zip(sec_names, h_rel)]
sections_rk = [{"name": f"{n} РК", "D": d_rk + (D_rk - d_rk) * h} for n, h in zip(sec_names, h_rel)]

print("\n--- 3.1. Профилирование СА и РК по сечениям ---")
c1ui_arr = []
ui_arr = []
w_1i_arr = []
print(f"\n{'Сечение':<15} | {'D, мм':<8} | {'u, м/с':<8} | {'c1, м/с':<9} | {'c1u, м/с':<9} | {'w1, м/с':<9} | {'β1, град':<9} | {'λc1':<6} | {'t, мм':<8}")
print("-" * 105)

for sec in sections_ca:
    D_i = sec["D"] / D_sr
    u_i = u_sr * D_i

    c_1i = c_1 / (D_i ** ((phi ** 2) * (math.cos(alpha_1_rad) ** 2)))
    c_1ai = c_1i * math.sin(alpha_1_rad)
    c_1ui = c_1i * math.cos(alpha_1_rad)
    beta_1i_rad = math.atan2(c_1ai, c_1ui - u_i)
    w_1i = math.sqrt((c_1ui - u_i) ** 2 + c_1ai ** 2)
    lambda_c1i = c_1i / math.sqrt(2 * k_g * R_g * T_in / (k_g + 1))
    c1ui_arr.append(c_1ui)
    ui_arr.append(u_i)
    w_1i_arr.append(w_1i)
    t_ca_i = math.pi * sec["D"] / Z_ca_bl
    print(
        f"{sec['name']:<15} | {sec['D'] * 1000:<8.1f} | {u_i:<8.1f} | {c_1i:<9.1f} | {c_1ui:<9.1f} | {w_1i:<9.1f} | {math.degrees(beta_1i_rad):<9.2f} | {lambda_c1i:<6.3f} | {t_ca_i * 1000:<8.2f}")

rho_i_arr = []
lambda_w2_arr = []
L_u_sr = u_sr * (c_1u - c_2u)  # Удельная работа на среднем радиусе (Батурин п. 106)

print(f"\n{'Сечение':<15} | {'D, мм':<8} | {'u, м/с':<8} | {'c2, м/с':<9} | {'c2u, м/с':<9} | {'w2, м/с':<9} | {'β2, град':<9} | {'ρ':<8} | {'λw2':<6} | {'t, мм':<8} | {'α2, град'}")
print("-" * 120)

for i, sec in enumerate(sections_rk):
    D2_i = sec["D"] / D_sr
    u2_i = u_sr * D2_i

    c2_ui = (c1ui_arr[i] * ui_arr[i] - L_u_sr) / u2_i
    w2_ui = u2_i - c2_ui
    # ДОПУЩЕНИЕ: c2a = const по высоте (как в расчёте ТВД)
    c2_ai = c_2a
    c2_i = math.sqrt(c2_ui ** 2 + c2_ai ** 2)
    w2_i = math.sqrt(w2_ui ** 2 + c2_ai ** 2)
    betta2_i = math.atan2(c2_ai, w2_ui)
    alpha2i = math.degrees(math.atan2(c2_ai, -c2_ui))

    t_rk_i = math.pi * sec["D"] / Z_rk_bl
    rho_i = ((w2_i / psi) ** 2 - w_1i_arr[i] ** 2 + ui_arr[i] ** 2 - u2_i ** 2) / (2 * L_tad)
    # λ_w2 по относительной температуре торможения, в каждом сечении своя (ротальпия: T_w ~ T* + (u²)/2cp)
    T2_wi = T2_w + (u2_i ** 2 - u_sr ** 2) / (2 * cp_prime)
    lambda_w2i = w2_i / math.sqrt(2 * k_g * R_g * T2_wi / (k_g + 1))
    rho_i_arr.append(rho_i)
    lambda_w2_arr.append(lambda_w2i)

    print(
        f"{sec['name']:<15} | {sec['D'] * 1000:<8.1f} | {u2_i:<8.1f} | {c2_i:<9.1f} | {-c2_ui:<9.1f} | {w2_i:<9.1f} | {math.degrees(betta2_i):<9.2f} | {rho_i:<8.3f} | {lambda_w2i:<6.3f} | {t_rk_i * 1000:<8.2f} | {alpha2i:<8.2f}")
print("c2u, β2, α2 на выходе — относительно направления, противоположного вращению (как у Батурина)")
if rho_i_arr[0] < 0.1:
    print(f"⚠️ Реактивность на втулке {rho_i_arr[0]:.3f} < 0,1 — увеличьте rho_st (Батурин п. 44)")
if max(lambda_w2_arr) > 0.95:
    print(f"⚠️ lambda_w2 на периферии {max(lambda_w2_arr):.3f} > 0,95 — близко к запиранию РК")

print("\n--- 3.3. Проверка углов раскрытия проточной части (Батурин п. 43) ---")
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
    section_name=f"СТ, среднее сечение (D = {D_sr * 1000:.0f} мм)"
)
