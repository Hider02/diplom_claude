import math
import plot_turb as p
from termo.termodynamiks import (k_g, R_g, cp_g, T_g, alpha, eta_tvd, g_gas, P3, N_k, P4, L_tvd, T4,
                                     eta_mex, l_agr)

# ==========================================
# 1. ИСХОДНЫЕ ДАННЫЕ (из расчета компрессора)
# ==========================================
n_tk = 45000  # Частота вращения ротора, об/мин
D_k = 0.23  # Наружный диаметр колеса компрессора, м

# Параметры газа (продукты сгорания)
k_ratio = k_g / (k_g - 1)  # k/(k-1) для газа (k_g согласован с cp_g, см. termodynamiks.py)
cp_prime = cp_g(T_g, alpha)  # Изобарная теплоемкость газа, Дж/(кг·К)

# Предварительно принимаемые коэффициенты
eta_i_guess = eta_tvd  # Начальное приближение внутреннего КПД турбины
phi = 0.98  # Коэффициент скорости соплового аппарата
psi = 0.96  # Коэффициент скорости рабочего колеса (Шестаков 0,93…0,97; Батурин 0,96…0,97)
alpha_1_deg = 20.0  # Угол выхода газа из сопел, град
alpha_1_rad = math.radians(alpha_1_deg)
rho_st = 0.5  # Степень реактивности
alpha_0_deg = 90  # Угол входа в СА
# Угол beta_2 не задаётся, а находится из уравнения неразрывности (Шестаков, п. 1.4)
# при высоте лопатки РК l_рк = l_c + (0…2 мм), чтобы канал не сужался по тракту
dl_rk = 0.001  # Превышение высоты лопатки РК над высотой лопатки СА, м

C_l_CA = 0.85  # Число цвайфеля для СА
C_l_RK = 0.85  # Число цвайфеля для РК
# ==========================================
# 2. ИТЕРАЦИОННЫЙ РАСЧЕТ ТУРБИНЫ (1D)
# ==========================================
print("--- НАЧАЛО РАСЧЕТА ОСЕВОЙ ТУРБИНЫ ---\n")

for iteration in range(5):
    L_tad = L_tvd / eta_i_guess
    print(f"=== Итерация {iteration + 1} (eta_i_guess = {eta_i_guess:.3f}) ===")
    print(f"Адиабатный теплоперепад L_tad* = {L_tad:.2f} Дж/кг")

    c_ad = math.sqrt(2 * L_tad)
    print(f"Условная скорость c_ad = {c_ad:.1f} м/с")
    print(f"Давление перед турбиной p_T* = {P3:.1f} Па")

    D_sr = 0.9 * D_k
    print(f"Средний диаметр D_cp = {D_sr * 1000:.0f} мм")

    u_sr = (math.pi * D_sr * n_tk) / 60
    print(f"Окружная скорость u_cp = {u_sr:.1f} м/с")
    print(f'Коэффициент нагрузки {L_tad / u_sr ** 2}')

    x_1 = u_sr / c_ad
    H = 1 / (x_1 ** 2)
    print(f"Отношение скоростей x_1 = {x_1:.2f}, Коэффициент напора H = {H:.2f}")

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
    S_sa = round(l_c / 1.3, 3)  # Осевая ширина СА: (h/S)_са = 1,3

    # Параметр Цвайфеля для СА (cot(90°) = 0)
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

    T1_w = T_3_4_stat + w_1 ** 2 / (2 * cp_prime)
    lambda_c1 = c_1 / math.sqrt(2 * k_g * R_g * T_g / (k_g + 1))
    lambda_w1 = w_1 / math.sqrt(2 * k_g * R_g * T1_w / (k_g + 1))
    if lambda_w1 >= 1 or lambda_c1 >= 1:
        print("Запирание потока!")
    # atan2 даёт верный угол и при u > c1u (beta_1 > 90°), в отличие от asin
    beta_1_rad = math.atan2(c_1a, c_1u - u_sr)
    beta_1_deg = math.degrees(beta_1_rad)

    w_2 = psi * math.sqrt(w_1 ** 2 + 2 * L_rk)
    # Статические параметры за РК — из расчёта самой ступени (Шестаков, п. 1.4):
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
    if error <= 0.02:
        print(f"\n>>> РАСЧЕТ СХОДИТСЯ! Погрешность {error * 100:.2f}% <= 2%.\n")
        break
    else:
        print(f"Погрешность {error * 100:.2f}% > 2%. Повторяем расчет...\n")
        eta_i_guess = eta_i
else:
    print(f"⚠️ Расчёт НЕ сошёлся за 5 итераций (погрешность {error * 100:.2f}%)\n")

N_t = L_i * g_gas * eta_mex      # мощность, передаваемая на вал компрессора (с учётом мех. КПД)
N_k_treb = N_k * (1 + l_agr)     # потребная мощность: компрессор + коробка приводов
print("--- ИТОГОВЫЕ РЕЗУЛЬТАТЫ ---")
print(f"Мощность турбины на валу N_t = {N_t:.1f} Вт (с учётом eta_mex = {eta_mex})")
print(f"Потребная мощность N_k*(1+l_agr) = {N_k_treb:.1f} Вт")
print(f"Отличие N_t от N_k*(1+l_agr): {abs(N_t - N_k_treb) / N_k_treb * 100:.2f}%")

# Сопоставление с термодинамическим расчётом.
# eta_i (по Шестакову) включает потерю с выходной скоростью, т.е. отнесён к статическому давлению p2.
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
# Числа Маха
M_c1 = c_1 / math.sqrt(k_g * R_g * T_3_4_stat)
M_c2 = c_2 / math.sqrt(k_g * R_g * T_4_stat)

# Относительная скорость w2
T2_w = T_4_stat + w_2**2 / (2 * cp_prime)
a_cr_w2 = math.sqrt(2 * k_g * R_g * T2_w / (k_g + 1))
lambda_w2 = w_2 / a_cr_w2

print(f"M_c1 = {M_c1:.3f}, M_c2 = {M_c2:.3f}, lambda_w2 = {lambda_w2:.3f}")

if M_c1 > 0.95 or M_c2 > 0.95 or lambda_w2 > 0.95:
    print("⚠️ ВНИМАНИЕ: скорость близка к звуковой!")
    print("   Увеличьте rho_st или T_g, уменьшите pik/phi.")

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
c1ui_arr = []
ui_arr = []
w_1i_arr = []
c_1i_arr = []
t_ca_i_arr = []
print(f"\n{'Сечение':<15} | {'D, мм':<8} | {'u, м/с':<8} | {'c1, м/с':<9} | {'c1u, м/с':<9} | {'w1, м/с':<9} | {'β1, град':<9} | {'t, мм':<8}")
print("-" * 95)


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
        a = math.atan2(c_1ai, w_1ui)
    betta1_angle_arr.append(math.degrees(beta_1i_rad))
    c1ui_arr.append(c_1ui)
    ui_arr.append(u_i)
    w_1i = math.sqrt(w_1ui ** 2 + c_1ai ** 2)
    w_1i_arr.append(w_1i)
    c_1i_arr.append(c_1i)
    t_ca_i = math.pi * sec["D"] / Z_ca_bl
    t_ca_i_arr.append(t_ca_i)
    print(
        f"{sec['name']:<15} | {sec['D'] * 1000:<8.1f} | {u_i:<8.1f} | {c_1i:<9.1f} | {c_1ui:<9.1f} | {w_1i:<9.1f} | {math.degrees(beta_1i_rad):<9.2f} | {t_ca_i * 1000:<8.2f}")


c2_ui_arr = []
c2_ai_arr = []
c2_i_arr = []
w2_i_arr = []
betta2_i_arr = []
t_rk_i_arr = []
rho_i_arr = []
alpha2_i_arr = []

print(f"\n{'Сечение':<15} | {'D, мм':<8} | {'u, м/с':<8} | {'c2, м/с':<9} | {'c2u, м/с':<9} | {'w2, м/с':<9} | {'β2, град':<9} | {'ρ':<8} | {'t, мм':<8} | {'α2, град'}")
print("-" * 105)

for i, sec in enumerate(sections_rk):
    D2_i = sec["D"] / D_sr
    u2_i = u_sr * D2_i

    # Постоянство работы по высоте (Батурин п. 106): u1*c1u - u2*c2u = L_u_ср
    L_u_sr = u_sr * (c_1u - c_2u)  # Удельная работа на среднем радиусе
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

    c2_ui_arr.append(c2_ui)
    c2_ai_arr.append(c2_ai)
    c2_i_arr.append(c2_i)
    w2_i_arr.append(w2_i)
    betta2_i_arr.append(math.degrees(betta2_i))
    t_rk_i_arr.append(t_rk_i)
    rho_i_arr.append(rho_i)

    print(
        f"{sec['name']:<15} | {sec['D'] * 1000:<8.1f} | {u2_i:<8.1f} | {c2_i:<9.1f} | {-c2_ui:<9.1f} | {w2_i:<9.1f} | {math.degrees(betta2_i):<9.2f} | {rho_i:<8.3f} | {t_rk_i * 1000:<8.2f}| {alpha2i:<8.2f}")
print("c2u, β2, α2 на выходе — относительно направления, противоположного вращению (как у Батурина)")
if rho_i_arr[0] < 0.1:
    print(f"⚠️ Реактивность на втулке {rho_i_arr[0]:.3f} < 0,1 — увеличьте rho_st (Батурин п. 44)")

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