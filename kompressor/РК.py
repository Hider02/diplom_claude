import numpy as np
from kompressor.plot import plot_comp
from gdf import *


def PK_calculation(G_v, pik, P0, T0, eta_k, H_t, n, D1_k, D1_vt, beta2l, S_coef, t_k, t_vt, verbose=True):
    # ================= ШАГИ 1-11: БАЗОВЫЕ ПАРАМЕТРЫ =================
    Lks = cp * T0 * (pik ** ((k - 1) / k) - 1.0)
    Lz = Lks / eta_k
    Pk = P0 * pik
    Tk = T0 + Lz / cp

    u2 = np.sqrt(Lz / H_t)
    D2 = 60.0 * u2 / (np.pi * n)

    D1 = D1_k * D2
    D_vt = D1_vt * D2
    D1_cp = (D1 + D_vt) / 2.0
    h1 = (D1 - D_vt) / 2.0
    S = S_coef * D2
    u1_cp = np.pi * D1_cp * n / 60.0

    rho0 = P0 / (R * T0)
    F_p = 4.0 * G_v / (np.pi * rho0 * D2 ** 2 * u2)

    # Определяем сплиттеры ДО изменения z_rk
    z_rk_raw = int(np.round(beta2l / 4.0 + (105.0 - beta2l) * (beta2l - 10.0) / 200.0))
    use_splitters = z_rk_raw > 15

    z_rk = z_rk_raw
    if z_rk < 8: z_rk = 12
    if z_rk > 15: z_rk = int(np.ceil(z_rk / 2))

    m_gdf = 0.0404  # Для воздуха k=1.4

    # ================= ШАГИ 12-22: ВХОД В РК (ИТЕРАЦИЯ) =================
    c1u = 0.0
    beta1_cp = 35.0
    for _ in range(10):
        F1a = np.pi / 4.0 * (D1 ** 2 - D_vt ** 2) - z_rk * 0.5 * h1 * (t_vt + t_k) / np.sin(np.radians(beta1_cp))
        if F1a <= 0: raise ValueError("F1a <= 0")  # Защита от деления на ноль/корня из минуса

        q_c1a = G_v * np.sqrt(T0) / (m_gdf * P0 * F1a)
        lambda_c1a = solve_lambda(q_c1a)
        if lambda_c1a < 0: raise ValueError("lambda_c1a < 0")

        c1a = lambda_c1a * np.sqrt(2.0 * k / (k + 1.0) * R * T0)
        c1 = np.sqrt(c1a ** 2 + c1u ** 2)
        lambda_c1 = c1 / np.sqrt(2.0 * k / (k + 1.0) * R * T0)

        p1 = P0 * gdf_pi(lambda_c1)
        T1 = T0 * gdf_tau(lambda_c1)
        rho1 = p1 / (R * T1)

        w1u = u1_cp - c1u if c1u >= 0 else u1_cp + abs(c1u)
        beta1_new = np.degrees(np.arctan2(c1a, w1u))
        if abs(beta1_new - beta1_cp) < 0.1: break
        beta1_cp = beta1_new

    # ================= ШАГИ 23-25: ЗАТОРМОЖЕННЫЕ ПАРАМЕТРЫ =================
    w1 = np.sqrt(w1u ** 2 + c1a ** 2)
    T1w_star = T1 + w1 ** 2 / (2.0 * cp)
    lambda_w1 = w1 / np.sqrt(2.0 * k / (k + 1.0) * R * T1w_star)
    p1w_star = p1 / gdf_pi(lambda_w1)

    # ================= ШАГИ 26-35: РАСЧЁТ ПО ВЫСОТЕ НА ВХОДЕ =================
    num_sections = 3
    D_sections = np.linspace(D_vt, D1, num_sections)

    beta1_tip_l = 0.0
    w1_array = np.zeros(num_sections)
    betta1l_array = np.zeros(num_sections)
    for i in range(num_sections):
        D_i = D_sections[i]
        u_i = np.pi * D_i * n / 60.0
        c1ai = c1a
        c1ui = c1u * (D1_cp / D_i)
        c1i = np.sqrt(c1ui ** 2 + c1ai ** 2)
        w1ui = u_i - c1ui if c1ui >= 0 else u_i + abs(c1ui)
        w1i = np.sqrt(w1ui ** 2 + c1ai ** 2)
        w1_array[i] = w1i

        T1wi_star = T1 + w1i ** 2 / (2.0 * cp)
        lambda_w1i = w1i / np.sqrt(2.0 * k / (k + 1.0) * R * T1wi_star)
        beta1i = np.degrees(np.arctan2(c1ai, w1ui))

        i_attack = 2.0
        beta1li = beta1i + i_attack
        betta1l_array[i] = beta1li
        if i == num_sections - 1:
            beta1_tip_l = beta1li

    # ================= ШАГИ 36-53: ВЫХОД ИЗ РК (ИТЕРАЦИЯ β_тр) =================
    beta_tr = 1.03
    for _ in range(5):
        c2u = (Lz / (1.0 + beta_tr) - c1u * u1_cp) / u2
        mu = 1.0 - np.sqrt(np.sin(np.radians(beta2l))) / (z_rk ** 0.7)
        if mu <= 0.1: mu = 0.1

        c2u_inf = c2u / mu
        w2u_inf = u2 - c2u_inf
        c2r = w2u_inf * np.tan(np.radians(beta2l))
        if c2r <= 0: raise ValueError("c2r <= 0")

        alpha2 = np.degrees(np.arctan2(c2r, c2u))
        w2u = u2 - c2u
        w2 = np.sqrt(w2u ** 2 + c2r ** 2)
        c2 = np.sqrt(c2u ** 2 + c2r ** 2)

        T2 = T1 + (Lz - (c2 ** 2 - c1 ** 2) / 2.0) / cp
        if T2 <= T1: raise ValueError("T2 <= T1")

        p2 = p1 * (T2 / T1) ** (k / (k - 1.0) * eta_k)
        rho2 = p2 / (R * T2)

        T2_star = T2 + c2 ** 2 / (2.0 * cp)
        lambda_c2 = c2 / np.sqrt(2.0 * k / (k + 1.0) * R * T2_star)
        p2_star = p2 / gdf_pi(lambda_c2)

        F2a = G_v / (c2r * rho2)
        sin_b2 = np.sin(np.radians(beta2l))
        if sin_b2 < 0.1: sin_b2 = 0.1

        if z_rk < 15:
            h2 = F2a / ((np.pi * D2) - z_rk * 2 * (t_k + t_vt) / (2 * sin_b2))
        else:
            h2 = F2a / ((np.pi * D2) - z_rk * (t_k + t_vt) / (2 * sin_b2))

        if h2 <= 0: raise ValueError("h2 <= 0")

        delta_beta = 0.172e-3 / (H_t * (c2r / u2) * (h2 / D2))
        beta_tr_new = 1.0 + delta_beta
        if abs(beta_tr_new - beta_tr) < 0.005: break
        beta_tr = beta_tr_new

    # ================= ШАГИ 54-57: ПОТЕРИ И КПД РК =================
    rho_st = 1.0 - (c2 ** 2 - c1 ** 2) / (2.0 * u2 * c2u)
    Ro = u2 / w1_array[2]

    k_E = 0.75 if use_splitters else 0.6

    term1 = 1.0 - w2 / w1_array[2]
    term2 = (k_E * H_t * Ro) / ((z_rk / np.pi) * (1.0 - D1 / D2) + 2.0 * D1 / D2)
    Df = term1 + term2

    xi_prof = 0.1 * Ro ** 2 * Df ** 2
    dh_prof = xi_prof * w1_array[1] ** 2 / 2.0

    dh_out = 0.5 * (c2r - F_p * u2 / (4.0 * rho2 / rho1 * h2 / D2)) ** 2

    Lu = c2u * u2 - c1u * u1_cp
    eta_rk = 1.0 - (dh_prof + dh_out + Lu * (beta_tr - 1)) / (Lu * (beta_tr))

    # Печать результатов только если запущен одиночный расчет (verbose=True)
    if verbose:
        print("=" * 65)
        print("ПОЛНЫЙ РАСЧЁТ СТУПЕНИ ЦБК (ШАГИ 1-90)")
        print("=" * 65)
        print(f"1-5  | u2={u2:.0f} м/с, D2={D2 * 1000:.1f} мм, Lz={Lz:.0f} Дж/кг")
        print(f"10   | F_p={F_p:.3f} (норма: 0.05-0.12)")
        print(f"11   | z_рк={z_rk}, μ={mu:.3f}")
        print(f"22   | β1_ср={beta1_cp:.1f}°, λw1_ср={lambda_w1:.3f}")
        print(f"35   | β1л_периф={beta1_tip_l:.1f}° (лимит ≤25°)")
        print(f"35   | β1л_вт={betta1l_array[0]:.1f}° ")
        print(f"41   | α2={alpha2:.1f}°, c2r={c2r:.0f} м/с")
        print(f"52   | h2={h2 * 1000:.1f} мм (лимит >5 мм) {h2 / D2:.3f}<0.15 ")
        print(f"54-57| η_рк={eta_rk:.3f}, ρ_ст={rho_st:.3f}")
        print(
            f"D1 ={D1 * 1000:.1f} мм,D_vt = {D_vt * 1000:.1f} мм, D2 ={D2 * 1000:.1f} мм, h1={h1 * 1000:.1f} мм, h2={h2 * 1000:.1f} мм, S = {S * 1000:.1f} мм")
        plot_comp(S, r_vt=D_vt/2, r1=D1/2, r2=D2/2, h1=h1, h2=h2)

    # Возвращаем ВСЕ параметры, которые могут понадобиться оптимизатору
    return {
        'h2': h2, 'D2': D2, 'c2r': c2r, 'c2u': c2u, 'T2_star': T2_star,
        'alpha2': alpha2, 'p2': p2, 'T2': T2, 'eta_k': eta_k, 'p2_star': p2_star,
        'c2': c2, 'dh_prof': dh_prof, 'dh_out': dh_out, 'Lu': Lu, 'beta_tr': beta_tr,
        # --- Новые параметры для оптимизации ---
        'eta_rk': eta_rk, 'u2': u2, 'rho_st': rho_st, 'D1': D1, 'D_vt': D_vt,
        'h1': h1, 'S': S, 'beta1_tip_l': beta1_tip_l, 'lambda_w1': lambda_w1,
        'n': n, 'G_v': G_v, 'z_rk': z_rk, 'use_splitters': use_splitters, 'mu': mu
    }