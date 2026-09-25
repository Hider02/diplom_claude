import numpy as np
from gdf import *
# ================= ШАГИ 58-70: БЕЗЛОПАТОЧНЫЙ ДИФФУЗОР =================
# 58
def diff_calculation(h2, D2, c2r, c2u, T2_star, alpha2, p2, T2, eta_k, p2_star, c2, P0, dh_prof, dh_out, Lu, beta_tr):
    b2_D2 = h2 / D2
    #if b2_D2 > 0.06: h3_h2 = 0.77
    #elif 0.04 <= b2_D2 <= 0.06: h3_h2 = 1.0
    #else: h3_h2 = 1.15
    h3_h2=1
    h3 = h2 * h3_h2
    # 59-60
    D2_pr = 1.02 * D2
    D3 = 1.25 * D2
    F3r = np.pi*D3*h3
    # 62
    c3r = c2r * D2 * h2 / (D3 * h3)
    # 63
    rho3_rho2 = 1.02
    alpha3 = np.degrees(np.arctan(np.tan(np.radians(alpha2)) / ((h3/h2) * rho3_rho2)))
    # 64
    c3u = c3r / np.tan(np.radians(alpha3))
    c3 = np.sqrt(c3u**2 + c3r**2)
    # 65-67
    T3_star = T2_star
    T3 = T3_star - c3**2 / (2.0 * cp)
    lambda_c3 = c3 / np.sqrt(2.0*k/(k+1.0)*R*T3_star)
    p3 = p2 * (T3 / T2) ** (k / (k - 1.0) * eta_k)
    p3_star = p3 / gdf_pi(lambda_c3)
    # 68-70
    sigma_bld = p3_star / p2_star
    rho3 = p3 / (R * T3)

    # ================= ШАГ 71: ПОТЕРИ В БЛД =================
    # 71
    tg_half_υ = 2.0 * np.sqrt(h3/D3) / (1.0 + np.sqrt(D2/D3)) * np.sin(np.radians(alpha3))
    υ_bld = 2.0 * np.degrees(np.arctan(tg_half_υ)) #!!!!!
    print(f"υ_bld = {υ_bld}")
    xi_bld = 0.147 + 0.0046 * (υ_bld - 12.0)**2
    dh_bld = xi_bld * c2**2 / 2.0

    # ================= ШАГИ 72-89: ЛОПАТОЧНЫЙ ДИФФУЗОР =================
    # 72-73
    alpha3l = (alpha2 + alpha3) / 2.0
    alpha4l = alpha3l + 12.5
    # 74-75
    b_t3 = 2
    delta_alpha4l = 0.346 * (alpha4l - alpha3l) / b_t3
    alpha4 = alpha4l - delta_alpha4l
    # 76-77
    D4 = 1.45 * D2
    h4 = h3
    F4r = np.pi*D4*h4
    # 79-81
    c4r = c3r*D3 * h3 / (D4 * h4)
    c4u = c4r / np.tan(np.radians(alpha4))
    c4 = np.sqrt(c4u**2 + c4r**2)
    # 82-84
    T4_star = T3_star
    T4 = T4_star - c4**2 / (2.0 * cp)
    lambda_c4 = c4 / np.sqrt(2.0*k/(k+1.0)*R*T4_star)
    p4 = p3 * (T4 / T3) ** (k / (k - 1.0) * eta_k)
    p4_star = p4 / gdf_pi(lambda_c4)
    # 85-87
    pi_k_star = p4_star / P0
    sigma_ld = p4_star / p3_star
    rho4 = p4 / (R * T4)
    # 88
    z_ld = (b_t3) * 2.0*np.pi *np.sin(np.radians((alpha3+alpha4)/2.0)) / np.log(D4/D3)
    z_ld2 = 360*((F4r/F3r)-1)*np.sin(np.radians(alpha3))**2/(10*(D4/D3-1))
    z_ld = int(np.round(z_ld))
    if z_ld < 8: z_ld = 12
    # 89
    F_l = F4r/F3r
    ld = (D4-D3)/(np.sin(np.radians((alpha3+alpha4)/2.0))*2)
    υ_ld = 2*np.degrees(np.arctan((np.sqrt(F4r)-np.sqrt(F3r))/(ld*np.sqrt(np.pi))))
    # u_ld2 = 2.0 * np.degrees(np.arctan(2.0 * np.sqrt(h4/D4) / (1.0 + np.sqrt(D4/D3)) * np.sin(np.radians(alpha4))))
    k_f = 1.7 +0.03*υ_ld
    xi_ld0 = 3.5*k_f*((np.tan(np.radians(υ_ld/2)))**1.25)*(1-1/F_l)**1.65
    xi_ld = xi_ld0*(1+4.3*(lambda_c3-0.8)**2)
    dh_ld = xi_ld * c3**2 / 2.0

    # ================= ШАГ 90: КПД СТУПЕНИ И СХОДИМОСТЬ =================
    eta_k_calc = 1.0 - (dh_prof + dh_out+dh_bld+dh_ld+ Lu*(beta_tr-1)) / (Lu*(beta_tr))
    print(f"КПД ступени = {eta_k_calc}")
    print(f"Длинна БЛД = {(D3-D2)/2*1000}")
    print(f"D4 = {D4}")
    print(f'z = {z_ld}')
