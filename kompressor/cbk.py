"""
Газодинамический расчёт центробежного компрессора (ЦБК) с безлопаточным (БЛД) и лопаточным (ЛД) диффузорами.
Методика: «Осевые и центробежные компрессоры двигателей летательных аппаратов», разд. 6.10 (с. 150…165);
μ и типы лопаток — разд. 6.5, 6.6 (с. 135, 141); углы атаки ВНА — разд. 6.11 (с. 165, 166).
Номера пунктов (п. N) — по разд. 6.10.

ВНА здесь — вращающийся направляющий аппарат (входная осевая часть рабочего колеса).
Углы β отсчитываются от окружного направления (фронта решётки), α — тоже от окружного направления.
"""
import math
from kompressor.gdf import k, R, cp, s_v, a_cr, gdf_pi, gdf_tau, gdf_eps, lambda_from_q

# Параметры, выбираемые в рекомендуемых книгой пределах
DEFAULTS = dict(
    G_bar=0.60,  # Коэффициент производительности Ḡк (п. 2: 0,5…0,7)
    k_G=0.95,  # Неравномерность поля скоростей на входе (п. 3: 0,93…0,98)
    H_bar0=0.70,  # Начальное приближение коэффициента адиабатического напора H̄ (п. 4: 0,65…0,75
                  # для радиальных лопаток); уточняется из условия π*_расч = π*_задан (с. 164, п. 7)
    d_vt=0.50,  # Относительный диаметр втулки D_вт/D1 (п. 7: 0,5…0,7)
    alpha1=90.0,  # Угол закрутки потока на входе в ВНА, град (п. 8: 90° — без закрутки)
    m_zakr=0.0,  # Показатель закона закрутки по радиусу ВНА (с. 155: −1 ≤ m ≤ 1)
    psi_zakr=1.0,  # Коэффициент неравномерности ψ в формуле затраченной работы (с. 157)
    xi1=0.25,  # Коэффициент потерь в ВНА ξ1 (п. 11: 0,2…0,3)
    xi2=0.125,  # Коэффициент потерь в рабочем колесе ξ2 (п. 11: 0,1…0,15)
    sigma_in=1.0,  # δ_вх входного патрубка (потери воздухозаборника уже учтены в термодинамике)
    z_vna=21,  # Число лопаток ВНА = числу лопаток РК (с. 156, п. 11: 19…25)
    b_t_mean=1.4,  # Густота решётки ВНА на среднем радиусе (с. 156, п. 13: 1,3…1,5)
    i_tip=0.0,  # Угол атаки ВНА на периферии, град (с. 166: 0)
    i_hub=2.5,  # Угол атаки ВНА у корня, град (с. 166: 2…3°)
    beta2l=90.0,  # Лопаточный угол на выходе из РК, град (90° — радиальные лопатки; с. 136: 50…90°)
    alpha_end0=0.06,  # Коэффициент концевых потерь α (с. 157: 0,04…0,08 для одностороннего колеса)
    beta_fr=0.85e-3,  # β в формуле проверки α (с. 159, п. 11: (0,7…1,0)·10⁻³)
    D3_D2=1.10,  # Наружный диаметр БЛД (с. 159: D3 = (1,05…1,15)·D2)
    n_bld=1.55,  # Показатель политропы в БЛД (с. 160: 1,5…1,6)
    D4_D3=1.30,  # Диаметр выхода ЛД (с. 161: D4 = (1,25…1,35)·D3)
    d_alpha4=15.0,  # Угол расширения потока в ЛД α4 − α3, град (с. 161: 12…18°)
    theta_ld=9.0,  # Средний угол расширения каналов ЛД, град (с. 161: ≤ 8…11°)
    c_k=120.0,  # Скорость на выходе из выходного устройства, м/с (с. 163: 100…140)
    xi_vu=0.5,  # Коэффициент потерь выходного устройства (с. 163)
    n_sections=5,  # Число расчётных сечений по высоте ВНА
    S_D2=0.22,  # Осевая ширина колеса S/D2 (только для рисунка)
)


def _coprime_blade_number(z_raw, z_rk, z_min=9, z_max=25):
    """Ближайшее к z_raw число лопаток ЛД в диапазоне 9…25, не кратное числу лопаток РК (с. 161)"""
    z0 = min(max(round(z_raw), z_min), z_max)
    for d in range(0, z_max):
        for z in (z0 + d, z0 - d):
            if z_min <= z <= z_max and math.gcd(z, z_rk) == 1:
                return z
    return z0


def _stage(pik, G, p_n, T_n, H_bar, p):
    """Однократный расчёт ступени при заданном H̄ (разд. 6.10)"""
    res = {"H_bar": H_bar}
    # ---------- Определение основных параметров ----------
    H = cp * T_n * (pik ** ((k - 1) / k) - 1)  # п. 1
    D1 = math.sqrt(4 * G * math.sqrt(T_n) / (math.pi * s_v * p_n * p["G_bar"] * p["k_G"]))  # п. 3
    u2 = math.sqrt(H / H_bar)  # п. 5
    a1 = math.radians(p["alpha1"])
    d = p["d_vt"]
    q1 = p["G_bar"] / ((1 - d ** 2) * math.sin(a1))  # п. 8
    lam1 = lambda_from_q(q1)  # п. 9
    c1 = lam1 * a_cr(T_n)  # п. 10
    ctg1 = math.cos(a1) / math.sin(a1)
    k_alpha = (2 / (k + 1)) ** (1 / (k - 1)) * (
        math.sqrt((ctg1 / 2) ** 2 + 2 * (1 + ctg1 ** 2 + p["xi2"] / p["xi1"])) - ctg1 / 2)  # п. 11
    lam_u2 = u2 / a_cr(T_n)
    D1_D2 = k_alpha * p["G_bar"] / (lam_u2 * gdf_eps(lam1) * (1 - d ** 2))  # п. 11
    D2 = D1 / D1_D2  # п. 12
    D_vt = D1 * d  # п. 13
    u1 = u2 * D1_D2  # п. 14
    n_rps = u1 / (math.pi * D1)  # п. 15
    c1a, c1u = c1 * math.sin(a1), c1 * math.cos(a1)
    w1 = math.sqrt(c1 ** 2 + u1 ** 2 - 2 * c1 * u1 * math.cos(a1))  # п. 16 (периферия)
    T1 = T_n - c1 ** 2 / (2 * cp)  # п. 17
    M_w1 = w1 / math.sqrt(k * R * T1)  # п. 18
    c1a_bar = c1a / u2  # п. 19
    res.update(H=H, D1=D1, u2=u2, lam1=lam1, c1=c1, k_alpha=k_alpha, D1_D2=D1_D2, D2=D2, D_vt=D_vt, u1=u1,
               n_rpm=n_rps * 60, c1a=c1a, c1u=c1u, w1=w1, T1=T1, M_w1=M_w1, c1a_bar=c1a_bar)

    # ---------- Входной патрубок и ВНА: параметры по высоте (с. 154…156) ----------
    p1_star = p["sigma_in"] * p_n
    T1_star = T_n
    t_mean = math.pi * (D1 + D_vt) / 2 / p["z_vna"]
    chord = p["b_t_mean"] * t_mean  # хорда ВНА постоянна по высоте (п. 13, 14)

    def inlet_section(D_i):
        Dbar = D_i / D1
        if abs(p["alpha1"] - 90.0) < 1e-9:
            a1i, c1ai, c1ui = math.pi / 2, c1a, 0.0
        else:
            a1i = math.atan(Dbar ** p["m_zakr"] * math.tan(a1))
            c1ai = c1 * math.sin(a1) / Dbar ** ((1 - p["m_zakr"]) * math.cos(a1) ** 2)
            c1ui = c1ai / math.tan(a1i)
        u1i = u1 * Dbar
        c1i = c1ai / math.sin(a1i)
        lam_i = c1i / a_cr(T1_star)
        p1i, T1i = p1_star * gdf_pi(lam_i), T1_star * gdf_tau(lam_i)
        w1i = math.hypot(c1ai, u1i - c1ui)
        return dict(D=D_i, u=u1i, c=c1i, ca=c1ai, cu=c1ui, alpha=math.degrees(a1i), lam=lam_i, p=p1i, T=T1i,
                    rho=p1i / (R * T1i), w=w1i, beta=math.degrees(math.atan2(c1ai, u1i - c1ui)),
                    M_w=w1i / math.sqrt(k * R * T1i))

    sections = []
    n_sec = p["n_sections"]
    for j in range(n_sec):
        h_rel = j / (n_sec - 1)
        s = inlet_section(D_vt + (D1 - D_vt) * h_rel)
        s["h_rel"] = h_rel
        s["i"] = p["i_hub"] + (p["i_tip"] - p["i_hub"]) * h_rel  # с. 166: 2…3° у корня, 0 на периферии
        s["beta_l"] = s["beta"] + s["i"]
        s["theta"] = 90 - s["beta_l"]  # угол изгиба профиля ВНА (с. 165)
        s["t"] = math.pi * s["D"] / p["z_vna"]
        s["b_t"] = chord / s["t"]
        s["R_ml"] = chord / (2 * math.sin(math.radians(s["theta"]) / 2)) if s["theta"] > 0 else float("inf")
        sections.append(s)
    sr = inlet_section((D1 + D_vt) / 2)  # средний радиус входа
    res.update(sections=sections, sr=sr, chord_vna=chord, p1_star=p1_star)

    # ---------- Рабочее колесо (с. 157…159) ----------
    z_min = 2 * math.pi * u2 / sr["ca"]  # п. 1
    z = max(p["z_vna"], math.ceil(z_min))
    b2l = math.radians(p["beta2l"])
    x = 0.6 * (1 + math.sin(b2l))
    mu = 1 / (1 + 2 * x / z / (1 - D1_D2 ** 2))  # п. 2
    c2r = sr["ca"]  # с. 157: c2r ≈ c1a,ср
    c2u = mu * (u2 - c2r / math.tan(b2l))  # п. 4; μ = c2u/c2u∞, c2u∞ = u2 − c2r·ctg β2л (с. 135, 141)
    alpha_end = p["alpha_end0"]
    for _ in range(50):
        H_z = u2 * c2u + alpha_end * u2 ** 2 - p["psi_zakr"] * c1u * u2 * D1_D2 ** 2  # п. 3
        c2 = math.hypot(c2r, c2u)  # п. 5
        T2 = sr["T"] + (sr["c"] ** 2 - c2 ** 2) / (2 * cp) + H_z / cp  # п. 6
        # п. 7; w1, c1a — на среднем радиусе входа, как и остальные параметры «ср» в расчёте РК (с. 158)
        sum_HR = alpha_end * u2 ** 2 + p["xi1"] * sr["w"] ** 2 / 2 + p["xi2"] * sr["ca"] ** 2 / 2
        dE_z = H_z - (c2 ** 2 - sr["c"] ** 2) / 2
        if T2 <= sr["T"] or dE_z <= sum_HR:
            raise ValueError("РК не сжимает воздух при данных параметрах")
        nk = k / (k - 1) * (1 - sum_HR / dE_z)  # n/(n−1)
        p2 = sr["p"] * (T2 / sr["T"]) ** nk  # п. 8
        rho2 = p2 / (R * T2)
        b2 = G / (math.pi * rho2 * D2 * c2r)  # п. 10
        alpha_new = p["beta_fr"] / (c1a_bar * b2 / D2)  # п. 11
        if abs(alpha_new - alpha_end) < 1e-5:
            break
        alpha_end = alpha_new
    alpha2 = math.degrees(math.atan2(c2r, c2u))  # п. 9
    w2u = u2 - c2u
    w2 = math.hypot(c2r, w2u)
    beta2 = math.degrees(math.atan2(c2r, w2u))
    T2_star = T2 + c2 ** 2 / (2 * cp)
    lam_c2 = c2 / a_cr(T2_star)
    p2_star = p2 / gdf_pi(lam_c2)
    H_ad_rk = cp * sr["T"] * ((p2 / sr["p"]) ** ((k - 1) / k) - 1)  # п. 12
    eta_rk = (H_ad_rk + (c2 ** 2 - sr["c"] ** 2) / 2) / H_z  # п. 13
    rho_st = 1 - (c2 ** 2 - sr["c"] ** 2) / 2 / (H_z - alpha_end * u2 ** 2)  # степень реактивности
    res.update(z=z, z_min=z_min, mu=mu, c2r=c2r, c2u=c2u, c2=c2, alpha_end=alpha_end, H_z=H_z, T2=T2, p2=p2,
               rho2=rho2, b2=b2, n_poly_rk=nk / (nk - 1), alpha2=alpha2, w2=w2, w2u=w2u, beta2=beta2,
               T2_star=T2_star, p2_star=p2_star, lam_c2=lam_c2, M_c2=c2 / math.sqrt(k * R * T2),
               H_ad_rk=H_ad_rk, eta_rk=eta_rk, rho_st=rho_st)

    # ---------- Безлопаточный диффузор (с. 159, 160) ----------
    D3 = p["D3_D2"] * D2
    b3 = b2
    c3 = c2 * D2 * b2 / (D3 * b3)  # п. 2
    T3 = T2 + (c2 ** 2 - c3 ** 2) / (2 * cp)  # п. 3
    T3_star = T3 + c3 ** 2 / (2 * cp)  # п. 4
    M3 = c3 / math.sqrt(k * R * T3)  # п. 5
    p3 = p2 * (T3 / T2) ** (p["n_bld"] / (p["n_bld"] - 1))  # п. 7
    rho3 = p3 / (R * T3)
    p3_star = p3 * (1 + (k - 1) / 2 * M3 ** 2) ** (k / (k - 1))  # п. 8
    alpha3 = math.degrees(math.atan(b2 * rho2 / (b3 * rho3) * math.tan(math.radians(alpha2))))  # п. 9
    res.update(D3=D3, b3=b3, c3=c3, T3=T3, T3_star=T3_star, M3=M3, p3=p3, rho3=rho3, p3_star=p3_star,
               alpha3=alpha3, sigma_bld=p3_star / p2_star)

    # ---------- Лопаточный диффузор (с. 161…163) ----------
    D4 = p["D4_D3"] * D3  # п. 1
    alpha4 = alpha3 + p["d_alpha4"]  # п. 2
    b4 = b3
    F4 = math.pi * D4 * b4 * math.sin(math.radians(alpha4))  # п. 3
    F3 = math.pi * D3 * b3 * math.sin(math.radians(alpha3))  # п. 4
    z_ld_raw = 360 / p["theta_ld"] * (F4 / F3 - 1) / (D4 / D3 - 1) * math.sin(math.radians(alpha3)) ** 2  # п. 5
    z_ld = _coprime_blade_number(z_ld_raw, z)
    xi_ld = 0.17 if M3 <= 0.7 else 0.22  # п. 11
    delta_ld = 1.0  # п. 6
    for _ in range(50):
        p4_star = delta_ld * p3_star
        lam4 = lambda_from_q(G * math.sqrt(T3_star) / (s_v * p4_star * F4))  # п. 7, 8
        c4 = lam4 * a_cr(T3_star)  # п. 9
        T4 = T3 + (c3 ** 2 - c4 ** 2) / (2 * cp)  # п. 10
        dE_ld = (c3 ** 2 - c4 ** 2) / 2
        if dE_ld <= 0:
            raise ValueError(f"c4 = {c4:.0f} ≥ c3 = {c3:.0f} м/с — ЛД не тормозит поток")
        n_ld = k / (k - 1) * (1 - xi_ld * c3 ** 2 / 2 / dE_ld)  # п. 11
        p4 = p3 * (T4 / T3) ** n_ld  # п. 12
        p4_star_new = p4 / gdf_pi(lam4)  # п. 13 (p4* = p4/π(λ4))
        delta_new = p4_star_new / p3_star  # п. 14
        if abs(delta_new - delta_ld) < 1e-6:
            break
        delta_ld = delta_new
    res.update(D4=D4, b4=b4, alpha4=alpha4, F3=F3, F4=F4, z_ld_raw=z_ld_raw, z_ld=z_ld, xi_ld=xi_ld, lam4=lam4,
               c4=c4, T4=T4, p4=p4, p4_star=p4_star_new, sigma_ld=delta_ld, rho4=p4 / (R * T4))

    # ---------- Выходное устройство (с. 163, 164) ----------
    c_k = p["c_k"]
    T_k = T4 + (c4 ** 2 - c_k ** 2) / (2 * cp)  # п. 2
    dE_vu = (c4 ** 2 - c_k ** 2) / 2
    if dE_vu <= 0:
        raise ValueError(f"c_к = {c_k:.0f} ≥ c4 = {c4:.0f} м/с (с. 163: должно быть c_к < c4)")
    n_vu = k / (k - 1) * (1 - p["xi_vu"] * c4 ** 2 / 2 / dE_vu)  # п. 3
    p_k = p4 * (T_k / T4) ** n_vu  # п. 4
    T_k_star = T3_star  # п. 5
    lam_k = c_k / a_cr(T_k_star)
    p_k_star = p_k / gdf_pi(lam_k)  # п. 6
    pi_calc = p_k_star / p_n  # п. 7
    eta_k = H / H_z  # п. 8
    res.update(c_k=c_k, T_k=T_k, p_k=p_k, T_k_star=T_k_star, p_k_star=p_k_star, pi_calc=pi_calc, eta_k=eta_k,
               N=H_z * G, G=G, pik=pik, p_n=p_n, T_n=T_n, S=p["S_D2"] * D2)
    return res


def compressor(pik, G, p_n, T_n, **params):
    """
    Расчёт ЦБК на заданные π*_к, расход G (кг/с) и полные параметры на входе p_н* (Па), T_н* (К).
    H̄ подбирается так, чтобы расчётная π*_к совпала с заданной (с. 164, п. 7): вместо ручной корректировки
    на 2 % решается уравнение π*_расч(H̄) = π*_задан методом секущих.
    """
    p = dict(DEFAULTS)
    unknown = set(params) - set(p)
    if unknown:
        raise KeyError(f"Неизвестные параметры компрессора: {unknown}")
    p.update(params)

    def f(H_bar):
        return _stage(pik, G, p_n, T_n, H_bar, p)["pi_calc"] - pik

    x0, x1 = p["H_bar0"], p["H_bar0"] - 0.02
    f0, f1 = f(x0), f(x1)
    for _ in range(50):
        x2 = x1 - f1 * (x1 - x0) / (f1 - f0)
        x2 = min(max(x2, 0.35), 0.98)
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
        if abs(f1) / pik < 1e-5:
            break
    res = _stage(pik, G, p_n, T_n, x1, p)
    res["params"] = p
    res["warnings"] = _checks(res, p)
    return res


def _material(T_k_star):
    """Материал диска рабочего колеса по температуре на выходе (с. 152, п. 6)"""
    if T_k_star < 520:
        return "алюминиевый сплав"
    if T_k_star <= 770:
        return "титановый сплав"
    return "жаропрочная сталь или сплав"


def _checks(r, p):
    """Проверки по рекомендациям книги. Возвращает список предупреждений."""
    w = []

    def rng(name, val, lo, hi, ref):
        if not lo <= val <= hi:
            w.append(f"{name} = {val:.3f} вне рекомендуемых {lo}…{hi} ({ref})")

    rng("u2, м/с", r["u2"], 400, 550, "с. 151, п. 5")
    rng("D1/D2", r["D1_D2"], 0.55, 0.70, "с. 153, п. 11")
    rng("c1, м/с", r["c1"], 120, 150, "с. 152, п. 10")
    rng("c̄1a = c1a/u2", r["c1a_bar"], 0.25, 0.35, "с. 154, п. 19")
    rng("α2, град", r["alpha2"], 14, 18, "с. 158, п. 9")
    rng("α концевых потерь", r["alpha_end"], 0.04, 0.08, "с. 157, п. 3")
    rng("η_к", r["eta_k"], 0.78, 0.83, "с. 164, п. 8")
    rng("H̄/η_к", r["H_bar"] / r["eta_k"], 0.85, 1.1, "с. 164, п. 9")
    rng("Число лопаток ЛД", r["z_ld_raw"], 9, 25, "с. 161, п. 5")
    if abs(p["beta2l"] - 90) < 1e-9:
        rng("H̄", r["H_bar"], 0.65, 0.75, "с. 151, п. 4, радиальные лопатки")
    if r["M_w1"] > 0.9:
        w.append(f"M_w1 = {r['M_w1']:.3f} > 0,9 — нужна предварительная закрутка (α1 < 90°) (с. 154, п. 18)")
    if r["z_min"] > p["z_vna"]:
        w.append(f"z_ВНА = {p['z_vna']} < z_min = {r['z_min']:.1f}: число лопаток РК увеличено до {r['z']} "
                 f"(с. 157, п. 1)")
    if r["z"] > 25:
        w.append(f"z_РК = {r['z']} > 25 — много лопаток на входе: рассмотрите промежуточные (укороченные) лопатки")
    if r["M3"] > 0.95:
        w.append(f"M3 = {r['M3']:.3f} > 0,95 (с. 160, п. 5)")
    if r["c_k"] >= r["c4"]:
        w.append(f"c_к = {r['c_k']:.0f} ≥ c4 = {r['c4']:.0f} м/с — в выходном устройстве должно быть c_к < c4 (с. 163)")
    if r["eta_k"] >= r["eta_rk"]:
        w.append(f"η_к = {r['eta_k']:.3f} ≥ η_РК = {r['eta_rk']:.3f} — должно быть η_к < η_РК (с. 164, п. 8)")
    if r["b2"] < 5e-3:
        w.append(f"b2 = {r['b2'] * 1000:.2f} мм < 5 мм — высота лопатки на выходе мала")
    return w


def print_report(r):
    p = r["params"]
    print("=" * 78)
    print("   ГАЗОДИНАМИЧЕСКИЙ РАСЧЁТ ЦЕНТРОБЕЖНОГО КОМПРЕССОРА (разд. 6.10)")
    print("=" * 78)
    print(f"Исходные данные: π*_к = {r['pik']:.3f}, G = {r['G']:.3f} кг/с, p_н* = {r['p_n'] / 1000:.2f} кПа, "
          f"T_н* = {r['T_n']:.2f} К")

    print("\n─── Основные параметры ───")
    print(f"  Изоэнтропная работа H = {r['H'] / 1000:.2f} кДж/кг")
    print(f"  H̄ = {r['H_bar']:.4f} (подобран по π*), u2 = {r['u2']:.1f} м/с (400…550)")
    print(f"  D1 = {r['D1'] * 1000:.1f} мм, D_вт = {r['D_vt'] * 1000:.1f} мм (d̄_вт = {p['d_vt']}), "
          f"D2 = {r['D2'] * 1000:.1f} мм")
    print(f"  (D1/D2)_опт = {r['D1_D2']:.3f} (0,55…0,70), k_α = {r['k_alpha']:.3f}")
    print(f"  n = {r['n_rpm']:.0f} об/мин, u1 = {r['u1']:.1f} м/с")
    print(f"  λ1 = {r['lam1']:.3f}, c1 = {r['c1']:.1f} м/с (120…150), c̄1a = {r['c1a_bar']:.3f} (0,25…0,35)")
    print(f"  w1 (периферия) = {r['w1']:.1f} м/с, M_w1 = {r['M_w1']:.3f} (≤ 0,8…0,9)")

    print("\n─── ВНА: параметры и углы по высоте (β от окружного направления) ───")
    print(f"  z_ВНА = {p['z_vna']}, хорда b = {r['chord_vna'] * 1000:.1f} мм")
    print(f"  {'h/h':<5} | {'D, мм':<7} | {'u1, м/с':<8} | {'c1a, м/с':<8} | {'w1, м/с':<8} | {'M_w1':<6} | "
          f"{'β1, °':<7} | {'i, °':<5} | {'β1л, °':<7} | {'θ, °':<6} | {'b/t':<5} | {'R_ср.л, мм':<9}")
    print("  " + "-" * 110)
    for s in r["sections"]:
        print(f"  {s['h_rel']:<5.2f} | {s['D'] * 1000:<7.1f} | {s['u']:<8.1f} | {s['ca']:<8.1f} | {s['w']:<8.1f} | "
              f"{s['M_w']:<6.3f} | {s['beta']:<7.2f} | {s['i']:<5.2f} | {s['beta_l']:<7.2f} | {s['theta']:<6.2f} | "
              f"{s['b_t']:<5.2f} | {s['R_ml'] * 1000:<9.1f}")
    print("  β1л = β1 + i; θ = 90° − β1л — угол изгиба лопатки ВНА (выход ВНА осевой-радиальный, β = 90°)")

    print("\n─── Рабочее колесо ───")
    print(f"  z_РК = {r['z']} (z_min = {r['z_min']:.1f}), β2л = {p['beta2l']:.1f}°, μ = {r['mu']:.4f}")
    print(f"  α концевых потерь = {r['alpha_end']:.4f} (0,04…0,08), H_z = {r['H_z'] / 1000:.2f} кДж/кг")
    print(f"  c2r = {r['c2r']:.1f} м/с, c2u = {r['c2u']:.1f} м/с, c2 = {r['c2']:.1f} м/с, λ_c2 = {r['lam_c2']:.3f}")
    print(f"  w2 = {r['w2']:.1f} м/с, β2 (поток) = {r['beta2']:.2f}°, α2 = {r['alpha2']:.2f}° (14…18°)")
    print(f"  T2 = {r['T2']:.1f} К, p2 = {r['p2'] / 1000:.1f} кПа, T2* = {r['T2_star']:.1f} К, "
          f"p2* = {r['p2_star'] / 1000:.1f} кПа")
    print(f"  Показатель политропы n = {r['n_poly_rk']:.3f}, ρ2 = {r['rho2']:.3f} кг/м³")
    print(f"  Ширина лопаток на выходе b2 = {r['b2'] * 1000:.2f} мм (b2/D2 = {r['b2'] / r['D2']:.4f})")
    print(f"  η_РК = {r['eta_rk']:.4f}, степень реактивности ρ = {r['rho_st']:.3f}")

    print("\n─── Безлопаточный диффузор ───")
    print(f"  D3 = {r['D3'] * 1000:.1f} мм, b3 = {r['b3'] * 1000:.2f} мм, c3 = {r['c3']:.1f} м/с, "
          f"M3 = {r['M3']:.3f} (≤ 0,95), α3 = {r['alpha3']:.2f}°")
    print(f"  p3* = {r['p3_star'] / 1000:.1f} кПа, σ_БЛД = {r['sigma_bld']:.4f}")

    print("\n─── Лопаточный диффузор ───")
    print(f"  D4 = {r['D4'] * 1000:.1f} мм, α3 = {r['alpha3']:.2f}° → α4 = {r['alpha4']:.2f}°, "
          f"z_ЛД = {r['z_ld']} (по формуле {r['z_ld_raw']:.1f})")
    print(f"  ξ_ЛД = {r['xi_ld']}, λ4 = {r['lam4']:.3f}, c4 = {r['c4']:.1f} м/с, "
          f"p4* = {r['p4_star'] / 1000:.1f} кПа, σ_ЛД = {r['sigma_ld']:.4f}")

    print("\n─── Выходное устройство и компрессор в целом ───")
    print(f"  c_к = {r['c_k']:.0f} м/с, T_к* = {r['T_k_star']:.1f} К, p_к* = {r['p_k_star'] / 1000:.1f} кПа")
    print(f"  π*_к расч = {r['pi_calc']:.4f} (задано {r['pik']:.4f})")
    print(f"  η_к = H/H_z = {r['eta_k']:.4f} (0,78…0,83), H̄/η_к = {r['H_bar'] / r['eta_k']:.3f} (0,85…1,1)")
    print(f"  Мощность привода N = H_z·G = {r['N'] / 1000:.1f} кВт")
    print(f"  Материал колеса при T_к* = {r['T_k_star']:.0f} К: {_material(r['T_k_star'])} (с. 152, по температуре)")
    if _material(r["T_k_star"]) == "алюминиевый сплав" and r["u2"] > 450:
        print(f"  ⚠️ При u2 = {r['u2']:.0f} м/с алюминиевое колесо обычно не проходит по прочности (запас k_m ≈ 2, "
              f"с. 151, п. 6) — принять титановый сплав")

    if r["warnings"]:
        print("\n⚠️ Замечания:")
        for msg in r["warnings"]:
            print("   - " + msg)
    else:
        print("\n✅ Все параметры в рекомендуемых пределах")
