"""
Газодинамический расчёт центробежного компрессора (ЦБК) по средней линии тока:
осевой вход без входного направляющего аппарата → рабочее колесо (РК) с загнутыми назад лопатками →
радиальный безлопаточный диффузор (БЛД) → переходное колено (поворот потока из радиального направления в осевое) →
осевой лопаточный диффузор (ЛД, спрямляющие венцы) → камера сгорания.

Источники:
[К]  «Осевые и центробежные компрессоры двигателей летательных аппаратов», разд. 6.10 (с. 150…165) —
     последовательность расчёта, расчёт входа по сечениям, рекомендуемые диапазоны параметров.
[Oh] Oh H.W., Yoon E.S., Chung M.K. An optimum set of loss models for performance prediction of centrifugal
     compressors. Proc. IMechE, Part A, 1997, 211, 331–338 — набор моделей потерь РК:
     скольжение — Wiesner (1967); нагрузка лопаток — Coppage (1956); трение и зазор — Jansen (1967);
     смешение — Johnston & Dean (1966); трение диска — Daily & Nece (1960); рециркуляция — Oh (1997);
     утечки — Aungier (1995); удар (incidence) — Conrad (1980). Сводка формул: Gong X., Chen R. Total pressure loss
     mechanism of centrifugal compressors. Mechanical Engineering Research, 2014, 4(2), 45–59.
[L]  Lieblein S. Loss and stall analysis of compressor cascades. ASME J. Basic Eng., 1959 (NASA SP-36) —
     эквивалентный диффузорный фактор и толщина потери импульса для венцов осевого ЛД.

Углы: β — от окружного направления (как в [К]); α2 и углы в осевом ЛД — см. комментарии у формул.
"""
import math
from kompressor.gdf import k, R, cp, s_v, a_cr, gdf_pi, gdf_tau, lambda_from_q

DEFAULTS = dict(
    n_rpm=40000,  # Частота вращения, об/мин
    nu_hub=0.40,  # Относительный диаметр втулки на входе D1вт/D1к
    k_D1=1.00,  # D1к / D1к_опт (оптимум — минимум w1 на периферии)
    k_G=0.97,  # Коэффициент загромождения входа ([К]: 0,93…0,98)
    beta2l=55.0,  # Лопаточный угол на выходе, град от окружного направления (загнуты назад на 90° − β2л)
    z_main=13,  # Число основных лопаток
    z_split=13,  # Число промежуточных (укороченных) лопаток
    split_weight=0.5,  # Вклад промежуточной лопатки в эффективное число лопаток
    cr_ratio=0.85,  # c2r / c1a,ср — отношение меридиональных скоростей (задаёт b2)
    i_tip=0.0, i_hub=2.5,  # Углы атаки на входе, град ([К], с. 166: 0 на периферии, 2…3° у втулки)
    Lz_D2=0.30,  # Осевая длина колеса L_z/D2 (для длины канала в модели трения Jansen)
    eps_tip=0.30e-3,  # Радиальный зазор, м
    Cf=0.005,  # Коэффициент трения в каналах РК [Oh]
    f_inc=0.6,  # Коэффициент потерь на удар (Conrad: 0,5…0,7)
    eps_wake=0.15,  # Относительная ширина следа (Johnston & Dean)
    # БЛД
    D3_D2=1.12,  # Наружный диаметр радиального БЛД
    b3_b2=0.90,  # Сужение канала БЛД (b3/b2 ≤ 1 стабилизирует поток)
    cf_vd=0.006,  # Коэффициент трения стенок в БЛД и колене
    # Переходное колено (радиальное → осевое)
    Rb_b3=2.0,  # Радиус средней линии колена / ширина канала
    AR_bend=1.10,  # Отношение площадей выход/вход колена
    zeta_bend=0.10,  # Коэффициент потерь от кривизны (по меридиональному скоростному напору)
    # Осевой ЛД
    n_rows=4,  # Число спрямляющих венцов (закрутку 55…65° один венец без срыва не снимает)
    sigma_ld=1.8,  # Густота решётки b/t
    alpha_out=0.0,  # Угол выхода из ЛД от оси, град
    AR_ld=1.0,  # Расширение кольцевого канала ЛД (площадь выход/вход)
    n_sections=5,
)


def _wiesner(beta2l_deg, z_eff, D1t_D2):
    """Коэффициент скольжения Wiesner с поправкой на D1/D2"""
    sb = math.sin(math.radians(beta2l_deg))  # cos угла, отсчитанного от радиуса
    sigma = 1 - math.sqrt(sb) / z_eff ** 0.7
    eps_lim = math.exp(-8.16 * sb / z_eff)
    if D1t_D2 > eps_lim:
        sigma *= 1 - ((D1t_D2 - eps_lim) / (1 - eps_lim)) ** 3
    return sigma


def _static(T0, p0, c):
    T = T0 - c ** 2 / (2 * cp)
    return T, p0 * (T / T0) ** (k / (k - 1))


def _inlet(G, p01, T01, D1t, D1h, kG):
    """Скорость на входе без закрутки по уравнению расхода (газодинамическая функция q(λ))"""
    A1 = math.pi / 4 * (D1t ** 2 - D1h ** 2)
    lam = lambda_from_q(G * math.sqrt(T01) / (s_v * p01 * A1 * kG))
    return lam * a_cr(T01)


def _optimal_D1t(G, p01, T01, n, nu, kG):
    """D1к, минимизирующий относительную скорость на периферии (метод золотого сечения)"""
    def w1t(D):
        try:
            c = _inlet(G, p01, T01, D, nu * D, kG)
        except ValueError:
            return 1e9
        return math.hypot(c, math.pi * D * n / 60)
    a, b = 0.03, 0.6
    g = (math.sqrt(5) - 1) / 2
    for _ in range(80):
        c1, c2 = b - g * (b - a), a + g * (b - a)
        if w1t(c1) < w1t(c2):
            b = c2
        else:
            a = c1
    return (a + b) / 2


def _friction_march(state, L, geom, cf, G, steps=200):
    """
    Течение с закруткой в безлопаточном канале с трением о стенки (БЛД, колено).
    geom(s) → (r, b): радиус и ширина канала по длине средней линии s.
    Уравнения: момент количества движения d(r·cu)/ds = −cf·c·cu·r/(b·cm);
    диссипация dq = cf·c³/(b·cm)·ds; ds_энтропии = dq/T; p0 = p0·exp(−Δs/R); T0 = const.
    """
    T0, p0, rcu = state["T0"], state["p0"], state["r"] * state["cu"]
    r, b = geom(0.0)
    cm = state["cm"]
    ds_len = L / steps
    for j in range(steps):
        s = (j + 1) * ds_len
        r, b = geom(s)
        cu = rcu / r
        for _ in range(30):  # плотность и cm из уравнения расхода
            c = math.hypot(cu, cm)
            T, p = _static(T0, p0, c)
            cm_new = G / (p / (R * T) * 2 * math.pi * r * b)
            if abs(cm_new - cm) < 1e-6:
                break
            cm = cm_new
        c = math.hypot(cu, cm)
        dq = cf * c ** 3 / (b * cm) * ds_len
        rcu -= cf * c * cu * r / (b * cm) * ds_len
        p0 *= math.exp(-dq / (T * R))
    cu = rcu / r
    c = math.hypot(cu, cm)
    T, p = _static(T0, p0, c)
    return dict(T0=T0, p0=p0, r=r, b=b, cu=cu, cm=cm, c=c, T=T, p=p, rho=p / (R * T))


def _stage(pik_target, G, p01, T01, u2, p):
    res = {}
    n = p["n_rpm"]
    # ---------------- Вход (станция 1) ----------------
    D1t = p["k_D1"] * _optimal_D1t(G, p01, T01, n, p["nu_hub"], p["k_G"])
    D1h = p["nu_hub"] * D1t
    c1a = _inlet(G, p01, T01, D1t, D1h, p["k_G"])
    T1, p1 = _static(T01, p01, c1a)
    rho1 = p1 / (R * T1)
    D1m = math.sqrt((D1t ** 2 + D1h ** 2) / 2)  # среднеквадратичный диаметр
    sections = []
    for j in range(p["n_sections"]):
        h = j / (p["n_sections"] - 1)
        D = D1h + (D1t - D1h) * h
        u = math.pi * D * n / 60
        w = math.hypot(c1a, u)
        beta = math.degrees(math.atan2(c1a, u))
        i = p["i_hub"] + (p["i_tip"] - p["i_hub"]) * h
        sections.append(dict(h_rel=h, D=D, u=u, ca=c1a, cu=0.0, alpha=90.0, w=w, beta=beta, i=i, beta_l=beta + i,
                             M_w=w / math.sqrt(k * R * T1)))
    u1t, u1h, u1m = (math.pi * d * n / 60 for d in (D1t, D1h, D1m))
    w1t, w1h, w1m = (math.hypot(c1a, u) for u in (u1t, u1h, u1m))
    beta1m = math.degrees(math.atan2(c1a, u1m))
    res.update(D1t=D1t, D1h=D1h, D1m=D1m, c1a=c1a, T1=T1, p1=p1, rho1=rho1, u1t=u1t, w1t=w1t, w1h=w1h, w1m=w1m,
               M_w1t=w1t / math.sqrt(k * R * T1), M_c1=c1a / math.sqrt(k * R * T1), sections=sections,
               beta1m=beta1m)

    # ---------------- Рабочее колесо (станция 2) ----------------
    D2 = 60 * u2 / (math.pi * n)
    z_eff = p["z_main"] + p["split_weight"] * p["z_split"]
    z_all = p["z_main"] + p["z_split"]
    sigma = _wiesner(p["beta2l"], z_eff, D1t / D2)
    c2r = p["cr_ratio"] * c1a
    b2l = math.radians(p["beta2l"])
    c2u = sigma * u2 - c2r / math.tan(b2l)
    dh_euler = u2 * c2u  # c1u = 0
    c2 = math.hypot(c2r, c2u)
    w2 = math.hypot(c2r, u2 - c2u)
    alpha2_rad = math.atan2(c2u, c2r)  # от радиального направления (для моделей потерь)
    # Нагрузка на лопатку (Coppage)
    Df = 1 - w2 / w1t + 0.75 * (dh_euler / u2 ** 2) / (
        (w1t / w2) * ((z_all / math.pi) * (1 - D1t / D2) + 2 * D1t / D2))
    dh_bl = 0.05 * Df ** 2 * u2 ** 2
    # Удар (Conrad): составляющая относительной скорости, теряемая при рассогласовании угла на среднем диаметре
    i_m = (p["i_hub"] + p["i_tip"]) / 2
    dh_inc = p["f_inc"] * (w1m * math.sin(math.radians(i_m))) ** 2 / 2
    # Трение (Jansen): углы от меридионального направления
    b1m = math.radians(90 - beta1m)
    b2m = math.radians(90 - p["beta2l"])
    # b2 — из уравнения расхода; плотность ρ2 уточняется итерациями
    rho2 = rho1 * 2.0
    for _ in range(30):
        b2 = G / (rho2 * math.pi * D2 * c2r)
        D_hyd = D2 * (math.cos(b2m) / (z_all / math.pi + D2 * math.cos(b2m) / b2)
                      + 0.5 * (D1t / D2 + D1h / D2) * math.cos(b1m) /
                      (z_all / math.pi + (D1t + D1h) / (D1t - D1h) * math.cos(b1m)))
        L_b = math.pi / 8 * (D2 - (D1t + D1h) / 2 - b2 + 2 * p["Lz_D2"] * D2) * 2 / (math.cos(b1m) + math.cos(b2m))
        W_mean = (c1a + c2 + w1t + 2 * w1h + 3 * w2) / 8
        dh_sf = 2 * p["Cf"] * L_b / D_hyd * W_mean ** 2
        # Зазор (Jansen)
        r1t, r1h, r2 = D1t / 2, D1h / 2, D2 / 2
        dh_cl = 0.6 * p["eps_tip"] / b2 * c2u * math.sqrt(
            4 * math.pi / (b2 * z_all) * (r1t ** 2 - r1h ** 2) / ((r2 - r1t) * (1 + rho2 / rho1)) * c2u * c1a)
        # Смешение за колесом (Johnston & Dean), b* = b3/b2
        b_star = p["b3_b2"]
        dh_mix = 1 / (1 + math.tan(alpha2_rad) ** 2) * ((1 - p["eps_wake"] - b_star) / (1 - p["eps_wake"])) ** 2 \
            * c2 ** 2 / 2
        dh_int = dh_bl + dh_inc + dh_sf + dh_cl + dh_mix
        # Паразитные потери: трение диска (Daily & Nece), рециркуляция (Oh), утечки (Aungier)
        T02 = T01 + dh_euler / cp  # первое приближение для вязкости и плотности
        T2 = T02 - c2 ** 2 / (2 * cp)
        mu2 = 1.716e-5 * (T2 / 273.15) ** 1.5 * (273.15 + 110.4) / (T2 + 110.4)
        Re_df = u2 * r2 * rho2 / mu2
        f_df = 2.67 / Re_df ** 0.5 if Re_df < 3e5 else 0.0622 / Re_df ** 0.2
        dh_df = f_df * (rho1 + rho2) / 2 * r2 ** 2 * u2 ** 3 / (4 * G)
        dh_rc = 8e-5 * math.sinh(3.5 * alpha2_rad ** 3) * Df ** 2 * u2 ** 2
        dP_cl = G * r2 * c2u / (z_all * (r1t + r2) / 2 * (b2 + (r1t - r1h)) / 2 * L_b)
        U_cl = 0.816 * math.sqrt(2 * dP_cl / rho2)
        m_cl = rho2 * z_all * p["eps_tip"] * L_b * U_cl
        dh_lk = m_cl * U_cl * u2 / (2 * G)
        dh_par = dh_df + dh_rc + dh_lk
        T02 = T01 + (dh_euler + dh_par) / cp
        p02 = p01 * (1 + (dh_euler - dh_int) / (cp * T01)) ** (k / (k - 1))
        T2, p2 = _static(T02, p02, c2)
        rho2_new = p2 / (R * T2)
        if abs(rho2_new - rho2) < 1e-6:
            break
        rho2 = rho2_new
    eta_rk = cp * T01 * ((p02 / p01) ** ((k - 1) / k) - 1) / (cp * (T02 - T01))
    res.update(D2=D2, u2=u2, z_eff=z_eff, z_all=z_all, slip=sigma, c2r=c2r, c2u=c2u, c2=c2, w2=w2,
               alpha2=math.degrees(math.atan2(c2r, c2u)), alpha2_rad_deg=math.degrees(alpha2_rad),
               beta2=math.degrees(math.atan2(c2r, u2 - c2u)), dh_euler=dh_euler, Df=Df, b2=b2, D_hyd=D_hyd, L_b=L_b,
               losses_int=dict(удар=dh_inc, нагрузка=dh_bl, трение=dh_sf, зазор=dh_cl, смешение=dh_mix),
               losses_par=dict(трение_диска=dh_df, рециркуляция=dh_rc, утечки=dh_lk),
               T02=T02, p02=p02, T2=T2, p2=p2, rho2=rho2, eta_rk=eta_rk, W2_W1t=w2 / w1t,
               M_c2=c2 / math.sqrt(k * R * T2), lam_c2=c2 / a_cr(T02), n_rpm=n)

    # ---------------- Радиальный БЛД (станция 3) ----------------
    r3 = p["D3_D2"] * r2
    b3 = p["b3_b2"] * b2
    st2 = dict(T0=T02, p0=p02, r=r2, cu=c2u, cm=c2r)
    st3 = _friction_march(st2, r3 - r2, lambda s: (r2 + s, b2 + (b3 - b2) * s / (r3 - r2)), p["cf_vd"], G)
    res.update(D3=2 * r3, b3=b3, st3=st3, alpha3=math.degrees(math.atan2(st3["cm"], st3["cu"])),
               M3=st3["c"] / math.sqrt(k * R * st3["T"]))

    # ---------------- Переходное колено: радиальное → осевое (станция 4) ----------------
    Rb = p["Rb_b3"] * b3
    r4 = r3 + Rb
    h4 = p["AR_bend"] * r3 * b3 / r4  # высота кольцевого канала на выходе из колена
    L_bend = math.pi / 2 * Rb

    def bend_geom(s):
        phi = s / Rb
        r = r3 + Rb * math.sin(phi)
        A = 2 * math.pi * r3 * b3 + (2 * math.pi * r4 * h4 - 2 * math.pi * r3 * b3) * s / L_bend
        return r, A / (2 * math.pi * r)

    st4 = _friction_march(st3, L_bend, bend_geom, p["cf_vd"], G)
    # Потери от кривизны: доля меридионального скоростного напора
    dq_b = p["zeta_bend"] * st4["cm"] ** 2 / 2
    st4["p0"] *= math.exp(-dq_b / (st4["T"] * R))
    st4["T"], st4["p"] = _static(st4["T0"], st4["p0"], st4["c"])
    alpha4 = math.degrees(math.atan2(st4["cu"], st4["cm"]))  # угол закрутки от оси
    res.update(D4m=2 * r4, h4=h4, Rb=Rb, st4=st4, alpha4=alpha4, M4=st4["c"] / math.sqrt(k * R * st4["T"]))

    # ---------------- Осевой ЛД: спрямляющие венцы (станция 5) ----------------
    # Углы от оси. Поворот распределён так, чтобы отношение cos α_вых/cos α_вх (торможение) во всех венцах
    # было одинаковым — тогда венцы нагружены равномерно по D_eq; c_a меняется за счёт расширения канала.
    rows = []
    T0, p0 = st4["T0"], st4["p0"]
    ca, cu = st4["cm"], st4["cu"]
    a_in = alpha4
    n_rows = p["n_rows"]
    A_in = 2 * math.pi * r4 * h4
    for j in range(n_rows):
        q = (math.cos(math.radians(p["alpha_out"])) / math.cos(math.radians(alpha4))) ** (1 / n_rows)
        a_out = math.degrees(math.acos(min(1.0, math.cos(math.radians(alpha4)) * q ** (j + 1))))
        A_out = A_in * p["AR_ld"] ** ((j + 1) / n_rows)
        c_in = math.hypot(ca, cu)
        T_in, p_in = _static(T0, p0, c_in)
        ca_out = ca
        for _ in range(30):  # осевая скорость на выходе венца из уравнения расхода
            c_out = ca_out / math.cos(math.radians(a_out))
            T_o, p_o = _static(T0, p0, c_out)
            ca_new = G / (p_o / (R * T_o) * A_out)
            if abs(ca_new - ca_out) < 1e-6:
                break
            ca_out = ca_new
        a1r, a2r = math.radians(a_in), math.radians(a_out)
        s_ = p["sigma_ld"]
        # [L]: эквивалентный диффузорный фактор (при c_a1 ≠ c_a2 — через отношение скоростей)
        D_eq = (c_in / c_out) * (1.12 + 0.61 * math.cos(a1r) ** 2 / s_ * (math.tan(a1r) - (ca_out / ca) * math.tan(a2r)))
        theta_c = 0.004 * (1 + 3.1 * (D_eq - 1) ** 2 + 0.4 * (D_eq - 1) ** 8)
        omega = 2 * theta_c * s_ / math.cos(a2r) * (math.cos(a1r) / math.cos(a2r)) ** 2
        D_f = 1 - c_out / c_in + abs(ca * math.tan(a1r) - ca_out * math.tan(a2r)) / (2 * s_ * c_in)
        if D_eq > 3.0:
            raise ValueError(f"Срыв в венце ЛД {j + 1}: D_eq = {D_eq:.2f}")
        p0_out = p0 - omega * (p0 - p_in)
        # Лопаточные углы (от оси): вход — без угла атаки, выход — с отставанием по правилу Картера
        theta_est = abs(a_in - a_out)
        for _ in range(20):
            dev = 0.23 * theta_est / math.sqrt(s_)
            theta_new = abs(a_in - (a_out - dev))
            if abs(theta_new - theta_est) < 1e-6:
                break
            theta_est = theta_new
        rows.append(dict(a_in=a_in, a_out=a_out, c_in=c_in, c_out=c_out, D_eq=D_eq, D_f=D_f, omega=omega,
                         M_in=c_in / math.sqrt(k * R * T_in), dev=dev, blade_in=a_in, blade_out=a_out - dev,
                         camber=theta_est, sigma=s_))
        p0 = p0_out
        ca, cu = ca_out, ca_out * math.tan(a2r)
        a_in = a_out
    c5 = math.hypot(ca, cu)
    T5, p5 = _static(T0, p0, c5)
    res.update(rows=rows, T05=T0, p05=p0, c5=c5, h5=A_in * p["AR_ld"] / (2 * math.pi * r4))

    pi_calc = p0 / p01
    H = cp * T01 * (pi_calc ** ((k - 1) / k) - 1)
    res.update(pi_calc=pi_calc, H=H, H_z=cp * (T0 - T01), eta_k=H / (cp * (T0 - T01)), N=G * cp * (T0 - T01),
               G=G, pik=pik_target, p01=p01, T01=T01, S=p["Lz_D2"] * D2)
    return res


def compressor(pik, G, p_n, T_n, **params):
    """Расчёт ЦБК на заданные π*_к и расход; u2 (а при заданной n — и D2) подбирается из условия π*_расч = π*_к"""
    p = dict(DEFAULTS)
    unknown = set(params) - set(p)
    if unknown:
        raise KeyError(f"Неизвестные параметры компрессора: {unknown}")
    p.update(params)

    def f(u2):
        return _stage(pik, G, p_n, T_n, u2, p)["pi_calc"] - pik

    x0, x1 = 520.0, 560.0
    f0, f1 = f(x0), f(x1)
    for _ in range(60):
        if f1 == f0:
            break
        x2 = min(max(x1 - f1 * (x1 - x0) / (f1 - f0), 300.0), 800.0)
        x0, f0, x1, f1 = x1, f1, x2, f(x2)
        if abs(f1) / pik < 1e-6:
            break
    res = _stage(pik, G, p_n, T_n, x1, p)
    res["params"] = p
    res["warnings"] = _checks(res, p)
    return res


def _checks(r, p):
    w = []
    if abs(r["pi_calc"] - r["pik"]) / r["pik"] > 1e-3:
        w.append(f"π*_расч = {r['pi_calc']:.3f} не сошлась с заданной {r['pik']}")
    if r["M_w1t"] > 1.3:
        w.append(f"M_w1 на периферии = {r['M_w1t']:.2f} > 1,3 — трансзвуковой вход с большими потерями")
    if r["u2"] > 620:
        w.append(f"u2 = {r['u2']:.0f} м/с > 620 м/с — предел для титанового колеса")
    if r["W2_W1t"] < 0.6:
        w.append(f"w2/w1к = {r['W2_W1t']:.3f} < 0,6 — чрезмерное торможение в РК (критерий де Галлера)")
    if r["b2"] < 5e-3:
        w.append(f"b2 = {r['b2'] * 1000:.2f} мм < 5 мм")
    if r["D1t"] / r["D2"] > 0.7:
        w.append(f"D1к/D2 = {r['D1t'] / r['D2']:.3f} > 0,7")
    if not 60 <= r["alpha2_rad_deg"] <= 75:
        w.append(f"Угол потока за РК от радиуса α2 = {r['alpha2_rad_deg']:.1f}° вне 60…75° "
                 f"(больше — неустойчивость БЛД, меньше — большая c2r)")
    if r["M4"] > 0.85:
        w.append(f"M на входе в осевой ЛД = {r['M4']:.2f} > 0,85")
    for j, row in enumerate(r["rows"], 1):
        if row["D_eq"] > 2.0:
            w.append(f"Венец ЛД {j}: D_eq = {row['D_eq']:.2f} > 2,0 — срыв (Lieblein)")
        if row["D_f"] > 0.6:
            w.append(f"Венец ЛД {j}: диффузорный фактор D = {row['D_f']:.2f} > 0,6")
    if not 90 <= r["c5"] <= 150:
        w.append(f"Скорость на выходе c_к = {r['c5']:.0f} м/с вне 100…140 м/с ([К], с. 163)")
    return w


def _material(T):
    if T < 520:
        return "алюминиевый сплав"
    if T <= 770:
        return "титановый сплав"
    return "жаропрочная сталь или сплав"


def print_report(r):
    p = r["params"]
    print("=" * 84)
    print("   ГАЗОДИНАМИЧЕСКИЙ РАСЧЁТ ЦБК: РК + радиальный БЛД + колено + осевой ЛД")
    print("=" * 84)
    print(f"π*_к = {r['pik']:.3f}, G = {r['G']:.3f} кг/с, p* = {r['p01'] / 1000:.2f} кПа, T* = {r['T01']:.2f} К, "
          f"n = {r['n_rpm']:.0f} об/мин")

    print("\n─── Вход в РК (без входного направляющего аппарата, c1u = 0) ───")
    print(f"  D1к = {r['D1t'] * 1000:.1f} мм, D1вт = {r['D1h'] * 1000:.1f} мм (ν = {p['nu_hub']}), "
          f"c1a = {r['c1a']:.1f} м/с, M_c1 = {r['M_c1']:.3f}")
    print(f"  w1к = {r['w1t']:.1f} м/с, M_w1к = {r['M_w1t']:.3f}")
    print(f"  {'h/h':<5} | {'D, мм':<7} | {'u1, м/с':<8} | {'w1, м/с':<8} | {'M_w1':<6} | {'β1, °':<7} | {'i, °':<5} | "
          f"{'β1л, °':<7}")
    print("  " + "-" * 70)
    for s in r["sections"]:
        print(f"  {s['h_rel']:<5.2f} | {s['D'] * 1000:<7.1f} | {s['u']:<8.1f} | {s['w']:<8.1f} | {s['M_w']:<6.3f} | "
              f"{s['beta']:<7.2f} | {s['i']:<5.2f} | {s['beta_l']:<7.2f}")

    print("\n─── Рабочее колесо ───")
    print(f"  u2 = {r['u2']:.1f} м/с, D2 = {r['D2'] * 1000:.1f} мм, D1к/D2 = {r['D1t'] / r['D2']:.3f}")
    print(f"  Лопаток {p['z_main']} + {p['z_split']} промежуточных, β2л = {p['beta2l']:.1f}° "
          f"(загиб назад {90 - p['beta2l']:.1f}°), скольжение σ = {r['slip']:.4f}")
    print(f"  c2r = {r['c2r']:.1f}, c2u = {r['c2u']:.1f}, c2 = {r['c2']:.1f} м/с, M_c2 = {r['M_c2']:.3f}, "
          f"w2 = {r['w2']:.1f} м/с, w2/w1к = {r['W2_W1t']:.3f}")
    print(f"  α2 = {r['alpha2']:.2f}° от окружного ({r['alpha2_rad_deg']:.1f}° от радиуса), β2 = {r['beta2']:.2f}°")
    print(f"  b2 = {r['b2'] * 1000:.2f} мм (b2/D2 = {r['b2'] / r['D2']:.4f}), Df = {r['Df']:.3f}")
    print(f"  Работа Эйлера {r['dh_euler'] / 1000:.2f} кДж/кг")
    print("  Внутренние потери, кДж/кг: " + ", ".join(f"{kk} {v / 1000:.2f}" for kk, v in r["losses_int"].items()))
    print("  Паразитные потери, кДж/кг: " + ", ".join(f"{kk} {v / 1000:.2f}" for kk, v in r["losses_par"].items()))
    print(f"  p2* = {r['p02'] / 1000:.1f} кПа, T2* = {r['T02']:.1f} К, η_РК = {r['eta_rk']:.4f}")

    s3, s4 = r["st3"], r["st4"]
    print("\n─── Радиальный БЛД ───")
    print(f"  D3 = {r['D3'] * 1000:.1f} мм, b3 = {r['b3'] * 1000:.2f} мм, c3 = {s3['c']:.1f} м/с, M3 = {r['M3']:.3f}, "
          f"α3 = {r['alpha3']:.1f}° от окружного, σ_БЛД = {s3['p0'] / r['p02']:.4f}")
    print("\n─── Переходное колено (радиальное → осевое) ───")
    print(f"  R_кол = {r['Rb'] * 1000:.1f} мм, D_ср = {r['D4m'] * 1000:.1f} мм, h = {r['h4'] * 1000:.2f} мм, "
          f"c4 = {s4['c']:.1f} м/с (c_a = {s4['cm']:.1f}, c_u = {s4['cu']:.1f}), M4 = {r['M4']:.3f}")
    print(f"  Угол закрутки на входе в ЛД α = {r['alpha4']:.1f}° от оси, σ_кол = {s4['p0'] / s3['p0']:.4f}")
    print("\n─── Осевой лопаточный диффузор (углы от оси) ───")
    print(f"  {'Венец':<5} | {'α_вх':<6} | {'α_вых':<6} | {'c_вх':<6} | {'c_вых':<6} | {'M_вх':<5} | {'D_eq':<5} | "
          f"{'D':<5} | {'ω':<6} | {'δ':<5} | {'лоп. вх':<7} | {'лоп. вых':<8} | {'θ':<5}")
    for j, row in enumerate(r["rows"], 1):
        print(f"  {j:<5} | {row['a_in']:<6.1f} | {row['a_out']:<6.1f} | {row['c_in']:<6.1f} | {row['c_out']:<6.1f} | "
              f"{row['M_in']:<5.2f} | {row['D_eq']:<5.2f} | {row['D_f']:<5.2f} | {row['omega']:<6.4f} | "
              f"{row['dev']:<5.2f} | {row['blade_in']:<7.1f} | {row['blade_out']:<8.1f} | {row['camber']:<5.1f}")
    print(f"  Выход: c_к = {r['c5']:.1f} м/с, высота канала {r['h5'] * 1000:.1f} мм, "
          f"σ_ЛД = {r['p05'] / s4['p0']:.4f}")

    print("\n─── Компрессор в целом ───")
    print(f"  π*_к = {r['pi_calc']:.4f}, η_к* = {r['eta_k']:.4f}, T_к* = {r['T05']:.1f} К, "
          f"N = {r['N'] / 1000:.1f} кВт")
    print(f"  Материал колеса по T_к* ([К], с. 152): {_material(r['T05'])}; при u2 = {r['u2']:.0f} м/с — титановый сплав")
    if r["eta_k"] > 0.83:
        print(f"  Примечание: η_к = {r['eta_k']:.3f} выше типичных для малоразмерных ЦБК 0,78…0,83 ([К], с. 164) — "
              f"модель потерь по средней линии не учитывает, например, неравномерность потока на входе в диффузор")
    if r["warnings"]:
        print("\n⚠️ Замечания:")
        for msg in r["warnings"]:
            print("   - " + msg)
    else:
        print("\n✅ Все проверки пройдены")
