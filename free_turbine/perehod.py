"""
Переходный канал ТВД → СТ (свободная турбина) и входная часть СА СТ.

Компоновка — по чертежу реального двигателя books/гг.png: за РК ТВД канал идёт назад с небольшим подъёмом
(на чертеже ≈ 5…6°) и почти постоянной высотой (не диффузор и не конфузор), задняя стенка КС скошена и
опирается на наружную стенку канала; увеличение высоты проточной части до лопаток СТ происходит в самой
ступени СТ. Поэтому здесь:
  - канал — прямой, постоянной радиальной высоты h = высоте лопатки РК ТВД на выходе, от выхода РК ТВД до
    входа в СА СТ;
  - СА СТ на входе имеет ту же высоту h, а к выходу раскрывается до расчётной высоты лопатки СТ l_c (углы
    раскрытия — по Батурину п. 43: на периферии ≤ 15…20°, на втулке ≤ 8…12°; раскрытие делится между
    периферией и втулкой пропорционально предельным tg углов);
  - длина канала — наименьшая, при которой скошенная задняя стенка КС (линия x + r = C, 45°) садится на
    наружную стенку канала не ближе MARGIN к СА СТ.
Угол подъёма средней линии канала сравнивается с чертежом (GAMMA_MAX); если он больше, уменьшается D_ср СТ
(free_turbine/gasdynamiks.py).
"""
import math

GAMMA_MAX = 6.0  # допустимый угол подъёма средней линии канала, град (чертёж books/гг.png: ≈ 5…6°)
FLARE_TIP, FLARE_HUB = 20.0, 12.0  # предельные углы раскрытия, град (Батурин п. 43)
MARGIN = 0.010  # от посадки задней стенки КС на канал до входа в СА СТ, м


def duct(x_rk, r_rk_hub, r_rk_tip, r_ct_hub, r_ct_tip, S_sa_ct, C_wall):
    """
    x_rk — осевая координата выхода РК ТВД; r_rk_hub, r_rk_tip — радиусы втулки и корпуса на выходе РК ТВД;
    r_ct_hub, r_ct_tip — радиусы втулки и корпуса на выходе СА СТ; S_sa_ct — осевая ширина СА СТ;
    C_wall — задняя стенка КС: x + r = C_wall (все величины — в одних единицах длины).
    """
    h = r_rk_tip - r_rk_hub
    dl = (r_ct_tip - r_ct_hub) - h  # раскрытие в СА СТ
    tt, th = math.tan(math.radians(FLARE_TIP)), math.tan(math.radians(FLARE_HUB))
    r_in_tip = r_ct_tip - dl * tt / (tt + th)
    r_in_hub = r_ct_hub + dl * th / (tt + th)
    rise = (r_in_hub + r_in_tip) / 2 - (r_rk_hub + r_rk_tip) / 2
    gamma = 0.0
    for _ in range(100):  # длина канала и угол подъёма согласуются через посадку стенки КС
        tg = math.tan(gamma)
        x_b = (C_wall - r_rk_tip + x_rk * tg) / (1 + tg)  # x + r_tip(x) = C_wall
        L = max(x_b - x_rk, 0.0) + MARGIN
        g_new = math.atan2(rise, L)
        if abs(g_new - gamma) < 1e-12:
            break
        gamma = g_new
    tg = math.tan(gamma)
    return dict(h=h, L=L, gamma=math.degrees(gamma), rise=rise, x_rk=x_rk, x_ct=x_rk + L,
                x_b=x_b, r_b=r_rk_tip + (x_b - x_rk) * tg,
                r_rk_hub=r_rk_hub, r_rk_tip=r_rk_tip, r_in_hub=r_in_hub, r_in_tip=r_in_tip,
                r_ct_hub=r_ct_hub, r_ct_tip=r_ct_tip, area_ratio=(r_in_hub + r_in_tip) / (r_rk_hub + r_rk_tip),
                flare_tip=math.degrees(math.atan((r_ct_tip - r_in_tip) / S_sa_ct)),
                flare_hub=math.degrees(math.atan((r_in_hub - r_ct_hub) / S_sa_ct)))


def from_models(ks, tvd, ct):
    """Канал по результатам расчёта КС (kamera_sgoraniya/mingazov.py), модулей ТВД и СТ, в координатах КС (м):
    СА ТВД начинается в плоскости выхода газосборника x_c, зазор СА–РК ТВД = 0,25·S_рк (как на схемах)."""
    g = ks["geo"]
    x_rk = g["x_c"] + tvd.S_sa + 0.25 * tvd.S_rk + tvd.S_rk
    return duct(x_rk, tvd.d_rk / 2, tvd.D_rk / 2, ct.d_sa / 2, ct.D_sa / 2, ct.S_sa, g["C_wall"])


def warnings(d):
    w = []
    if d["gamma"] > GAMMA_MAX:
        w.append(f"угол подъёма канала {d['gamma']:.1f}° > {GAMMA_MAX}° — уменьшите D_ср СТ")
    if d["flare_tip"] > FLARE_TIP or d["flare_hub"] > FLARE_HUB:
        w.append(f"раскрытие СА СТ {d['flare_tip']:.1f}° (периферия) / {d['flare_hub']:.1f}° (втулка) больше "
                 f"{FLARE_TIP:.0f}° / {FLARE_HUB:.0f}° (Батурин п. 43) — удлините СА СТ")
    return w


def print_report(d, mm=1000.0):
    print("\n─── Переходный канал ТВД → СТ (по чертежу books/гг.png) ───")
    print(f"  Постоянная высота h = {d['h'] * mm:.1f} мм (= лопатка РК ТВД на выходе), длина L = {d['L'] * mm:.1f} мм, "
          f"подъём средней линии {d['rise'] * mm:.1f} мм → угол {d['gamma']:.1f}° (≤ {GAMMA_MAX}°)")
    print(f"  Вход: ⌀{2 * d['r_rk_hub'] * mm:.1f}/⌀{2 * d['r_rk_tip'] * mm:.1f} мм → выход (вход в СА СТ): "
          f"⌀{2 * d['r_in_hub'] * mm:.1f}/⌀{2 * d['r_in_tip'] * mm:.1f} мм; отношение площадей "
          f"{d['area_ratio']:.2f} (только за счёт роста радиуса)")
    print(f"  СА СТ раскрывается до ⌀{2 * d['r_ct_hub'] * mm:.1f}/⌀{2 * d['r_ct_tip'] * mm:.1f} мм: "
          f"{d['flare_tip']:.1f}° на периферии (≤ {FLARE_TIP:.0f}°), {d['flare_hub']:.1f}° на втулке (≤ {FLARE_HUB:.0f}°)")
    for msg in warnings(d):
        print("  ⚠️ " + msg)
