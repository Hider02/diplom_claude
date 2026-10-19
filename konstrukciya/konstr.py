"""
Схематичная конструкция роторов, опор и выходного устройства в координатах общей схемы (мм, x — по оси,
y — радиус). Компоновка — по аналогу ВК-800С (books/гг.png):
  - газогенератор: колесо ЦБК и диск ТВД на одном полом валу ⌀36, внутри — центральная стяжка (труба) с гайкой
    спереди; задняя цапфа колеса и передний фланец диска ТВД — болтовое соединение на r = 30 мм;
    опора 1 — шариковый радиально-упорный подшипник в 57 мм перед входом в колесо (во входном устройстве),
    опора 2 — роликовый в ~37 мм за диском ТВД, корпус опоры — под переходным каналом;
  - СТ: диск сплошной, без отверстия; передняя цапфа (короткий конус) и задний полый вал ⌀50 (⌀40 у аналога увеличен по критической частоте) крепятся к диску
    12 болтами М6 на r = 37 мм; опора 3 — роликовый в 70 мм за диском, опора 4 — шариковый в 190 мм; вал СТ
    проходит внутри кожуха выходного устройства и шлицами входит в редуктор сзади;
  - выходное устройство — по расчёту konstrukciya/vyhodnoe_ustr.py.
Профили дисков и ступиц назначены конструктивно (в расчёте проточной части их нет) и используются
также для масс в расчёте критических частот (konstrukciya/raschet.py).
"""
import math

from konstrukciya.geometriya import tract
from konstrukciya.vyhodnoe_ustr import calc as calc_vu

# ─── Назначенные размеры, мм ───
GG = dict(r_sh=18.0, r_bore=11.0, r_tie=7.0, r_tie_in=4.0, x_front=-78.0, x_rear=216.0,
          brg1_dx=-57.0, brg2_dx=37.0,  # опора 1 — от входа в колесо; опора 2 — от задней грани диска ТВД
          brg_d=40.0, brg_D=80.0, brg_B=18.0, bolt_r=30.0, bolt_d=8.0, bolt_n=8)
CT = dict(r_sh=25.0, r_bore=19.0, brg3_dx=70.0, brg4_dx=190.0,  # от средней плоскости диска СТ
          brg3=(55.0, 90.0, 18.0), brg4=(55.0, 100.0, 21.0), bolt_r=37.0, bolt_d=6.0, bolt_n=12, x_end=615.0)
INLET = dict(x0=-115.0, r_hub=52.0, r_tip=92.0, x_bend=-35.0)


def _poly_area_moments(pts):
    """Площадь и статический момент относительно оси (для объёма тела вращения по Паппу–Гульдину)."""
    A = Sy = Sx = 0.0
    n = len(pts)
    for i in range(n):
        x1, y1 = pts[i]
        x2, y2 = pts[(i + 1) % n]
        cr = x1 * y2 - x2 * y1
        A += cr
        Sy += (y1 + y2) * cr
        Sx += (x1 + x2) * cr
    A /= 2
    return abs(A), abs(Sy / 6), (Sx / 6) / A if A else 0.0  # A, ∫y dA, x_центра


def revolve_mass(pts, rho):
    """Масса тела вращения профиля pts (мм) при плотности rho (кг/м³), и осевая координата центра масс."""
    A, Sy, xc = _poly_area_moments(pts)
    return 2 * math.pi * Sy * 1e-9 * rho, xc


def build():
    T = tract()
    vu = calc_vu()
    c, tv, ct, du = T["cbk"], T["tvd"], T["ct"], T["duct"]
    P = dict(main=[], thin=[], axis=[], bearings=[], bolts=[], labels=[], dims=[], masses={})

    def poly(pts, closed=True, key="main"):
        pts = list(pts)
        P[key].append(pts + [pts[0]] if closed else pts)

    # ─── Входное устройство (схематично) и опора 1 ───
    xb1 = c["x0"] + GG["brg1_dx"]
    poly([(INLET["x0"], INLET["r_hub"]), (INLET["x_bend"], INLET["r_hub"]), (-12.0, c["r_vt"] + 6), (-2.0, c["r_vt"] + 1)],
         closed=False)
    poly([(INLET["x0"], INLET["r_tip"]), (INLET["x_bend"], INLET["r_tip"]), (c["x0"], c["r1"])], closed=False)
    # стойка входного устройства (корпус опоры 1 — к наружному корпусу)
    poly([(-92.0, INLET["r_hub"]), (-88.0, INLET["r_tip"]), (-72.0, INLET["r_tip"]), (-76.0, INLET["r_hub"])])
    # корпус опоры 1
    rD1 = GG["brg_D"] / 2
    poly([(xb1 - 16, rD1), (xb1 + 16, rD1), (xb1 + 16, rD1 + 4), (INLET["x_bend"] + 4, INLET["r_hub"]),
          (xb1 - 16, INLET["r_hub"])])
    P["bearings"].append(dict(x=xb1, d=GG["brg_d"], D=GG["brg_D"], B=GG["brg_B"], kind="ball", name="Опора 1"))

    # ─── Колесо ЦБК: тыльная стенка и задняя цапфа с фланцем ───
    S, r2 = c["S"], c["r2"]
    x_fl = 112.0  # плоскость стыка фланцев колесо ЦБК — диск ТВД
    poly([(S, r2), (S, 46.0), (x_fl - 6, 46.0), (x_fl - 6, 48.0), (x_fl, 48.0), (x_fl, GG["r_sh"])], closed=False)
    wheel = [(c["x0"] + S * math.sin(t / 30 * math.pi / 2),
              c["r_vt"] + (r2 - c["r_vt"]) * (1 - math.cos(t / 30 * math.pi / 2))) for t in range(31)]
    wheel += [(S, 46.0), (x_fl, 46.0), (x_fl, GG["r_sh"]), (c["x0"], GG["r_sh"])]
    m_w, xc_w = revolve_mass(wheel, 4500.0)  # титан ВТ8
    P["masses"]["ЦБК"] = (m_w * 1.10, xc_w)  # +10 % на лопатки

    # ─── Диск ТВД ───
    z = tv["z"]
    x_rk0, x_rk1 = z[2], z[3]
    xm = (x_rk0 + x_rk1) / 2
    r_rim = tv["hub"][3]
    disk_t = [(x_rk0, r_rim), (x_rk1, r_rim), (x_rk1, 82.0), (xm + 4, 78.0), (xm + 4, 48.0), (xm + 13, 44.0),
              (xm + 13, GG["r_sh"]), (xm - 13, GG["r_sh"]), (xm - 13, 44.0), (xm - 4, 48.0), (xm - 4, 78.0),
              (x_rk0, 82.0)]
    poly(disk_t)
    # передний фланец диска ТВД к цапфе колеса
    poly([(x_fl, GG["r_sh"]), (x_fl, 48.0), (x_fl + 6, 48.0), (x_fl + 6, 40.0), (xm - 13, 40.0)], closed=False)
    m_t, xc_t = revolve_mass(disk_t, 8100.0)
    m_bl = 8100.0 * math.pi * (tv["tip"][3] ** 2 - r_rim ** 2) * 1e-9 * tv["S_rk"] * 0.25
    P["masses"]["ТВД"] = (m_t + m_bl, xm)
    P["bolts"].append(dict(x0=x_fl - 6 - 4, x1=x_fl + 6 + 4, r=GG["bolt_r"], d=GG["bolt_d"],
                           note=f"{GG['bolt_n']} болтов М{GG['bolt_d']:.0f}"))

    # ─── Вал газогенератора, стяжка, опора 2 ───
    xb2 = x_rk1 + GG["brg2_dx"]  # 37 мм за задней кромкой обода диска ТВД
    rj = GG["brg_d"] / 2
    sh = [(GG["x_front"], GG["r_sh"]), (xb1 - 12, GG["r_sh"]), (xb1 - 12, rj), (xb1 + 12, rj), (xb1 + 12, GG["r_sh"]),
          (xb2 - 12, GG["r_sh"]), (xb2 - 12, rj), (xb2 + 12, rj), (xb2 + 12, GG["r_sh"]), (GG["x_rear"], GG["r_sh"]),
          (GG["x_rear"], GG["r_bore"]), (GG["x_front"], GG["r_bore"])]
    poly(sh)
    # стяжка (труба) с гайкой спереди и головкой сзади
    poly([(GG["x_front"] - 10, GG["r_tie"]), (GG["x_rear"] + 8, GG["r_tie"])], closed=False)
    poly([(GG["x_front"] - 10, GG["r_tie_in"]), (GG["x_rear"] + 8, GG["r_tie_in"])], closed=False)
    poly([(GG["x_front"] - 10, GG["r_tie_in"]), (GG["x_front"], GG["r_tie_in"]), (GG["x_front"], 16.0),
          (GG["x_front"] - 10, 16.0)])
    poly([(GG["x_rear"], GG["r_tie_in"]), (GG["x_rear"] + 8, GG["r_tie_in"]), (GG["x_rear"] + 8, 15.0),
          (GG["x_rear"], 15.0)])
    P["bearings"].append(dict(x=xb2, d=GG["brg_d"], D=GG["brg_D"], B=GG["brg_B"], kind="roller", name="Опора 2"))
    # корпус опоры 2: от втулки переходного канала вниз, с конической диафрагмой
    hub_duct = lambda x: du["hub"][0] + (du["hub"][1] - du["hub"][0]) * (x - du["x0"]) / (du["x1"] - du["x0"])
    poly([(xb2 - 16, rD1), (xb2 + 16, rD1), (xb2 + 16, rD1 + 4), (xb2 + 30, hub_duct(xb2 + 30)),
          (xb2 + 22, hub_duct(xb2 + 22)), (xb2 - 16, rD1 + 4)])
    # стойка в переходном канале (подвеска корпуса опоры 2)
    tip_duct = lambda x: du["tip"][0] + (du["tip"][1] - du["tip"][0]) * (x - du["x0"]) / (du["x1"] - du["x0"])
    poly([(xb2 + 12, hub_duct(xb2 + 12)), (xb2 + 16, tip_duct(xb2 + 16)), (xb2 + 34, tip_duct(xb2 + 34)),
          (xb2 + 30, hub_duct(xb2 + 30))])
    P["masses"]["опоры ГГ"] = (xb1, xb2)

    # ─── Диск СТ (сплошной), цапфа, вал СТ ───
    z = ct["z"]
    x0, x1 = z[2], z[3]
    xm = (x0 + x1) / 2
    r_rim = ct["hub"][3]
    disk_c = [(x0, r_rim), (x1, r_rim), (x1, 80.0), (xm + 5, 76.0), (xm + 5, 50.0), (xm + 15, 46.0), (xm + 15, 0.0),
              (xm - 15, 0.0), (xm - 15, 46.0), (xm - 5, 50.0), (xm - 5, 76.0), (x0, 80.0)]
    poly(disk_c)
    m_c, _ = revolve_mass(disk_c, 8100.0)
    m_bl = 8100.0 * math.pi * (ct["tip"][3] ** 2 - r_rim ** 2) * 1e-9 * ct["S_rk"] * 0.25
    P["masses"]["СТ"] = (m_c + m_bl, xm)
    xf, xr = xm - 15, xm + 15
    poly([(xf, 30.0), (xf - 6, 30.0), (xf - 6, 44.0), (xf, 44.0)])  # фланец передней цапфы
    poly([(xf - 6, 30.0), (xf - 22, 22.0), (xf - 22, 16.0), (xf - 6, 16.0)], closed=False)  # короткий конус
    xb3, xb4 = xm + CT["brg3_dx"], xm + CT["brg4_dx"]
    j3, j4 = CT["brg3"][0] / 2, CT["brg4"][0] / 2
    shc = [(xr, 44.0), (xr + 6, 44.0), (xr + 6, CT["r_sh"]), (xb3 - 10, CT["r_sh"]), (xb3 - 10, j3), (xb3 + 10, j3),
           (xb3 + 10, CT["r_sh"]), (xb4 - 12, CT["r_sh"]), (xb4 - 12, j4), (xb4 + 12, j4), (xb4 + 12, CT["r_sh"]),
           (CT["x_end"], CT["r_sh"]), (CT["x_end"], CT["r_bore"]), (xr, CT["r_bore"])]
    poly(shc)
    # шлицы под редуктор — тонкой линией
    poly([(CT["x_end"] - 30, CT["r_sh"] - 1.5), (CT["x_end"], CT["r_sh"] - 1.5)], closed=False, key="thin")
    P["bolts"].append(dict(x0=xf - 6 - 4, x1=xr + 6 + 4, r=CT["bolt_r"], d=CT["bolt_d"],
                           note=f"{CT['bolt_n']} болтов М{CT['bolt_d']:.0f}"))
    P["bearings"].append(dict(x=xb3, d=CT["brg3"][0], D=CT["brg3"][1], B=CT["brg3"][2], kind="roller", name="Опора 3"))
    P["bearings"].append(dict(x=xb4, d=CT["brg4"][0], D=CT["brg4"][1], B=CT["brg4"][2], kind="ball", name="Опора 4"))
    P["masses"]["опоры СТ"] = (xb3, xb4)

    # ─── Выходное устройство ───
    xe0 = x1
    xe1 = xe0 + vu["L_ax"]
    r_h, r_o0, r_o1 = vu["d_h"] / 2, vu["D_in"] / 2, vu["D_1"] / 2
    Rt = vu["R_turn"]
    x_f = xe1 + Rt
    x_b = x_f + vu["b_rad"]
    Ri = 25.0
    arc_o = [(xe1 + Rt * math.cos(a), r_o1 + Rt + Rt * math.sin(a))
             for a in [math.radians(-90 + 90 * i / 12) for i in range(13)]]
    arc_i = [(x_b - Ri + Ri * math.cos(a), r_h + Ri + Ri * math.sin(a))
             for a in [math.radians(-90 + 90 * i / 12) for i in range(13)]]
    r_v = 205.0  # начало раскрытия патрубка к фланцу
    a_f = vu["a_flange"]
    xc_f = (x_f + x_b) / 2
    tilt = math.tan(math.radians(vu["tilt"]))
    fl0 = (xc_f - a_f / 2, 280.0)
    fl1 = (xc_f + a_f / 2, 280.0 + a_f * tilt)
    front = [(xe0, r_o0), (xe1, r_o1)] + arc_o + [(x_f, r_v), fl0]
    rear = [(xe0, r_h)] + [(x_b - Ri, r_h)] + arc_i + [(x_b, r_v), fl1]
    poly(front, closed=False)
    poly(rear, closed=False)
    # выходной фланец патрубка
    nrm = (-math.sin(math.radians(vu["tilt"])), math.cos(math.radians(vu["tilt"])))
    ex = (math.cos(math.radians(vu["tilt"])), math.sin(math.radians(vu["tilt"])))
    for p, s in ((fl0, -1), (fl1, 1)):
        q = (p[0] + s * 16 * ex[0], p[1] + s * 16 * ex[1])
        poly([p, q, (q[0] - 8 * nrm[0], q[1] - 8 * nrm[1]), (p[0] - 8 * nrm[0], p[1] - 8 * nrm[1])])
    # стойка в кольцевом диффузоре (маслоподвод к опоре 3, подвеска кожуха)
    r_o = lambda x: r_o0 + (r_o1 - r_o0) * (x - xe0) / (xe1 - xe0)
    poly([(352.0, r_h), (356.0, r_o(356.0)), (380.0, r_o(380.0)), (376.0, r_h)])
    # корпуса опор 3 и 4 внутри кожуха
    for xb, (d, D, B) in ((xb3, CT["brg3"]), (xb4, CT["brg4"])):
        poly([(xb - B / 2 - 4, D / 2), (xb + B / 2 + 4, D / 2), (xb + B / 2 + 4, D / 2 + 4), (xb + 18, r_h),
              (xb - 18, r_h), (xb - B / 2 - 4, D / 2 + 4)])
    # стенка кожуха в задней стенке патрубка до редуктора
    poly([(x_b, r_h + Ri), (x_b, CT["r_sh"] + 26), (CT["x_end"] - 15, CT["r_sh"] + 26), (CT["x_end"] - 15, 60.0),
          (CT["x_end"] - 5, 60.0), (CT["x_end"] - 5, CT["r_sh"] + 4)], closed=False)

    # ─── Корпусные фланцевые болтовые соединения ───
    for xj, r0, note in ((c["x_out"], c["r_out_tip"], "24 болта М6"), (xe0, r_o0, "16 болтов М6")):
        poly([(xj - 6, r0), (xj - 6, r0 + 18), (xj, r0 + 18), (xj, r0)], closed=False)
        poly([(xj, r0 + 18), (xj + 6, r0 + 18), (xj + 6, r0)], closed=False)
        P["bolts"].append(dict(x0=xj - 10, x1=xj + 10, r=r0 + 11, d=6.0, note=note))

    P["labels"] += [
        ((INLET["x0"] + INLET["x_bend"]) / 2, INLET["r_tip"] + 6, "Входное устройство", 3.5),
        (xb1, -12.0, "Опора 1 (шар.)", 3.5), (xb2, -12.0, "Опора 2 (рол.)", 3.5),
        (xb3, -12.0, "Опора 3 (рол.)", 3.5), (xb4, -12.0, "Опора 4 (шар.)", 3.5),
        ((xe0 + xe1) / 2, r_o1 + 52, "Выходное устройство", 5.0),
        (xc_f, 240.0, "Выхлопной патрубок", 3.5), (xc_f, 233.0, "(выход вверх)", 3.5),
        (CT["x_end"] - 10, CT["r_sh"] + 34, "к редуктору", 3.5),
        (GG["x_front"] + 40, -24.0, "Стяжка", 3.5),
    ]
    P["geo"] = dict(xb1=xb1, xb2=xb2, xb3=xb3, xb4=xb4, x_tvd=(x_rk0 + x_rk1) / 2, x_ct=xm, x_rk_tvd=x_rk1,
                    xe0=xe0, xe1=xe1, r_h=r_h, r_o0=r_o0, r_o1=r_o1, x_f=x_f, x_b=x_b, fl0=fl0, fl1=fl1,
                    r_rad=vu["r_rad"], x_out_c=c["x_out"], r_out_tip=c["r_out_tip"])
    P["vu"] = vu
    return P


if __name__ == "__main__":
    P = build()
    for k, v in P["masses"].items():
        print(k, v)
    print(P["geo"])
