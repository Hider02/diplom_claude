"""
Конструкция на листе КОМПАС поверх проточной части (konstrukciya/konstr.py): валы, стяжка, диски ТВД и СТ,
опоры (подшипники — упрощённо по ГОСТ 2.420: контур и тело качения), болтовые соединения (ось — штрихпунктиром,
стержень, головка и гайка — схематично), выходное устройство по расчёту konstrukciya/vyhodnoe_ustr.py,
размеры привязки опор и основные размеры выходного устройства.
"""
import math

from draw_meridional import KC
from konstrukciya import konstr


def bearing(k, b):
    x, r_in, r_out, B = b["x"], b["d"] / 2, b["D"] / 2, b["B"]
    h = r_out - r_in
    x0, x1 = x - B / 2, x + B / 2
    k.polyline([x0, x1, x1, x0, x0], [r_in, r_in, r_out, r_out, r_in])
    # кольца — тонкими линиями, тело качения — основной
    k.line(x0, r_in + 0.22 * h, x1, r_in + 0.22 * h, KC.ksCSThin)
    k.line(x0, r_out - 0.22 * h, x1, r_out - 0.22 * h, KC.ksCSThin)
    rc = (r_in + r_out) / 2
    if b["kind"] == "ball":
        k.circle(x, rc, 0.28 * h)
    else:
        w, hh = 0.55 * B, 0.28 * h
        k.polyline([x - w / 2, x + w / 2, x + w / 2, x - w / 2, x - w / 2], [rc - hh, rc - hh, rc + hh, rc + hh, rc - hh])


def bolt(k, b):
    x0, x1, r, d = b["x0"], b["x1"], b["r"], b["d"]
    hd, nut = 0.7 * d, 0.8 * d
    k.line(x0 - 3, r, x1 + 3, r, KC.ksCSAxial)
    k.polyline([x1 - hd, x1, x1, x1 - hd, x1 - hd], [r - 0.9 * d, r - 0.9 * d, r + 0.9 * d, r + 0.9 * d, r - 0.9 * d])
    k.polyline([x0, x0 + nut, x0 + nut, x0, x0], [r - 0.9 * d, r - 0.9 * d, r + 0.9 * d, r + 0.9 * d, r - 0.9 * d])
    k.line(x0 + nut, r - d / 2, x1 - hd, r - d / 2, KC.ksCSThin)
    k.line(x0 + nut, r + d / 2, x1 - hd, r + d / 2, KC.ksCSThin)


def dim_value(v):
    v = round(v, 1)
    return float(int(v)) if abs(v - round(v)) < 1e-9 else v


def draw_construction(k):
    P = konstr.build()
    g, vu = P["geo"], P["vu"]
    for pl in P["main"]:
        k.polyline([p[0] for p in pl], [p[1] for p in pl])
    for pl in P["thin"]:
        k.polyline([p[0] for p in pl], [p[1] for p in pl], KC.ksCSThin)
    for b in P["bearings"]:
        bearing(k, b)
    for b in P["bolts"]:
        bolt(k, b)
        k.text((b["x0"] + b["x1"]) / 2, b["r"] + 1.2 * b["d"] + 2, b["note"], 2.5)
    for x, y, s, h in P["labels"]:
        k.text(x, y, s, h)

    # Ось двигателя по всей длине (вместе с конструкцией)
    k.line(-135.0, 0.0, konstr.CT["x_end"] + 25, 0.0, KC.ksCSAxial)

    H = KC.ksLinDHorizontal
    # Привязка опор: ГГ — от входа в колесо ЦБК и от задней кромки диска ТВД; СТ — от средней плоскости диска
    rb = lambda b: (b["d"] + b["D"]) / 4
    b1, b2, b3, b4 = P["bearings"]
    k._dim((b1["x"], rb(b1)), (0.0, 0.0), ((b1["x"]) / 2, -36.0), H, value=dim_value(-b1["x"]))
    k._dim((g["x_rk_tvd"], 0.0), (b2["x"], rb(b2)), ((g["x_rk_tvd"] + b2["x"]) / 2, -36.0), H,
           value=dim_value(b2["x"] - g["x_rk_tvd"]))
    k._dim((g["x_ct"], 0.0), (b3["x"], rb(b3)), ((g["x_ct"] + b3["x"]) / 2, -36.0), H,
           value=dim_value(b3["x"] - g["x_ct"]))
    k._dim((g["x_ct"], 0.0), (b4["x"], rb(b4)), ((g["x_ct"] + b4["x"]) / 2, -50.0), H,
           value=dim_value(b4["x"] - g["x_ct"]))

    # Выходное устройство: длина кольцевого участка, ⌀ на его выходе, ширина радиального канала, фланец
    k._dim((g["xe0"], g["r_o0"]), (g["xe1"], g["r_o1"]), ((g["xe0"] + g["xe1"]) / 2, g["r_o1"] + 24), H,
           value=dim_value(g["xe1"] - g["xe0"]))
    k.dim_diam(g["xe1"], g["r_o1"], g["xe1"] - 8)
    k._dim((g["x_f"], 198.0), (g["x_b"], 198.0), ((g["x_f"] + g["x_b"]) / 2, 198.0), H,
           value=dim_value(g["x_b"] - g["x_f"]))
    fl0, fl1 = g["fl0"], g["fl1"]
    a = math.radians(vu["tilt"])
    off = 26.0
    k._dim(fl0, fl1, ((fl0[0] + fl1[0]) / 2 - off * math.sin(a), (fl0[1] + fl1[1]) / 2 + off * math.cos(a)),
           KC.ksLinDParallel, value=dim_value(vu["a_flange"]))
    xm, ym = (fl0[0] + fl1[0]) / 2, (fl0[1] + fl1[1]) / 2
    k.text(xm, ym + off + 18, f"Фланец выхода {vu['a_flange']:.0f} x {vu['c_flange']:.0f} мм, наклон {vu['tilt']:.0f}° к оси",
           3.5)
    k.text(g["x_f"] - 22, g["r_rad"] - 30, f"R{vu['R_turn']:.0f}", 3.5)
    return P
