"""
Меридиональное сечение проточной части турбины (ТВД → переходный канал → СТ) — чертёж в КОМПАС-3D через API7.

Геометрия берётся из расчётов проекта без их изменения (layout/meridional.py:_get_geometry →
turbine/gasdynamiks.py, free_turbine/gasdynamiks.py, free_turbine/perehod.py); построение ступеней
повторяет turbine/plot_turb.py:plot_meridional_section и layout/meridional.py:draw_axial_turbine/draw_duct.

Запуск из корня проекта (venv проекта, нужен pywin32 и установленный КОМПАС-3D):
    python kompas_export/draw_meridional.py
Результат: kompas_export/meridional_turbine.cdw и meridional_turbine.dxf, КОМПАС остаётся открытым.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

import pythoncom  # noqa: E402
from win32com.client import Dispatch, gencache  # noqa: E402

from layout.meridional import MM, _get_geometry  # noqa: E402

OUT_DIR = os.path.join(ROOT, "kompas_export")
CDW = os.path.join(OUT_DIR, "meridional_turbine.cdw")
DXF = os.path.join(OUT_DIR, "meridional_turbine.dxf")

K5 = gencache.EnsureModule("{0422828C-F174-495E-AC5D-D31014DBBE87}", 0, 1, 0)
K7 = gencache.EnsureModule("{69AC2981-37C0-4379-84FD-5DD2F3C0A520}", 0, 1, 0)
KC = gencache.EnsureModule("{75C9F5D0-B5B8-4526-8681-9903C567D2ED}", 0, 1, 0).constants


# ─── Геометрия (мм), x = 0 — вход в СА ТВД, y = 0 — ось двигателя ───
def stage(x0, g, sa_in=None):
    """Точки ступени СА + зазор + РК (как plot_meridional_section)."""
    gap = 0.25 * g["S_rk"]
    z = [x0, x0 + g["S_sa"], x0 + g["S_sa"] + gap, x0 + g["S_sa"] + gap + g["S_rk"]]
    r_hub, r_tip = g["d_sa"] / 2, g["D_sa"] / 2
    r0_hub, r0_tip = sa_in if sa_in else (r_hub, r_tip)
    hub = [r0_hub, r_hub, r_hub, g["d_rk"] / 2]
    tip = [r0_tip, r_tip, r_tip, g["D_rk"] / 2]
    return dict(z=z, hub=hub, tip=tip, gap=gap, S_sa=g["S_sa"], S_rk=g["S_rk"])


def geometry():
    _, _, tvd, ct, duct = _get_geometry()
    s1 = stage(0.0, tvd)
    x_ct = s1["z"][3] + (duct["x_ct"] - duct["x_rk"]) * MM
    s2 = stage(x_ct, ct, sa_in=(duct["r_in_hub"] * MM, duct["r_in_tip"] * MM))
    d = dict(x0=s1["z"][3], x1=x_ct, hub=(s1["hub"][3], s2["hub"][0]), tip=(s1["tip"][3], s2["tip"][0]),
             h=duct["h"] * MM, gamma=duct["gamma"])
    return s1, s2, d


# ─── КОМПАС ───
class Sheet:
    def __init__(self, fmt=KC.ksFormatA3, origin=(95.0, 95.0)):
        self.app = K7.IApplication(Dispatch("Kompas.Application.7")._oleobj_.QueryInterface(
            K7.IApplication.CLSID, pythoncom.IID_IDispatch))
        self.k5 = K5.KompasObject(Dispatch("Kompas.Application.5")._oleobj_.QueryInterface(
            K5.KompasObject.CLSID, pythoncom.IID_IDispatch))
        self.app.Visible = True
        self.doc = self.app.Documents.Add(KC.ksDocumentDrawing, True)
        doc2d = K7.IKompasDocument2D(self.doc)

        ls = self.doc.LayoutSheets.Item(0)
        ls.Format.Format = fmt
        ls.Format.VerticalOrientation = False
        ls.Update()
        self.layout = ls

        view = doc2d.ViewsAndLayersManager.Views.Add(KC.vt_Normal)
        view.X, view.Y = origin
        view.Scale, view.Name = 1.0, "Меридиональное сечение"
        view.Update()
        view.Current = True
        self.dc = K7.IDrawingContainer(view)
        self.sc = K7.ISymbols2DContainer(view)

    def line(self, x1, y1, x2, y2, style=KC.ksCSNormal):
        s = self.dc.LineSegments.Add()
        s.X1, s.Y1, s.X2, s.Y2, s.Style = x1, y1, x2, y2, style
        s.Update()

    def polyline(self, xs, ys, style=KC.ksCSNormal):
        for i in range(len(xs) - 1):
            self.line(xs[i], ys[i], xs[i + 1], ys[i + 1], style)

    def text(self, x, y, s, height=3.5, angle=0.0, left=False):
        t = self.dc.DrawingTexts.Add()
        t.X, t.Y, t.Angle, t.Allocation = x, y, angle, KC.ksAlLeft if left else KC.ksAlCentre
        item = K7.IText(t).Add().Add()
        item.Str = s
        K7.ITextFont(item).Height = height
        t.Update()

    def _dim(self, p1, p2, p3, orient, value=None, diam=False):
        d = self.sc.LineDimensions.Add()
        d.X1, d.Y1 = p1
        d.X2, d.Y2 = p2
        d.X3, d.Y3 = p3
        d.Orientation = orient
        prm = K7.IDimensionParams(d)
        prm.ArrowPos = KC.ksDimArrowAuto
        dt = K7.IDimensionText(d)
        if diam:
            dt.Sign = 1  # знак ⌀
        if value is not None:
            dt.AutoNominalValue = False
            dt.NominalValue = value
        d.Update()
        return d

    def dim_h(self, x1, y1, x2, y2, y_dim):
        """Осевой размер (горизонтальный)."""
        return self._dim((x1, y1), (x2, y2), ((x1 + x2) / 2, y_dim), KC.ksLinDHorizontal)

    def dim_diam(self, x_surf, r, x_dim):
        """Диаметр тела вращения на полусечении: от оси до поверхности, надпись ⌀2r."""
        return self._dim((x_surf, 0.0), (x_surf, r), (x_dim, r / 2), KC.ksLinDVertical, value=round(2 * r), diam=True)

    def stamp(self, cells):
        st = self.layout.Stamp
        for n, s in cells.items():
            st.Text(n).Str = s
        st.Update()


def draw():
    s1, s2, d = geometry()
    k = Sheet()

    # Контуры ступеней: втулка, периферия, кромки СА и РК
    for s in (s1, s2):
        z, hub, tip = s["z"], s["hub"], s["tip"]
        k.polyline(z, hub)
        k.polyline(z, tip)
        for i in range(4):
            k.line(z[i], hub[i], z[i], tip[i])
        # Осевой зазор между СА и РК — тонкой штриховой серединой, как зазор на графике
        zm = (z[1] + z[2]) / 2
        k.line(zm, hub[1], zm, tip[1], KC.ksCSThin)

    # Переходный канал ТВД → СТ
    k.line(d["x0"], d["hub"][0], d["x1"], d["hub"][1])
    k.line(d["x0"], d["tip"][0], d["x1"], d["tip"][1])

    # Ось двигателя
    x_end = s2["z"][3]
    k.line(-60.0, 0.0, x_end + 65.0, 0.0, KC.ksCSAxial)

    # Подписи венцов и узлов
    for s, name, h_lab in ((s1, "ТВД", 16.0), (s2, "СТ (свободная турбина)", 22.0)):
        z, tip = s["z"], s["tip"]
        k.text((z[0] + z[1]) / 2, max(tip[0], tip[1]) + 4, "СА", 5.0)
        k.text((z[2] + z[3]) / 2, tip[3] + 4, "РК", 5.0)
        k.text((z[0] + z[3]) / 2, max(tip) + h_lab, name, 7.0)
    k.text((d["x0"] + d["x1"]) / 2, d["hub"][0] - 12,
           f"переходный канал: h = {d['h']:.0f} мм, подъём {d['gamma']:.1f}°", 3.5)

    # Диаметры: вход в СА ТВД (слева), выход из РК СТ (справа) — как на графике;
    # выход РК ТВД и вход СА СТ — у кромок, внутри переходного канала
    z, hub, tip = s1["z"], s1["hub"], s1["tip"]
    k.dim_diam(z[0], hub[0], z[0] - 0.16 * tip[0])
    k.dim_diam(z[0], tip[0], z[0] - 0.34 * tip[0])
    k.dim_diam(z[3], hub[3], z[3] + 5)
    k.dim_diam(z[3], tip[3], z[3] + 12)
    z, hub, tip = s2["z"], s2["hub"], s2["tip"]
    k.dim_diam(z[0], hub[0], z[0] - 12)
    k.dim_diam(z[0], tip[0], z[0] - 5)
    k.dim_diam(z[3], hub[3], z[3] + 0.16 * tip[1])
    k.dim_diam(z[3], tip[3], z[3] + 0.34 * tip[1])

    # Осевые длины: СА, РК, зазор (под втулкой), переходный канал и общая длина (над периферией)
    for s in (s1, s2):
        z, hub = s["z"], s["hub"]
        off = 0.14 * s["tip"][1]
        k.dim_h(z[0], hub[0], z[1], hub[1], hub[1] - off)
        k.dim_h(z[2], hub[2], z[3], hub[3], hub[3] - off)
        k.dim_h(z[1], hub[1], z[2], hub[2], hub[1] - 2.1 * off)
    y_top = max(s2["tip"])
    k.dim_h(d["x0"], d["tip"][0], d["x1"], d["tip"][1], y_top + 8)
    k.dim_h(s1["z"][0], s1["tip"][0], x_end, s2["tip"][3], y_top + 42)

    k.stamp({1: "Проточная часть турбины.\nМеридиональное сечение", 2: "ГТД-800.ТУРБ.001", 3: "", 6: "1:1"})

    k.doc.SaveAs(CDW)  # возвращаемое значение у КОМПАС v23 ненадёжно — проверяем файл
    if not os.path.isfile(CDW):
        raise RuntimeError(f"КОМПАС не сохранил {CDW} (файл открыт или занят?)")
    doc5 = k.k5.ActiveDocument2D()
    ok_dxf = doc5.ksSaveToDXF(DXF)
    print(f"Сохранено: {CDW}")
    print(f"DXF: {DXF} ({'ok' if ok_dxf else 'ошибка'})")


if __name__ == "__main__":
    draw()
