"""
Продольный разрез (меридиональная схема) всего тракта двигателя ЦБК → КС → ТВД → переходный канал → СТ —
чертёж в КОМПАС-3D, тот же, что строит layout/meridional.py в matplotlib.

Схема строится самим layout/meridional.py:main() (без изменения файлов проекта: plt.show и savefig
подавлены), затем с осей matplotlib снимается геометрия в мм и переносится в КОМПАС:
  - чёрные линии и контуры узлов — основной линией; отверстия ЖТ (белые линии на графике) вырезаются из стенок;
  - границы окрашенных областей (кромки колеса, БЛД, колена, венцы ЛД), не совпадающие с контуром, — тоже основной;
  - направления потоков, струи отверстий, пояса охлаждения, форсунка — тонкими линиями со стрелками;
  - ось двигателя — штрихпунктирной;
  - размеры (⌀ и осевые длины) — размерами КОМПАС в тех же местах и с теми же значениями, что на графике;
  - подписи узлов и легенда отверстий ЖТ — текстом.

Запуск из корня проекта: python kompas_export/draw_layout.py
Результат: kompas_export/meridional_section.cdw (А2, 1:1) и meridional_section.dxf, КОМПАС остаётся открытым.
"""
import math
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.figure  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from draw_meridional import KC, Sheet  # noqa: E402

OUT_DIR = os.path.join(ROOT, "kompas_export")
CDW = os.path.join(OUT_DIR, "meridional_section.cdw")
DXF = os.path.join(OUT_DIR, "meridional_section.dxf")

FLOW_COLORS = {"royalblue", "firebrick", "navy", "teal"}
TOL = 0.35  # мм — допуск совпадения линий


# ─── Снятие геометрии с графика matplotlib ───
def capture():
    import layout.meridional as L

    dims = []
    orig = {n: getattr(L, n) for n in ("_dim_diam", "_dim_width", "_dim_narrow")}

    def rec(name):
        def f(*a, **k):
            dims.append((name, a[1:]))
            return orig[name](*a, **k)
        return f

    for n in orig:
        setattr(L, n, rec(n))
    show, savefig = plt.show, matplotlib.figure.Figure.savefig
    plt.show = lambda *a, **k: None
    matplotlib.figure.Figure.savefig = lambda *a, **k: None
    try:
        L.main()
    finally:
        plt.show, matplotlib.figure.Figure.savefig = show, savefig
        for n, f in orig.items():
            setattr(L, n, f)
    fig = plt.gcf()
    fig.canvas.draw()
    return fig.axes[0], dims


def color_name(c):
    if isinstance(c, str) and c in ("k", "black"):
        return "black"
    if isinstance(c, str):
        return c
    rgba = mcolors.to_rgba(c)
    for name in ("black", "white", "dimgray", "gray") + tuple(FLOW_COLORS):
        if np.allclose(mcolors.to_rgba(name)[:3], rgba[:3], atol=0.01):
            return name
    return "other" if rgba[3] > 0 else "none"


def segs(xs, ys):
    return [((xs[i], ys[i]), (xs[i + 1], ys[i + 1])) for i in range(len(xs) - 1)
            if math.hypot(xs[i + 1] - xs[i], ys[i + 1] - ys[i]) > 1e-6]


def to_data(ax, path_display):
    inv = ax.transData.inverted()
    return [inv.transform(p) for p in path_display.to_polygons(closed_only=False)]


def point_seg(p, a, b):
    """Расстояние от точки до отрезка и параметр проекции."""
    a, b, p = np.asarray(a), np.asarray(b), np.asarray(p)
    ab = b - a
    t = float(np.dot(p - a, ab) / np.dot(ab, ab))
    q = a + min(max(t, 0.0), 1.0) * ab
    return float(np.linalg.norm(p - q)), t


def covered(s, kept):
    """Отрезок s целиком лежит на одном из уже взятых отрезков."""
    for k in kept:
        if all(point_seg(p, *k)[0] < TOL for p in (s[0], s[1], ((s[0][0] + s[1][0]) / 2, (s[0][1] + s[1][1]) / 2))):
            return True
    return False


def cut(segments, holes):
    """Вырезает из отрезков участки, перекрытые белыми линиями (отверстия и окна поясов охлаждения ЖТ)."""
    out = []
    for a, b in segments:
        removed = []
        for p, q in holes:
            (dp, tp), (dq, tq) = point_seg(p, a, b), point_seg(q, a, b)
            # расстояния до прямой AB, а не до отрезка
            ab = np.subtract(b, a)
            nrm = np.array([-ab[1], ab[0]]) / np.linalg.norm(ab)
            if abs(np.dot(np.subtract(p, a), nrm)) < TOL and abs(np.dot(np.subtract(q, a), nrm)) < TOL:
                lo, hi = sorted((tp, tq))
                if hi > 0 and lo < 1:
                    removed.append((max(lo, 0.0), min(hi, 1.0)))
        t0 = 0.0
        for lo, hi in sorted(removed):
            if lo > t0 + 1e-6:
                out.append((lerp(a, b, t0), lerp(a, b, lo)))
            t0 = max(t0, hi)
        if t0 < 1 - 1e-6:
            out.append((lerp(a, b, t0), b))
    return out


def lerp(a, b, t):
    return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t)


def collect(ax):
    main, thin, dashed, axis, holes, flow, nozzle = [], [], [], [], [], [], []
    for ln in ax.lines:
        xs, ys = (np.asarray(v, float) for v in ln.get_data())
        c, ls = color_name(ln.get_color()), ln.get_linestyle()
        s = segs(xs, ys)
        if c == "black" and ls == "-.":
            axis += s
        elif c == "black" and ls == "-":
            # пояса охлаждения (lw 1.3) — тонкие, остальное — контур
            (thin if abs(ln.get_linewidth() - 1.3) < 1e-6 else main).extend(s)
        elif c == "white":
            holes += s
        elif c == "gray" and ls == "--":
            dashed += s
        elif c == "dimgray" and ln.get_linewidth() > 2:
            nozzle += s
        # dimgray ':' и '-' — размерные/выносные линии графика, заменяются размерами КОМПАС

    fills = []
    for p in ax.patches:
        ec = p.get_edgecolor()
        verts = p.get_transform().transform(p.get_path().vertices)
        verts = ax.transData.inverted().transform(verts)
        v = list(verts)
        if not np.allclose(v[0], v[-1]):
            v.append(v[0])
        s = segs([q[0] for q in v], [q[1] for q in v])
        if ec[3] == 0:
            continue
        if color_name(ec) == "black":
            main += s
        elif p.get_linestyle() == "--":
            continue  # осевые зазоры СА–РК: их границы — кромки СА и РК
        else:
            fills += s

    for a in ax.texts:
        if hasattr(a, "arrow_patch") and a.arrow_patch is not None:
            c = color_name(a.arrow_patch.get_edgecolor())
            if c in FLOW_COLORS or c == "firebrick":
                for poly in to_data(ax, a.arrow_patch.get_path()):
                    flow += segs(poly[:, 0], poly[:, 1])

    kept = []
    for s in main:
        if not covered(s, kept):
            kept.append(s)
    main = cut(kept, holes)
    fills = [s for s in fills if not covered(s, kept + thin)]
    fills_kept = []
    for s in fills:
        if not covered(s, fills_kept):
            fills_kept.append(s)
    return dict(main=main + fills_kept, thin=thin, dashed=dashed, axis=axis, flow=flow, nozzle=nozzle)


def texts(ax):
    out = []
    for t in ax.texts:
        s = t.get_text().strip()
        if not s or hasattr(t, "arrow_patch") or color_name(t.get_color()) == "dimgray":
            continue
        x, y = t.get_position()
        fs = t.get_fontsize()
        h = 7.0 if fs >= 12 else 5.0 if fs >= 10 else 3.5
        va = t.get_verticalalignment()
        y = y - h if va == "top" else y - h / 2 if va == "center" else y
        out.append((float(x), float(y), s.replace("подъём 5.5", "подъём 5,5"), h))
    return out


def legend_lines(ax):
    lines = []
    for leg in [c for c in ax.get_children() if isinstance(c, matplotlib.legend.Legend)]:
        title = leg.get_title().get_text()
        if title:
            lines.append(title)
        lines += [t.get_text().replace("\n", " ") for t in leg.get_texts()]
    return list(dict.fromkeys(lines))  # легенда отверстий на графике добавлена дважды


def num(label):
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)", label)
    return float(m.group(1)) if m else None


# ─── Перенос в КОМПАС ───
def draw():
    ax, dims = capture()
    g = collect(ax)
    x_min, x_max = ax.get_xlim()

    k = Sheet(fmt=KC.ksFormatA2, origin=(85.0, 125.0))

    for (a, b) in g["main"]:
        k.line(*a, *b)
    for (a, b) in g["thin"] + g["flow"]:
        k.line(*a, *b, KC.ksCSThin)
    for (a, b) in g["dashed"]:
        k.line(*a, *b, KC.ksCSDashed)
    for (a, b) in g["axis"]:
        k.line(*a, *b, KC.ksCSAxial)
    # форсунка (на графике — толстая линия 3 pt) — полоса шириной 3 мм
    for (a, b) in g["nozzle"]:
        d = np.subtract(b, a)
        n = np.array([-d[1], d[0]]) / np.linalg.norm(d) * 1.5
        pa, pb, ma, mb = a + n, b + n, a - n, b - n
        k.line(*pa, *pb), k.line(*ma, *mb), k.line(*pa, *ma), k.line(*pb, *mb)

    for x, y, s, h in texts(ax):
        k.text(x, y, s, h)

    # Размеры — в тех же местах и с теми же значениями, что на графике
    for name, a in dims:
        if name == "_dim_diam":
            x_block, x_dim, r, label = a[0], a[1], a[2], a[3]
            d = k._dim((x_dim, 0.0), (x_block, r), (x_dim, r / 2), KC.ksLinDVertical, value=num(label), diam=True)
        elif name == "_dim_width":
            x0, x1, y0, y1, y_dim, label = a[:6]
            k._dim((x0, y0), (x1, y1), ((x0 + x1) / 2, y_dim), KC.ksLinDHorizontal, value=dim_value(x1 - x0))
        elif name == "_dim_narrow":
            x0, x1, y, gap, label = a[:5]
            k._dim((x0, y), (x1, y), ((x0 + x1) / 2, y), KC.ksLinDHorizontal, value=dim_value(x1 - x0))

    # Легенда отверстий ЖТ — текстом под схемой
    y = -32.0
    for i, s in enumerate(legend_lines(ax)):
        k.text(x_min + 5, y - 6 * i, re.sub(r"(\d)\.(\d)", r"\1,\2", s.replace("⌀", "диам. ")), 3.5, left=True)

    k.stamp({1: "Двигатель ГТД 800 л.с.\nПродольный разрез проточной части", 2: "ГТД-800.00.001", 6: "1:1"})

    k.doc.SaveAs(CDW)  # возвращаемое значение у КОМПАС v23 ненадёжно — проверяем файл
    if not os.path.isfile(CDW):
        raise RuntimeError(f"КОМПАС не сохранил {CDW} (файл открыт или занят?)")
    ok = k.k5.ActiveDocument2D().ksSaveToDXF(DXF)
    print(f"Сохранено: {CDW}\nDXF: {DXF} ({'ok' if ok else 'ошибка'})")
    print(f"линий: основных {len(g['main'])}, тонких {len(g['thin'])}, потоков {len(g['flow'])}, размеров {len(dims)}")


def dim_value(v):
    v = round(v, 1)
    return float(int(v)) if abs(v - round(v)) < 1e-9 else v


if __name__ == "__main__":
    draw()
