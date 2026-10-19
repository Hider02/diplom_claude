"""
Координаты узлов проточной части в системе общей схемы layout/meridional.py (мм; x = 0 — вход в колесо ЦБК,
y — радиус от оси двигателя). Геометрия берётся теми же функциями, что строят график (на невидимой оси),
поэтому конструкция (валы, опоры, выходное устройство) стыкуется с трактом точно.
"""
import contextlib
import io
import os
import sys
from functools import lru_cache

import matplotlib

matplotlib.use("Agg")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib.pyplot as plt  # noqa: E402


@lru_cache(maxsize=1)
def tract():
    import layout.meridional as L
    from kamera_sgoraniya.plot import draw_combustor

    with contextlib.redirect_stdout(io.StringIO()):
        r, ks, tvd, ct, duct = L._get_geometry()
        fig, ax = plt.subplots()
        x_in_c, hub_in_c, tip_in_c, x_out_c, hub_out_c, tip_out_c = L.draw_compressor(ax, r, 0.0, out_dims=False)
        c = draw_combustor(ax, ks, x_out_c, dims=False, duct=duct)
        x_ct = L.draw_duct(ax, duct, x_out_c)
        plt.close(fig)

    def stage(x0, g, sa_in=None):
        gap = 0.25 * g["S_rk"]
        z = [x0, x0 + g["S_sa"], x0 + g["S_sa"] + gap, x0 + g["S_sa"] + gap + g["S_rk"]]
        hub0, tip0 = sa_in if sa_in else (g["d_sa"] / 2, g["D_sa"] / 2)
        return dict(z=z, hub=[hub0, g["d_sa"] / 2, g["d_sa"] / 2, g["d_rk"] / 2],
                    tip=[tip0, g["D_sa"] / 2, g["D_sa"] / 2, g["D_rk"] / 2], S_rk=g["S_rk"])

    mm = 1000.0
    return dict(
        cbk=dict(x0=0.0, S=r["S"], r_vt=r["D1h"] / 2, r1=r["D1t"] / 2, r2=r["D2"] / 2, b2=r["b2"],
                 x_out=x_out_c, r_out_in=hub_out_c, r_out_tip=tip_out_c),
        ks=dict(x_back=c["x_back"], x_end=c["x_end"], x_wall=c["x_wall"], r_n=ks["D_n"] * mm / 2,
                r_floor=duct["r_rk_tip"] * mm),
        tvd=stage(c["xc"], tvd),
        duct=dict(x0=x_out_c + duct["x_rk"] * mm, x1=x_ct,
                  hub=(duct["r_rk_hub"] * mm, duct["r_in_hub"] * mm), tip=(duct["r_rk_tip"] * mm, duct["r_in_tip"] * mm)),
        ct=stage(x_ct, ct, sa_in=(duct["r_in_hub"] * mm, duct["r_in_tip"] * mm)),
    )


if __name__ == "__main__":
    import pprint
    pprint.pprint(tract())
