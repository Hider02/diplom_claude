"""
Расчёт центробежного компрессора, связанный с термодинамическим расчётом двигателя.
π*_к, расход и параметры ЦБК задаются в termo/termodynamiks.py (pik, g_v_pr, KOMP);
подобрать их под заданную мощность можно скриптом kompressor/podbor.py.

Запуск из корня проекта: python kompressor/main.py
"""
import termo.termodynamiks as t
from kompressor.cbk import compressor, print_report
from kompressor import plot

komp = t.komp or compressor(t.pik, t.g_v, t.P1, t.T1, **t.KOMP)
print_report(komp)

print("\n─── Согласование с термодинамическим расчётом ───")
print(f"  η_к в термодинамике = {t.eta_k:.4f} ({'из ЦБК' if t.komp else 'оценочный'}), по газодинамике = {komp['eta_k']:.4f}")
print(f"  N_к термодинамика = {t.N_k / 1000:.1f} кВт, газодинамика = {komp['N'] / 1000:.1f} кВт")
print(f"  T_к* термодинамика = {t.T2:.1f} К, газодинамика = {komp['T_k_star']:.1f} К")
print(f"  Мощность двигателя N_ст = {t.N_ct / t.HP:.1f} л.с. (задано {t.N_target_hp} л.с.)")
print(f"  n = {komp['n_rpm']:.0f} об/мин и D2 = {komp['D2'] * 1000:.1f} мм передаются в расчёт ТВД (общий вал)")

plot.plot_meridional(komp)
plot.plot_velocity_triangles(komp)
