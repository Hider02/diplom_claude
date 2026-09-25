import numpy as np
from scipy.optimize import differential_evolution
from kompressor.РК import PK_calculation  # Импортируем готовый расчет
from termo.termodynamiks import *

# ================= ЖЕСТКО ЗАФИКСИРОВАННЫЕ ПАРАМЕТРЫ (ТЗ) =================
# k = 1.4
# R = 287.0
# cp = k * R / (k - 1)
#
# pik = 7
# P0 = 101325.0
# T0 = 288.0
# eta_k = 0.8

h2_target_mm = 10
h2_min = h2_target_mm / 1000.0

t_k = 1.0e-3
t_vt = 1.5e-3
m_gdf = 0.0404
S_coef = 0.2


# ================= 1. ОБОЛОЧКА ДЛЯ ВЫЗОВА РАСЧЕТА =================
def compute_stage(x):
    H_t, beta2l, D1_k, D1_vt, n, G_v = x

    # Базовая геометрическая проверка
    if D1_k <= D1_vt + 0.05:
        return None

    try:
        # Вызываем готовую функцию.
        # ВАЖНО: verbose=False, чтобы не замусоривать консоль при 8000+ итерациях!
        res = PK_calculation(
            G_v=G_v, pik=pik, P0=P0, T0=T0, eta_k=eta_k,
            H_t=H_t, n=n, D1_k=D1_k, D1_vt=D1_vt, beta2l=beta2l,
            S_coef=S_coef, t_k=t_k, t_vt=t_vt, verbose=False
        )

        # Дополнительная проверка на адекватность результатов (защита от NaN)
        if res['h2'] <= 0 or res['eta_rk'] <= 0 or np.isnan(res['eta_rk']):
            return None

        return res

    except Exception:
        # Если внутри PK_calculation возникла математическая ошибка (корень из минуса и т.д.),
        # мы просто возвращаем None, и оптимизатор отбрасывает эту точку.
        return None


# ================= 2. ЦЕЛЕВАЯ ФУНКЦИЯ ДЛЯ ОПТИМИЗАТОРА =================
def evaluate_compressor(x):
    res = compute_stage(x)
    if res is None:
        return 1000.0  # Большое значение штрафа для невалидных точек

    penalty = 0.0

    if res['h2'] < h2_min:
        penalty += 10000.0 * (h2_min - res['h2'])

    if res['u2'] > 700: penalty += 50.0 * (res['u2'] - 600)
    if res['beta1_tip_l'] > 25.0: penalty += 50.0 * (res['beta1_tip_l'] - 25.0)
    if res['lambda_w1'] > 1.15: penalty += 100.0 * (res['lambda_w1'] - 1.15)
    if (res['D1'] - res['D_vt']) < 0.005: penalty += 100.0 * (0.005 - (res['D1'] - res['D_vt']))

    if res['rho_st'] < 0.3: penalty += 20.0 * (0.4 - res['rho_st'])
    if res['rho_st'] > 0.7: penalty += 20.0 * (res['rho_st'] - 0.75)

    return -res['eta_rk'] + penalty


# ================= ЗАПУСК ОПТИМИЗАЦИИ =================
print(f"🚀 Запуск оптимизации с ЦЕЛЕВЫМ ПРИОРИТЕТОМ: h₂ ≥ {h2_target_mm} мм...")

bounds = [
    (0.6, 0.8),  # H_t
    (40.0, 90.0),  # beta2l
    (0.40, 0.90),  # D1_k
    (0.25, 0.50),  # D1_vt
    (35000, 45000),  # n
    (2.5, 3.3),  # G_v
]

# maxiter и popsize можно увеличить для более точного поиска
result = differential_evolution(evaluate_compressor, bounds, seed=42,
                                maxiter=200, popsize=40, tol=1e-4)

# ================= КРАСИВЫЙ ВЫВОД РЕЗУЛЬТАТОВ =================
best_params = result.x
# Для финального вывода включаем verbose=True, чтобы увидеть полную красоту расчета
final_res = PK_calculation(
    G_v=best_params[5], pik=pik, P0=P0, T0=T0, eta_k=eta_k,
    H_t=best_params[0], n=best_params[4], D1_k=best_params[2], D1_vt=best_params[3],
    beta2l=best_params[1], S_coef=S_coef, t_k=t_k, t_vt=t_vt, verbose=True
)

if final_res is not None:
    h2_status = f"✅ ≥ {h2_target_mm} мм" if final_res['h2'] >= h2_min else f"❌ < {h2_target_mm} мм (Физ. предел)"
    splitter_status = "✅ ДА (z > 15)" if final_res['use_splitters'] else "❌ НЕТ (однорядное)"

    print("\n" + "=" * 65)
    print(f"🏆 РЕЗУЛЬТАТ ОПТИМИЗАЦИИ (Приоритет: h₂ ≥ {h2_target_mm} мм)")
    print("=" * 65)
    print("Фиксированные параметры ТЗ:")
    print(f"  n = {final_res['n']} об/мин | G_v = {final_res['G_v']} кг/с | πк = {pik}")
    print("-" * 65)
    print("Найденные оптимальные параметры (Табл. 3.1):")
    print(f"  H_t (Коэфф. напора)   = {best_params[0]:.3f}")
    print(f"  β2л (Угол выхода)     = {best_params[1]:.1f}°")
    print(f"  D1_k (Отн. D периф.)  = {best_params[2]:.3f}")
    print(f"  D1_vt (Отн. D втул.)  = {best_params[3]:.3f}")
    print("-" * 65)
    print("📊 ГАЗОДИНАМИКА И ГЕОМЕТРИЯ (Точный расчет):")
    print(f"  u₂ = {final_res['u2']:.0f} м/с {'✅ <550' if final_res['u2'] < 550 else '❌ >550'}")
    print(f"  D₂ = {final_res['D2'] * 1000:.1f} мм")
    print(f"  D₁ = {final_res['D1'] * 1000:.1f} мм, D_вт = {final_res['D_vt'] * 1000:.1f} мм")
    print(f"  h₁ = {final_res['h1'] * 1000:.1f} мм")
    print(f"  h₂ = {final_res['h2'] * 1000:.2f} мм {h2_status}")
    print(f"  ρ_ст = {final_res['rho_st']:.3f} {'✅' if 0.4 <= final_res['rho_st'] <= 0.75 else '⚠️'}")
    print(f"  Число лопаток z₂ = {final_res['z_rk']} | Сплиттеры: {splitter_status}")
    print("-" * 65)
    print("🔥 ЭФФЕКТИВНОСТЬ:")
    print(f"  КПД рабочего колеса (η_рк) = {final_res['eta_rk']:.3f}  (или {final_res['eta_rk'] * 100:.1f}%)")
    print(f"  Ожидаемый КПД ступени (η_к)= {eta_k:.3f}  (заложено в ТЗ)")
    print("=" * 65)
else:
    print("❌ Оптимизатор не нашел валидного решения.")