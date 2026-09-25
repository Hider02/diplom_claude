from plot import plot_comp
from РК import *
from diffuzor import *
from termo.termodynamiks import g_v_pr, pik, P0, T0

eta_k = 0.8      # Ожидаемый КПД ступени
H_t = 0.8        # Коэфф. теоретического напора
n = 45000   # Частота вращения, об/мин
D1_k = 0.769        # Относительный диаметр входа
D1_vt = 0.414        # Относительный втулочный диаметр
beta2l = 40  # Угол лопатки на выходе РК, град
S_coef = 0.21       # Относительная ширина РК
t_k = 1.0e-3    # Толщина выходной кромки, м (1 мм)
t_vt = 1.0e-3   # Толщина входной кромки, м (1 мм)
resultPK = PK_calculation(pik =pik, P0 = P0, T0 = T0, eta_k = eta_k, G_v = g_v_pr, H_t = H_t, n = n, D1_k = D1_k, D1_vt = D1_vt, beta2l = beta2l, S_coef = S_coef, t_k = t_k, t_vt = t_vt)
diff_calculation(h2=resultPK['h2'], D2=resultPK['D2'], c2r=resultPK['c2r'], c2u=resultPK['c2u'], T2_star=resultPK['T2_star'], alpha2=resultPK['alpha2'], p2=resultPK['p2'],T2=resultPK['T2'], eta_k=eta_k, p2_star=resultPK['p2_star'], c2= resultPK['c2'], P0=P0, dh_prof=resultPK['dh_prof'], dh_out=resultPK['dh_out'], Lu= resultPK['Lu'], beta_tr=resultPK['beta_tr'])