import numpy as np

c0 = 1.527                              #концентрация NaO2 в жидкой фазе шламовой пульпы сгустителя, г/с
consumption_of_sudge_h2o = 2.51         # объем воды, необходимый на промывку 1 т шлама, м^3 - 7-8 м3 на 1 т глинозема '''
consumption_of_circus_h2o = 5.53        #расход оборотной воды, м3/т'''
liquid_solid = 2.5                      #'отношение ж:т'''
C_sudge_h2o = 2.86                      #концентрация натрия в подшламовой воде г/л'''
C_circus_h2o = 0.67                     #концентрация натрия в оборотной воде г/л'''
p_Al = 1.282                            #плотность алюминатного раствора т/м3'''

V_liq = liquid_solid/p_Al
#объем жидкой фазы в шламе сгустителей'''
summ_of_h2o = consumption_of_sudge_h2o + consumption_of_circus_h2o

'''c1,c2,c3,c4,c5'''
left_side = np.array([[summ_of_h2o+V_liq,summ_of_h2o*(-1),0,0,0],
                      [V_liq,(summ_of_h2o + V_liq)*(-1), summ_of_h2o,0,0],
                      [0,V_liq,(summ_of_h2o + V_liq)*(-1),summ_of_h2o,0],
                      [0,0,V_liq,(summ_of_h2o + V_liq)*(-1),summ_of_h2o],
                      [0,0,0,V_liq,(summ_of_h2o+V_liq)*(-1)]])
right_side = np.array([c0*100*V_liq, 0,0,0,(consumption_of_sudge_h2o * C_sudge_h2o + consumption_of_circus_h2o * C_circus_h2o)*(-1)])

rez = np.linalg.inv(left_side).dot(right_side) #вычисляет обратную матрицу, а затем результат
# умножается на вектор правой части, чтобы получить решения.

mass=['c1', 'c2', 'c3', 'c4', 'c5']
for i in range(4):
    print(mass[i],' = ', round(rez[i],3), 'г/л')

#определим степень отмывки шлама
Nu = summ_of_h2o*rez[0]*100/(V_liq*c0*100+consumption_of_sudge_h2o * C_sudge_h2o + consumption_of_circus_h2o * C_circus_h2o)
print(Nu//0.01/100, ' %')
