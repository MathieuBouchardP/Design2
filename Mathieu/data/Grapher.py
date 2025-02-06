from getcsvData import getcsvData_dict
import matplotlib.pyplot as plt

T_1 = getcsvData_dict("Mathieu\data\data_logT1.csv")
temps, resistance, température = T_1.keys()
T1_temps = T_1[temps]
T1_resistance = T_1[resistance]
T1_température = T_1[température]


T_2 = getcsvData_dict("Mathieu\data\data_logT2.csv")
temps, resistance, température = T_2.keys()
T2_temps = T_2[temps]
T2_resistance = T_2[resistance]
T2_température = T_2[température]

T_3 = getcsvData_dict("Mathieu\data\data_logT3.csv")
temps, resistance, température = T_3.keys()
T3_temps = T_3[temps]
T3_resistance = T_3[resistance]
T3_température = T_3[température]


plt.plot(T1_temps, T1_température)
plt.plot(T2_temps, T2_température)
plt.plot(T3_temps, T3_température)
plt.show()

