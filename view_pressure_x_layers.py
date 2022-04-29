import scipy.interpolate as interp
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA



def linear_scale(serie_or_np_arr, min_val, max_val):
    # https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html
    #   X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    # X_scaled = X_std * (max - min) + min
    min, max = 0, 1
    y = lambda x: ((x - min_val) / (max_val - min_val)) * (max - min) + min
    return y(serie_or_np_arr)


def delinear_scale(serie_or_np_arr, min_val, max_val):
    # https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.MinMaxScaler.html
    #   X_std = (X - X.min(axis=0)) / (X.max(axis=0) - X.min(axis=0))
    # X_scaled = X_std * (max - min) + min
    min, max = 0, 1
    y = lambda x: min_val + (max_val - min_val) * (x - min) / (max - min)
    return y(serie_or_np_arr)


year = '2014'

f = open('/home/denis/_COM_BACKUP/NN_BAM1d/bam1d/work/model/datain/GOAMAZON_{}/SOND_IN_original'.format(year))

arrPress = []
arrTemp = []
for line in f:
    spl = line.split()
    arrPress.append(float(spl[0]))
    arrTemp.append(float(spl[1]))
f.close()

# primeiro maior
arr_levs = [42, 28]
arrPress_rescaled = []
arrTemp_rescaled = []
y=[]
for k, i in zip(arr_levs, range(len(arr_levs))):
    levs = np.array(range(k, 0, -1))
    levs_scaled = linear_scale(levs, 1, k)
    arrPress_scaled = linear_scale(arrPress, np.min(arrPress), np.max(arrPress))
    arrPress_interp = interp.interp1d(np.arange(arrPress_scaled.size), arrPress_scaled)
    arrPress_compress = arrPress_interp(np.linspace(0, arrPress_scaled.size - 1, levs_scaled.size))
    arrPress_rescaled.append(delinear_scale(arrPress_compress, np.min(arrPress), np.max(arrPress)))

    arrTemp_scaled = linear_scale(arrTemp, np.min(arrTemp), np.max(arrTemp))
    arrTemp_interp = interp.interp1d(np.arange(arrTemp_scaled.size), arrTemp_scaled)
    arrTemp_compress = arrTemp_interp(np.linspace(0, arrTemp_scaled.size - 1, levs_scaled.size))
    arrTemp_rescaled.append(delinear_scale(arrTemp_compress, np.min(arrTemp), np.max(arrTemp)))

    y.append(np.linspace(1, k, k))

    print(y)
    print(arrPress_rescaled)
    print(arrTemp_rescaled)

max_press = np.max(arrPress)
min_press = np.min(arrPress)
max_temp = np.max(arrTemp_rescaled[0])
min_temp = np.min(arrTemp_rescaled[0])
ticks_temp = [max_temp] + list(range(25, -80, -5)) + [min_temp]
ticks_press = [max_press] + list(range(950, 0, -50)) + [min_press]

barPress = host_subplot(111, axes_class=AA.Axes)
barPress.figure.set_figwidth(40)
barPress.figure.set_figheight(30)
plt.subplots_adjust(bottom=0.20)
plt.subplots_adjust(top=0.80)

# plt.subplots_adjust(top=0.85)
# plt.subplots_adjust(right=0.75)

plotPress_x = barPress.twinx()
plotPress_y = barPress.twiny()
plotPress = plotPress_x.twiny()
plotTemp28 = barPress.twiny()
plotTemp18 = barPress.twiny()

offset = 40

# new_fixed_axis = plotPress_x.get_grid_helper().new_fixed_axis
# plotPress_x.axis["top"] = new_fixed_axis(loc="top", axes=plotPress_x, offset=(0, 0))
# plotPress_x.set_xlabel(f'Column Level for {arr_levs[1]} Levels')

# plotPress_x.axis["top"].label = f'Column Level for {arr_levs[1]} Levels'
# TODO - nao funcionando
plotPress_x.axis["top"].axis.label = f'Column Level for {arr_levs[1]} Levels'


new_fixed_axis = plotTemp28.get_grid_helper().new_fixed_axis
plotTemp28.axis["bottom"] = new_fixed_axis(loc="bottom", axes=plotTemp28, offset=(0, -offset))
# plotTemp28.axis["bottom"].toggle(all=True)


barPress.set_xlabel(f'Column Level for {arr_levs[0]} Levels')
barPress.set_ylabel("Pressure (mb)")

# plotPress_x.set_ylabel("Pressure (mb) for 18 Levels")

# plotTemp18.set_xlabel("Temperature (C) for 18 levels")
plotTemp28.set_xlabel("Temperature (C)")

p1 = barPress.bar(y[0], arrPress_rescaled[0], label=f'Press {arr_levs[0]} Levels', color='orange', alpha=1.0)
barPress.set_xticks(y[0])
# barPress.set_yticks(ticks_temp)
barPress.invert_yaxis()
barPress.set_ylim(max_press, min_press)
barPress.grid(axis='y')

p2 = plotPress.scatter(y[1], arrPress_rescaled[1],  color='b', alpha=1.0)
plotPress.set_xticks(y[1])
plotPress.invert_yaxis()
plotPress.set_ylim(max_press, min_press)

p2y = plotPress_y.scatter(y[1], arrPress_rescaled[1],  label=f'Press {arr_levs[1]} Levels',  color='b', alpha=1.0)
plotPress_y.set_xticks(y[1])
plotPress_y.invert_yaxis()
plotPress_y.set_ylim(max_press, min_press)

# p3 = plotTemp28.plot(arrTemp_rescaled[0], arrPress_rescaled[0], '-o', color='red', label='Temp 28 Levels', alpha=1.0)
p3 = plotTemp28.plot(arrTemp_rescaled[0], arrPress_rescaled[0], '-o', color='red', label=f'Temp {arr_levs[0]} Levels', alpha=1.0)
# plotTemp28.set_xticks(arrTemp_rescaled[0])
# plotTemp28.set_yticks(arrPress_rescaled[0])
plotTemp28.invert_xaxis()
plotTemp28.set_xticks(ticks_temp)
plotTemp28.set_yticks(ticks_press)

# p4 = plotTemp18.plot(arrTemp_rescaled[1], arrPress_rescaled[1], '-o', color='g', label='Temp 18 Levels', alpha=1.0)
p4 = plotTemp18.plot(arrTemp_rescaled[1], arrPress_rescaled[1], '-o', color='g', label=f'Temp {arr_levs[1]} Levels', alpha=1.0)
plotTemp18.invert_xaxis()
# plotTemp18.set_xticks(arrTemp_rescaled[1])
# plotTemp18.set_yticks(arrPress_rescaled[1])

# plotPress.legend(loc="lower center")
barPress.legend(loc="lower center")

plt.draw()
plt.show()


# fig, ax1 = plt.subplots()
# # press
# ax1.bar(y[0], arrPress_rescaled[0], label='Press 28 Levels', color='orange', alpha=1.0)
# ax1.set_xticks(y[0])
# ax1.set_yticks(arrPress_rescaled[0])
# ax1.invert_yaxis()
# ax1.set_xlabel('Column Level')
# ax1.set_ylabel('Pressure (mb) for 28 Levels')
# ax1.set_ylim(max_press, min_press)
#
# ax2x = ax1.twinx()
# ax2 = ax2x.twiny()
# ax2.scatter(y[1], arrPress_rescaled[1], label='Press 18 Levels', color='b', alpha=0.5)
# ax2.set_xticks(y[1])
# ax2.set_yticks(arrPress_rescaled[1])
# ax2.invert_yaxis()
# ax2x.set_ylabel('Pressure (mb) for 18 Levels')
# ax2.set_ylim(max_press, min_press)
#
# # temp
#
# ax3x = ax2.twinx()
# ax3 = ax2x.twiny()
# ax3.scatter(y[1], arrTemp_rescaled[1], label='Temp 18 Levels', color='g', alpha=1.0)
# ax3.set_xticks(y[1])
# ax3.set_yticks(arrTemp_rescaled[1])
# ax3.invert_yaxis()
# ax3x.set_ylabel('Temperature (C) for 18 Levels')
# ax3.set_ylim(max_temp, min_temp)
#
#
# ax1.grid(axis='y', color='orange')
# # ax2x.grid(axis='y', color='b')
# ax1.legend()
# ax2.legend()
# ax3.legend()
#
# fig.tight_layout()
#
# plt.show()



