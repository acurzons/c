import numpy as np
import matplotlib.pyplot as plt
from fit_str8_line import str8_line

x = [0, 3, 21, 41, 61, 81, 101, 121, 171, 221, 271, 321, 371, 571, 871, 1171, 1471, 1771, 2071, 2371, 2671, 2886, 2971, 3371, 3671, 4071, 4471, 4871, 5156, 5371, 5771, 6071, 6371]
y = [981, 982, 983, 983, 984, 984, 984, 985, 987, 989, 991, 993, 994, 999, 997, 992, 991, 994, 1002, 1017, 1042, 1069, 1050, 953, 874, 760, 641, 517, 427, 355, 218, 122, 0]

x2 = x[21:]
y2 = y[21:]

line_x = np.arange(2800, 7000, 1)

beta, var, var_y, r, p = str8_line(x2, y2, line_x = line_x)
line_y = beta[0] + beta[1] * line_x
print r, p

plt.figure()
plt.scatter(x,y, color = 'k')
plt.plot(line_x, line_y)
plt.xlim((0,6500))
plt.ylim((0,1100))
plt.xlabel("g (cm/s^2)")
plt.ylabel("Depth (km)")
plt.show()
