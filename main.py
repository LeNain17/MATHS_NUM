from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib.patches import Circle

euleur_position1T = np.loadtxt(open("TAPEI_2000euleur_positionT.txt", 'rt').readlines())
euleur_2000position_t   = euleur_position1T[:, 0]
euleur_2000position_x1   = euleur_position1T[:, 1]
euleur_2000position_x2   = euleur_position1T[:, 2]

euleur_position2T = np.loadtxt(open("TAPEI_20000euleur_positionT.txt", 'rt').readlines())
euleur_20000position_2t   = euleur_position2T[:, 0]
euleur_20000position_2x1   = euleur_position2T[:, 1]
euleur_20000position_2x2   = euleur_position2T[:, 2]

euleur_position3T = np.loadtxt(open("TAPEI_200000euleur_positionT.txt", 'rt').readlines())
euleur_200000position_3t   = euleur_position3T[:, 0]
euleur_200000position_3x1   = euleur_position3T[:, 1]
euleur_200000position_3x2   = euleur_position3T[:, 2]

fig1, axes1 = plt.subplots(nrows=3, ncols=1, figsize=(8, 10))

axes1[0].grid()
axes1[1].grid()
axes1[2].grid()

axes1[0].plot(euleur_2000position_t, euleur_2000position_x1, color='green', label='position1 euleur 2000')
axes1[0].plot(euleur_2000position_t, euleur_2000position_x2, color='pink', label='position2 euler 2000')

axes1[1].plot(euleur_20000position_2t, euleur_20000position_2x1, color='green', label='position1 euleur 20000')
axes1[1].plot(euleur_20000position_2t, euleur_20000position_2x2, color='pink', label='position2 euler 20000')

axes1[2].plot(euleur_200000position_3t, euleur_200000position_3x1, color='green', label='position1 euleur 200000')
axes1[2].plot(euleur_200000position_3t, euleur_200000position_3x2, color='pink', label='position2 euler 200000')

axes1[0].legend(loc='upper right')
axes1[1].legend(loc='upper right')
axes1[2].legend(loc='upper right')

#FIGURE 2 RK4

RK4_position1T = np.loadtxt(open("TAPEI_2000RK4_positionT.txt", 'rt').readlines())
RK4_2000position_t   = RK4_position1T[:, 0]
RK4_2000position_x1   = RK4_position1T[:, 1]
RK4_2000position_x2   = RK4_position1T[:, 2]

RK4_position2T = np.loadtxt(open("TAPEI_20000RK4_positionT.txt", 'rt').readlines())
RK4_20000position_2t   = RK4_position2T[:, 0]
RK4_20000position_2x1   = RK4_position2T[:, 1]
RK4_20000position_2x2   = RK4_position2T[:, 2]

RK4_position3T = np.loadtxt(open("TAPEI_200000RK4_positionT.txt", 'rt').readlines())
RK4_200000position_3t   = RK4_position3T[:, 0]
RK4_200000position_3x1   = RK4_position3T[:, 1]
RK4_200000position_3x2   = RK4_position3T[:, 2]

fig2, axes2 = plt.subplots(nrows=3, ncols=1, figsize=(8, 10))

axes2[0].grid()
axes2[1].grid()
axes2[2].grid()

axes2[0].plot(RK4_2000position_t, RK4_2000position_x1, color='indianred', label='position1 RK4 2000')
axes2[0].plot(RK4_2000position_t, RK4_2000position_x2, color='cyan', label='position2 RK4 2000')

axes2[1].plot(RK4_20000position_2t, RK4_20000position_2x1, color='indianred', label='position1 RK4 20000')
axes2[1].plot(RK4_20000position_2t, RK4_20000position_2x2, color='cyan', label='position2 RK4 20000')

axes2[2].plot(RK4_200000position_3t, RK4_200000position_3x1, color='indianred', label='position1 RK4 200000')
axes2[2].plot(RK4_200000position_3t, RK4_200000position_3x2, color='cyan', label='position2 RK4 200000')

axes2[0].legend(loc='upper right')
axes2[1].legend(loc='upper right')
axes2[2].legend(loc='upper right')

plt.show()
