from parameters import Parameters
from repressilator_s_ODE import repressilator_S_ODE
import math
import numpy as np
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

from repressilator_s_ode_obj import Repressilator


def shift_right(arr):
    result = np.empty_like(arr)
    result[:, 1:] = arr[:, :-1]
    result[:, 0] = arr[:, -1]
    return result


def shift_left(arr):
    result = np.empty_like(arr)
    result[:, :-1] = arr[:, 1:]
    result[:, -1] = arr[:, 0]
    return result


def shift_up(arr):
    result = np.empty_like(arr)
    result[:-1, :] = arr[1:, :]
    result[-1, :] = arr[0, :]
    return result


def shift_down(arr):
    result = np.empty_like(arr)
    result[1:, :] = arr[:-1, :]
    result[0, :] = arr[-1, :]
    return result


def simulate(showPlots = False):

    start_time = time.time()

    # fix random seed for testing purposes
    # np.random.seed(1)

    # ali rob predstavlja konec prostora ali so meje neskončne?
    # periodic_bounds = 1
    # nalaganje shranjene konfiguracije?
    # load_conf = 0
    # shranjevanje končne konfiguracije?
    # save_conf = 0
    # fiksni robovi ali spremenljivi
    # borderfixed = 0
    # snemanje videa - časovno potratno
    # movie_on = 0

    # nalaganje vrednosti parametrov
    p = params = Parameters()

    # #### fiksirani parametri #### #
    # velikost stranice kvadrata, ki predstavlja svet
    size = p.size

    # gostota "zivih" celic
    density = p.density

    # number of cells to be added to CELLS matrix
    n_cells = math.ceil(density * size**2)

    # diffusion rate (hitrost sirjenja)
    D1 = p.D1

    # # environment # #
    # stevilo casovnih ciklov
    t_end = p.t_end

    # velikost casovnega cikla
    dt = p.dt

    # Grid size: in micro meters - E coli size 1 um x 2 um (volume = 1 um^3)
    h = p.h
    h2 = p.h2
    # ######## #

    # #### te parametre je treba pogruntat z optimizac. algoritmi #### #
    alpha = p.alpha
    alpha0 = p.alpha0
    Kd = p.Kd
    delta_m = p.delta_m
    delta_p = p.delta_p
    n = p.n
    beta = p.beta
    kappa = p.kappa
    kS0 = p.kS0
    kS1 = p.kS1
    kSe = p.kSe
    eta = p.eta
    # ######## #

    # #### matrike in ostale spremenljivke za delovanje modela #### #

    # S=zeros(size,size) #Initialise species S
    S_e = np.random.rand(size,size)
    S_i = np.zeros((size,size))
    A = np.zeros((size,size))
    B = np.zeros((size,size))
    C = np.zeros((size,size))
    mA = np.zeros((size,size))
    mB = np.zeros((size,size))
    mC = np.zeros((size,size))

    # this creates a matrix (CELLS) with randomly located cells
    # number of cells = n_cells
    # 1 represents a cell, 0 represents an empty place
    # indices of cells are in cell_idx (1D)
    # and cell_matrix_idx (2D)

    # flattened 10*10 matrix = array with 100 elements
    nums = np.zeros(100)

    # some places contain cells (number of cells = n_cells)
    nums[:n_cells] = 1

    # shuffle the array
    np.random.shuffle(nums)

    # reshape array into the matrix
    CELLS = np.reshape(nums, (size, size))

    # locations of cells in 1D array
    cell_idx = np.argwhere(nums == 1)

    # locations of cells in 2D matrix, already sorted
    cell_matrix_idx = np.argwhere(CELLS == 1)

    first_idx = cell_idx[0]
    first_matrix_idx = cell_matrix_idx[0]

    for (x, y) in cell_matrix_idx:
        A[x, y] = 100*np.random.rand()
        mA[x, y] = 100*np.random.rand()
        B[x, y] = 100*np.random.rand()
        mB[x, y] = 100*np.random.rand()
        C[x, y] = 100*np.random.rand()
        mC[x, y] = 100*np.random.rand()

    A_series = np.zeros((1, int(t_end/dt)))
    S_e_series = np.zeros((1, int(t_end/dt)))
    A_full = np.zeros((int(t_end/dt), n_cells))

    A_series[0, 0] = A[first_matrix_idx[0], first_matrix_idx[1]]
    S_e_series[0, 0] = S_e[first_matrix_idx[0], first_matrix_idx[1]]
    A_full[0, :] = A[cell_matrix_idx[:, 0], cell_matrix_idx[:, 1]]

    # TODO: LOAD CONFIGURATION

    t = 0
    k = 0
    step = 0

    # constants
    D1_div_h2 = D1 / h2

    setup_time_sum = 0
    integr_time_sum = 0
    mul_time_sum = 0
    other_time_sum = 0

    r = Repressilator(CELLS, alpha, alpha0, Kd, beta, delta_m, delta_p, n, kS0, kS1, kSe, kappa, eta)

    while t < t_end:
        # print("t:", t)
        # print("step:", step)

        setup_start = time.time()

        two_times_Se = 2 * S_e
        S_e_xx = (shift_right(S_e) + shift_left(S_e) - two_times_Se) * D1_div_h2
        S_e_yy = (shift_down(S_e) + shift_up(S_e) - two_times_Se) * D1_div_h2

        setup_time_sum += time.time() - setup_start

        D2S_e = S_e_xx + S_e_yy

        # Calculate dx/dt
        integr_start = time.time()
        [dmA, dmB, dmC, dA, dB, dC, dS_i, dS_e] = r.s_ode(mA, mB, mC, A, B, C, S_i, S_e)
        integr_time_sum += time.time() - integr_start

        other_start = time.time()
        dS_e = dS_e + D2S_e
        other_time_sum += time.time() - other_start

        mul_start = time.time()
        mA = mA + (dt * dmA)
        mB = mB + (dt * dmB)
        mC = mC + (dt * dmC)
        A = A + (dt * dA)
        B = B + (dt * dB)
        C = C + (dt * dC)
        S_i = S_i + (dt * dS_i)
        S_e = S_e + (dt * dS_e)
        mul_time_sum += time.time() - mul_start

        other_start = time.time()
        A_series[0, step] = A[first_matrix_idx[0], first_matrix_idx[1]]
        S_e_series[0, step] = S_e[first_matrix_idx[0], first_matrix_idx[1]]
        A_full[step,:] = A[cell_matrix_idx[:,0], cell_matrix_idx[:,1]]

        # increments AFTER ... see if it actually works
        t += dt
        step += 1
        other_time_sum += time.time() - other_start

    # TODO: SAVE CONFIGURATION

    print("--- %s seconds ---" % (time.time() - start_time))
    print("--- %s seconds for setup ---" % setup_time_sum)
    print("--- %s seconds for integration ---" % integr_time_sum)
    print("--- %s seconds for mul ---" % mul_time_sum)
    print("--- %s seconds for other ---" % other_time_sum)

    # GRAPHS

    if showPlots:
        print("preparing plots")

        T = np.arange(0, t_end, dt)[np.newaxis]

        # doesn't show data
        # plt.figure(1)
        # plt.plot(T, S_e_series)
        # plt.xlabel('Time [min]')
        # plt.ylabel('Concentration [nM]')

        # doesn't show data
        # plt.figure(2)
        # plt.plot(T,A_series, T, S_e_series)
        # plt.xlabel('Time [min]')
        # plt.ylabel('Concentration [nM]')

        TT = T.T
        TMat = np.repeat(TT, n_cells, 1)

        y = np.arange(0, n_cells)[np.newaxis]
        yT = y.T
        yMat = np.repeat(yT, TT.size, 1)
        yMat = yMat.T

        # graf prikazuje koncentracijo molekule v celicah skozi cas
        fig, ax1 = plt.subplots(subplot_kw={'projection': '3d'})
        ax1.plot_wireframe(TMat, yMat, A_full)
        ax1.set_title("concentration in cells through time")

        # graf prikazuje razporeditev celic
        plt.figure(4)
        plt.imshow(CELLS, cmap='binary')
        plt.xlabel(r'$\mu m$')
        plt.xticks(np.arange(0, size), np.arange(0, size/2, step=0.5))
        plt.ylabel(r'$\mu m$')
        plt.yticks(np.arange(0, size), np.arange(0, size/2, step=0.5))

        # za vse plote prikazat
        plt.show()


if __name__ == "__main__":
    print("starting simulation")
    showPlots = True
    simulate(showPlots)
    print("simulation ended")
