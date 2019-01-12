from parameters import Parameters
import math
import numpy as np
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time

from ode_funs import dmAdt, dmBdt, dmCdt, dAdt, dBdt, dCdt, dS_idt, dS_edt
from scipy.integrate import odeint


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

    # flattened size*size matrix = array with size**2 elements
    nums = np.zeros(size**2)

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

    t = np.linspace(0, t_end, t_end/dt)
    print(t)

    print(CELLS.flatten().shape, mA.flatten().shape, A.flatten().shape)

    # funkcije za klic v scipy ... ODE
    dmA = odeint(dmAdt, mA.flatten(), t, args=(CELLS.flatten(), alpha, C.flatten(), Kd, n, alpha0, delta_m))
    dmB = odeint(dmBdt, mB.flatten(), t, args=(CELLS.flatten(), alpha, A.flatten(), Kd, n, alpha0, delta_m))
    dmC = odeint(dmCdt, mC.flatten(), t, args=(CELLS.flatten(), alpha, B.flatten(), Kd, n, alpha0, delta_m, kappa, S_i.flatten()))

    dA = odeint(dAdt, A.flatten(), t, args=(CELLS.flatten(), beta, mA.flatten(), delta_p))
    dB = odeint(dBdt, B.flatten(), t, args=(CELLS.flatten(), beta, mB.flatten(), delta_p))
    dC = odeint(dCdt, C.flatten(), t, args=(CELLS.flatten(), beta, mC.flatten(), delta_p))

    dS_i = odeint(dS_idt, S_i.flatten(), t, args=(CELLS.flatten(), kS0, kS1, A.flatten(), eta, S_e.flatten()))

    dS_e = odeint(dS_edt, S_e.flatten(), t, args=(kSe, CELLS.flatten(), eta, S_i.flatten(), D1_div_h2))




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
    showPlots = False
    simulate(showPlots)
    print("simulation ended")
