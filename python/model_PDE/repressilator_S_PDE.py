from parameters import Parameters
import math
import numpy as np
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from scipy.optimize import differential_evolution, dual_annealing

import os
from multiprocessing import Process

from repressilator_s_ode_obj import Repressilator
from optimize_params import get_sync_index


def simulate(inputs):

    print('starting process with PID {}'.format(os.getpid()))

    # fix random seed for testing purposes
    # np.random.seed(1)

    # ali rob predstavlja konec prostora ali so meje neskoncne?
    # periodic_bounds = 1
    # nalaganje shranjene konfiguracije?
    # load_conf = 0
    # shranjevanje koncne konfiguracije?
    # save_conf = 0
    # fiksni robovi ali spremenljivi
    # borderfixed = 0
    # snemanje videa - casovno potratno
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
    alpha = inputs[0]
    alpha0 = inputs[1]
    Kd = inputs[2]
    delta_m = inputs[3]
    delta_p = inputs[4]
    n = inputs[5]
    beta = inputs[6]
    kappa = inputs[7]
    kS0 = inputs[8]
    kS1 = inputs[9]
    kSe = inputs[10]
    eta = inputs[11]
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

    integr_time_sum = 0
    mul_time_sum = 0
    other_time_sum = 0

    r = Repressilator(CELLS, alpha, alpha0, Kd, beta, delta_m, delta_p, n, kS0, kS1, kSe, kappa, eta, D1_div_h2)

    start_time = time.time()

    try:

        while t < t_end:
            # print("t:", t)
            # print("step:", step)

            # Calculate dx/dt
            integr_start = time.time()
            with np.errstate(invalid='raise'):
                [dmA, dmB, dmC, dA, dB, dC, dS_i, dS_e] = r.s_ode(mA, mB, mC, A, B, C, S_i, S_e)
            integr_time_sum += time.time() - integr_start

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

        printTimes = False

        if printTimes:
            print("--- %s seconds ---" % (time.time() - start_time))
            print("--- %s seconds for integration ---" % integr_time_sum)
            print("--- %s seconds for mul ---" % mul_time_sum)
            print("--- %s seconds for other ---" % other_time_sum)
            print()

        print("racunam primernost parametrov")

        # tabela,
        # razlika med min in max,
        # stevilo zaporednih casovnih korakov, ki zadostujejo razliki med min in max
        first_synced_index = get_sync_index(A_full, 2, 3)

        print("first synced", first_synced_index)


        # GRAPHS
        showPlots = False
        if (first_synced_index < math.inf) and showPlots:
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
            # fig, ax1 = plt.subplots(subplot_kw={'projection': '3d'})
            # ax1.plot_wireframe(TMat, yMat, A_full)
            # ax1.set_title("concentration in cells through time")

            # graf prikazuje razporeditev celic
            # plt.figure(4)
            # plt.imshow(CELLS, cmap='binary')
            # plt.xlabel(r'$\mu m$')
            # plt.xticks(np.arange(0, size), np.arange(0, size/2, step=0.5))
            # plt.ylabel(r'$\mu m$')
            # plt.yticks(np.arange(0, size), np.arange(0, size/2, step=0.5))

            plt.figure(5)
            plt.plot(TT, A_full)
            plt.xlabel('time')
            plt.ylabel('concentration')

            # za vse plote prikazat
            plt.show()

        saveConfig = False
        if (first_synced_index < math.inf) and saveConfig:
            name = './reports/' + str(os.getpid()) + '_' + str(time.time()) + 'score_' + str(first_synced_index) + '.txt'
            with open(name, 'w') as file:
                string = '\n'.join(str(e) for e in inputs)
                file.write(string + '\n-----\n' + str(first_synced_index))

        print('ending process with PID {}\n'.format(os.getpid()))

        # return first synced index so you know in scipy.optimize.anneal what is better
        return first_synced_index

    except:
        return math.inf


if __name__ == "__main__":
    # num_of_processes = 4
    # print("starting {} simulations".format(num_of_processes))
    # showPlots = False
    # processes = [Process(target=simulate, args=(showPlots, )) for _ in range(num_of_processes)]
    #
    # for p in processes:
    #     p.start()
    #     time.sleep(1)
    #
    # for p in processes:
    #     p.join()

    print("starting simulation with optimization algorithms")

    b_alpha = (0.001, 10)
    b_alpha0 = (0.001, 10)
    b_Kd = (0.01, 100)
    b_delta_m = (0.001, 10)
    b_delta_p = (0.001, 10)
    b_n = (1, 4)
    b_beta = (0.001, 10)
    b_kappa = (0.001, 10)
    b_kS0 = (0.001, 10)
    b_kS1 = (0.001, 10)
    b_kSe = (0.001, 10) # problems - overflow
    b_eta = (0.01, 100) # problems - overflow

    bounds = [b_alpha, b_alpha0, b_Kd, b_delta_m, b_delta_p, b_n, b_beta, b_kappa, b_kS0, b_kS1, b_kSe, b_eta]

    alpha = 5
    alpha0 = 0.001 * alpha
    Kd = 10
    delta_m = 0.1
    delta_p = 0.1
    n = 2
    beta = 1
    kappa = 0.2
    kS0 = 1
    kS1 = 0.01
    kSe = 0.01
    eta = 2

    init = [[alpha, alpha0, Kd, delta_m, delta_p, n, beta, kappa, kS0, kS1, kSe, eta],
            [alpha, alpha0, Kd, delta_m, delta_p, n, beta, kappa, kS0, kS1, kSe, eta],
            [alpha, alpha0, Kd, delta_m, delta_p, n, beta, kappa, kS0, kS1, kSe, eta],
            [alpha, alpha0, Kd, delta_m, delta_p, n, beta, kappa, kS0, kS1, kSe, eta],
            [alpha, alpha0, Kd, delta_m, delta_p, n, beta, kappa, kS0, kS1, kSe, eta]]

    x0 = [alpha, alpha0, Kd, delta_m, delta_p, n, beta, kappa, kS0, kS1, kSe, eta]

    # diff evol
    # result = differential_evolution(simulate, bounds, init=init, updating='deferred', workers=4, strategy='currenttobest1bin', maxiter=250)

    # dual anneal, CHANGE maxiter to get BETTER results
    result = dual_annealing(simulate, bounds=bounds, seed=1234, maxiter=1, x0=x0)


    print(result.x)
    print(result.fun)
    print("simulation ended")

