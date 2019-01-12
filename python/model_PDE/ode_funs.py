import math
import numpy as np
import scipy


def shift_right(arr):
    arr = arr.reshape(10, 10)
    result = np.empty_like(arr)
    result[:, 1:] = arr[:, :-1]
    result[:, 0] = arr[:, -1]
    return result.flatten()


def shift_left(arr):
    arr = arr.reshape(10, 10)
    result = np.empty_like(arr)
    result[:, :-1] = arr[:, 1:]
    result[:, -1] = arr[:, 0]
    return result.flatten()


def shift_up(arr):
    arr = arr.reshape(10, 10)
    result = np.empty_like(arr)
    result[:-1, :] = arr[1:, :]
    result[-1, :] = arr[0, :]
    return result.flatten()


def shift_down(arr):
    arr = arr.reshape(10, 10)
    result = np.empty_like(arr)
    result[1:, :] = arr[:-1, :]
    result[0, :] = arr[-1, :]
    return result.flatten()


def dmAdt(mA, t, CELLS, alpha, C, Kd, n, alpha0, delta_m):
    return np.multiply(CELLS, (np.true_divide(alpha, (1 + np.power((np.true_divide(C, Kd)), n))) + alpha0 - delta_m * mA))


def dmBdt(mB, t, CELLS, alpha, A, Kd, n, alpha0, delta_m):
    return np.multiply(CELLS, (np.true_divide(alpha, (1 + np.power((np.true_divide(A, Kd)), n))) + alpha0 - delta_m * mB))


def dmCdt(mC, t, CELLS, alpha, B, Kd, n, alpha0, delta_m, kappa, S_i):
    return np.multiply(CELLS, (np.true_divide(alpha, (
                1 + np.power((np.true_divide(B, Kd)), n))) + alpha0 - delta_m * mC + np.true_divide((kappa * S_i),(1 + S_i))))


def dAdt(A, t, CELLS, beta, mA, delta_p):
    print(CELLS.shape, A.shape, mA.shape)
    return np.multiply(CELLS, (beta * mA - delta_p * A))


def dBdt(B, t, CELLS, beta, mB, delta_p):
    return np.multiply(CELLS, (beta * mB - delta_p * B))


def dCdt(C, t, CELLS, beta, mC, delta_p):
    return np.multiply(CELLS, (beta * mC - delta_p * C))


def dS_idt(S_i, t, CELLS, kS0, kS1, A, eta, S_e):
    return np.multiply(CELLS, (- kS0 * S_i + kS1 * A - eta * (S_i - S_e)))


def dS_edt(S_e, t, kSe, CELLS, eta, S_i, D1_div_h2):

    two_times_Se = 2 * S_e
    S_e_xx = (shift_right(S_e) + shift_left(S_e) - two_times_Se) * D1_div_h2
    S_e_yy = (shift_down(S_e) + shift_up(S_e) - two_times_Se) * D1_div_h2

    D2S_e = S_e_xx + S_e_yy

    dS_e = - kSe * S_e + np.multiply(CELLS, (eta * (S_i - S_e)))
    dS_e = dS_e + D2S_e

    return dS_e