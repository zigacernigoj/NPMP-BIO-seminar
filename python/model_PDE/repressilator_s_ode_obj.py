import numpy as np

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

class Repressilator:
    def __init__(self, cells, alpha, alpha0, kd, beta, delta_m, delta_p, n, k_s0, k_s1, k_se, kappa, eta, D1_div_h2):
        self.cells = cells
        self.alpha = alpha
        self.alpha0 = alpha0
        self.kd = kd
        self.beta = beta
        self.delta_m = delta_m
        self.delta_p = delta_p
        self.n = n
        self.k_s0 = k_s0
        self.k_s1 = k_s1
        self.k_se = k_se
        self.kappa = kappa
        self.eta = eta
        self.D1_div_h2 = D1_div_h2

    def s_ode(self, mA, mB, mC, A, B, C, S_i, S_e):
        def calcDm(A, mA, plus=0):
            return (self.cells * (
                    np.true_divide(self.alpha, (1 + np.power((np.true_divide(A, self.kd)), self.n))) + self.alpha0 - self.delta_m * mA + plus))

        def calcD(A, mA):
            return self.cells * (self.beta * mA - self.delta_p * A)

        two_times_Se = 2 * S_e
        S_e_xx = (shift_right(S_e) + shift_left(S_e) - two_times_Se) * self.D1_div_h2
        S_e_yy = (shift_down(S_e) + shift_up(S_e) - two_times_Se) * self.D1_div_h2

        D2S_e = S_e_xx + S_e_yy

        dS_i = (self.cells * (- self.k_s0 * S_i + self.k_s1 * A - self.eta * (S_i - S_e)))
        dS_e = - self.k_se * S_e + (self.cells * (self.eta * (S_i - S_e)))

        dS_e = dS_e + D2S_e

        return [calcDm(C, mA), calcDm(A, mB), calcDm(B, mC, np.true_divide((self.kappa * S_i), (1 + S_i))),
                calcD(A, mA), calcD(B, mB), calcD(C, mC), dS_i, dS_e]
