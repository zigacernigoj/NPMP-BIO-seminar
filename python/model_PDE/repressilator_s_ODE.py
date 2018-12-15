# function  [dmA, dmB, dmC, dA, dB, dC, dS_i, dS_e] 
# =
# repressilator(
#    CELLS, mA, mB, mC, A, B, C, S_i, S_e, 
#    alpha, alpha0, Kd, beta, delta_m, delta_p, 
#    n, kS0, kS1, kSe, kappa, eta)
import numpy as np

#TODO checck
def repressilator_S_ODE(CELLS, mA, mB, mC, A, B, C, S_i, S_e, alpha, alpha0, Kd, beta, delta_m, delta_p, n, kS0, kS1, kSe, kappa, eta):
    #rewrite from matlab file
    print("repressilator_S_ODE")
    dmA = np.multiply(CELLS, (np.true_divide(alpha, (1 + np.power((np.true_divide(C, Kd)), n))) + alpha0 - delta_m * mA))
    dmB = np.multiply(CELLS, (np.true_divide(alpha, (1 + np.power((np.true_divide(A, Kd)), n))) + alpha0 - delta_m * mB))
    dmC = np.multiply(CELLS, (np.true_divide(alpha, (1 + np.power((np.true_divide(B, Kd)), n))) + alpha0 - delta_m * mC + np.true_divide((kappa * S_i), (1 + S_i))))

    dA = np.multiply(CELLS, (beta * mA - delta_p * A))
    dB = np.multiply(CELLS, (beta * mB - delta_p * B))
    dC = np.multiply(CELLS, (beta * mC - delta_p * C))

    dS_i = np.multiply(CELLS, (- kS0 * S_i + kS1 * A - eta * (S_i - S_e)))
    dS_e = - kSe * S_e + CELLS. * (eta * (S_i - S_e))

    return [dmA, dmB, dmC, dA, dB, dC, dS_i, dS_e]


# test data:
# TODO: define the values for variables
# if __name__ == "__main__":
#     repressilator_S_ODE()