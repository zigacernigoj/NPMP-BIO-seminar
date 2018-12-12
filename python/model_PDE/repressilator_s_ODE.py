# function  [dmA, dmB, dmC, dA, dB, dC, dS_i, dS_e] 
# =
# repressilator(
#    CELLS, mA, mB, mC, A, B, C, S_i, S_e, 
#    alpha, alpha0, Kd, beta, delta_m, delta_p, 
#    n, kS0, kS1, kSe, kappa, eta)

def repressilator_S_ODE(CELLS, mA, mB, mC, A, B, C, S_i, S_e, alpha, alpha0, Kd, beta, delta_m, delta_p, n, kS0, kS1, kSe, kappa, eta):
    # TODO: rewrite from matlab file
    print("repressilator_S_ODE")


# test data:
# TODO: define the values for variables
# if __name__ == "__main__":
#     repressilator_S_ODE()