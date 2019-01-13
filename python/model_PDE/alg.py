import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import time
import math
from parameters import Parameters
from model_ode import Repressilator


start_time = time.time()


# initial condition
# nalaganje vrednosti parametrov
p = params = Parameters()

# #### fiksirani parametri #### #
# velikost stranice kvadrata, ki predstavlja svet
size = p.size

# gostota "zivih" celic
density = p.density

# number of cells to be added to CELLS matrix
n_cells = math.ceil(density * size ** 2)

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
S_e = np.random.rand(size, size)
S_i = np.zeros((size, size))
A = np.zeros((size, size))
B = np.zeros((size, size))
C = np.zeros((size, size))
mA = np.zeros((size, size))
mB = np.zeros((size, size))
mC = np.zeros((size, size))

# this creates a matrix (CELLS) with randomly located cells
# number of cells = n_cells
# 1 represents a cell, 0 represents an empty place
# indices of cells are in cell_idx (1D)
# and cell_matrix_idx (2D)

# flattened size*size matrix = array with size**2 elements
nums = np.zeros(size ** 2)

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
    A[x, y] = 100 * np.random.rand()
    mA[x, y] = 100 * np.random.rand()
    B[x, y] = 100 * np.random.rand()
    mB[x, y] = 100 * np.random.rand()
    C[x, y] = 100 * np.random.rand()
    mC[x, y] = 100 * np.random.rand()


# time points
t = np.linspace(0, t_end, t_end/dt)

mA = mA.flatten()
mB = mB.flatten()
mC = mC.flatten()
A = A.flatten()
B = B.flatten()
C = C.flatten()
S_i = S_i.flatten()
S_e = S_e.flatten()

flattened = np.concatenate((mA, mB, mC, A, B, C, S_i, S_e))
D1_div_h2 = D1 / h2

r = Repressilator(CELLS, alpha, alpha0, Kd, beta, delta_m, delta_p, n, kS0, kS1, kSe, kappa, eta, D1_div_h2)

print("zacenjam ode")
# solve ODE
y = odeint(r.rep_ode, flattened, t)
print("y sh", y.shape)

# store solution
split_stuff = np.split(y, 8, axis=1)
print("a sh", split_stuff[0].shape)

print("vse -- %s seconds ---" % (time.time() - start_time))

# plot results
plt.plot(t,split_stuff[0])
plt.xlabel('time')
plt.ylabel('z(t)')
plt.show()