from parameters import Parameters
import math
import numpy as np

# ali rob predstavlja konec prostora ali so meje neskončne?
periodic_bounds = 1
# nalaganje shranjene konfiguracije?
load_conf = 0
# shranjevanje končne konfiguracije?
save_conf = 0
# fiksni robovi ali spremenljivi
borderfixed = 0
# snemanje videa - časovno potratno
movie_on = 0

# nalaganje vrednosti parametrov
p = params = Parameters()
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
D1 = p.D1
eta = p.eta

size = p.size                               
density = p.density
# number of cells to be added to CELLS matrix
n_cells = math.ceil(density * size**2)

t_end = p.t_end
dt = p.dt
h = p.h 
h2 = p.h2

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
cell_idx = np.argwhere(nums==1)

# locations of cells in 2D matrix
cell_matrix_idx = np.argwhere(CELLS==1)

print("cell idx 1D", cell_idx, cell_idx.shape)
print("cell idx 2D", cell_matrix_idx, cell_matrix_idx.shape)

print(CELLS)
print(CELLS.shape)

if __name__ == "__main__":
    # params = Parameters()
    # print(params.size)
    print("hello world")
