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

# locations of cells in 2D matrix, already sorted
cell_matrix_idx = np.argwhere(CELLS==1)


# TODO: convert from matlab code below
# cell_idx = ceil(size^2 * rand(1, n_cells));
# CELLS = zeros(size,size);
# CELLS(cell_idx) = 1;

first_idx = cell_idx[0]
first_matrix_idx = cell_matrix_idx[0]

for (x,y) in cell_matrix_idx:
    A[x, y] = 100*np.random.rand()
    mA[x, y] = 100*np.random.rand()
    B[x, y] = 100*np.random.rand()
    mB[x, y] = 100*np.random.rand()
    C[x, y] = 100*np.random.rand()
    mC[x, y] = 100*np.random.rand()

A_series = np.zeros((1, int(t_end/dt)))
S_e_series = np.zeros((1, int(t_end/dt)))
A_full = np.zeros((int(t_end/dt), n_cells))

A_series[0] = A[first_matrix_idx[0], first_matrix_idx[1]]
S_e_series[0] = S_e[first_matrix_idx[0], first_matrix_idx[1]]
A_full[0,:] = A[cell_matrix_idx[:,0], cell_matrix_idx[:,1]]

# TODO: convert from matlab code below
# if (load_conf)
#     cc = load('final_state');
#     %cc = load('end_conf.mat');
#     %a = cc.a;
#     %i = cc.i;
# end;

# TODO: convert from matlab code below
# if movie_on == 1
#     % Setup image
#     ih=imagesc(S_e); set(ih,'cdatamapping','direct')
#     %colormap(jet); 
#     axis image off; th=title('');
#     set(gcf,'position',[100 200 768 768],'color',[1 1 1],'menubar','none')
#
#     % Create 'Quit' pushbutton in figure window
#     uicontrol('units','normal','position',[.45 .02 .13 .07], ...
#         'callback','set(gcf,''userdata'',1)',...
#         'fontsize',10,'string','Quit');
#     max_val = 0.1;
#     min_val = 0;
# end; 

t = 0
k = 0
step = 0
while t <= t_end:
    print(t)
    if periodic_bounds:
        S_e_xx = D1 * (np.roll(S_e, 1, axis=1) + np.roll(S_e, -1, axis=1) - 2 * S_e) / h2
        S_e_yy = D1 * (np.roll(S_e, 1, axis=0) + np.roll(S_e, -1, axis=0) - 2 * S_e) / h2                
    else:
        # Create padded matrix to incorporate Neumann boundary conditions

        onlyZero = np.matrix(0)
        rowBefore = np.concatenate((onlyZero, S_e[1,:], onlyZero), axis=1)
        rowAfter =  np.concatenate((onlyZero, S_e[-2,:], onlyZero), axis=1)
        
        columnBefore = S_e[:,1]
        columnAfter = S_e[:, -2]
        between = np.concatenate((columnBefore, S_e, columnAfter), axis=1)

        SS_e = np.concatenate((rowBefore, between, rowAfter), axis=0)

        # Calculate diffusion part of the equations
        # S_e_xx
        # S_e_yy


if __name__ == "__main__":
    # params = Parameters()
    # print(params.size)
    print("hello world")
