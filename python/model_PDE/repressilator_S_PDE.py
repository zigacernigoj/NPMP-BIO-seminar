from parameters import Parameters
#from repressilator_s_ODE import repressilator_S_ODE
import math
import numpy as np
import scipy
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from multiprocessing import Pool

import time


def calcDm(args):
    A, mA, cells, alpha, alpha0, Kd, delta_m, n, plus = args
    return np.multiply(cells, (
            np.true_divide(alpha, (1 + np.power((np.true_divide(A, Kd)), n))) + alpha0 - delta_m * mA + plus))


def calcD(args):
    A, mA, beta, delta_p = args
    return np.multiply(CELLS, (beta * mA - delta_p * A))


def repressilator_S_ODE(CELLS, mA, mB, mC, A, B, C, S_i, S_e, alpha, alpha0, Kd, beta, delta_m, delta_p, n, kS0, kS1,
                        kSe, kappa, eta, process_pool):
    dS_i = np.multiply(CELLS, (- kS0 * S_i + kS1 * A - eta * (S_i - S_e)))
    dS_e = - kSe * S_e + np.multiply(CELLS, (eta * (S_i - S_e)))

    dm_results = process_pool.map_async(calcDm, [(C, mA, CELLS, alpha, alpha0, Kd, delta_m, n, 0),
                                                 (A, mB, CELLS, alpha, alpha0, Kd, delta_m, n, 0),
                                                 (B, mC, CELLS, alpha, alpha0, Kd, delta_m, n, np.true_divide((kappa * S_i), (1 + S_i)))])

    [d1, d2, d3] = process_pool.map(calcD, [(A, mA, beta, delta_p),
                                                (B, mB, beta, delta_p),
                                                (C, mC, beta, delta_p)])

    [dm1, dm2, dm3] = dm_results.get()

    return [dm1, dm2, dm3, d1, d2, d3, dS_i, dS_e]


start_time = time.time()

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

# I DON'T KNOW
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

A_series[0,0] = A[first_matrix_idx[0], first_matrix_idx[1]]
S_e_series[0,0] = S_e[first_matrix_idx[0], first_matrix_idx[1]]
A_full[0,:] = A[cell_matrix_idx[:,0], cell_matrix_idx[:,1]]

# LOAD CONF
# TODO: convert from matlab code below
# if (load_conf)
#     cc = load('final_state');
#     %cc = load('end_conf.mat');
#     %a = cc.a;
#     %i = cc.i;
# end;

# RECORD THE ACTION
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

# timing sums
integr_time_sum = 0
mul_time_sum = 0
setup_time_sum = 0
other_time_sum = 0

# extract as a constant
onlyZero = np.matrix(0)
D1_div_h2 = D1 / h2

process_pool = Pool(10)

while t < t_end:
    # print("t:", t)
    # print("step:", step)

    if periodic_bounds:
        setup_start = time.time()

        S_e_xx = (np.roll(S_e, 1, axis=1) + np.roll(S_e, -1, axis=1) - 2 * S_e) * D1_div_h2
        S_e_yy = (np.roll(S_e, 1, axis=0) + np.roll(S_e, -1, axis=0) - 2 * S_e) * D1_div_h2

        setup_time_sum += time.time() - setup_start
    else:
        # Create padded matrix to incorporate Neumann boundary conditions

        rowBefore = np.concatenate((onlyZero, S_e[1, :], onlyZero), axis=1)
        rowAfter =  np.concatenate((onlyZero, S_e[-2, :], onlyZero), axis=1)
        
        columnBefore = S_e[:, 1]
        columnAfter = S_e[:, -2]
        between = np.concatenate((columnBefore, S_e, columnAfter), axis=1)

        SS_e = np.concatenate((rowBefore, between, rowAfter), axis=0)

        # Calculate diffusion part of the equations

        leftPart = SS_e[1:-1, 0:-2]
        rightPart = SS_e[1:-1, 2:]
        S_e_xx = D1 * (leftPart + rightPart - 2 * S_e) / h2

        upPart = SS_e[0:-2, 1:-1]
        downPart = SS_e[2:, 1:-1]
        S_e_yy = D1 * (upPart + downPart - 2 * S_e) / h2

    D2S_e = S_e_xx + S_e_yy

    integr_start = time.time()
    # Calculate dx/dt
    [dmA, dmB, dmC, dA, dB, dC, dS_i, dS_e] = repressilator_S_ODE(CELLS, mA, mB, mC, A, B, C, S_i, S_e, alpha, alpha0, Kd, beta, delta_m, delta_p, n, kS0, kS1, kSe, kappa, eta, process_pool)
    integr_time_sum += time.time() - integr_start
    
    dS_e = dS_e + D2S_e

    if borderfixed == 1:
        width = len(dS_e)
        dS_e[0:width-1, 0] = 0
        dS_e[0:width-1, width-1] = 0
         
        dS_e[0, 0:width-1] = 0
        dS_e[width-1, 0:width-1] = 0

    mul_start = time.time()
    mA = mA + np.multiply(dt, dmA)
    mB = mB + np.multiply(dt, dmB)
    mC = mC + np.multiply(dt, dmC)
    A = A + np.multiply(dt, dA)
    B = B + np.multiply(dt, dB)
    C = C + np.multiply(dt, dC)
    S_i = S_i + np.multiply(dt, dS_i)
    S_e = S_e + np.multiply(dt, dS_e)
    mul_time_sum += time.time() - mul_start

    other_start = time.time()
    A_series[0, step] = A[first_matrix_idx[0], first_matrix_idx[1]]
    S_e_series[0, step] = S_e[first_matrix_idx[0], first_matrix_idx[1]]
    A_full[step,:] = A[cell_matrix_idx[:,0], cell_matrix_idx[:,1]]
    other_time_sum += time.time() - other_start

    # increments AFTER ... see if it actually works
    t += dt
    step += 1


    # RECORD THE ACTION
    # TODO: convert from matlab code below
    # if movie_on == 1
    #     obs = S_e;
    #   
    #     % Map v to image grayscale value
    #     m = 1+round(62*(obs - min_val)/max_val); m=max(m,1); m=min(m,63);
    #     %m=1+round(63*v); m=max(m,1); m=min(m,64);
    #
    #      % Update image and text 
    #     set(ih,'cdata',m);
    #     set(th,'string',sprintf('%d  %0.2f   %0.2f',t,obs(first_idx)))
    #     drawnow
    #
    #     % Write every 500th frame to movie 
    #     %if rem(n,500)==1
    #     if rem(step,100)==1
    #         k=k+1;
    #         mov(k)=getframe;
    #     end
    #
    # 	if ~isempty(get(gcf,'userdata')), break; end % Quit if user clicks on 'Quit' button.
    # end;


    # print()

# SAVE CONF
# TODO: convert from matlab code below
# if (save_conf)
#     save('final_state.mat','A');
# end

# RECORD THE ACTION
# TODO: convert from matlab code below
# if movie_on == 1
#     %Write movie as AVI
#     if isunix, sep='/'; else sep='\'; end
#     [fn,pn]=uiputfile([pwd sep 'mov.avi'],'Save movie as:');
#     if ischar(fn)
#         movie2avi(mov,[pn fn],'quality',75)
#     else
#         disp('User pressed cancel')
#     end
#
#     close(gcf)
# end;


print("--- %s seconds ---" % (time.time() - start_time))
print("--- %s seconds for setup ---" % setup_time_sum)
print("--- %s seconds for integration ---" % integr_time_sum)
print("--- %s seconds for other ---" % other_time_sum)
print("--- %s seconds for mul ---" % mul_time_sum)


# GRAPHS
# TODO: convert from matlab code below
# T=0:dt:t_end-dt;

T = np.arange(0, t_end, dt)[np.newaxis]

print("t shape", T.shape)

# % figure(1)
# % plot(T, S_e_series, 'red')
# % hold on
# % xlabel('Time [min]')
# % ylabel('Concentration [nM]')
# % legend('S_e')
# % hold off
# % 
# % figure(2)
# % plot(T,A_series, T, S_e_series)
# % hold on
# % xlabel('Time [min]')
# % ylabel('Concentration [nM]')
# % legend('A','S_e')
# % hold off


TT = T.T
TMat = np.repeat(TT, n_cells, 1)

y = np.arange(0, n_cells)[np.newaxis]
yT = y.T
yMat = np.repeat(yT, TT.size, 1)
yMat = yMat.T

# graf prikazuje koncentracijo molekule v celicah skozi cas
fig, ax1 = plt.subplots(subplot_kw={'projection': '3d'})
ax1.plot_wireframe(TMat, yMat, A_full, rstride=10, cstride=0)
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
    # params = Parameters()
    # print(params.size)
    print("hello world")
