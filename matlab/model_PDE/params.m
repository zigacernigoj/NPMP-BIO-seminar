% repressilator parameters
%alpha = 216;
%alpha0 = 0.001 * alpha;
%n = 2;
%beta = 5; 

alpha = 5; % min**(-1) 
alpha0 = 0.001 * alpha; % min**(-1) 
beta = 1; % min**(-1) 
Kd = 10; % nM
delta_p = 0.1;   % min**(-1) 
delta_m = 0.1;   % min**(-1) 
n = 2;

kappa = 0.2;
kS0 = 1;
kS1 = 0.01;
kSe = 0.01;
eta = 2;

% diffusion rate
D1=0.5;

% environment
size = 10; % size of the space equals [size x size]
density = 0.8;


% simulation parameters 
t_end = 9999;
dt = 0.1;
h=0.5; % Grid size: in micro meters - E coli size 1 um x 2 um (volume = 1 um^3)
h2=h^2;

save('params')

%repressilator_PDE