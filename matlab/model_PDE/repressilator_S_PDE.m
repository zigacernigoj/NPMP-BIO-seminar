% ali rob predstavlja konec prostora ali so meje neskonène?
periodic_bounds = 1;
% nalaganje shranjene konfiguracije?
load_conf = 0;
% shranjevanje konène konfiguracije?
save_conf = 0;
% fiksni robovi ali spremenljivi
borderfixed = 0;
% snemanje videa - èasovno potratno
movie_on = 0;

% nalaganje vrednosti parametrov
p = load('params.mat');
alpha = p.alpha;
alpha0 = p.alpha0;
Kd = p.Kd;
delta_m = p.delta_m;
delta_p = p.delta_p;
n = p.n;
beta = p.beta; 
kappa = p.kappa;
kS0 = p.kS0;
kS1 = p.kS1;
kSe = p.kSe;
D1 = p.D1;
eta = p.eta;

size = p.size;                               
density = p.density;
n_cells = ceil(density * size^2);

t_end = p.t_end;
dt = p.dt;
h = p.h; 
h2 = p.h2;

%S=zeros(size,size);                    % Initialise species S
S_e = rand(size,size);
S_i = zeros(size,size);
A = zeros(size,size);
B = zeros(size,size);
C = zeros(size,size);
mA = zeros(size,size);
mB = zeros(size,size);
mC = zeros(size,size);


CELLS = zeros(size,size);
cell_idx = zeros(1,n_cells);
for i = 1:n_cells
    idx = randi(size^2);
    while CELLS(idx) == 1
        idx = randi(size^2);
    end;
    CELLS(idx) = 1;
    cell_idx(i) = idx;
end;


cell_idx = sort(cell_idx);

cell_idx = ceil(size^2 * rand(1, n_cells));
CELLS = zeros(size,size);
CELLS(cell_idx) = 1;

first_idx = cell_idx(1);


A(cell_idx) = 100 * rand(length(cell_idx), 1);
mA(cell_idx) = 100 * rand(length(cell_idx), 1);
B(cell_idx) = 100 * rand(length(cell_idx), 1);
mB(cell_idx) = 100 * rand(length(cell_idx), 1);
C(cell_idx) = 100 * rand(length(cell_idx), 1);
mC(cell_idx) = 100 * rand(length(cell_idx), 1);


A_series = zeros(1,t_end/dt);
S_e_series = zeros(1,t_end/dt);
A_full = zeros(t_end/dt, n_cells);

A_series(1) = A(first_idx);
S_e_series(1) = S_e(first_idx);
A_full(1,:) = A(cell_idx);



if (load_conf)
    cc = load('final_state');
    %cc = load('end_conf.mat');
    %a = cc.a;
    %i = cc.i;
end;

if movie_on == 1
    % Setup image
    ih=imagesc(S_e); set(ih,'cdatamapping','direct')
    %colormap(jet); 
    axis image off; th=title('');
    set(gcf,'position',[100 200 768 768],'color',[1 1 1],'menubar','none')

    % Create 'Quit' pushbutton in figure window
    uicontrol('units','normal','position',[.45 .02 .13 .07], ...
        'callback','set(gcf,''userdata'',1)',...
        'fontsize',10,'string','Quit');
    max_val = 0.1;
    min_val = 0;
end; 


t = 0;
k = 0;
step = 0;
while t <= t_end
      
    if (periodic_bounds)
        S_e_xx= D1 * ([S_e(:,end),S_e(:,1:end-1)] + [S_e(:,2:end),S_e(:,1)] -2*S_e)/h2; 
        S_e_yy= D1 * ([S_e(end,:);S_e(1:end-1,:)] + [S_e(2:end,:);S_e(1,:)] -2*S_e)/h2;                 
    else
        % Create padded matrix to incorporate Neumann boundary conditions 
        SS_e=[[0 S_e(2,:) 0];[S_e(:,2) S_e S_e(:,end-1)];[0 S_e(end-1,:) 0]];
            
        % Calculate diffusion part of the equations
        S_e_xx= D1 * (SS_e(2:end-1,1:end-2) + SS_e(2:end-1,3:end) -2*S_e)/h2; 
        S_e_yy= D1 * (SS_e(1:end-2,2:end-1) + SS_e(3:end,2:end-1) -2*S_e)/h2;        
    end;
    
    D2S_e = S_e_xx + S_e_yy;
    
	
    
	% Calculate dx/dt
    [dmA, dmB, dmC, dA, dB, dC, dS_i, dS_e] = repressilator_S_ODE(CELLS, mA, mB, mC, A, B, C, S_i, S_e, alpha, alpha0, Kd, beta, delta_m, delta_p, n, kS0, kS1, kSe, kappa, eta);
    
    dS_e = dS_e + D2S_e;
    
    
    if (borderfixed == 1)
        % leave border as distrotion centers
        width = length(dS_e);
        dS_e(1:width,[1 width])=0;
        dS_e([1 width],1:width)=0;
        
    end;
        
	mA = mA + dt .* dmA;
    mB = mB + dt .* dmB;
    mC = mC + dt .* dmC;
    A = A + dt .* dA;
    B = B + dt .* dB;
    C = C + dt .* dC;
    S_i = S_i + dt .* dS_i;
    S_e = S_e + dt .* dS_e;
        
	t = t + dt;
	step = step + 1;
    
    A_series(step) = A(first_idx);
    S_e_series(step) = S_e(first_idx);
    A_full(step,:) = A(cell_idx);
	
    if movie_on == 1
        obs = S_e;
        
        % Map v to image grayscale value
        m = 1+round(62*(obs - min_val)/max_val); m=max(m,1); m=min(m,63);
        %m=1+round(63*v); m=max(m,1); m=min(m,64);


         % Update image and text 
        set(ih,'cdata',m);
        set(th,'string',sprintf('%d  %0.2f   %0.2f',t,obs(first_idx)))
        drawnow

        % Write every 500th frame to movie 
        %if rem(n,500)==1
        if rem(step,100)==1
            k=k+1;
            mov(k)=getframe;
        end

    	if ~isempty(get(gcf,'userdata')), break; end % Quit if user clicks on 'Quit' button.
    end;
end;

if (save_conf)
    save('final_state.mat','A');
end
if movie_on == 1
    %Write movie as AVI
    if isunix, sep='/'; else sep='\'; end
    [fn,pn]=uiputfile([pwd sep 'mov.avi'],'Save movie as:');
    if ischar(fn)
        movie2avi(mov,[pn fn],'quality',75)
    else
        disp('User pressed cancel')
    end


    close(gcf)
end;

T=0:dt:t_end-dt;
% figure(1)
% plot(T, S_e_series, 'red')
% hold on
% xlabel('Time [min]')
% ylabel('Concentration [nM]')
% legend('S_e')
% hold off
% 
% figure(2)
% plot(T,A_series, T, S_e_series)
% hold on
% xlabel('Time [min]')
% ylabel('Concentration [nM]')
% legend('A','S_e')
% hold off

figure(3)
TT = T.';
TMat = repmat(TT,1,n_cells);
y = 1:n_cells;
yMat = repmat(y, numel(TT), 1); %//For plot3

plot3(TMat, yMat, A_full,'b')
xlabel('Time [min]'); 
zlabel('Concentration [nM]');
view(0,100);
set(gca, 'XDir','reverse')

grid;

%%%


pos = (0 : size-1)*h;
hm = HeatMap(CELLS, 'RowLabels', pos, 'ColumnLabels', pos);
addXLabel(hm, '\mu m');
addYLabel(hm, '\mu m');
hold off




