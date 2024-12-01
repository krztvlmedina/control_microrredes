run('Tarea4_init.m')
 
Ny = 15;
Nu = 10;

lambda1 = 7;
lambda2 = 800;
lambda3 = 0.001;

parametros_mpc = [Ny, Nu, lambda1, lambda2, lambda3];

a_adj = ones(3);
a1 = a_adj(1, :);
a2= a_adj(2, :);
a3 = a_adj(3, :);

mpc_sample_time = 3; %segundos
tau = 1*mpc_sample_time; % retardo comunicaciones


%% test
% [u, x_pred] = dpmc(1, V0_dc,lambda1, lambda2, lambda3, Ny, Nu,Mdc_1, V0_dc, V0_dc, 1200, Pmax_1, Pmax_1-300, Vmin_dc, Vmax_dc, a_Costo_1, b_Costo_1, a1, 1/Rl_1_DC, {zeros(4,6), zeros(4,6), zeros(4,6)});

%%
options = simset('SrcWorkspace','current');
sim('Tarea4_base.slx',[],options)