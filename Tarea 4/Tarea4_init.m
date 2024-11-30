clear
clc

%% PARAMETROS comparison
    Tsim=110;

%% Generation units 
    % abilitation
    e_DG1 = 0;
    e_DG2 = 0;
    e_DG3 = 0;

    %turn on
    %turn off 

%% connect loads
    Ton_L1   =1;
    Ton_L2   =1;
    Ton_L3   =20;
    Ton_L4   =20;
    Ton_L1_DC=1;
    Ton_L2_DC=3;
    Ton_L3_DC=50;
    Ton_L4_DC=50;
    
    %% Disconnect loads
    Toff_L1   =1005;
    Toff_L2   =70;
    Toff_L3   =90;
    Toff_L4   =70;
    Toff_L1_DC=260;
    Toff_L2_DC=250;
    Toff_L3_DC=300;
    Toff_L4_DC=350;

    %% loads (no tocar)
    Rl_1_DC=100;
    Rl_2_DC=100;
    Rl_3_DC=125;
    Rl_4_DC=125;

%% Control activation
e_droop = 0;
e_secondary = 100;    
    
%%    
%%%%%%GENERALES DEL SISTEMA______________________________________________
format short;
V0_dc=400; %Tensión nominal de la micro-red DC
Vmin_dc=V0_dc-5;
Vmax_dc=V0_dc+5;
Ts_sim=1e-3;
Ts=1/16e3;

%% 
%Capacidades máximas de cada unidad
Pmax_1=2500;
Pmax_2=Pmax_1;
Pmax_3=Pmax_1;

Mdc_1=-30/Pmax_1;
Mdc_2=-30/Pmax_2;
Mdc_3=-30/Pmax_3;

%%%%%% ELECTRICOS______________________________________________
%Inductancias de salida
Lout=2.5e-3;
Lout_1=2.5e-3;
Lout_2=2.5e-3;
Lout_3=2.5e-3;
Lout_4=2.5e-3;
Lout_5=2.5e-3;

Lload=2.5e-1;
Cload=200e-6;

%% CONTROL PRIMARIO______________________________________________________
%%%%%%%%%% droop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--> Calculo del filtro  de estimación de P
fc=1;
N_filter=1; %orden de filtro
% Wc=0.5*2*pi;% frecuencia de corte%0.4
Wc=1*2*pi;% frecuencia de corte%0.4
Ws=2*pi/Ts;% frecuencia de muestreo
Wc_norm= Wc/(0.5*Ws);% normalización de frecuencia
[num_f_droop,den_f_droop]=butter(N_filter,Wc_norm,'low'); % calculo de filtro tipo butterworth ([b(1)+b(2)Z^-1+...+b(n+1)Z^-n]/[1+a(2)Z^-1+...+a(n+1)Z^-n])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Control Primario 
Cf=70e-6; %60uF + 10uF
Lf=0.85e-3;

zeta = 0.8;
wnv = 2*pi*50;
Kpv = 2*zeta*wnv*Cf*10;
Kiv = wnv*wnv*Cf*10;

wni = 2*pi*500;
Kpi = 2*zeta*wni*Lf;
Kii = Lf*wni*wni;

%% Despacho Economico
a_Costo_1=0.444;
a_Costo_2=0.264;
a_Costo_3=0.5;

b_Costo_1=0.111;
b_Costo_2=0.067;
b_Costo_3=0.125;
   