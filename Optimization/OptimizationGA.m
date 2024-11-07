%{ 
    Author: Luis C. Timaná E.
    email : lctimanae@unal.edu.co
    The values ​​of V_crit, K_C and K_G are tuned to find the minimum error. However, initial values ​​for V_crit, K_C and K_G are set, found by trial and error
%}

clear
close all

addpath('..');%Add the previous folder to use the simulation functions
Tmax=12.5e-6;% maximum simulation time
delta_t=1e-8;% simulation step time
N=115;% Number of line segments
L=2222; % total line length (m)
f=100e3;% Frequency of electrical parameters (Hz)
w=2*pi*f;% parameter frequency (rad/s) 

Test=2;% 1: diameter=50.8 mm ACSR 
        % 2: diameter=41.9 mm copper 
        % 3: diameter=23.54 mm ACSR 
Model=2;% 1: Equal capacitance at both ends of the line. Based on average line terminal voltage (VDLM)
        % 2: Different capacitance at both ends of the line (AVDLM)

switch Test
    case {3,8}
        R_atp=3.092278e1; % resistance ohms/km     % Earth resistivity = 56ohm.m, f=100 kHz
        X_atp=1.054232e3; % reactance ohms/km 
        B_atp=4.329376e-3; % reactance s/km 
        % R_atp=2.136361e1; % resistance ohms/km   % Earth resistivity = 20ohm.m, f=100 kHz
        % X_atp=1.039569e3; % reactance ohms/km 
        % B_atp=4.329376e-3; % reactance s/km 
        V_crit=220e3; % Corona inception voltage (V) (value to be optimized). V_crit=470e3 according to Thassio Pereira
        sigma_C=22.5;% Corona model parameters for Tidd Line
        sigma_G=0.35e7;% Corona model parameters for Tidd Line
        r=11.775e-3;% Corona model parameters for Tidd Line (m)
        h=18.89;% Corona model parameters for Tidd Line (m)
        distan_exp = [622  1280  2222];%Distances on the line (from the source) where the measurements were taken (in m)
        load("../Data_Wagner_3.mat")%Experimental data: diameter=23.54 mm ACSR - 1600 kV
        v_exp_1 = v_exp;
        t_exp_1 = t_exp_; 
        load("../Data_Wagner_8.mat")%Experimental data: diameter=23.54 mm ACSR - 1300 kV
        v_exp_2 = v_exp;
        t_exp_2 = t_exp_;
        % t_exp_: sampling times for the voltages at different points on the line. All times start at zero. The initial time is subsequently included
        % v_exp: voltage sampled at different distances on the line
    case {2,6}
        R_atp=2.987880e1; % resistance ohms/km     % Earth resistivity = 56ohm.m, f=100 kHz 
        X_atp=9.808055e2; % reactance ohms/km 
        B_atp=4.662065e-3; % reactance s/km 
        % R_atp=2.031962e1; % resistance ohms/km   % Earth resistivity = 20ohm.m, f=100 kHz
        % X_atp=9.661424e2; % reactance ohms/km 
        % B_atp=4.662065e-3; % reactance s/km 
        V_crit=410e3; % Corona inception voltage (V) (value to be optimized). V_crit=470e3 according to Thassio Pereira
        sigma_C=21;% Corona model parameters for Tidd Line
        sigma_G=1.25e7;% Corona model parameters for Tidd Line
        r=20.95e-3;% Corona model parameters for Tidd Line (m)
        h=18.89;% Corona model parameters for Tidd Line (m)
        distan_exp = [622  1280  2222];%Distances on the line (from the source) where the measurements were taken (in m)
        load("../Data_Wagner_2.mat")%Experimental data: diameter=41.9 mm copper - 1600 kV
        v_exp_1 = v_exp;
        t_exp_1 = t_exp_; 
        load("../Data_Wagner_6.mat")%Experimental data: diameter=41.9 mm ACSR - 1300 kV
        v_exp_2 = v_exp;
        t_exp_2 = t_exp_; 
        % t_exp_: sampling times for the voltages at different points on the line. All times start at zero. The initial time is subsequently included
        % v_exp: voltage sampled at different distances on the line
    case {1,7}
        % R0=0.02; % series resistance (Ohm/m) % Data obtained from (Pereira, 2021)
        % L0=1.49e-6; % series inductance (H/m)
        % C0=7.61e-12; % shunt capacitance (F/m)
        % ZL=443.65; % resistive load connected at the line terminal
        R_atp=2.966895e1; % resistance ohms/km   % Earth resistivity = 56 ohm.m, f=100 kHz
        X_atp=9.563933e2; % reactance ohms/km 
        B_atp=4.784987e-3; % reactance s/km 
        ZL=432; % Data obtained from (Wagner, 1954)
        % R_atp=2.010978e1; % resistance ohms/km     % Earth resistivity = 20 ohm.m, f=100 kHz 
        % X_atp=9.417303e2; % reactance ohms/km 
        % B_atp=4.784987e-3; % reactance s/km 
        V_crit=410e3; % Corona inception voltage (V) (value to be optimized). V_crit=470e3 according to Thassio Pereira
        sigma_C=22;% Corona model parameters for Tidd Line
        sigma_G=0.65e7;% Corona model parameters for Tidd Line
        r=25.4e-3;% Corona model parameters for Tidd Line (m)
        h=18.89;% Corona model parameters for Tidd Line (m)
        distan_exp = [658  1298  2185];%Distances on the line (from the source) where the measurements were taken (in m)
        load("../Data_Wagner.mat")%Experimental data: diameter=50.8 mm ACSR - 1600 kV
        v_exp_1 = v_exp;
        t_exp_1 = t_exp_;
        load("../Data_Wagner_7.mat")%Experimental data: diameter=50.8 mm ACSR - 1300 kV
        v_exp_2 = v_exp;
        t_exp_2 = t_exp_;
        % t_exp_: sampling times for the voltages at different points on the line. All times start at zero. The initial time is subsequently included
        % v_exp: voltage sampled at different distances on the line
end

if ~exist('R0','var')% If there is no R0 variable
    R0=R_atp/1000; % series resistance (Ohm/m)
    L0=X_atp/(1000*w); % series inductance (H/m)
    C0=B_atp/(w*1000); % shunt capacitance (F/m)
end
if ~exist('ZL','var')% If there is no ZL variable
    ZL=(L0/C0)^0.5; % resistive load connected at the line terminal = (L0/C0)^0.5 
end

%%

% Preparing data for optimization
tt=(L0*C0)^0.5*distan_exp;%propagation times from the voltage source to the distances given by vector "d"
Ne=length(distan_exp);

energ_1 = [t_exp_1(:,1)  v_exp_1(:,1)]; % time (t must start at zero) and energizing voltage
v_exp_1 = v_exp_1(:,2:end);% Experimental measurements voltages at values ​​of "t_exp_". "v_exp" can have NaN values
t_exp_1 = t_exp_1(:,2:end);% Sampling times of the measured voltages v_exp
for k=1:Ne
    t_exp_1(:,k)=t_exp_1(:,k)+tt(k);%Initial times are updated for each measurement
end

energ_2 = [t_exp_2(:,1)  v_exp_2(:,1)]; % time (t must start at zero) and energizing voltage
v_exp_2 = v_exp_2(:,2:end);% Experimental measurements voltages at values ​​of "t_exp_". "v_exp" can have NaN values
t_exp_2 = t_exp_2(:,2:end);% Sampling times of the measured voltages v_exp
for k=1:Ne
    t_exp_2(:,k)=t_exp_2(:,k)+tt(k);%Initial times are updated for each measurement
end
%%
% Optimization of parameters V_crit, K_C, K_G using genetic algorithms
% x=[V_crit sigma_C  sigma_G ]; variable to be optimized

%Genetic algorithm optimizes simulations with experimental results for two voltage pulses
objectiveFunction = @(x) exper_simul_comparison(energ_1, energ_2, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, x , r , h , t_exp_1 , v_exp_1, t_exp_2 , v_exp_2, Model);
options = optimoptions('ga','PlotFcn','gaplotbestf','InitialPopulationMatrix',[V_crit sigma_C sigma_G],'MaxGenerations',20);% MaxGenerations=100*numberOfVariables
lb=[0 0 0];% Lower bounds for V_crit, K_C and K_G
ub=[1e6 100 1e8];% Upper bounds for V_crit, K_C and K_G
tic
[x_opt, fval] = ga(objectiveFunction, 3, [], [], [], [], lb, ub, [], options);
timeElapsed = toc;% execution time of the GA
V_crit_op = x_opt(1); % Corona inception voltage (V) - optimal value
sigma_C_op = x_opt(2);% Corona model parameters for Tidd Line - optimal value
sigma_G_op = x_opt(3);% Corona model parameters for Tidd Line - optimal value

K_C_op=sigma_C*(r/(2*h))^0.5*1e-11; %Constant associated with shunt capacitance (optimized)
K_G_op=sigma_G*(r/(2*h))^0.5*1e-11; %Constant associated with shunt resistance (optimized)
%%

switch Model
    case 1
        [t_sim_dist_1,v_sim_dist_1] = exper_simul_comparison_ini(energ_1, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit_op, K_C_op, K_G_op);
        [t_sim_dist_2,v_sim_dist_2] = exper_simul_comparison_ini(energ_2, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit_op, K_C_op, K_G_op);
    case 2
        [t_sim_dist_1,v_sim_dist_1] = exper_simul_comparison_2(energ_1, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit_op, K_C_op, K_G_op);
        [t_sim_dist_2,v_sim_dist_2] = exper_simul_comparison_2(energ_2, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit_op, K_C_op, K_G_op);
end
graphComp(tt, energ_1, t_sim_dist_1, v_sim_dist_1, t_exp_1, v_exp_1)
graphComp(tt, energ_2, t_sim_dist_2, v_sim_dist_2, t_exp_2, v_exp_2)