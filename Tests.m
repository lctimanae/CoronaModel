%{ 
    Author: Luis C. Timaná E.
    email : lctimanae@unal.edu.co
    Example: energization (at t=0) of a line with a resistive load with resistance value equal to the characteristic impedance - N line segments in series
%}

clear
close all

Tmax=12.5e-6;% maximum simulation time
delta_t=1e-8;% simulation step time
N=115;% Number of line segments
L=2222; % total line length (m)
f=100e3;% Frequency of electrical parameters (Hz)
w=2*pi*f;% parameter frequency (rad/s) 

Test=6;% 1: diameter=50.8 mm ACSR - 1600 kV
        % 2: diameter=41.9 mm copper - 1600 kV
        % 3: diameter=23.54 mm ACSR - 1600 kV
        % 6: diameter=41.9 mm ACSR - 1300 kV
        % 7: diameter=50.8 mm ACSR - 1300 kV
        % 8: diameter=23.54 mm ACSR - 1300 kV
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
        V_crit=220.0e3; % Corona inception voltage (V) (value to be optimized). V_crit=470e3 according to Thassio Pereira
        sigma_C=22.44;% Corona model parameters for Tidd Line
        sigma_G=0.8382e7;% Corona model parameters for Tidd Line
        r=11.775e-3;% Corona model parameters for Tidd Line (m)
        h=18.89;% Corona model parameters for Tidd Line (m)
        distan_exp = [622  1280  2222];%Distances on the line (from the source) where the measurements were taken (in m)
        if Test==3
            load("Data_Wagner_3.mat")%Experimental data: diameter=23.54 mm ACSR - 1600 kV
        else 
            load("Data_Wagner_8.mat")%Experimental data: diameter=23.54 mm ACSR - 1300 kV
        end
        % t_exp_: sampling times for the voltages at different points on the line. All times start at zero. The initial time is subsequently included
        % v_exp: voltage sampled at different distances on the line
    case {2,6}
        R_atp=2.987880e1; % resistance ohms/km     % Earth resistivity = 56ohm.m, f=100 kHz 
        X_atp=9.808055e2; % reactance ohms/km 
        B_atp=4.662065e-3; % reactance s/km 
        % R_atp=2.031962e1; % resistance ohms/km   % Earth resistivity = 20ohm.m, f=100 kHz
        % X_atp=9.661424e2; % reactance ohms/km 
        % B_atp=4.662065e-3; % reactance s/km 
        V_crit=410.0e3; % Corona inception voltage (V) (value to be optimized). V_crit=470e3 according to Thassio Pereira
        sigma_C=21.06;% Corona model parameters for Tidd Line
        sigma_G=1.361e7;% Corona model parameters for Tidd Line
        r=20.95e-3;% Corona model parameters for Tidd Line (m)
        h=18.89;% Corona model parameters for Tidd Line (m)
        distan_exp = [622  1280  2222];%Distances on the line (from the source) where the measurements were taken (in m)
        if Test==2
            load("Data_Wagner_2.mat")%Experimental data: diameter=41.9 mm copper - 1600 kV
        else
            load("Data_Wagner_6.mat")%Experimental data: diameter=41.9 mm ACSR - 1300 kV
        end
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
        V_crit=421.8e3; % Corona inception voltage (V) (value to be optimized). V_crit=470e3 according to Thassio Pereira
        sigma_C=21.00;% Corona model parameters for Tidd Line
        sigma_G=0.8382e7;% Corona model parameters for Tidd Line
        r=25.4e-3;% Corona model parameters for Tidd Line (m)
        h=18.89;% Corona model parameters for Tidd Line (m)
        distan_exp = [658  1298  2185];%Distances on the line (from the source) where the measurements were taken (in m)
        if Test==1
            load("Data_Wagner.mat")%Experimental data: diameter=50.8 mm ACSR - 1600 kV
        else
            load("Data_Wagner_7.mat")%Experimental data: diameter=50.8 mm ACSR - 1300 kV
        end
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
K_C=sigma_C*(r/(2*h))^0.5*1e-11; %Constant associated with shunt capacitance (value to be optimized)
K_G=sigma_G*(r/(2*h))^0.5*1e-11; %Constant associated with shunt resistance (value to be optimized)

energ = [t_exp_(:,1)  v_exp(:,1)]; % time (t must start at zero) and energizing voltage
v_exp = v_exp(:,2:end);% Experimental measurements voltages at values ​​of "t_exp_". "v_exp" can have NaN values
t_exp_ = t_exp_(:,2:end);% Sampling times of the measured voltages v_exp
tt=(L0*C0)^0.5*distan_exp;%propagation times from the voltage source to the distances given by vector "d"
Ne=length(distan_exp);
for k=1:Ne
    t_exp_(:,k)=t_exp_(:,k)+tt(k);%Initial times are updated for each measurement
end

switch Model
    case 1
        [t_sim_dist,v_sim_dist] = exper_simul_comparison_ini(energ, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit, K_C, K_G);
    case 2
        [t_sim_dist,v_sim_dist] = exper_simul_comparison_2(energ, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit, K_C, K_G);
end
tm=tt+4.9e-6;%maximum time of available experimental data. 4.9e-6 is the maximum time span recorded for each voltage (see vector t_exp_)
figure(1)
plot(energ(1:end,1),energ(1:end,2),'k')% experimental data - energizing voltage
hold on
for k=1:Ne
    t_temp=t_sim_dist(t_sim_dist<tm(k));%Simulated data are plotted over a time interval similar to the experimental data
    N_temp=length(t_temp);
    plot(t_temp,v_sim_dist(k,1:N_temp),'r')%simulated data
    plot(t_exp_(:,k),v_exp(:,k),'b')%experimental data
end
legend('Energizing voltage','Simulated data','Experimental data')
grid on
xlabel('Time [s]')
ylabel('Voltage [V]')

error = error_calculation(t_sim_dist , v_sim_dist ,  t_exp_ , v_exp);
X = ['Percentage error between simulated and experimental values: ',num2str(error*100),'%'];
disp(X)