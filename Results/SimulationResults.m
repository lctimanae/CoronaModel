clear
close all

Test=1; % 1: diameter=50.8 mm ACSR - 1600 kV
        % 2: diameter=41.9 mm copper - 1600 kV
        % 3: diameter=23.54 mm ACSR - 1600 kV

Voltage=1; % 1: 1600 kV
           % 2: 1300 kV

f=100e3;% Frequency of electrical parameters (Hz)
w=2*pi*f;% parameter frequency (rad/s) 

if (Test==1 && Voltage==1)
    load("../Data_Wagner.mat");
    load("ACSR_50mm_1600kV_Model1.mat");
    t_sim_dist_1 = t_sim_dist;
    v_sim_dist_1 = v_sim_dist;
    load("ACSR_50mm_1600kV_Model2.mat");
    t_sim_dist_2 = t_sim_dist;
    v_sim_dist_2 = v_sim_dist;
    distan_exp = [658  1298  2185];%Distances on the line (from the source) where the measurements were taken (in m)
    R_atp=2.966895e1; % resistance ohms/km   % Earth resistivity = 56 ohm.m, f=100 kHz
    X_atp=9.563933e2; % reactance ohms/km 
    B_atp=4.784987e-3; % reactance s/km 
    ZL=432; % Data obtained from (Wagner, 1954)
elseif (Test==1 && Voltage==2)
    load("../Data_Wagner_7.mat");
    load("ACSR_50mm_1300kV_Model1.mat");
    t_sim_dist_1 = t_sim_dist;
    v_sim_dist_1 = v_sim_dist;
    load("ACSR_50mm_1300kV_Model2.mat");
    t_sim_dist_2 = t_sim_dist;
    v_sim_dist_2 = v_sim_dist;
    distan_exp = [658  1298  2185];%Distances on the line (from the source) where the measurements were taken (in m)
    R_atp=2.966895e1; % resistance ohms/km   % Earth resistivity = 56 ohm.m, f=100 kHz
    X_atp=9.563933e2; % reactance ohms/km 
    B_atp=4.784987e-3; % reactance s/km 
    ZL=432; % Data obtained from (Wagner, 1954)
elseif (Test==2 && Voltage==1)
    load("../Data_Wagner_2.mat");
    load("ACSR_41mm_1600kV_Model1.mat");
    t_sim_dist_1 = t_sim_dist;
    v_sim_dist_1 = v_sim_dist;
    load("ACSR_41mm_1600kV_Model2.mat");
    t_sim_dist_2 = t_sim_dist;
    v_sim_dist_2 = v_sim_dist;
    distan_exp = [622  1280  2222];%Distances on the line (from the source) where the measurements were taken (in m)
    R_atp=2.987880e1; % resistance ohms/km     % Earth resistivity = 56ohm.m, f=100 kHz 
    X_atp=9.808055e2; % reactance ohms/km 
    B_atp=4.662065e-3; % reactance s/km 
elseif (Test==2 && Voltage==2)
    load("../Data_Wagner_6.mat");
    load("ACSR_41mm_1300kV_Model1.mat");
    t_sim_dist_1 = t_sim_dist;
    v_sim_dist_1 = v_sim_dist;
    load("ACSR_41mm_1300kV_Model2.mat");
    t_sim_dist_2 = t_sim_dist;
    v_sim_dist_2 = v_sim_dist;
    distan_exp = [622  1280  2222];%Distances on the line (from the source) where the measurements were taken (in m)
    R_atp=2.987880e1; % resistance ohms/km     % Earth resistivity = 56ohm.m, f=100 kHz 
    X_atp=9.808055e2; % reactance ohms/km 
    B_atp=4.662065e-3; % reactance s/km 
elseif (Test==3 && Voltage==1)
    load("../Data_Wagner_3.mat");
    load("ACSR_23mm_1600kV_Model1.mat");
    t_sim_dist_1 = t_sim_dist;
    v_sim_dist_1 = v_sim_dist;
    load("ACSR_23mm_1600kV_Model2.mat");
    t_sim_dist_2 = t_sim_dist;
    v_sim_dist_2 = v_sim_dist;
    distan_exp = [622  1280  2222];%Distances on the line (from the source) where the measurements were taken (in m)
    R_atp=3.092278e1; % resistance ohms/km     % Earth resistivity = 56ohm.m, f=100 kHz
    X_atp=1.054232e3; % reactance ohms/km 
    B_atp=4.329376e-3; % reactance s/km 
elseif (Test==3 && Voltage==2)
    load("../Data_Wagner_8.mat");
    load("ACSR_23mm_1300kV_Model1.mat");
    t_sim_dist_1 = t_sim_dist;
    v_sim_dist_1 = v_sim_dist;
    load("ACSR_23mm_1300kV_Model2.mat");
    t_sim_dist_2 = t_sim_dist;
    v_sim_dist_2 = v_sim_dist;
    distan_exp = [622  1280  2222];%Distances on the line (from the source) where the measurements were taken (in m)
    R_atp=3.092278e1; % resistance ohms/km     % Earth resistivity = 56ohm.m, f=100 kHz
    X_atp=1.054232e3; % reactance ohms/km 
    B_atp=4.329376e-3; % reactance s/km 
end

if ~exist('R0','var')% If there is no R0 variable
    R0=R_atp/1000; % series resistance (Ohm/m)
    L0=X_atp/(1000*w); % series inductance (H/m)
    C0=B_atp/(w*1000); % shunt capacitance (F/m)
end
if ~exist('ZL','var')% If there is no ZL variable
    ZL=(L0/C0)^0.5; % resistive load connected at the line terminal = (L0/C0)^0.5 
end

energ = [t_exp_(:,1)  v_exp(:,1)]; % time (t must start at zero) and energizing voltage
v_exp = v_exp(:,2:end);% Experimental measurements voltages at values ​​of "t_exp_". "v_exp" can have NaN values
t_exp_ = t_exp_(:,2:end);% Sampling times of the measured voltages v_exp
tt=(L0*C0)^0.5*distan_exp;%propagation times from the voltage source to the distances given by vector "d"
Ne=length(distan_exp);
for k=1:Ne
    t_exp_(:,k)=t_exp_(:,k)+tt(k);%Initial times are updated for each measurement
end

tm=tt+4.9e-6;%maximum time of available experimental data. 4.9e-6 is the maximum time span recorded for each voltage (see vector t_exp_)
N_exp=size(energ,1);%Number of experimental data
figure(1)
plot(1e6*energ(1:end,1),1e-3*energ(1:end,2),'Color',[0 0.4470 0.7410],'Marker','*','MarkerIndices',1:10:N_exp)% experimental data - energizing voltage
hold on
for k=1:Ne
    t_temp=t_sim_dist(t_sim_dist_1<tm(k));%Simulated data (model VDLM) are plotted over a time interval similar to the experimental data
    N_temp=length(t_temp);
    plot(1e6*t_temp,1e-3*v_sim_dist_1(k,1:N_temp),'Color',[0.4660 0.6740 0.1880],'Marker','square','MarkerIndices',1:100:N_temp)%simulated data - VDLM

    t_temp=t_sim_dist(t_sim_dist_2<tm(k));%Simulated data (model AVDLM) are plotted over a time interval similar to the experimental data
    N_temp=length(t_temp);
    plot(1e6*t_temp,1e-3*v_sim_dist_2(k,1:N_temp),'Color',[0.8500 0.3250 0.0980],'Marker','o','MarkerIndices',50:100:N_temp)%simulated data - AVDLM

    plot(1e6*t_exp_(:,k),1e-3*v_exp(:,k),'Color',[0 0.4470 0.7410],'Marker','*','MarkerIndices',1:10:N_exp)%experimental data
end
legend('Measured','VDLM','AVDLM')
grid on
xlabel('Time [\mus]')
ylabel('Voltage [kV]')