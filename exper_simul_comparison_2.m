function [t,v_sim_dist] = exper_simul_comparison_2(energ, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit, K_C, K_G)
%{ 
    Author: Luis C. Timaná E.
    exper_simul_comparison: Calculates the error, given by the sum of the squares of the differences between the simulated and measured data
    Simulation: energization (at t=0) of a line with a resistive load - N line segments in series
%}

% energ: Experimental data. The first column corresponds to the sampling time. The second column corresponds to the energization voltage                               
% v_exp: measured voltages
% distan_exp: (horizontal vector) Distances on the line (from the source) where the measurements were taken 
% distant dimension = Number of columns of matrix exper_data - 1
% delta_t: simulation step time

% w: parameter frequency (rad/s) <=> 60 Hz
% R0: series resistance (Ohm/m)
% L0: series inductance (H/m)
% C0: shunt capacitance (F/m)
% N: Number of line segments
% L: total line length (m)
% ZL: resistive load connected at the line terminal
% model     1: Series Resistance Lumped at Each End  (Pereira and Tavares, 2021)
          % 2: Corrected model proposed by Pereira and Tavares, 2021
% V_crit: Corona inception voltage (V) (value to be optimized)
% K_C: Constant associated with shunt capacitance (value to be optimized)
% K_G: Constant associated with shunt resistance (value to be optimized)

r = isnan(energ(:,1));
energ(r,:)=[];%All non-numeric values ​​are removed from matriz energ

Ls=L/N; % length of line segments (m)
t_exp=energ(:,1); % time of the experimental voltage to energize the line
t=0:delta_t:Tmax;%simulation time
vf_exp=energ(:,2);% Experimental energization voltage (V)
vf = interp1(t_exp, vf_exp, t,'linear','extrap');% Energization voltage (V)

line = BergeronModel([R0 , L0 , C0 , Ls , V_crit , K_C , K_G]);% Model proposed by Timaná, 2024

% The "historicalCurrents" function depends on the past voltages (t-tao) at the ends of each segment
tao=ones(1,N)*line.taoLine(0);% travel time for each segment, only for the first iteration
% Initially "tao" and "Z" are equal for all line segments

Nt=length(t);
v_k=zeros(N+1,Nt);% voltages at the terminals of the line segments
i_k_kp1=zeros(N,Nt);% currents between segments: i_k,k+1    k=1,...,N   -->
i_kp1_k=zeros(N,Nt);% currents between segments: i_k+1,k    k=1,...,N   <--
I_k_ab=zeros(N,2);% historical currents => first row: I1a, I2a, ..., INa
                    %                    => second row: I1b, I2b, ..., INb

Ne=length(distan_exp);%number of points on the line where the tensions were measured
distan_sim = 0:Ls:L;%distances from the voltage source, at which the voltages are simulated
v_sim_dist=zeros(Ne,Nt);%Simulated tensions at the points where the measured tensions were taken at user-defined times

for n=1:Nt
    
    % Calculation of Norton impedances: Z1, Z2, Z3, ..., ZN+1
    if n==1
        Z=ones(1,N+1)*line.ZInitialLine();
    else
        for k=1:N+1
            Z(k)=line.CalculateZLine(v_k(k,n-1));
        end
    end

    % Calculation of historical current sources
    for k=1:N
        t_tao_k = t(n)-tao(k);
        t_taok_deltat = t(n)-tao(k)-delta_t;

        if t_tao_k<0
            i_k_kp1_t_taok = 0;
            v_k_t_taok = 0;

            i_kp1_k_t_taok = 0;
            v_kp1_t_taok = 0;
        else 
            i_k_kp1_t_taok = interp1( t(1:n-1) , i_k_kp1(k,1:n-1) , t_tao_k );
            v_k_t_taok = interp1( t(1:n-1) , v_k(k,1:n-1) , t_tao_k );

            i_kp1_k_t_taok = interp1( t(1:n-1) , i_kp1_k(k,1:n-1) , t_tao_k );
            v_kp1_t_taok = interp1( t(1:n-1) , v_k(k+1,1:n-1) , t_tao_k );
        end

        if t_taok_deltat<0
            v_k_t_taok_deltat = 0;

            v_kp1_t_taok_deltat = 0;
        else
            v_k_t_taok_deltat = interp1( t(1:n-1) , v_k(k,1:n-1) , t_taok_deltat );

            v_kp1_t_taok_deltat = interp1( t(1:n-1) , v_k(k+1,1:n-1) , t_taok_deltat );
        end

        if n~=1
            v_k_t_deltat = v_k(k,n-1);
            v_kp1_t_deltat = v_k(k+1,n-1);
        else
            v_k_t_deltat = 0;
            v_kp1_t_deltat = 0;
        end

        [I_k_ab(k,1),I_k_ab(k,2)]=line.historicalCurrents(  i_k_kp1_t_taok , i_kp1_k_t_taok  ,  v_k_t_taok , v_kp1_t_taok ,  v_k_t_deltat , v_kp1_t_deltat , v_k_t_taok_deltat , v_kp1_t_taok_deltat);
    end

    v_k(1,n)=vf(n);
    for k=2:N
        v_k(k,n) = -(Z(k)/2)*( I_k_ab(k-1,2) + I_k_ab(k,1) );
    end
    v_k(N+1,n) = -((Z(N+1)*ZL)/(Z(N+1)+ZL))*I_k_ab(N,2);

    for k=1:N-1
        i_kp1_k(k,n) = v_k(k+1,n)/Z(k+1) + I_k_ab(k,2);
    end
    i_kp1_k(N,n)=I_k_ab(N,2)*(Z(N+1)/(Z(N+1)+ZL));

    i_k_kp1(1,n) = v_k(1,n)/Z(1) + I_k_ab(1,1);
    for k=2:N
        k2=k-1;
        i_k_kp1(k,n) = -i_kp1_k(k2,n);
    end

    v_sim_dist(:,n) = interp1(distan_sim, v_k(:,n), distan_exp');%The vector is transposed to convert it into a column

    for k=1:N % tao values ​​are updated for the next iteration
        tao(k) = line.taoLineSection( v_k(k,n) , v_k(k+1,n) );
    end
 
end



