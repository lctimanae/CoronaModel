function error = error_calculation(t , v_sim_dist ,  t_exp_ , v_exp)

Ne=size(v_exp,2);
error=0;
for k=1:Ne
    v_exp2 = v_exp(:,k);
    v_exp2(isnan(v_exp2)) = [];%Removing non-numeric data

    t_exp_2 = t_exp_(:,k);
    t_exp_2(isnan(t_exp_2)) = [];%Removing non-numeric data

    t_max = max(t_exp_(:,k));% Maximum experimental time recorded for each distance
    t_min = min(t_exp_(:,k));% Minimum experimental time recorded for each distance
    t2 = t(t<=t_max & t>=t_min);% Time to interpolate must be greater than the minimum experimental time and less than the maximum experimental time
    N_min=find(t>=t_min); N_min=N_min(1);%Find the minimum and maximum simulation range to compare with the experimental results
    N_max=find(t<=t_max); N_max=N_max(end);

    v_exp_inter = interp1( t_exp_2 ,  v_exp2 , t2' ,'linear','extrap' );% Interpolation of experimental data
    per = abs(  (v_sim_dist(k,N_min:N_max)' - v_exp_inter)./v_exp_inter  );% Percentage error for each distance
    per = per(10:end);%The first ten values ​​are not taken into account because there are usually significant differences in the first moments.
    error = error + mean( per );
end
error = error/Ne;
