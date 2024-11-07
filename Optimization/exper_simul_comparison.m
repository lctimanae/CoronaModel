function error=exper_simul_comparison(energ_1, energ_2, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, x ,r , h , t_exp_1 , v_exp_1,  t_exp_2 , v_exp_2, Model)

addpath('..');%Add the previous folder to use the simulation functions

V_crit=x(1);
K_C=x(2)*(r/(2*h))^0.5*1e-11; %Constant associated with shunt capacitance (initial condiciton - value to be optimized)
K_G=x(3)*(r/(2*h))^0.5*1e-11; %Constant associated with shunt resistance (initial condiciton - value to be optimized)

switch Model
    case 1 %(VDLM)
        [t_sim_dist_1,v_sim_dist_1] = exper_simul_comparison_ini(energ_1, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit, K_C, K_G);
        [t_sim_dist_2,v_sim_dist_2] = exper_simul_comparison_ini(energ_2, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit, K_C, K_G);
    case 2 %(AVDLM)
        [t_sim_dist_1,v_sim_dist_1] = exper_simul_comparison_2(energ_1, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit, K_C, K_G);
        [t_sim_dist_2,v_sim_dist_2] = exper_simul_comparison_2(energ_2, distan_exp, delta_t, Tmax, R0, L0, C0, N, L, ZL, V_crit, K_C, K_G);
end

error_1 = error_calculation(t_sim_dist_1 , v_sim_dist_1 ,  t_exp_1 , v_exp_1);
error_2 = error_calculation(t_sim_dist_2 , v_sim_dist_2 ,  t_exp_2 , v_exp_2);
error = error_1 + error_2;