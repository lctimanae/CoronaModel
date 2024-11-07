function graphComp(tt, energ, t_sim_dist, v_sim_dist, t_exp_, v_exp)

Ne=length(tt);
tm=tt+4.9e-6;%maximum time of available experimental data. 4.9e-6 is the maximum time span recorded for each voltage (see vector t_exp_)
figure
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