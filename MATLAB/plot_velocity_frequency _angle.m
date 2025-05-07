% Plot figure 4b in "The acoustic signature of fluid flow in complex porous media"
% (Jakobsen et al. 2003)



function   plot_velocity_frequency_angle( out )




figure(1)
%subplot(3,1,1);
semilogx(out.frequency, out.Vp,'-k');
hold on;
ylabel('Vp (m/s)');
hold on;
figure(2)
%subplot(3,1,2);
semilogx(out.frequency, out.Vsv,'-k');
ylabel('Vsv (m/s');xlabel('Frequency (Hz)');
hold on;
figure(3)
%subplot(3,1,3);
semilogx(out.frequency, out.Vsh ,'-k');
ylabel('Vsh (m/s');xlabel('Frequency (Hz)');
hold on;

