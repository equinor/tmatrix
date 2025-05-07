% Plot figure 4b in "The acoustic signature of fluid flow in complex porous media"
% (Jakobsen et al. 2003)



function ans = plot_velocity_frequency(out)



figure(1)
subplot(2,1,1);
semilogx(out.frequency, out.Vp,'-k');
hold on;
ylabel('Vp (m/s)');
hold on;
subplot(2,1,2);
semilogx(out.frequency, out.invQp,'-k');
ylabel('1/Qp');xlabel('Frequency (Hz)');
hold on;

figure(2)
subplot(2,1,1);
semilogx(out.frequency, out.Vs,'-k');
hold on;
ylabel('Vs (m/s)');
hold on;
subplot(2,1,2);
semilogx(out.frequency, out.invQs,'-k');
ylabel('1/Qs');xlabel('Frequency (Hz)');
hold on;

