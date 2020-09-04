% =========================================================================
% Title       : nML detector for massive-MU-MIMO-OFDM using BOX1OFDM
% -------------------------------------------------------------------------
% Date: 06/29/2019
% -------------------------------------------------------------------------
% Author: Seyed Hadi Mirfarshbafan
% =========================================================================
fullpath = mfilename('fullpath');
my_path = strrep(fullpath, mfilename(), '');
t_n = -4;
t_p = 4;
N = 160;
t = linspace(-5,5, N);
tp = linspace(-5,5, N/4);

asymp = t;
asymp(t >= 0) = 0;
asymp(t < 0) = -t(t < 0);

omega = exp(-t.^2/2)./(sqrt(2*pi)*normcdf(t));
omegatilde = exp(-tp.^2/2)./(sqrt(2*pi)*normcdf(tp));
omegatilde(tp > t_p) = 0;
omegatilde(tp < t_n) = -tp(tp < t_n);

figure
%plot(t, omega, t, asymp, tp, omegatilde, 'x', 'LineWidth', 1.5);
plot(t, omega, 'r', 'LineWidth', 2);
hold on;
plot(t, asymp, 'b-', 'LineWidth', 2);
hold on;
plot(tp, omegatilde, 'k:', 'LineWidth', 4);
hold off;
xlabel('t');
ylabel('$\omega(t)$, $\tilde{\omega}(t)$, asymptotes', 'Interpreter','latex');
grid on;
ylim([0 5])
legend({'exact value, $\omega(t)$','asymptotes', 'approx. value, $\tilde{\omega}(t)$'}, 'Interpreter','latex', 'FontSize',14);
% h = annotation('textarrow');
% set(h,'Parent',gca)
% set(h,'String',{'t_p'})
% set(h,'FontSize',18)
% set(h,'X',[-4 0]);
% set(h,'Y',[-3 1]);



matlab2tikz('standalone',true,'showInfo', false,[my_path '..\results\figures\tikz\omegafunction.tex']);