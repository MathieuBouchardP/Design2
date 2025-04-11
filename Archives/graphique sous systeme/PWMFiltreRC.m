data = evalin('base','out');
info = data.PWMsortie;
temps = info(:,1);

%pwm = info.signals(1).values(:,1);
pwm = info(:,2);
subplot(2,1,1);

% 1. Tronquer axeY à 1000 échantillons
%axeY = axeY(1:1000);
 
% 2. Tronquer filtre à 250 échantillons
%pwm = pwm(1:250);
%temps = temps(1:250);



% 
hold on
%plot(temps, pwm,'LineWidth',1)
plot(temps, pwm,'--', 'LineWidth',3)
ylabel('Tension [V]')
xlabel('Temps [s]')
legend('PWM simulé','Location','southeast')
axeX = moku.data(86800:86900, 1)-12.4;
axeY = moku.data(86800:86900, 2) -1;

subplot(2,1,2);
plot(axeX, axeY,'LineWidth',2, 'Color','r');
xlim([0.00160799,5e-03]);
ylabel('Tension [V]')
xlabel('Temps [s]')
legend('PWM mesuré','Location','southeast')