data = evalin('base','out');
info = data.filtrePWM;
temps = info.time(:,1);

%pwm = info.signals(1).values(:,1);
filtre = info.signals.values(:,1);
subplot(3,1,[1 2]);
hold on
%plot(temps, pwm,'LineWidth',1)
plot(temps, filtre,'--', 'LineWidth',3)
xlabel('Temps [s]')
ylabel('Tension [V]')
%figure ;
%hold on;
axeX = moku.data(86800:87800, 1)-12.4;
axeY = moku.data(86800:87800, 3) -1;
plot(axeX, axeY,'LineWidth',3, 'Color','r');
ylim([0.25 3.2])
legend('filtre simulé','filtre mesuré','Location','southeast')

