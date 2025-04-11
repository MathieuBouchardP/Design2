data = evalin('base','out');
info = data.filtrePWM;
temps = info.time(:,1);

pwm = info.signals(1).values(:,1);
filtre = info.signals(2).values(:,1);

hold on
%plot(temps, pwm,'LineWidth',1)
plot(temps, pwm, 'LineWidth',3)
xlabel('Temps [s]')
ylabel('Tension [V]')