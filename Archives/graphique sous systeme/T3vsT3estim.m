data = evalin('base','estimationT3');
time = data.time(: ,1);
T3 = data.signals(1).values(:,1);
T3estim = data.signals(2).values(:,1);

for i = 1:size(T3,1)
    T3(i) = T3(i) + 25;
end

for i = 1: size(T3estim)
    T3estim(i) = T3estim(i) + 25;
end

disp(T3);
subplot(3,1,[1 2]);
 plot(time, T3, '-','LineWidth', 3);
 hold on
 plot(time, T3estim,'--','LineWidth', 3)
    legend('T3 mesurée','T3 estimée', 'Location', 'southeast');
    ylim([24,31])
    xlim([0,650])
    ylabel('Température [^oC]');
    colororder("gem")
    fontsize(gcf,scale=1.3)

subplot(3,1,3);
plot(time,T3-T3estim);
xlabel('Temps [s]');
ylabel('Différence [^oC]');
ylim([0,0.25])
xlim([0, 650])
fontsize(gcf,scale=1.3)

