data = evalin('base','PertubationSIM');
time = data.time(:,1);

T3 = data.signals.values(:,3);
T2 = data.signals.values(:,2);
T1 = data.signals.values(:,1);

for i = 1:size(T1,1)
    T1(i) = T1(i) + 24;
end

for i = 1: size(T2,1)
    T2(i) = T2(i) + 24;
end

for i = 1: size(T3,1)
    T3(i) = T3(i) + 24;
end

hold on
subplot(3,1,[1 2])
plot(time,T1,'-','LineWidth', 3);
plot(time,T2,'--','LineWidth', 3);
plot(time,T3,'-.','LineWidth', 3);
ylabel('Température [^oC]');
xlabel('Temps [s]')

plot(perturbationsOpen(2:end, 1), perturbationsOpen(2:end, 4),'o', 'LineWidth', 1.5);
plot(perturbationsOpen(2:end, 1), perturbationsOpen(2:end, 5),'s','LineWidth', 1.5);
plot(perturbationsOpen(2:end, 1), perturbationsOpen(2:end, 6),'^','LineWidth', 1.5);
legend('T1 (simulation)','T2 (simulation)','T3 (simulation)','T1 (prototype)','T2 (prototype)','T3 (prototype)','Location','southeast')
subplot(3,1,3)
plot()