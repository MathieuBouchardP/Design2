data = evalin('base','PertubationSIM');
time = data.time(:,1);

T3 = data.signals.values(:,3);
T2 = data.signals.values(:,2);
T1 = data.signals.values(:,1);

for i = 1:size(T1,1)
    T1(i) = T1(i) + 25;
end

for i = 1: size(T2,1)
    T2(i) = T2(i) + 25;
end

for i = 1: size(T3,1)
    T3(i) = T3(i) + 25;
end

hold on
plot(time,T1,'-','LineWidth', 3);
plot(time,T2,'--','LineWidth', 3);
plot(time,T3,'-.','LineWidth', 3);
ylabel('Température [^oC]');
xlabel('Temps [s]')
legend('T1','T2','T3','Location','southeast')