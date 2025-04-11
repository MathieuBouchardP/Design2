addpath("simuTempCourbe");
% ignorer la permiere rangée
T = csvread('simulation_data_1.2W.csv',1,0);
% première colonnes c'est le temps
temps = T(:,1);
T3_12 = T(:,4);

A = csvread('simulation_data_0.8W.csv',1,0);
T3_08 = A(:,4);

B = csvread('simulation_data_0.4W.csv',1,0);
T3_04 = B(:,4);

D = csvread('simulation_data_0.2W.csv',1,0);
T3_02 = D(:,4);

C = csvread('simulation_data-0.2W.csv',1,0);
T3_002 = C(:,4);

F = csvread('simulation_data-0.4W.csv',1,0);
T3_004 = F(:,4);

G = csvread('simulation_data-0.8W.csv',1,0);
T3_008 = G(:,4);

H = csvread('simulation_data-1.2W.csv',1,0);
T3_012 = H(:,4);
%% fabrication des figures
plot(temps,T3_012,'-',LineWidth=3);
%title("Réponses à l'échelon pour T2")
hold on
plot(temps,T3_008,'--',LineWidth=3);
plot(temps,T3_004,'-.',LineWidth=3);
plot(temps,T3_002, ':',LineWidth=3);
ylim([9 30]);
xlim([0 1050]);
ylabel('Température [{\circ}C]')
xlabel('Temps [s]')
legend('T3 à -1.2W','T3 à -0.8W','T3 à -0.4W','T3 à -0.2W')
colororder("dye")
fontsize(gcf,scale=1.3)
hold off
