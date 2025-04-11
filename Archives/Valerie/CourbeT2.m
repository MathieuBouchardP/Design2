addpath("simuTempCourbe");
% ignorer la permiere rangée
T = csvread('simulation_data_1.2W.csv',1,0);
% première colonnes c'est le temps
temps = T(:,1);
T2_12 = T(:,3);

A = csvread('simulation_data_0.8W.csv',1,0);
T2_08 = A(:,3);

B = csvread('simulation_data_0.4W.csv',1,0);
T2_04 = B(:,3);

D = csvread('simulation_data_0.2W.csv',1,0);
T2_02 = D(:,3);

C = csvread('simulation_data_-0.2W.csv',1,0);
T2_002 = C(:,3);

F = csvread('simulation_data_-0.4W.csv',1,0);
T2_004 = F(:,3);

G = csvread('simulation_data_-0.8W.csv',1,0);
T2_008 = G(:,3);

H = csvread('simulation_data_-1.2W.csv',1,0);
T2_012 = H(:,3);
%% fabrication des figures
%% fabrication des figures
plot(temps,T2_012,'-',LineWidth=3);
%title("Réponses à l'échelon pour T2")
hold on
plot(temps,T2_008,'--',LineWidth=3);
plot(temps,T2_004,'-.',LineWidth=3);
plot(temps, T2_002, ':',LineWidth=3);
ylim([5 28]);
xlim([0 1050]);
ylabel('Température [{\circ}C]')
xlabel('Temps [s]')
legend('T2 à -1.2W','T2 à -0.8W','T2 à -0.4W','T2 à -0.2W')
colororder("reef")
fontsize(gcf,scale=1.3)
hold off
