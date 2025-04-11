%% faire les courbes pour aller dans le chapitre 9
addpath("simuTempCourbe");
% ignorer la permiere rangée
T = csvread('simulation_data_1.2W.csv',1,0);
% première colonnes c'est le temps
temps = T(:,1);
%prendre les données du csv
T1_12 = T(:,2);

A = csvread('simulation_data_0.8W.csv',1,0);
T1_08 = A(:,2);

B = csvread('simulation_data_0.4W.csv',1,0);
T1_04 = B(:,2);

D = csvread('simulation_data_0.2W.csv',1,0);
T1_02 = D(:,2);

C = csvread('simulation_data-0.2W.csv',1,0);
T1_002 = C(:,2);

F = csvread('simulation_data-0.4W.csv',1,0);
T1_004 = F(:,2);

G = csvread('simulation_data-0.8W.csv',1,0);
T1_008 = G(:,2);

H = csvread('simulation_data-1.2W.csv',1,0);
T1_012 = H(:,2);
% set les couleurs

%% fabrication des figures
%plot(temps,T1_12,'-',LineWidth=3);
%title("Réponses à l'échelon pour T1")
%hold on
%plot(temps,T1_08,'--',LineWidth=3);
%plot(temps,T1_04,'-.',LineWidth=3);
%plot(temps, T1_02, ':',LineWidth=3);
plot(temps,T1_012,'-',LineWidth=3);
hold on
plot(temps,T1_008,'--',LineWidth=3);
plot(temps,T1_004,'-.',LineWidth=3);
plot(temps, T1_002, ':',LineWidth=3);
ylim([5 28]);
xlim([0 1050]);
ylabel('Température [{\circ}C]')
xlabel('Temps [s]')
%legend('T1 à 1.2W','T1 à 0.8W','T1 à 0.4W','T1 à 0.2W''Location','northeast')
legend('T1 à -0.2W','T1 à -0.4W','T1 à -0.8W','T1 à -1.2W','Location','northeast')
colororder("gem")
fontsize(gcf,scale=1.3)
hold off
