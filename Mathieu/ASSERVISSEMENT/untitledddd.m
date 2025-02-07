load data_t123.mat
test = table2array(datalog);

cut = 21;
temps = test(cut:end, 1) - test(cut, 1);
t1 = test(cut:end, 2) - test(cut, 2);

t2 = test(cut:end, 3) - test(cut, 3);
t3 = test(cut:end, 4) - test(cut, 4);

chose = test(:, 1);
N = size(t3, 1);
consigne = ones(N, 1);

%plot(temps,t1,'o');
%plot(temps,t2,'o');
%plot(temps,t3,'o');

load 1.mat
load 2.mat
load 3.mat

T_s = 0.001481707317073;
stop = round(max(temps)/T_s);
NS = numel(thermistance1);
%NS = stop;
temps_s = (0:NS-1)' * T_s;

thermistance1_0 = thermistance1(:) - thermistance1(1);
thermistance2_0 = thermistance2(:) - thermistance2(1);
thermistance3_0 = thermistance3(:) - thermistance3(1);


figure;
hold on;
plot(temps, t1, 'o-', 'DisplayName', 'T1')
plot(temps_s(:), thermistance1_0(:)) % Assure que les deux sont des vecteurs colonne

plot(temps, t2, 'o-', 'DisplayName', 'T1')
plot(temps_s(:), thermistance2_0(:)) % Assure que les deux sont des vecteurs colonne

plot(temps, t3, 'o-', 'DisplayName', 'T1')
plot(temps_s(:), thermistance3_0(:)) % Assure que les deux sont des vecteurs colonne

xlabel('Temps (s)');
ylabel('Température');
title('Évolution de la thermistance');