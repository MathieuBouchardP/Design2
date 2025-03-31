%plotter
%file = "Essai1_0.2V.csv";
file = "Essai1_0.6V.csv";
%file = "Essai1_-0.2V.csv";
%file = "Essai1_-0.8V.csv";
%file = "Essai1_0.4V.csv";
%file = "Essai1_-0.4V.csv";
%file = "Essai1_-0.6V.csv";


addpath("support")
res = read_csv(file);
disp(res_exp)
y3 = res.T3(:,1) - res.T3(1,1);
y2 = res.T2(:,1) - res.T2(1,1);
y1 = res.T1(:,1) - res.T1(1,1);
 
t = res.Temps(:,1);
figure(1)
hold on;
plot(t, y3);
plot(t, y2);
plot(t, y1);
grid on;

len = length(t);
factor_stuff = 0.2;
cuttof_index = round((1-factor_stuff)*len);

f_t_1 = mean(y1(cuttof_index:end,1));
disp(f_t_1);
f_t_2 = mean(y2(cuttof_index:end,1));
disp(f_t_2);
f_t_3 = mean(y3(cuttof_index:end,1));
disp(f_t_3);
title(['Valeurs moyennes : t_1 = ', num2str(f_t_1), ...
       ', t_2 = ', num2str(f_t_2), ...
       ', t_3 = ', num2str(f_t_3)], 'FontSize', 14);
