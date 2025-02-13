clear
% Lire les données du csv
res = read_csv();

% Récupérer les vecteur
t1 = res.Temp_0__C_;
t2 = res.Temp_1__C_;
t3 = res.Temp_2__C_;
temps =  res.Temps_s_;
u = res.Echelon_V_;

hold on
plot(temps, t1);
plot(temps, t2);
plot(temps, t3);