addpath("support")

res = read_csv("data_logs.csv");

tension = res.DeltaV_V_(:, 1);
courant = res.Courant_V_(:, 1);
temps = res.Temps_s_(:, 1);
echelon = res.Echelon_V_(:,1);
k_c = 1.34375;
courant = courant/ k_c;

P = times(courant, tension);
hold on;
plot(temps, P);
plot(temps, echelon);
