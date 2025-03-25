clear
% Lire les données du csv
addpath("support")
res = read_csv("data_logs");

% Récupérer les vecteur
t3 = res.Temp_0__C_;
t2 = res.Temp_1__C_;
t1 = res.Temp_2__C_;
t1 = t1(10:end, 1);
temps =  res.Temps_s_;
temps = temps(10:end, 1);
u = res.Echelon_V_;

hold on
plot(temps, t1);
%plot(temps, t2);
%plot(temps, t3);


n_fit = 10;                     % Nombre de points pour la régression
rho = 2582.33;                 % Densité de l'aluminium (kg/m³)
cp = 897;                      % Capacité thermique massique en J/(kg·K)
z = 5*10^-3;
x = 14.84*10^-3;
y = 15.8*10^-3;

function [slope, power_total] = est_P(t, T, n_fit, rho, cp, A, h)
% estimatePowerTotal : Estime la puissance totale appliquée sur un petit élément de surface.
%
% Syntaxe:
%   [slope, power_total] = estimatePowerTotal(t, T, n_fit, rho, cp, A, h)
%
% Entrées:
%   t    - Vecteur temps (en s)
%   T    - Vecteur température (en °C)
%   n_fit- Nombre de points à utiliser pour la régression linéaire (par défaut 5)
%   rho  - Masse volumique (kg/m^3)
%   cp   - Capacité thermique massique (J/(kg·K))
%   A    - Surface de l'élément sur lequel la chaleur est appliquée (en m^2)
%   h    - Épaisseur de l'élément considéré (en m)
%
% Sorties:
%   slope       - Estimation de la dérivée de T par rapport au temps (°C/s)
%   power_total - Puissance totale appliquée sur l'élément (W)
%
% Le modèle utilisé est :
%   puissance volumique (W/m^3) = rho * cp * (dT/dt)
%   volume de l'élément = A * h
%   puissance totale = puissance volumique * volume
%
% La fonction affiche également les résultats et trace la courbe de température
% avec le fit linéaire sur les n_fit premiers points.

    % Vérification du nombre de points disponibles
    if length(t) < n_fit
        error('Le vecteur t doit contenir au moins %d points.', n_fit);
    end

    % Régression linéaire sur les n_fit premiers points pour estimer dT/dt
    coeff = polyfit(t(1:n_fit), T(1:n_fit), 1);
    slope = coeff(1); % dT/dt en °C/s

    % Calcul de la puissance volumique (W/m^3)
    power_density = rho * cp * slope;

    % Calcul du volume de l'élément
    volume = A * h;
    
    % Calcul de la puissance totale appliquée sur l'élément (W)
    power_total = power_density * volume;

    % Affichage des résultats
    fprintf('Dérivée de T au début (dT/dt) = %.4f °C/s\n', slope);
    fprintf('Puissance volumique (rho*cp*dT/dt) = %.4f W/m^3\n', power_density);
    fprintf('Puissance totale sur l''élément (A*h*(rho*cp*dT/dt)) = %.4f W\n', power_total);

    % Visualisation : tracé de la courbe de température et du fit linéaire
    figure;
    hold on;
    grid on;
    % Tracé des données expérimentales
    plot(t, T, 'b.-', 'DisplayName', 'Données expérimentales');
    % Tracé du fit linéaire sur les n_fit premiers points
    t_fit = t(1:n_fit);
    T_fit = polyval(coeff, t_fit);
    plot(t_fit, T_fit, 'r-', 'LineWidth', 2, 'DisplayName', 'Fit linéaire');
    xlabel('Temps (s)');
    ylabel('Température (°C)');
    title('Estimation de dT/dt et de la puissance totale sur un élément');
    legend('Location','Best');
    hold off;
end
[dTdt, p_surf] = est_P(temps, t1, n_fit, rho, cp, x*y, z);

%surface_TEC = y*x;
%power_total = p_surf * surface_TEC;
%fprinf(power_total)