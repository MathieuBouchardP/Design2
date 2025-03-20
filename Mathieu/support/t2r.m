function RT = t2r(T)
% Définition des coefficients
A = -14.65719769;
B = 4798.84200000;
C = -115334.00000000;
D = -3730535.00000000;
R25 = 10000; % Exemple de R25 (à adapter selon ton besoin)

% Calcul de la résistance
RT = R25 * exp(A + B./(T+273.15) + C./(T+273.15).^2 + D./(T+273.15).^3);
end