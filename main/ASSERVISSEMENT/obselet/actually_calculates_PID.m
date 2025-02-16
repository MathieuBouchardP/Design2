function [gc_tf] = actually_calculates_PID(modele, gp)
%% Calculateur de régulateur FONCTIONNE JUSTE SI YA PAS DE ZÉROS
    if nargin == 2
        num = gp(1, :);
        den = gp(2, :);
        retard = gp(3, 1);
        modele = tf(num, den, 'InputDelay', retard); 
    end

    [C, ~] = pidtune(modele, 'pidf'); % Génère directement le PID
    gc_tf = tf(C); 
    
    %num = Kc*conv([Ti 1], [Td 1]);
    %den = Ti*[Tf, 1, 0];

    %gc = zeros(2, max(size(num, 2), size(den, 2)));
    %gc(1, :) = num;
    %gc(2, :) = den;
    %Gc_ideal = pid(Kc, Kc/Ti, Kc*Td, Tf, 'IFormula', 'Ideal');
    %disp(Gc_ideal);
end