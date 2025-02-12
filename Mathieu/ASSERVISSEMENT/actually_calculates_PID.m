function [gc] = actually_calculates_PID(gp)
%% Calculateur de régulateur FONCTIONNE JUSTE SI YA PAS DE ZÉROS
    num = gp(1, :);
    den = gp(2, :);
    retard = gp(3, 1);

    Gp = tf(num, den, 'InputDelay', retard); % Retard de 2 secondes
    [C, ~] = pidtune(Gp, 'pidf'); % Génère directement le PID

% Extraction des paramètres sous la forme idéale
    Kc = C.Kp;
    Ti = Kc / C.Ki;
    Td = C.Kd / Kc;
    Tf = C.Tf;
    if Tf == 0
        Tf = 0.01; % Petite valeur pour éviter les erreurs numériques
    end
    % convertir pour que simulink comprenne dequoi parce qu'yé trop épais
    % pour qu'on mettre plusieurs []
    num = Kc*conv([Ti 1], [Td 1]);
    den = Ti*[Tf, 1, 0];

    gc = zeros(2, max(size(num, 2), size(den, 2)));
    gc(1, :) = num;
    gc(2, :) = den;
    %Gc_ideal = pid(Kc, Kc/Ti, Kc*Td, Tf, 'IFormula', 'Ideal');
    %disp(Gc_ideal);
end