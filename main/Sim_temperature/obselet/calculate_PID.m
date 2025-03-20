function [gc] = calculate_PID(gp)
%% Calculateur de régulateur FONCTIONNE JUSTE SI YA PAS DE ZÉROS
    num = gp(1, :);
        K_p = num(end); % Aller chercher gain
    den = gp(2, :);
    retard = gp(3, 1);

    %identifier l'ordre du procédé
    ordre = what_order(gp);
    if ordre == 2
        % Conversion de pôlynome à (1 + T1)(1 + T2)
        roots_denom = roots(den);
        T1 = -1/max(roots_denom);
        T2 = -1/min(roots_denom);
    end
    condition_de_andre_debiens = retard/T1;

    if condition_de_andre_debiens <= 2
        % un calcul de Ti
        Ti = (1 + 0.175*(retard / T1) + 0.3*(T2/T1)^2 + 0.2*(T2/T1)) * T1;
    elseif condition_de_andre_debiens > 2
        % un autre calcul de Ti because andre said so
        Ti = (0.65 + 0.35*(retard/T1) + 0.3*(T2/T1)^2 + 0.2*(T2/t1)) * T1;
    end
    T0 = 0; % Ici je considère qu'il y a pas de zéro, j'ai pas envie (pour l'instant dumoins)

    % Calculer w_0. C'est pas une fonction linéaire, donc on va demander a
    % matlab de fairede de la magie
    % Définition de l'équation: 
    eq = @(w0) 1.015 - ( -pi/2 + atan(w0*Ti) - atan(w0*T0) - atan(w0*T1) - atan(w0*T2) - w0*retard + pi );
    w0_initial = 3/T1*4; % Valeur initiale pour la recherche
    w0 = fsolve(eq, w0_initial);
    
    % Trouver le gain du régulateur
    K_c = Ti/K_p * ( ((T1*T2)^2*w0^6 + (T1^2 + T2^2)*w0^4 + w0^2) / ((Ti*T0)^2*w0^4 + (Ti^2 + T0^2)*w0^2 +1) );

    % Retrouver la forme souhaité et séparer en num, dénum
    num = [Ti, 1]*K_c;
    denum = [Ti, 0];
    %denom_coeffs = poly([-1/T1, -1/T2]);
    %disp(denom_coeffs);

    % mettre des coefficients du polynôme B*s^2 + C*s + D
    gc = zeros(2, max(size(num, 2), size(denum, 2)));
    gc(1, :) = num;
    gc(2, :) = denum;
end



function [ordre] = what_order(modele)
    ordre = size(modele, 1) - 1;
end