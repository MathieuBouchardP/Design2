function [res] = identify(y, u, t, ordre, nbr_zero, comparer)
% La fonction retourne un vecteur donc la première rangé est le numérateur
% (en puissance décroissante de s), la seconde le dénominateur et la
% dernière le retard. On a donc dequoi comme ça (sans les s explicitement)
%       a*s^2       b*s       c
%       d*s^2       3*s       5
%       retard      0         0
%
    % Si le pas d'échantillonage n'est pas uniforme, il faut le rendre
    % uniforme et adapter les données en conséquence (interpolation)
    if t(2) - t(1) == t(4) - t(3)         
        t_uniform = linspace(min(t), max(t), length(t)); % uniformiser le pas
        y = interp1(t, y, t_uniform(1, :), 'linear'); % interpoler les données
        y = y(:);                                     % s'assurer qu'on a le bon format de matrice
        sample_rate = t_uniform(2) - t_uniform(1);             % Pas d'échantillonage 
    else
        sample_rate = t(2) - t(1);
    end
    data = iddata(y, u, sample_rate);           % Initialiser le data d'identification
    iodelay = NaN;                              % Mettre le delay a NaN pour que tfest l'identifie
    model = tfest(data, ordre, nbr_zero, iodelay); % Identifier automatiquement la ft

    [num, den] = tfdata(model, 'v');  % récupérer les données du modèle
    retard = model.IODelay; % aller chercher la valeur du retard
    res = [];
    res(1,:) = num;
    res(2,:) = den;
    res(3,1) = retard;
    if comparer
       compare(data, model);
    end

    num_str = poly_to_string(num);
    den_str = poly_to_string(den);

    if retard > 0
        fprintf('G(s) = (%s) * exp(-%.4f s) / (%s)\n', num_str, retard, den_str);
    else
        fprintf('G(s) = (%s) / (%s)\n', num_str, den_str);
    end
end

function str = poly_to_string(poly)
    % Convertit un polynôme en chaîne de caractères sous forme lisible
    order = length(poly) - 1;
    terms = arrayfun(@(c, p) sprintf('%.4f s^%d', c, order - p), poly, 0:order, 'UniformOutput', false);
    
    % Supprimer les termes nuls
    terms = terms(poly ~= 0);
    
    % Gérer le cas où il ne reste plus rien (ex: numérateur nul)
    if isempty(terms)
        str = '0';
    else
        str = strjoin(terms, ' + ');
    end
end