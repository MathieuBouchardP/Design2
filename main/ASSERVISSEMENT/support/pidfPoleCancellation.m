function [Gc, T_i, T_d, T_f, Kc] = pidfPoleCancellation(Gp, push_gain)
    if nargin < 2
        push_gain = 1;
    end
    % Valeur par défaut du damping ratio (marge de phase) :
    zeta = 0.707*1.01;
    % Facteur pour le choix de T_f : T_f = T_d / factor (par défaut factor = 4)
    factor = 2;
    
    % Récupération des pôles du procédé
    p = pole(Gp);
    if length(p) ~= 2
        error('Le procédé doit être du 2e ordre avec exactement 2 pôles réels.');
    end
    
    % On s'assure que les pôles sont réels et négatifs
    if any(imag(p) ~= 0) || any(real(p) >= 0)
        error('Les pôles du procédé doivent être réels et dans le plan de gauche.');
    end
    
    % Tri en ordre croissant (p1 est le pôle lent, p2 le pôle rapide)
    p = sort(real(p), 'ascend');
    % Notons p1 = -a et p2 = -b (avec a < b)
    a = -p(1);   % pôle rapide
    b = -p(2);   % pôle lent
    
    % Calcul de Tᵢ et T_d
    T_i = 1 / b *0.89;
    %T_i = 1 / b;
    T_d = 1 / a *1.09;
    
    % Choix de T_f : on prend T_f = T_d / factor (par défaut factor = 4)
    T_f = T_d / factor;
    
    % Calcul du gain statique du procédé k_p :
    % Gp(0) = k_p/(a*b)  →  k_p = dcgain(Gp) * a * b.
    k_p = dcgain(Gp) * a * b;
    
    % Calcul du gain du régulateur Kc selon :
    % Kc = b/(4*zeta^2*k_p*T_f)
    Kc = b / (4 * zeta^2 * k_p * T_f)*push_gain;
    
    % Construction du PIDF en forme "classique" (série)
    % Gc(s) = Kc*(1+Tᵢ s)*(1+T_d s) / [Tᵢ s*(1+T_f s)]
    numGc = Kc * conv([T_i, 1], [T_d, 1]); % (T_i s + 1)(T_d s + 1)
    denGc = conv([T_i, 0], [T_f, 1]);       % T_i s*(T_f s + 1)
    %%%%%%%%%%%%%%% à enlever
    %denGc = [T_i, 0];
    %%%%%%%%%%%%%%%%%%%%%

    Gc_raw = tf(numGc, denGc);
    Gc = minreal(Gc_raw);
    
    % Affichage des paramètres dérivés (pour vérification)
    fprintf('Paramètres dérivés :\n');
    fprintf('  Tᵢ = %.4g,  T_d = %.4g,  T_f = %.4g\n', T_i, T_d, T_f);
    %fprintf('  Gain du procédé, k_p = %.4g\n', k_p);
    fprintf('  Gain du régulateur, Kc = %.4g\n', Kc);
end