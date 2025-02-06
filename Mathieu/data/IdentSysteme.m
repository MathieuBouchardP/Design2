% CHARGEMENT DES DONNÉES
%
% Chargement et préparation des données
load data_t123.mat 
N = size(datalog,1); % Nombre de lignes de datalog
T_s = 2; % Période d'échantillonnage en secondes

% Création du vecteur temps
t = (0:N-1)' * T_s; 
dataid = [t table2array(datalog)];

% Extraction des colonnes
Ncol = size(dataid,2);
t = dataid(:,1);
u = dataid(:,2);
y = dataid(:,Ncol);

% Visualisation des données brutes
figure(1);
subplot(211);
plot(t,u,'o');
ylabel('u'); xlabel('t');
title('Données brutes');
subplot(212);
plot(t,y,'o');
ylabel('y'); xlabel('t');

% Soustraction du point d'opération
u_id = u - u(1);
y_id = y - y(1);

% Visualisation des données transformées
figure(2);
subplot(211);
plot(t,u_id,'o');
ylabel('u_{id}'); xlabel('t');
title('Données d''identification');
subplot(212);
plot(t,y_id,'o');
ylabel('y_{id}'); xlabel('t');

% CRÉATION DES DONNÉES D'IDENTIFICATION
id_data = iddata(y_id, u_id, T_s);

% IDENTIFICATION AUTOMATIQUE DE LA FONCTION DE TRANSFERT OPTIMALE
% Estimation d'une fonction de transfert de type (num, den)
ordres = [1 2; 2 2; 2 3; 3 3]; % Essai avec plusieurs ordres
meilleur_modele = [];
meilleur_critere = Inf;

for i = 1:size(ordres,1)
    modele_test = tfest(id_data, ordres(i,1), ordres(i,2));
    critere = modele_test.Report.Fit.FitPercent;
    
    if critere < meilleur_critere
        meilleur_critere = critere;
        meilleur_modele = modele_test;
    end
end

% Affichage du meilleur modèle
disp('Meilleure fonction de transfert identifiée:');
display(meilleur_modele);

% Visualisation de la réponse du modèle
figure(3);
compare(id_data, meilleur_modele);
