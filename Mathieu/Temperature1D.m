close all
clear all
clc

F(length(zoom)) = struct('cdata',[],'colormap',[]);
writerObj = VideoWriter('TemperatureDistribution1D.avi');
open(writerObj);

%% Paramètres 

TempsTotal = 500;           % Durée de la simulation

Lx = 120e-3;                % Longueur [m]
thickness = 1.5e-3;         % Épaisseur de la plaque [m]

Nx = 120;                   % Nombre d'élements en x

% Propriétes du matériau, ici on a pour l'Aluminium
k = 205;                    % Conductivité Thermique [W/m·K]
rho = 2700;                 % Densité [kg/m^3]
cp = 900;                   % Chaleur spécifique  [J/kg·K]
alpha = k / (rho * cp);     % Diffusivité Thermique [m^2/s]

h_conv = 20;                % Coeff. de convection[W/m^2·K]
h_conv = 0;                 % commenter pour enlever la convection
%% Paramètres calculés

dx = Lx / Nx;               % Pas de discrétisation en x
dy = thickness;             % épaisseur en y
dz = thickness;             % épaisseur en z

dt = dx^2/(8*alpha);        % Pas en temps [s]
                            % Choisi pour assurer la stabilité 

Nt = round(TempsTotal/dt);  % Nombre ditérations temporelles

                            % Pour la convection
aire_bouts = dy*dz;         % aire des bouts de la barre
aire_sides = dx*dz;         % aire des côtés , pour chaque éléement
aire_top = dx*dy;           % aire du dessus et du dessous

volume = dx*dy*dz;          % Volume dun élément

Temps = (0:Nt-1)*dt;        % Vecteur de temps pour la simulation
Position = (0:Nx-1).*dx;    % Vecteur de position 

%% Puissance déposée par lactuacteur

% Variable power input (example: 5 W localized input)
Pin = 1.5;                  % Puissance [W]
%Pin = 0;                   % commenter pour mettre de la puissance 

Pin_loc_x =  Lx/4;          % Localisation sur la barre [m]
Pin_loc_x =round(Pin_loc_x/dx); %Élément qui reçoit la puissance de lactuateur

P = zeros(Nx, 1);           % Matrice de puissance à ajouter à chaque elt
P(Pin_loc_x) = Pin;         % La puissance est mise sur un seul élément
                            % En pratique, la puissance déposée sera probablement 
                            % répartie sur plusieurs éléments. 
%% Conditions initiales
T_piece = 273+25;           % Tempéraure pièce [K]

T = T_piece * ones(Nx, 1);  % Temperature de tous les éléments

T_loc_x = Lx/2;             % Localisation dun dirac en Température [m]
T_loc_x = round(T_loc_x/dx);% Element qui est plus chaud

%T(T_loc_x) = 273+45;        % Un élement plus chaud

Therm_loc_x = 3*Lx/4;       %Position de la mesure de température
Therm_loc_x = round(Therm_loc_x/dx);

%% Création et positionnement de la figure

f1 = figure(1)
f1.Position = [1200 1000 1400 460];
f1.Position =[600 1000 2000 460]

%% Préallocation des vecteurs qui seront utilisés dans la boucle

energy_added = zeros(1,Nt);
energy_loss = zeros(1,Nt);

thermistance = zeros(1,Nt);
Tnew = 0*T;


%% Itération
%
% En pratique on va vouloir précalculer le maximum de coefficients pour accélérer la boucle
% On aura aussi avantage à vectoriser le code pour Matlab
% ici on privilégie la clarté 

for t = 1:Nt            % Itération sur chaque temps


    for i = 1:Nx        % Itération sur les points en x
                        % Sans la normalisation par dt/(rho cp), l'équation est en [W/m^3]

        T_new(i) = T(i);                                                                            % Température précédente 
        
        if(i == 1)      % Le 1er elt n'a pas de voisin à gauche
            T_new(i) = T_new(i) + dt / (rho * cp) * k * (T(i+1) - T(i) ) / dx^2;                    % Conduction de l'unique voisin
            T_new(i) = T_new(i) + dt / (rho * cp) * h_conv * (T_piece - T(i))*aire_bouts/volume;    % convection 1 surface du bout


        elseif(i==Nx)   % Le dernier elt n'a pas de voisin à droite
            T_new(i) = T_new(i) + dt / (rho * cp) * k * (-T(i) + T(i-1) ) / dx^2;                   % Conduction de l'unique voisin
            T_new(i) = T_new(i) + dt / (rho * cp) * h_conv * (T_piece - T(i))*aire_bouts/volume;    % convection 1 surface du bout
        else
            T_new(i) = T_new(i) + dt / (rho * cp) * k * (T(i+1) - 2*T(i) + T(i-1)) / dx^2;          % Conduction des voisins
        end

        T_new(i) = T_new(i) + dt / (rho * cp) * P(i)/volume;                                        % Puissance actuateur / perturbation
        T_new(i) = T_new(i) + dt / (rho * cp) * h_conv * (T_piece - T(i))*2*aire_sides/volume;      % convection 2 surfaces de côté
        T_new(i) = T_new(i) + dt / (rho * cp) * h_conv * (T_piece - T(i))*2*aire_top/volume;        % convection 2 surfaces dessus / dessous


    end
    
    % On a terminé la mise à jour de la température pour ce tour;
    T = T_new;

    thermistance(t) = T(Therm_loc_x);

    % Pour vérification, bilan d'énergie

    % Énergie ajoutée dans le systeme par l'actuateur
    energy_added(t) = sum(P(:)) * dt;

    % Énergie dissipée par convection
    energy_loss_sides  = h_conv * sum(T - T_piece)   * 2*aire_sides * dt;
    energy_loss_top    = h_conv * sum(T - T_piece)   * 2*aire_top * dt;
    energy_loss_bout_1 = h_conv * (T(1) - T_piece)   * aire_bouts * dt;
    energy_loss_bout_2 = h_conv * (T(end) - T_piece) * aire_bouts * dt;

    % Énergie totale dissipée
    energy_loss(t) = energy_loss_sides + energy_loss_top + energy_loss_bout_1 + energy_loss_bout_2;

    % Affichage
    if mod(t, round(Nt/1000)) == 0 || t==1
        clf
        subplot(131)
        plot(Position,T-273);
        title(['Température à t = ', num2str(t*dt), ' s'],'FontSize',16);
        xlabel('Position sur la barre [mm]','FontSize',16);
        ylabel('Température (C)','FontSize',16);
        %ylim([15 35])
        %xlim([0 120])
        grid on
        ax = gca; % Get current axes
        ax.FontSize = 16; % Set font size for tick labels
        
        subplot(132)
        plot(Temps(1:t),thermistance(1:t)-273)
        hold all
        grid on
        %ylim([20 35])

        ax = gca; % Get current axes
        ax.FontSize = 16; % Set font size for tick labels
        xlabel('Temps [s]','FontSize',16)
        ylabel('Température [C]','FontSize',16)
        title('Température à la thermistance','FontSize',16)
          
        subplot(133)
        
        hold all
        plot(Temps(1:t),energy_added(1:t))
        plot(Temps(1:t),energy_loss(1:t))
        xlabel('Time [s]','FontSize',16)
        ylabel('Énergie dnas l''itération','FontSize',16)
        legend('Energie déposée','Energie dissipée par convection','FontSize',16,'Location','southeast')
        grid on
        drawnow;

        F(i) = getframe(gcf);
        writeVideo(writerObj,F(i)); 
    end
end


close(writerObj);