import numpy as np
import matplotlib as plt
"""
#Définition des paramètres: 

# Durée de la simulation
temps_tot = 500


# Dimension de la plaque
Lx = 117*10**(-3)               # Longueur [m]
Nx = 117                        # Nbr éléments en x
ly = 61*10**(-3)                # Largeur [m]
ny = 61                         # Nbr éléments en y
thickness = 1.5*10**(-3)        # Épaisseur [m] (1 élément)


# Propriété du matériau
k = 205                         # Conductivité thermique [W/m*K]
rho = 2700                      # Densité [kg/m^3]
cp = 900                        # Chaleur spécifique [J/kg*K]  (Capacité calorigique massique)
alpha = k / (rho * cp)          # Diffusivité thermique [m^2/s]


# Propriété convection
h_conv = 20                     # Coefficient de convection [W/m^2*K]
h_conv = 0                      # Aucune convection


# Paramètre calculé
dx = Lx / Nx                    # Pas de discrétisation en x
dy = ly / ny                    # Pas de discrétisation en y
dz = thickness                  # Épaisseur en z

dt = dx**2/(8*alpha)            # Pas en temps [s]    --->   Attention, il faut le chosiir pour assurer la stablilité
Nt = round(temps_tot/dt)         # Nombre d'itération temporelles

air_ends = dz * dy              # Aire des bouts de la plaque, pour chaque élément
air_top_et_bot = dx * dy        # Aire du dessus ou du dessous, pour chaque élément
air_sides = dz * dx             # Aire des côtés de la plaque, pour chaque élément

volume = dx*dy*dz               # Volume d'un élément

temps = np.arange(Nt) * dt
position_x = np.arange(Nx) * dx
position_y = np.arange(ny) * dy


# Puissance déposé par l'actuateur
P_in = 1.5                      # Puissance [W]
#P_in = 0                       # Commenter pour mettre la puissance

P_in_loc_x = Lx/4               # Localisation en x de l'actuateur sur la plaque
P_in_loc_x = round(P_in_loc_x/dx)   # Élément qui reçoit la puissance
P_in_loc_y = ly/2               # Localisation en x de l'actuateur sur la plaque
P_in_loc_y = round(P_in_loc_y/dy)   # Élément qui reçoit la puissance


Power = np.zeros((Nx, ny))                # Matrice de puissance à ajouter à chaque élément
Power[P_in_loc_x, P_in_loc_y] = P_in    # Élément sur lequel la puissance est ajouté


# Condition initiales
T_pièce = 273.15 + 25           # Température de la pièce [K]
T = T_pièce * np.ones((Nx,ny))   # Temmpérature de tous les éléments

T_loc_x = Lx/2;             # Localisation en x dun dirac en Température [m]
T_loc_x = round(T_loc_x/dx) # Element qui est plus chaud
T_loc_y = ly/2;             # Localisation en y dun dirac en Température [m]
T_loc_y = round(T_loc_y/dy) # Element qui est plus chaud
T[T_loc_x, T_loc_y] = 273.15+45;        # Un élement plus chaud

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

# Paramètres
TEMPS_TOTAL = 2000  # Durée totale de la simulation [s]
Lx, ly, thickness = 117e-3, 61e-3, 1.5e-3  # Dimensions de la plaque [m]
Nx, Ny = 117, 61  # Discrétisation spatiale

# Propriétés du matériau
k, rho, cp = 205, 2700, 900  # Aluminium
alpha = k / (rho * cp)

# Convection
h_conv = 0  # Pas de convection

# Discrétisation
Dx, Dy, Dz = Lx / Nx, ly / Ny, thickness
Dt = min(Dx**2, Dy**2) / (8 * alpha)  # Assurer la stabilité
Nt = round(TEMPS_TOTAL / Dt)

# Matrices de calcul
T = np.ones((Nx, Ny)) * (273.15 + 25)  # Température initiale
T_new = np.copy(T)
Power = np.zeros((Nx, Ny))

# Puissance appliquée
P_in, P_x, P_y = 1.5, round((Lx/4) / Dx), round((ly/2) / Dy)
Power[P_x, P_y] = P_in

# Indices de mesure
T_x, T_y = round((Lx/2) / Dx), round((ly/2) / Dy)
Therm_x, Therm_y = round((3*Lx/4) / Dx), round((ly/2) / Dy)
T[T_x, T_y] = 273.15 + 45  # Point initialement plus chaud

# Simulation
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
x, y = np.meshgrid(np.linspace(0, Lx * 1e3, Nx), np.linspace(0, ly * 1e3, Ny))
def update(frame):
    global T, T_new
    for _ in range(Nt // 50):  # Accélérer la simulation
        T_new[:, :] = T  # Copier la température précédente
        T_new[1:-1, 1:-1] += Dt * alpha * (
            (T[:-2, 1:-1] - 2*T[1:-1, 1:-1] + T[2:, 1:-1]) / Dx**2 +
            (T[1:-1, :-2] - 2*T[1:-1, 1:-1] + T[1:-1, 2:]) / Dy**2
        )
        T_new += Dt / (rho * cp) * Power / (Dx * Dy * Dz)
        T = np.copy(T_new)
    
    ax.clear()
    ax.plot_surface(x, y, T.T - 273.15, cmap=cm.inferno)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("Température (°C)")
    ax.set_title(f"Évolution thermique - t = {frame * Dt * (Nt // 500):.1f} s")
    #ax.set_zlim(20, 60)

ani = animation.FuncAnimation(fig, update, frames=2000, interval=50)
plt.show()

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

# Paramètres
TEMPS_TOTAL = 2000  # Durée totale de la simulation [s]
Lx, ly, thickness = 117e-3, 61e-3, 1.5e-3  # Dimensions de la plaque [m]
Nx, Ny = 117, 61  # Discrétisation spatiale

# Propriétés du matériau
k, rho, cp = 205, 2700, 900  # Aluminium
alpha = k / (rho * cp)  # Diffusivité thermique [m²/s]

# Convection
h_conv = 0  # Pas de convection

# Discrétisation
Dx, Dy, Dz = Lx / Nx, ly / Ny, thickness
Dt = min(Dx**2, Dy**2) / (16 * alpha)  # Pas de temps réduit pour une meilleure précision
Nt = round(TEMPS_TOTAL / Dt)  # Nombre total d'itérations

# Matrices de calcul
T = np.ones((Nx, Ny)) * (273.15 + 25)  # Température initiale [K]
T_new = np.copy(T)
Power = np.zeros((Nx, Ny))  # Matrice de puissance appliquée

# Puissance appliquée
P_in, P_x, P_y = 1.5, round((Lx / 4) / Dx), round((ly / 2) / Dy)
Power[P_x, P_y] = P_in  # Appliquer la puissance à un point spécifique

# Point initialement plus chaud
T_x, T_y = round((Lx / 2) / Dx), round((ly / 2) / Dy)
T[T_x, T_y] = 273.15 + 45  # Température initiale plus élevée

# Fonction pour mettre à jour la température
def update_temp(T):
    T_new = np.copy(T)
    T_new[1:-1, 1:-1] += Dt * alpha * (
        (T[:-2, 1:-1] - 2 * T[1:-1, 1:-1] + T[2:, 1:-1]) / Dx**2 +
        (T[1:-1, :-2] - 2 * T[1:-1, 1:-1] + T[1:-1, 2:]) / Dy**2
    )
    T_new += Dt / (rho * cp) * Power / (Dx * Dy * Dz)  # Ajouter la puissance
    return T_new

# Initialisation de la figure
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Grille spatiale
x = np.linspace(0, Lx * 1e3, Nx)  # Axe x en mm
y = np.linspace(0, ly * 1e3, Ny)  # Axe y en mm
x, y = np.meshgrid(x, y)  # Grille 2D pour x et y

# Fonction d'animation
def update(frame):
    global T, T_new
    for _ in range(Nt // 500):  # Accélérer la simulation
        T_new[:, :] = T  # Copier la température précédente
        T_new[1:-1, 1:-1] += Dt * alpha * (
            (T[:-2, 1:-1] - 2*T[1:-1, 1:-1] + T[2:, 1:-1]) / Dx**2 +
            (T[1:-1, :-2] - 2*T[1:-1, 1:-1] + T[1:-1, 2:]) / Dy**2
        )
        T_new += Dt / (rho * cp) * Power / (Dx * Dy * Dz)
        T = np.copy(T_new)
    
    ax.clear()
    ax.plot_surface(x, y, T.T - 273.15, cmap=cm.inferno)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("Température (°C)")
    ax.set_title(f"Évolution thermique - t = {frame * Dt * (Nt // 500):.1f} s")
    ax.set_zlim(20, 60)

ani = animation.FuncAnimation(fig, update, frames=200, interval=50)
plt.show()
"""