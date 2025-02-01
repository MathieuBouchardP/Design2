import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
import matplotlib.colors as mcolors

"""Non fonctionnel pour l'instance parce que j'ai pas encore fini d'implémenter la propagation
 de la chaleur entre les élément du bords et le reste
"""
#Définition des paramètres: 

# Durée de la simulation
temps_tot = 500


# Dimension de la plaque
facteur = 0.2
Lx = 117*10**(-3)               # Longueur [m]
Nx = round(facteur*117)         # Nbr éléments en x
ly = 61*10**(-3)                # Largeur [m]
ny = round(facteur*61)          # Nbr éléments en y
thickness = 1.5*10**(-3)        # Épaisseur [m] (1 élément)


# Propriété du matériau
k = 205                         # Conductivité thermique [W/m*K]
rho = 2700                      # Densité [kg/m^3]
cp = 900                        # Chaleur spécifique [J/kg*K]  (Capacité calorigique massique)
alpha = k / (rho * cp)          # Diffusivité thermique [m^2/s]


# Propriété convection
h_conv = 0.5                     # Coefficient de convection [W/m^2*K]
#h_conv = 0                      # Aucune convection


# Paramètre calculé
dx = Lx / Nx                    # Pas de discrétisation en x
dy = ly / ny                    # Pas de discrétisation en y
dz = thickness                  # Épaisseur en z

dt = dx**2/(8*alpha)            # Pas en temps [s]    --->   Attention, il faut le chosiir pour assurer la stablilité
Nt = round(temps_tot/dt)         # Nombre d'itération temporelles
print(dt)
air_ends = dz * dy              # Aire des bouts de la plaque, pour chaque élément
air_top_et_bot = dx * dy        # Aire du dessus ou du dessous, pour chaque élément
air_sides = dz * dx             # Aire des côtés de la plaque, pour chaque élément

volume = dx*dy*dz               # Volume d'un élément

temps = np.arange(Nt) * dt
position_x = np.arange(Nx) * dx
position_y = np.arange(ny) * dy


# Puissance déposé par l'actuateur
P_in = 50                       # Puissance [W]
#P_in = 0                       # Commenter pour mettre la puissance

P_in_loc_x = Lx/4               # Localisation en x de l'actuateur sur la plaque
P_in_loc_x = round(P_in_loc_x/dx)   # Élément qui reçoit la puissance
P_in_loc_y = ly/2               # Localisation en x de l'actuateur sur la plaque
P_in_loc_y = round(P_in_loc_y/dy)   # Élément qui reçoit la puissance


Power = np.zeros((Nx, ny))                # Matrice de puissance à ajouter à chaque élément
Power[P_in_loc_x, P_in_loc_y] = P_in    # Élément sur lequel la puissance est ajouté


# Condition initiales
T_piece = 273.15 + 25           # Température de la pièce [K]
T_plaque = 273.15 + 25
T_piece_mat = T_piece * np.ones((Nx,ny))   # Temmpérature de tous les éléments
T = T_plaque * np.ones((Nx,ny))

T_loc_x = Lx/2;             # Localisation en x dun dirac en Température [m]
T_loc_x = round(T_loc_x/dx) # Element qui est plus chaud
T_loc_y = ly/2;             # Localisation en y dun dirac en Température [m]
T_loc_y = round(T_loc_y/dy) # Element qui est plus chaud
T[T_loc_x, T_loc_y] = 273.15+45;        # Un élement plus chaud

energy_added = np.zeros(Nt)
energy_loss = np.zeros(Nt)
thermistance = np.zeros(Nt)
T_new = np.zeros_like(T)

"""
for t in range(Nt):
    for i, j in T:
        T_new[i] = T[i]  # Température précédente
        
        if i == 0:
            T_new[i] += dt / (rho * cp) * k * (T[i+1] - T[i]) / dx**2
            T_new[i] += dt / (rho * cp) * h_conv * (T_piece - T[i]) * aire_bouts / volume
        elif i == Nx - 1:
            T_new[i] += dt / (rho * cp) * k * (T[i-1] - T[i]) / dx**2
            T_new[i] += dt / (rho * cp) * h_conv * (T_piece - T[i]) * aire_bouts / volume
        else:
            T_new[i] += dt / (rho * cp) * k * (T[i+1] - 2*T[i] + T[i-1]) / dx**2
        
        T_new[i] += dt / (rho * cp) * P[i] / volume
        T_new[i] += dt / (rho * cp) * h_conv * (T_piece - T[i]) * 2 * aire_sides / volume
        T_new[i] += dt / (rho * cp) * h_conv * (T_piece - T[i]) * 2 * aire_top / volume
"""
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
x, y = np.meshgrid(np.linspace(0, Lx * 1e3, Nx), np.linspace(0, ly * 1e3, ny), indexing="ij")
t=0
norm = mcolors.Normalize(vmin=10, vmax=40)
def update(frame):
    global T, T_new, T_piece_mat
    T += dt / (rho * cp) * Power / volume
    T_new = T.copy()

    # Décalages pour obtenir les voisins
    haut = np.roll(T_new, -1, axis=0)
    bas = np.roll(T_new, 1, axis=0)
    gauche = np.roll(T_new, -1, axis=1)
    droite = np.roll(T_new, 1, axis=1)

    # Moyenne des voisins
    T_new[1:-1, 1:-1] += ((haut[1:-1, 1:-1]-2*T_new[1:-1, 1:-1]+bas[1:-1, 1:-1])/dx**2 +  # Échanger en x
                          (gauche[1:-1, 1:-1]- 2*T_new[1:-1, 1:-1] + droite[1:-1, 1:-1])/dy**2 +   # Échanger en y
                          ((T_piece_mat[1:-1, 1:-1] - T_new[1:-1, 1:-1]) *h_conv *  2*air_top_et_bot / volume))* alpha*dt   # Convection top et bot
    #T_new[1:-1, 1:-1] += ((T_piece_mat[1:-1, 1:-1] - T_new[1:-1, 1:-1]) *h_conv *  2*air_top_et_bot / volume)*alpha*dt
    
    # Gestion des bords : à faire
    T_new[0, 1:-1] += alpha* h_conv * ((T_piece_mat[0, 1:-1] - T_new[0, 1:-1]) * (2*air_top_et_bot + air_ends)*dt / volume     # Bord haut
                        (haut[1:-1, 1:-1]-T_new[1:-1, 1:-1])/dx**2 +  # Échanger en x sur le haut de la plaque
                          (gauche[1:-1, 1:-1]- 2*T_new[1:-1, 1:-1] + droite[1:-1, 1:-1])/dy**2)     #Échanger en 
    # Figure out les matrice décalé si c'est en haut en bas et implémenter la propagation entre les éléments sur les côté et coin et le reste, :)
    T_new[-1, 1:-1] += alpha* h_conv * (T_piece_mat[-1, 1:-1] - T_new[-1, 1:-1]) * (2*air_top_et_bot + air_ends)*dt / volume   # Bord bas
    T_new[1:-1, 0] += alpha* h_conv * (T_piece_mat[1:-1,0] - T_new[1:-1, 0]) * (2*air_top_et_bot + air_sides)*dt / volume     # Bord gauche
    T_new[1:-1, -1] += alpha* h_conv * (T_piece_mat[1:-1, -1] - T_new[1:-1, -1]) * (2*air_top_et_bot + air_sides)*dt / volume  # Bord droit
    

    T_new[0, 0] += alpha * h_conv * (T_piece_mat[0, 0] - T_new[0, 0]) * (2 * air_top_et_bot + air_ends + air_sides) * dt / volume  # Coin haut-gauche
    T_new[0, -1] += alpha * h_conv * (T_piece_mat[0, -1] - T_new[0, -1]) * (2 * air_top_et_bot + air_ends + air_sides) * dt / volume  # Coin haut-droit
    T_new[-1, 0] += alpha * h_conv * (T_piece_mat[-1, 0] - T_new[-1, 0]) * (2 * air_top_et_bot + air_ends + air_sides) * dt / volume  # Coin bas-gauche
    T_new[-1, -1] += alpha * h_conv * (T_piece_mat[-1, -1] - T_new[-1, -1]) * (2 * air_top_et_bot + air_ends + air_sides) * dt / volume  # Coin bas-droit
    T = np.copy(T_new)
    #T += dt / (rho * cp) * Power / volume
    
    ax.clear()
    ax.plot_surface(x, y, T - 273.15, cmap=cm.magma, norm=norm)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_zlabel("Température (°C)")
    #ax.set_title(f"Évolution thermique - t = {frame * dt * (Nt // 500):.1f} s")
    ax.set_zlim(10, 75)
    t=1

def penis(yes):
    if t == 0:
        update(1)
    for i in range(10):
        update(1)


ani = animation.FuncAnimation(fig, penis, frames=5000, interval=10)
plt.show()

