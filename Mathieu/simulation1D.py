import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Paramètres
temps_total = 500  # Durée de la simulation

Lx = 120e-3  # Longueur [m]
thickness = 1.5e-3  # Épaisseur de la plaque [m]

Nx = 120  # Nombre d'éléments en x

# Propriétés du matériau (Aluminium)
k = 205  # Conductivité Thermique [W/m·K]
rho = 2700  # Densité [kg/m^3]
cp = 900  # Chaleur spécifique [J/kg·K]
alpha = k / (rho * cp)  # Diffusivité Thermique [m^2/s]

h_conv = 20  # Coeff. de convection [W/m^2·K]
h_conv = 0  # Commenter pour enlever la convection

# Paramètres calculés
dx = Lx / Nx
dy = thickness
dz = thickness
dt = dx**2 / (8 * alpha)  # Pas en temps [s]
Nt = round(temps_total / dt)  # Nombre d'itérations temporelles

# Aires pour la convection
aire_bouts = dy * dz
aire_sides = dx * dz
aire_top = dx * dy

volume = dx * dy * dz

Temps = np.arange(Nt) * dt  # Vecteur de temps
Position = np.arange(Nx) * dx  # Vecteur de position

# Puissance déposée
Pin = 1.5  # Puissance [W]
Pin_loc_x = round((Lx / 4) / dx)  # Élément qui reçoit la puissance

P = np.zeros(Nx)
P[Pin_loc_x] = Pin

# Conditions initiales
T_piece = 273 + 25  # Température pièce [K]
T = np.ones(Nx) * T_piece

T_loc_x = round((Lx / 2) / dx)  # Élément plus chaud
# T[T_loc_x] = 273 + 45  # Un élément plus chaud

Therm_loc_x = round((3 * Lx / 4) / dx)

# Préallocation des vecteurs pour stockage
energy_added = np.zeros(Nt)
energy_loss = np.zeros(Nt)
thermistance = np.zeros(Nt)
T_new = np.zeros_like(T)

# Création de la figure
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Itération temporelle
for t in range(Nt):
    for i in range(Nx):
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
    
    T[:] = T_new[:]
    thermistance[t] = T[Therm_loc_x]
    energy_added[t] = np.sum(P) * dt
    energy_loss[t] = h_conv * np.sum(T - T_piece) * (2 * aire_sides + 2 * aire_top + aire_bouts) * dt
    
    # Affichage à certains intervalles
    if t % (Nt // 1000) == 0 or t == 0:
        axs[0].cla()
        axs[0].plot(Position, T - 273)
        axs[0].set_title(f'Température à t = {t*dt:.2f} s')
        axs[0].set_xlabel('Position [m]')
        axs[0].set_ylabel('Température [°C]')
        axs[0].grid()
        
        axs[1].cla()
        axs[1].plot(Temps[:t], thermistance[:t] - 273)
        axs[1].set_xlabel('Temps [s]')
        axs[1].set_ylabel('Température [°C]')
        axs[1].set_title('Température à la thermistance')
        axs[1].grid()
        
        axs[2].cla()
        axs[2].plot(Temps[:t], energy_added[:t], label='Énergie déposée')
        axs[2].plot(Temps[:t], energy_loss[:t], label='Énergie dissipée')
        axs[2].set_xlabel('Temps [s]')
        axs[2].set_ylabel('Énergie [J]')
        axs[2].legend()
        axs[2].grid()
        
        plt.pause(0.01)

plt.show()
