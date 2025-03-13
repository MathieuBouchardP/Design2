import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm

# Conversion en millimètres
mm = 10**(-3)

# Dimensions
longueur_Z = 12 * mm
hauteur_R = 3 * mm
pas = 0.1 * mm

dimension_en_z = int(longueur_Z / pas) + 1
dimension_en_r = int(hauteur_R / pas) + 1

# Grille de calcul
matrice_pot = np.zeros((dimension_en_r, dimension_en_z))
x, y = np.meshgrid(np.linspace(0, longueur_Z / mm, dimension_en_z),
                   np.linspace(-hauteur_R / mm, hauteur_R / mm, 2 * dimension_en_r - 1))

# Conditions frontières
for i in range(31):
    matrice_pot[i, i] = -300  # Mur en angle

matrice_pot[30, 30:] = -300  # Mur du haut

# Copie pour conservation
matrice_initiale = matrice_pot.copy()

# Paramètre de sur-relaxation
w = 0.878

def diffusion():
    """Applique une itération de l'algorithme de sur-relaxation à la matrice."""
    global matrice_pot
    potentiel = matrice_pot.copy()

    for r in range(29, -1, -1):
        if r == 0:
            for z in range(1, 46):
                potentiel[r, z] = (1 + w) * (4 * potentiel[1, z] + potentiel[0, z + 1] + potentiel[0, z - 1]) / 6 - w * potentiel[r, z]
        else:
            for z in range(119, -1, -1):
                if potentiel[r, z] == -300:
                    break
                R = r * 10
                potentiel[r, z] = (1 + w) * ((potentiel[r + 1, z] + potentiel[r - 1, z] +
                                              potentiel[r, z + 1] + potentiel[r, z - 1]) / 4 +
                                             (pas / (8 * R)) * (potentiel[r + 1, z] - potentiel[r - 1, z])) - w * potentiel[r, z]

    matrice_pot = potentiel  # Mise à jour de la matrice principale

# Animation
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

def update(frame):
    """Met à jour l'animation à chaque itération."""
    diffusion()
    ax.clear()

    matrice_pot_inv = np.flip(matrice_pot, axis=0)
    potentiel_affiché = np.concatenate((matrice_pot_inv, matrice_pot[1:, :]), axis=0)

    surf = ax.plot_surface(x, y, potentiel_affiché, cmap=cm.inferno)
    ax.set_xlabel("Z (mm)")
    ax.set_ylabel("Rayon (mm)")
    ax.set_zlabel("Potentiel (V)")
    ax.set_title("Propagation du Potentiel")
    ax.set_zlim(-350, 50)  # Échelle fixe pour éviter les sauts d'échelle

ani = animation.FuncAnimation(fig, update, frames=1000, interval=10)
plt.show()
