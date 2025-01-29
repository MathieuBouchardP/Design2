import serial
import csv
import time
import matplotlib.pyplot as plt
 
# Configuration du port série
PORT = 'COM5'  # Remplacez par votre port série
BAUDRATE = 19200  # Ajustez selon la configuration de votre périphérique
OUTPUT_FILE = 'data_log.csv'
 
# Fonction pour lire les données du port série
def lire_donnees(port, baudrate, duree=10):
    """
    Lit les données du port série pendant une durée spécifiée.
    """
    try:
        ser = serial.Serial(port, baudrate, timeout=1)
        print(f"Connexion au port {port} réussie.")
        debut = time.time()
        donnees = []
       
        while time.time() - debut < duree:  # Lire pendant "duree" secondes
            line = ser.readline().decode('utf-8').strip()
            if line:
                temps_relatif = time.time() - debut  # Temps écoulé depuis le début
                print(f"{temps_relatif:.2f} sec - {line}")
               
                # Extraction des données
                try:
                    parts = line.split('|')
                    resistance = float(parts[0].split(':')[1].strip())  # Résistance après "Résistance NTC (Ohm):"
                    temperature = float(parts[1].split(':')[1].strip())  # Température après "Température (°C):"
 
                    donnees.append((temps_relatif, resistance, temperature))
                except (IndexError, ValueError):
                    print(f"Ligne ignorée (format incorrect) : {line}")
        ser.close()
        return donnees
    except serial.SerialException as e:
        print(f"Erreur de connexion au port série : {e}")
        return []
 
# Fonction pour consigner les données dans un fichier CSV
def consigner_donnees_csv(fichier, donnees):
    """
    Écrit les données dans un fichier CSV.
    """
    try:
        with open(fichier, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['Temps (s)', 'Résistance (Ohms)', 'Température (°C)'])  # En-têtes
            writer.writerows(donnees)  # Écrire les lignes
        print(f"Les données ont été consignées dans {fichier}.")
    except Exception as e:
        print(f"Erreur lors de l'écriture du fichier CSV : {e}")
 
# Fonction pour tracer les données
def tracer_donnees(donnees):
    """
    Trace les résistances et températures en fonction du temps.
    """
    try:
        temps = [t[0] for t in donnees]  # Temps relatif
        resistances = [t[1] for t in donnees]  # Résistances
        temperatures = [t[2] for t in donnees]  # Températures
       
        plt.figure(figsize=(12, 6))
       
        # Tracé des résistances
        plt.plot(temps, resistances, label='Résistance (Ohms)', marker='o', linestyle='-', color='b')
       
        # Tracé des températures
        plt.plot(temps, temperatures, label='Température (°C)', marker='s', linestyle='--', color='r')
       
        plt.title("Résistances et températures en fonction du temps")
        plt.xlabel("Temps écoulé (s)")
        plt.ylabel("Valeurs")
        plt.legend()
        plt.grid(True)
        plt.show()
    except Exception as e:
        print(f"Erreur lors du traçage des données : {e}")
 
# Programme principal
if __name__ == "__main__":
    # Étape 1 : Lire les données
    donnees = lire_donnees(PORT, BAUDRATE, duree=600)  # Durée de 10 secondes (modifiable)
    if donnees:
        # Étape 2 : Consigner les données dans un fichier CSV
        consigner_donnees_csv(OUTPUT_FILE, donnees)
       
        # Étape 3 : Tracer les données
        tracer_donnees(donnees)