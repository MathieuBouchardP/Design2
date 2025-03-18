#include <avr/interrupt.h>
#include <math.h>
#define actuateur 5 // Broche de sortie 

volatile uint8_t currentChannel = 0;
volatile uint16_t adcValues[3];
volatile uint8_t interruptCount = 0;
char buffer[50];

/*Fonctions utilitaires : Décalage vers la droite ; Produit scalaire dot*/
//-----------------------------------------------------------------------------------------//
void rotate(float *array, int size) {
    if (size <= 1) return;
    for (int i = size - 1; i > 0; i--) {
        array[i] = array[i - 1];
    }
}

float dot(const float* A, const float* B, int taille) {
    float somme = 0.0f;
    for (int i = 0; i < taille; i++) {
        somme += A[i] * B[i];
    }
    return somme;
}

void split(String input, float* floatT, int maxValues) {
    int index = 0;
    int pos = 0;
    int lastPos = 0;
    // Trouver la position de la première virgule
    pos = input.indexOf(',');
    // Si aucune virgule n'est trouvée, convertir toute la chaîne en float
    if (pos == -1) {
        floatT[0] = input.toFloat();  // Convertir la chaîne en float
        return;  // Sortir de la fonction
    }
    // Si des virgules sont trouvées, diviser la chaîne
    while (pos != -1 && index < maxValues - 1) {
        // Extraire la sous-chaîne entre la dernière position et la position actuelle
        String substring = input.substring(lastPos, pos);
        // Convertir la sous-chaîne en float et la stocker dans le tableau
        floatT[index++] = substring.toFloat();
        // Mettre à jour la dernière position
        lastPos = pos + 1;
        // Trouver la position de la prochaine virgule
        pos = input.indexOf(',', lastPos);
    }
    // Ajouter la dernière sous-chaîne après la dernière virgule (si on a encore de la place)
    if (index < maxValues) {
        String substring = input.substring(lastPos);
        floatT[index] = substring.toFloat();
    }
}

/*Variables globales*/
//-----------------------------------------------------------------------------------------//
float f = 0.4; // moitié de la fréquence d'échantillonage
float R25 = 10000;    // Résistance à 25 dégré Celsius
//Coefficients de SteinHart-Hart
float A = 0.00335401643468053;
float B = 0.000256523550896126;
float C = 0.00000260597012072052;
float D = 0.000000063292612648746;
//Banque d'offsets et de gains
float Rfixe[3] = {9948, 9950, 9900};
float gain[3] = {2.40735, 3.13733,4.64232};
float offset[3] = {1.52, 1.74, 1.98};

/*Régulateur*/
//-----------------------------------------------------------------------------------------//

float error[2] = {0, 0}; // Tableau des erreurs (Consigne - Valeur mesurée)

float u[2] = {0,0}; // Commande (sortie du régulateur, entrée du variateur)
float y[2] = {0,0}; // 

// Coefficients du régulateur
float Kc = 0.05818; //Coefficient proportionnel
// Coefficient du PI
float coeff_u = 0.035427;  // numérateur 
float coeff_y[2] = {0, 0.96457}; //Dénominateur
//Action dérivative et filtre (DF)
float coeff_ym[2] = {1.67107357564755,  -1.01332596656011}; // Numérateur
float coeff_ym_a[2] = {0, 0.342252390912559}; // Dénominateur

//Variables d'entrée et de sortie du Derivative Kick
float T3_m[2] = {25, 25}; // Avant, T3 mesurée / estimé
float T3[2] = {25, 25}; // Après , T3 servant à calculer l'erreur

/*Estimation de T3*/
//-----------------------------------------------------------------------------------------//
//Coefficients pour estimer T3 à partir de T2
float num3_2[5] = {0, 0,  0, 0, 0.457504805944544} ; //numérateur
float den3_2[2] = {0, 0.474132406960294};  // Dénominateur
//Variables d'entrée et de sortie de l'estimation de T3 à partir de T2
float T2[5] = {25,25,25,25,25}; 
float T3_2[2] = {25,25} ;
 //Coefficients pour estimer T3 à partir de T1
float num3_1[6] = {0, 0,  0,  0,  0.0258353250191554, 0.0192717560628011} ;
float den3_1[3] = {0, 1.34773025690786, -0.414201051310975}; 
//Variables d'entrée et de sortie de l'estimation de T3 à partir de T1
float T1[6] = {25,25,25,25,25,25}; 
float T3_1[3] = {25,25,25}; 

/* Variables d'asservissement reçue de l'interface*/
//-----------------------------------------------------------------------------------------//
volatile float lastTemperature = 0.0;    // Dernière température calculée
float consigne = 25.0;
bool asservissementActif = false;
float Consigne = 25;
bool Ti_set = false, Kc_set = false, Consigne_set = false, Td_set = false, Tf_set = false;
String inputString = ""; // Stocker la ligne reçue
bool stringComplete = false;

/*Nettoyage du buffer d'entrée ou de sortie*/
//-----------------------------------------------------------------------------------------//
void flushSerial() {
    // Vider le tampon de réception
    while (Serial.available()) {
        Serial.read();  // Lire et jeter toutes les données
    }
}

/*Gestion du flux de données sur le port série*/
//-----------------------------------------------------------------------------------------//
void handleSerial() {
    String input = Serial.readString();  // Lire la chaîne envoyée par MATLAB
    input.trim();  // Enlever les espaces en début et en fin de la chaîne
    if (input.startsWith("START"))
    {
      Serial.println("Asservissement debute.");
      startAsservissement();

    }
        if (input.startsWith("STOP"))
    {
      stopAsservissement();
      Serial.println("Asservissement arrete."); 
    }
    
    if (input.startsWith("PIDF")) {
        // Enlever le préfixe "PIDF-"
        input = input.substring(5);  // La partie après "PIDF"
 
        // Extraire les valeurs séparées par des barres '|'
        Kc = input.substring(0, input.indexOf('|')).toFloat();
        input = input.substring(input.indexOf('|') + 1);
 
        String AA = input.substring(0, input.indexOf('|'));
        split(AA,coeff_ym,2);
        input = input.substring(input.indexOf('|') + 1);
 
       String BB = input.substring(0, input.indexOf('|'));
       split(BB,coeff_ym_a,2);
        input = input.substring(input.indexOf('|') + 1);
 
        coeff_u = input.substring(0, input.indexOf('|')).toFloat();
        input = input.substring(input.indexOf('|') + 1);
 
        String DD = input.substring(0, input.indexOf('|'));
        split(DD,coeff_y,2);
        input = input.substring(input.indexOf('|') + 1);
 
 
        Consigne = input.toFloat();  // Le dernier paramètre, pas besoin de chercher un autre tiret
         Serial.println("Commande valide."); 
    } else {
        Serial.println("Commande invalide.");
    }
}

/*Arrêt de l'asservissement et début de l'asservissement*/
//-----------------------------------------------------------------------------------------//

void  stopAsservissement()
{
  asservissementActif = false;
    OCR3A = 8000;
}
void startAsservissement()
{
 asservissementActif = true; 
}

/*Fonction pour calculer la température*/
//-----------------------------------------------------------------------------------------//

float calculerTemperature(float lect, float gain, float offset, float Rfixe) {
    float lecture_analogique = (lect / 1023.0) * 5.0;
    lecture_analogique = (lecture_analogique/gain) + offset;
    float Rntc = Rfixe / (4.9850/lecture_analogique - 1.0);
    float lnR = log(Rntc / R25);
    float T = 1.0 / (A + B * lnR + C * lnR * lnR + D * lnR * lnR * lnR);
    return T - 273.15;
}

/* Fonction qui estime T3 en fonction de T1 et T2*/
//-----------------------------------------------------------------------------------------//
float estimer_T3 (float temperature_0_mes,float temperature_1_mes){
  rotate(T2, 5) ;
  T2[0] = temperature_1_mes;
  rotate(T3_2, 2) ;
  rotate(T1, 6) ;
  T1[0] = temperature_0_mes;
  rotate(T3_1, 3);
  T3_1[0] = dot(num3_1,T1,6) + dot(T3_1, den3_1,3) ;
  T3_2[0] = dot(T2,num3_2,5) + dot(T3_1, den3_2,2) ;
  return 0.3*T3_1[0] + 0.7*T3_2[0];
}

/*Fonction pour calculer la commande */
//-----------------------------------------------------------------------------------------//
int calculerCommande(float temperature_cible) {
    rotate(error, 2);
    error[0] = temperature_cible - T3[0];
    rotate(u, 2);
    rotate(y,2);
    float uop = 0.001 ;
    y[0] =  coeff_u * (u[1] - uop) + dot(coeff_y, y ,2);
    u[0] = Kc * error[0] + y[0] + uop;
    int tension_commande = round(((u[0] + 1) / 2) * 15999);
    return constrain(tension_commande, 0, 15999);
}

/*Fonction pour initialiser les compteurs */
//-----------------------------------------------------------------------------------------//
void initTimers() {
    TCCR1A = 0;
    TCCR1B = (1 << WGM12) | (1 << CS12) | (1 << CS10);
    OCR1A = (16000000 / (1024 * f)) - 1;
    TIMSK1 |= (1 << OCIE1A);

    TCCR3A = (1 << COM3A1) | (1 << WGM31);
    TCCR3B = (1 << WGM33) | (1 << WGM32) | (1 << CS30);
    ICR3 = 15999;
    OCR3A = 8000;
}
/*Routine d'interruption pour l'échantillonnage */
//-----------------------------------------------------------------------------------------//
ISR(TIMER1_COMPA_vect) {
    interruptCount++;
    if (interruptCount == 2) {
        interruptCount = 0;
        currentChannel = 0;
        ADMUX = (ADMUX & 0xF8) | currentChannel;
        ADCSRA |= (1 << ADSC);
}}

/*Routine d'interruption après la demande conversion par l'ADC*/
//-----------------------------------------------------------------------------------------//
ISR(ADC_vect) {
    adcValues[currentChannel] = ADC;
    currentChannel++;
    if (currentChannel < 3) {
        ADMUX = (ADMUX & 0xF8) | currentChannel;
        ADCSRA |= (1 << ADSC);
    } else {
        float temperature_0_mes = calculerTemperature(adcValues[0], gain[0], offset[0], Rfixe[0]);
        float temperature_1_mes = calculerTemperature(adcValues[1], gain[1], offset[1], Rfixe[1]);
        float temperature_2_mes = calculerTemperature(adcValues[2], gain[2], offset[2], Rfixe[2]);
        rotate(T3_m, 2);
        T3_m[0] = temperature_2_mes;
        rotate(T3,2);
        T3[0] = dot(coeff_ym, T3_m, 2) + dot(T3, coeff_ym_a, 2); 

        if(asservissementActif){
        OCR3A = calculerCommande(Consigne);
        }
        sprintf(buffer, "Data%.2f,%.2f,%.2f,%.2f", temperature_0_mes, temperature_1_mes, temperature_2_mes, T3[0]); 
        Serial.println(buffer);
        }
}

/*Setup (première fonction exécutée dès l'alimentation) et Loop ( Fonction exéctée en boucle) */
//-----------------------------------------------------------------------------------------//
void setup() {
    Serial.begin(115200);
    pinMode(actuateur, OUTPUT);
    ADMUX = (1 << REFS0) | currentChannel;
    ADCSRA = (1 << ADEN) | (1 << ADIE) | (1 << ADPS2) | (1 << ADPS1);
    initTimers();
    sei();
}

void loop() {
  //flushSerial();
    if (Serial.available()) {
        handleSerial();
    }
}