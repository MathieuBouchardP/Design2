#include <avr/interrupt.h>  
#include <math.h>
 
// Déclaration des variables globales
volatile uint8_t currentChannel = 0;         // Canal ADC actuel (pour échantillonner les températures)
volatile uint16_t adcValues[3];              // Tableau pour stocker les valeurs des températures mesurées
volatile uint8_t interruptCount = 0;
#define actuateur 5                          // Pin pour le contrôle du PWM (par exemple, la sortie du régulateur)
 
// Paramètres pour la conversion des résistances en température (formule de Steinhart-Hart)
float R25 = 10000;                           // Résistance à 25°C en ohms
float A = 0.00335401643468053;               // Coefficients de la formule Steinhart-Hart
float B = 0.000256523550896126;
float C = 0.00000260597012072052;
float D = 0.000000063292612648746;
float Rfixe[3] = {9948, 9950, 9900};                         // Résistance fixe dans le pont diviseur (mesure de la résistance NTC) T1, T2, T3
float gain[3] = {2.40735, 3.13733,4.64232};                             // gains des circuits de conditionnement pour T1, T2, T3
float offset[3] = {1.52, 1.74, 1.98};                             // Offsets du circuit de conditionnement pour T1, T2, T3
 
// Variables pour le régulateur PIDF
float f =  0.4;                                //fréquence d'échantillonnage
float error[2] = {0, 0};                         // Erreur 0 -1 -2
float u[2] = {0,0};                                // commande -1 -2
float y[2] = {0,0};
float Kc = 0.05818;
float coeff_u = 0.035427; //coefficients du PI (équation récurrente) u(k)
float coeff_y[2] = {0, 0.96457};  //coefficients du PIDF (équation récurrente) y(k)- y(k-1) y(k-2)
//coefficients du DF
 
float coeff_ym[2] = {1.67107357564755,  -1.01332596656011};
float coeff_ym_a[2] = {0, 0.342252390912559}; // coefficients multipliés par-1
float T3_m[2] = {25, 25};
float T3[2] = {25, 25};
 
//Variables pour estimer T3
 
float num3_2[5] = {0, 0,  0, 0, 0.457504805944544} ;
float den3_2[2] = {0, 0.474132406960294}; //coefficients multipliés par -1
 
float T2[5] = {25,25,25,25,25}; //T2(k) T2(k-1) T2(k-2) T2(k-3) T2(k-4)
float T3_2[2] = {25,25} ;//T3_2(k) T3_2(k-1)
 
float num3_1[6] = {0, 0,  0,  0,  0.0258353250191554, 0.0192717560628011} ;
float den3_1[3] = {0, 1.34773025690786, -0.414201051310975}; //coefficients multipliés par -1
 
float T1[6] = {25,25,25,25,25,25}; //T1(k) T1(k-1) T1(k-2) T1(k-3) T1(k-4) T1(k-5)
float T3_1[3] = {25,25,25}; //T3_1(k) T3_1(k-1) T3_1(k-2)
 
 
// Fonction de calcul de la température à partir de la résistance du capteur NTC
float calculerTemperature(float lect, float gain, float offset, float Rfixe) {
  float lecture_analogique = (lect / 1023.0) * 5.0;
  lecture_analogique = (lecture_analogique/gain) + offset;
  float Rntc = Rfixe / (4.9850/lecture_analogique - 1.0);
  float lnR = log(Rntc / R25);               // Calcul du logarithme de la résistance par rapport à la résistance à 25°C
  float T = 1.0 / (A + B * lnR + C * lnR * lnR + D * lnR * lnR * lnR);  // Formule de Steinhart-Hart pour obtenir la température en Kelvin
  return T - 273.15;                         // Conversion de Kelvin en Celsius
}
 
// Fonction pour roter le tableau
void rotate(float *array, int size) {
    if (size <= 1) return; // Rien à faire si le tableau a 1 élément ou moins
    for (int i = size - 1; i > 0; i--) {
        array[i] = array[i - 1]; // Décalage des éléments vers la droite
    }
}
//Fonction pour faire un produit scalaire
float dot(const float* A, const float* B, int taille) {
    float somme = 0.0f;
    for (int i = 0; i < taille; i++) {
        somme += A[i] * B[i];
    }
    return somme;
}
 
 
// Fonction qui estime T3 en fonction de T1 et T2
float estimer_T3 (float temperature_0_mes,float temperature_1_mes){
  rotate(T2, 5) ;
  T2[0] = temperature_1_mes;
  rotate(T3_2, 2) ;
  rotate(T1, 6) ;
  T1[0] = temperature_0_mes;
  rotate(T3_1, 3);
  T3_1[0] = dot(num3_1,T1,6) + dot(T3_1, den3_1,3) ;
  T3_2[0] = dot(T2,num3_2,5) + dot(T3_1, den3_2,2) ;
  return (T3_1[0] + T3_2[0])/2;
}
 
// Fonction qui calcule la commande
int calculerCommande(float temperature_cible) {
    rotate(error, 2);  // Décale les erreurs dans le tableau
    error[0] = temperature_cible - T3[0];  // Calcul de l'erreur actuelle
 
    rotate(u, 2);  // Décale le tableau des commandes
    rotate(y,2);
    y[0] =  coeff_u * (u[1] - 0.05) + dot(coeff_y, y ,2);
    // Calcul du signal de commande en utilisant un filtre récursif
 
    u[0] = Kc*error[0] +y[0] + 0.05;
    Serial.println(u[0]);
    // Conversion du signal en une valeur PWM (0 à 15999)
    int tension_commande = round(((u[0] + 1) / 2) * 15999);
    return constrain(tension_commande, 0, 15999);
}
 
// Initialisation des timers pour l'échantillonnage des températures et pour le contrôle PWM
void initTimers() {
  // Timer1 pour l'échantillonnage des températures à une fréquence de 10 Hz
  // Timer1 pour l'échantillonnage des températures à une fréquence de 10 Hz
  TCCR1A = 0;                                 // Désactiver Timer1
  TCCR1B = 0;                                 // Désactiver Timer1
  TCNT1 = 0;                                  // Réinitialiser le compteur
  TCCR1B = (1 << WGM12) | (1 << CS12) | (1 << CS10); // Mode CTC (Clear Timer on Compare Match), prescaler 1024
  OCR1A = (16000000 / (1024 * f)) - 1;       // Valeur de comparaison pour obtenir une fréquence d'échantillonnage de 10 Hz (une mesure par seconde)
  TIMSK1 |= (1 << OCIE1A);                    // Activation de l'interruption pour Timer1 (chaque fois que le compteur atteint OCR1A)
 
  // Timer3 pour contrôler la sortie PWM
  // Configurer Timer3 en mode Fast PWM (mode 14), avec ICR3 comme TOP
  TCCR3A = (1 << COM3A1) | (1 << WGM31);  // Active PWM sur OC3A, Fast PWM mode 14
  TCCR3B = (1 << WGM33) | (1 << WGM32) | (1 << CS30);  // Fast PWM, prescaler = 1
  ICR3 = 15999;  // Définir le TOP à 15999 pour avoir 1KHz
  OCR3A = 8000;            // Valeur de comparaison pour obtenir une fréquence de 1 KHz pour le PWM
}
 
//Fonction Setup
void setup() {
  Serial.begin(115200);                       // Initialisation de la communication série pour le débogage
  pinMode(actuateur, OUTPUT);                // Initialisation de la pin pour le PWM comme sortie
  // Conffiguration de l'ADC pour lire les températures
  ADMUX = (1 << REFS0)| currentChannel;                     // Utilise la référence de tension à 5V (REFS0)
  ADCSRA = (1 << ADEN) | (1 << ADIE) | (1 << ADPS2) | (1 << ADPS1);  // Active l'ADC et les interruptions, avec un prescaler de 64
  ADCSRB = 0;                               // Mode de fonctionnement par défaut (pas de multiplexer externe)
  initTimers();                             // Appel de la fonction pour initialiser les timers
  sei();                                    // Activation des interruptions globales (nécessaire pour que les ISR fonctionnent) permet  à l'arduino de répondre aux interruptions genréres par les timers
}
 
// Interruption Timer1 : Démarre la conversion ADC sur A0
ISR(TIMER1_COMPA_vect) {
  interruptCount++;
  if (interruptCount == 2) {
  interruptCount = 0;
  currentChannel = 0;                           // Réinitialise le compteur de canal
  ADMUX = (ADMUX & 0xF8) | currentChannel;      // Sélectionne le canal A0
  ADCSRA |= (1 << ADSC);  
  }                      // Démarre la conversion ADC
 
}
 
// Interruption ADC : Stocke la valeur et passe au canal suivant
ISR(ADC_vect) {
  adcValues[currentChannel] = ADC;               // Stocke la valeur convertie
  currentChannel++;                              // Passe au canal suivant
  if (currentChannel < 3) {              // Encore des canaux à lire ?
    ADMUX = (ADMUX & 0xF8) | currentChannel;     // Sélectionne le prochain canal
    ADCSRA |= (1 << ADSC);                       // Lance la conversion suivante
  }
  else {
    // Calcul de la température mesurée à partir des valeurs de résistance
    float temperature_0_mes = calculerTemperature(adcValues[0], gain[0], offset[0], Rfixe[0]);
    float temperature_1_mes = calculerTemperature(adcValues[1], gain[1], offset[1], Rfixe[1]);
    float temperature_2_mes = calculerTemperature(adcValues[2], gain[2], offset[2],  Rfixe[2]);
    Serial.println(temperature_2_mes);// Affiche la température T3
    rotate(T3_m, 2);
    //T3_m[0] = estimer_T3 (temperature_0_mes, temperature_1_mes); // Asservissement avec T1, T2
    T3_m[0] = temperature_2_mes ; // Asservissement avec T3
    Serial.println(T3_m[0]);
    rotate(T3,2);
    T3[0] = dot(coeff_ym, T3_m, 2) + dot(T3, coeff_ym_a, 2);
    //Serial.println(calculerCommande(26, T3, error, u, coeff_u,coeff_y));
    // Calcul du PID pour ajuster le PWM en fonction de la température mesurée
    OCR3A = calculerCommande(25);  // Applique le signal PWM au pin actuateur
  }
}
void loop() {
}
