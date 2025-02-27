#include <avr/interrupt.h>  // Inclut la bibliothèque pour gérer les interruptions
#include <math.h>            // Inclut la bibliothèque mathématique pour les calculs (log, etc.)

// Déclaration des variables globales
volatile uint8_t currentChannel = 0;         // Canal ADC actuel (pour échantillonner les températures)
volatile uint16_t adcValues[3];              // Tableau pour stocker les valeurs des températures mesurées
#define actuateur 9                          // Pin pour le contrôle du PWM (par exemple, la sortie du régulateur)


// Paramètres pour la conversion des résistances en température (formule de Steinhart-Hart)
float R25 = 10000;                           // Résistance à 25°C en ohms
float A = 0.00335401643468053;               // Coefficients de la formule Steinhart-Hart
float B = 0.000256523550896126;
float C = 0.00000260597012072052;
float D = 0.000000063292612648746;
float Rfixe = 10000;                         // Résistance fixe dans le pont diviseur (mesure de la résistance NTC)
float temperature_cible = 25;                // Température cible pour le régulateur (en °C)
float tension_commande = 0;                  // Commande du régulateur (signal PWM)
float gain_0 = ;                             // gain du circuit de conditionnement pour T1
float gain_1 = 0;                             // gain du circuit de conditionnement pour T2
float gain_2 = 0;                             // gain du circuit de conditionnement pour T3
float offset_0 = 1.52;                             // Offset du circuit de conditionnement pour T1
float offset_1 = 1.74;                             // Offset du circuit de conditionnement pour T2
float offset_2 = 1.98;                             // Offset du circuit de conditionnement pour T3

// Variables pour le régulateur PIDF (Proportionnel, Intégral, Dérivé, filtre)
float f =  10;                                //fréquence d'échantillonnage
float error[3] = {0, 0, 0};                         // Erreur 0 -1 -2
float u[2] = {0, 0};                                // commande -1 -2
float coef[5] = {11.9124, 23.3739, 11.7876, -0.3696, 0.6304}; //coefficients du PIDF


// Fonction de calcul de la température à partir de la résistance du capteur NTC
float calculerTemperature(float Rntc) {
  float lnR = log(Rntc / R25);               // Calcul du logarithme de la résistance par rapport à la résistance à 25°C
  float T = 1.0 / (A + B * lnR + C * lnR * lnR + D * lnR * lnR * lnR);  // Formule de Steinhart-Hart pour obtenir la température en Kelvin
  return T - 273.15;                         // Conversion de Kelvin en Celsius
}



//Vu que le prof a dit qu'il est preferable d'utilsier un timer de 16 bits , j'ai pris Timer1 et Timer3 qui genenrent des valeurs jusqu'a 2**(16)
// Initialisation des timers pour l'échantillonnage des températures et pour le contrôle PWM
void initTimers() {
  // Timer1 pour l'échantillonnage des températures à une fréquence de 10 Hz
  TCCR1B = (1 << WGM12) | (1 << CS12);      // Mode CTC (Clear Timer on Compare Match), prescaler 256
  OCR1A = (16000000/(256*f)) - 1;                             // Valeur de comparaison pour obtenir une fréquence d'échantillonnage de 10 Hz (une mesure par seconde)
  TIMSK1 |= (1 << OCIE1A);                   // Activation de l'interruption pour Timer1 (chaque fois que le compteur atteint OCR1A)

  // Timer3 pour contrôler la sortie PWM 
  // Configurer Timer3 en mode Fast PWM (mode 14), avec ICR3 comme TOP
  TCCR3A = (1 << COM3A1) | (1 << WGM31);  // Active PWM sur OC3A, Fast PWM mode 14
  TCCR3B = (1 << WGM33) | (1 << WGM32) | (1 << CS30);  // Fast PWM, prescaler = 1

  ICR3 = 15999;  // Définir le TOP à 15999 pour avoir 1KHz
  OCR3A = 7999;            // Valeur de comparaison pour obtenir une fréquence de 1 KHz pour le PWM
  //TIMSK3 |= (1 << OCIE3A);                   // Activation de l'interruption pour Timer3
}


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
  currentChannel = 0;                           // Réinitialise le compteur de canal
  ADMUX = (ADMUX & 0xF8) | currentChannel;      // Sélectionne le canal A0
  ADCSRA |= (1 << ADSC);                        // Démarre la conversion ADC
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
    // Convertir les valeurs numériques 0 - 1023 en tensions 
    // On tient compte du conditionnement de signal , raison pour laquelle la formule est un peu complexe
    float lecture_analogique_0 = (adcValues[0] / 1023.0) * 5.0; 
    float lecture_analogique_1 = (adcValues[1] / 1023.0) * 5.0; 
    float lecture_analogique_2 = (adcValues[2] / 1023.0) * 5.0; 

    //Déconditionnement de signal
    lecture_analogique_0 = (lecture_analogique_0/gain_0) - offset_0;
    lecture_analogique_1 = (lecture_analogique_1/gain_1) - offset_1;
    lecture_analogique_2 = (lecture_analogique_2/gain_2) - offset_2;

    //Conversion en résistance
    float Rntc_0 = Rfixe * (5.0 /lecture_analogique_0 - 1.0);
    float Rntc_1 = Rfixe * (5.0 /lecture_analogique_1 - 1.0);
    float Rntc_2 = Rfixe * (5.0 / lecture_analogique_2 - 1.0);

    // Calcul de la température mesurée à partir des valeurs de résistance
    float temperature_0_mes = calculerTemperature(Rntc_0);
    float temperature_1_mes = calculerTemperature(Rntc_1);
    float temperature_2_mes = calculerTemperature(Rntc_2);


    // Calcul du PID pour ajuster le PWM en fonction de la température mesurée
    error[2] = error[1]; // décalage du tableau
    error[1] = error[0];
    error[0] = temperature_cible - temperature_2_mes;  // Erreur entre la température cible et la température mesurée
    // Calcul du signal de commande (PWM)
    u[1] = u[0];   // décalage du tableau
    u[0] = coeff[0]*error[0] + coeff[1]*error[1] + coeff[2]*error[2] + coeff[3]*u[0] + coeff[4]*u[1];
    tension_commande = round(((u[0] + 1) / 2) * 15999);
    OCR3A = tension_commande;  // Applique le signal PWM au pin actuateur
  }
}

void loop() {
  Serial.print("Température T1 : ");
  Serial.println(calculerTemperature((float)adcValues[0]));  // Affiche la température T1
  Serial.print("Température T2 : ");
  Serial.println(calculerTemperature((float)adcValues[1]));// Affiche la température T2
  Serial.print("Température T3 : ");
  Serial.println(calculerTemperature((float)adcValues[2]));// Affiche la température T3
}
