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
float Rfixe = 9820;                          // Résistance fixe dans le pont diviseur (mesure de la résistance NTC)
float temperature_cible = 25;                // Température cible pour le régulateur (en °C)
float tension_commande = 0;                  // Commande du régulateur (signal PWM)

// Variables pour le régulateur PID (Proportionnel, Intégral, Dérivé)
float error = 0.0;                           // Erreur actuelle (différence entre la consigne et la température mesurée)
float prevError = 0.0;                       // Erreur précédente pour le calcul dérivé
float integral = 0.0;                        // Somme des erreurs pour le calcul intégral
float Kp = 1.0;                              // Gain proportionnel
float Ki = 0.1;                              // Gain intégral
float Kd = 0.5;                              // Gain dérivé
float kf = 0.5

// Fonction de calcul de la température à partir de la résistance du capteur NTC
float calculerTemperature(float Rntc) {
  float lnR = log(Rntc / R25);               // Calcul du logarithme de la résistance par rapport à la résistance à 25°C
  float T = 1.0 / (A + B * lnR + C * lnR * lnR + D * lnR * lnR * lnR);  // Formule de Steinhart-Hart pour obtenir la température en Kelvin
  return T - 273.15;                         // Conversion de Kelvin en Celsius
}
//Vu que e prof a dit qu'il est preferable d'utilsier un timer de 16 bits , j'ai pris Timer1 et Timer3 qui genenrent des valeurs jusqu'a 2**(16)
// Initialisation des timers pour l'échantillonnage des températures et pour le contrôle PWM
void initTimers() {
  // Timer1 pour l'échantillonnage des températures à une fréquence de 10 Hz
  TCCR1B = (1 << WGM12) | (1 << CS12);      // Mode CTC (Clear Timer on Compare Match), prescaler 256
  OCR1A = 6249;                            // Valeur de comparaison pour obtenir une fréquence d'échantillonnage de 10 Hz (une mesure par seconde)
  TIMSK1 |= (1 << OCIE1A);                   // Activation de l'interruption pour Timer1 (chaque fois que le compteur atteint OCR1A)

  // Timer3 pour contrôler la sortie PWM (1 kHz,l'on poura se conserter pour la fréquence du PWM)
  TCCR3B = (1 << WGM32) | (1 << CS30);      // Mode CTC revient à 0 lorsque le Timer attient la valeur ORC3A
             // Valeur de comparaison pour obtenir une fréquence de 1 kHz pour le PWM
  ICR1 =15999;
  TIMSK3 |= (1 << ICR1);                   // Activation de l'interruption pour Timer3
}

void setup() {
  Serial.begin(9600);                       // Initialisation de la communication série pour le débogage
  pinMode(actuateur, OUTPUT);                // Initialisation de la pin pour le PWM comme sortie

  // Conffiguration de l'ADC pour lire les températures
  ADMUX = (1 << REFS0)| currentChannel;                     // Utilise la référence de tension à 5V (REFS0)
  ADCSRA = (1 << ADEN) | (1 << ADIE) | (1 << ADPS2) | (1 << ADPS1);  // Active l'ADC et les interruptions, avec un prescaler de 64
  ADCSRB = 0;                               // Mode de fonctionnement par défaut (pas de multiplexer externe)

  initTimers();                             // Appel de la fonction pour initialiser les timers
  sei();                                    // Activation des interruptions globales (nécessaire pour que les ISR fonctionnent) permet  à l'arduino de répondre aux 
  //interruptions genréres par les timers
}

// Interruption Timer1 pour démarrer la conversion ADC .L'interruption commence  à chaque fois que le compteur atteint la valeur ORC1A
ISR(TIMER1_COMPA_vect) {
  ADCSRA |= (1 << ADSC);                    // Démarre la conversion de l'ADC (l'ADC commencera à mesurer la température) et lorsque ca termine  , une autre interruption est soulevé pour traiter ses données cad ISR(ADC_vect)
}

// Interruption ADC pour récupérer la valeur mesurée et passer au canal suivant
ISR(ADC_vect) {
  adcValues[currentChannel] = ADC;          // Stocke la valeur mesurée de l'ADC dans le tableau adcValues c'est la  valur obtenu durnat la ^récedente interruption qui est stockée ici
  currentChannel = (currentChannel + 1) % 3;  // Passe au prochain canal pour la prochaine lecture
  ADMUX = (ADMUX & 0xF8) | currentChannel;  // Modifie le registre ADMUX pour sélectionner le prochain canal
  //Revoir la logique de la convertion ADC
  float lecture_analogique_0 = (adcValues[0] / 1023.0) * 5.0; 
  float lecture_analogique_1 = (adcValues[1] / 1023.0) * 5.0; 
  float lecture_analogique_2 = (adcValues[2] / 1023.0) * 5.0; 
  float Rntc_0 = Rfixe / (5.0 /  lecture_analogique_0 - 1.0);
  float Rntc_1 = Rfixe / (5.0 /lecture_analogique_1 - 1.0);
  float Rntc_2 = Rfixe / (5.0 / lecture_analogique_2 - 1.0);
  // Calcul de la température mesurée à partir des valeurs de résistance
  float temperature_0_mesu = calculerTemperature(Rntc_0);
  float temperature_1_mesu = calculerTemperature(Rntc_1);
  float temperature_2_mesu = calculerTemperature(Rntc_2);

  // Calcul du PID pour ajuster le PWM en fonction de la température mesurée
  error = temperature_cible - temperature_2_mesu;  // Erreur entre la température cible et la température mesurée
  integral += error;                               // Somme des erreurs pour l'intégrale
  float derivative = error - prevError;            // Calcul de la dérivée (différence entre l'erreur actuelle et précédente)
  prevError = error;                               // Mise à jour de l'erreur précédente

  // Calcul du signal de commande (PWM)//Nous supposons qu'on a un PID
  tension_commande = Kp * error + Ki * integral + Kd * derivative; 
  tension_pwm = (tension_commande * 65535)/5//convertir en valeur numérique
  OCR3A = tension_pwm
}


void loop() {
  Serial.print("Température T1 : ");
  Serial.println(calculerTemperature((float)adcValues[0]));  // Affiche la température T1
  Serial.print("Température T2 : ");
  Serial.println(calculerTemperature((float)adcValues[1]));// Affiche la température T2
  Serial.print("Température T3 : ");
  Serial.println(calculerTemperature((float)adcValues[2]));// Affiche la température T3
}
