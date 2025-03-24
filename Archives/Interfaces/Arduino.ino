#include <avr/interrupt.h>
#include <math.h>
volatile uint8_t currentChannel = 0;
volatile uint16_t adcValues[3];
volatile uint8_t interruptCount = 0;
#define actuateur 5

float R25 = 10000;
float A = 0.00335401643468053;
float B = 0.000256523550896126;
float C = 0.00000260597012072052;
float D = 0.000000063292612648746;
float Rfixe[3] = {9948, 9950, 9900};
float gain[3] = {2.40735, 3.13733,4.64232};
float offset[3] = {1.52, 1.74, 1.98};

float f = 0.4;
float error[2] = {0, 0};
float u[2] = {0,0};
float y[2] = {0,0};
float Kc = 0.05818;

float coeff_u = 0.035427;
float coeff_y[2] = {0, 0.96457};
float coeff_ym[2] = {1.67107357564755,  -1.01332596656011};
float coeff_ym_a[2] = {0, 0.342252390912559};
float T3_m[2] = {25, 25};
float T3[2] = {25, 25};

volatile float lastTemperature = 0.0;    // Dernière température calculée
float consigne = 25.0;
bool asservissementActif = false;
float coefv[3] = {0,0,0};
float coef_[3] = {0,0,0};
float Consigne = 25;
float Ti_1, Kc_1, Td_1, Tf_1;
bool Ti_set = false, Kc_set = false, Consigne_set = false, Td_set = false, Tf_set = false;
String inputString = ""; // Stocker la ligne reçue
bool stringComplete = false;
void flushSerial() {
    // Vider le tampon de réception
    while (Serial.available()) {
        Serial.read();  // Lire et jeter toutes les données
    }
}

void handleSerial() {
    String input = Serial.readString();  // Lire la chaîne envoyée par MATLAB
    input.trim();  // Enlever les espaces en début et en fin de la chaîne
    if (input.startsWith("START"))
    {
      Serial.println("Asservissement commence.");
      startAsservissement();

    }
        if (input.startsWith("STOP"))
    {
      stopAsservissement();
      Serial.println("Asservissement STOPPEE."); 
    }
    
    if (input.startsWith("PID-SET")) {
        // Enlever le préfixe "PID-SET-"
        input = input.substring(8);  // La partie après "PID-SET-"

        // Extraire les valeurs séparées par des tirets '-'
        Kc_1 = input.substring(0, input.indexOf('-')).toFloat();
        input = input.substring(input.indexOf('-') + 1);

        Ti_1 = input.substring(0, input.indexOf('-')).toFloat();
        input = input.substring(input.indexOf('-') + 1);

        Td_1 = input.substring(0, input.indexOf('-')).toFloat();
        input = input.substring(input.indexOf('-') + 1);

        Tf_1 = input.substring(0, input.indexOf('-')).toFloat();
        input = input.substring(input.indexOf('-') + 1);

        Consigne = input.toFloat();  // Le dernier paramètre, pas besoin de chercher un autre tiret
       Serial.println("Commande valide."); 
    } else {
        Serial.println("Commande invalide.");
    }
}



void  stopAsservissement()
{
  asservissementActif = false;
    OCR3A = 8000;
}
void startAsservissement()
{
 asservissementActif = true; 
}
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

float calculerTemperature(float lect, float gain, float offset, float Rfixe) {
    float lecture_analogique = (lect / 1023.0) * 5.0;
    lecture_analogique = (lecture_analogique/gain) + offset;
    float Rntc = Rfixe / (4.9850/lecture_analogique - 1.0);
    float lnR = log(Rntc / R25);
    float T = 1.0 / (A + B * lnR + C * lnR * lnR + D * lnR * lnR * lnR);
    return T - 273.15;
}

int calculerCommande(float temperature_cible) {
    rotate(error, 2);
    error[0] = temperature_cible - T3[0];
    rotate(u, 2);
    rotate(y,2);
    y[0] =  coeff_u * (u[1] - 0.05) + dot(coeff_y, y ,2);
    u[0] = Kc * error[0] + y[0] + 0.05;
    int tension_commande = round(((u[0] + 1) / 2) * 15999);
    //int tension_commande = 8000;
    return constrain(tension_commande, 0, 15999);
}

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

void setup() {
    Serial.begin(115200);
    pinMode(actuateur, OUTPUT);
    ADMUX = (1 << REFS0) | currentChannel;
    ADCSRA = (1 << ADEN) | (1 << ADIE) | (1 << ADPS2) | (1 << ADPS1);
    initTimers();
    sei();
}

ISR(TIMER1_COMPA_vect) {
  if(asservissementActif){
    interruptCount++;
    if (interruptCount == 2) {
        interruptCount = 0;
        currentChannel = 0;
        ADMUX = (ADMUX & 0xF8) | currentChannel;
        ADCSRA |= (1 << ADSC);
    }
}}

ISR(ADC_vect) {
  if(asservissementActif){
    adcValues[currentChannel] = ADC;
    currentChannel++;
    if (currentChannel < 3) {
        ADMUX = (ADMUX & 0xF8) | currentChannel;
        ADCSRA |= (1 << ADSC);
    } else {
        float temperature_0_mes = calculerTemperature(adcValues[0], gain[0], offset[0], Rfixe[0]);
        float temperature_1_mes = calculerTemperature(adcValues[1], gain[1], offset[1], Rfixe[1]);
        float temperature_2_mes = calculerTemperature(adcValues[2], gain[2], offset[2], Rfixe[2]);
       // lastTemperature = temperature_2_mes;
        //Serial.println("336");
        rotate(T3_m, 2);
        T3_m[0] = temperature_2_mes;
        rotate(T3,2);
        T3[0] = dot(coeff_ym, T3_m, 2) + dot(T3, coeff_ym_a, 2);
            OCR3A = calculerCommande(Consigne);
            Serial.println(T3[0]);
        }
  }
}

void loop() {
  //flushSerial();
    if (Serial.available()) {
        handleSerial();
    }
}
  
  // Une autre logique ou action peut être ajoutée ici, mais le code ne fait rien tant qu'il n'y a pas de donnée série
