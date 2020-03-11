#include <SoftwareSerial.h>
char rawBufer;
String bufer;
String number_bufer = "0";
int point;
int left_speed;
int right_speed;
int delta = 60;
int max_velocity = 120;
// Motor A
int ENA = 6;
int IN1 = 8;
int IN2 = 7;

// Motor B
int ENB = 5;
int IN3 = 4;
int IN4 = 2;

SoftwareSerial mySerial(10, 11); // RX, TX

void setup() {
  // Open serial communications and wait for port to open:
  Serial.begin(9600);
  Serial.println("Goodnight moon!");
  mySerial.begin(38400);
  pinMode (ENA, OUTPUT);
  pinMode (ENB, OUTPUT);
  pinMode (IN1, OUTPUT);
  pinMode (IN2, OUTPUT);
  pinMode (IN3, OUTPUT);
  pinMode (IN4, OUTPUT);
}

void loop() { // run over and over
  if (mySerial.available()) {
    rawBufer = mySerial.read();
    bufer = String(rawBufer);
    if(bufer == "."){
      point = number_bufer.toInt();
      left_speed = max_velocity + 2*(((float)point/100)-0.5)*delta;
      right_speed = max_velocity - 2*(((float)point/100)-0.5)*delta;
      number_bufer = ""; 
    }else{
      number_bufer += bufer;
    }
  }
  //Direccion motor A
  digitalWrite (IN1, LOW);
  digitalWrite (IN2, HIGH);
  analogWrite (ENA, left_speed); //Velocidad motor A
  //Direccion motor B
  digitalWrite (IN3, HIGH);
  digitalWrite (IN4, LOW);
  analogWrite (ENB, right_speed); //Velocidad motor B
}
