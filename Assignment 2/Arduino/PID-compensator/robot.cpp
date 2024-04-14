
#include "robot.h"

bool Robot::init() {
  MECOtron::init(); // Initialize the MECOtron

  // Initializing the robot's specific variables
  for(int k=0; k<2; k++){
    u_A[k] = 0.0; 
    u_B[k] = 0.0;   
    e_B[k] = 0.0;   
    e_B[k] = 0.0;
    
    numA[0] = 0.7386;
    numA[1] = -0.6741;
    denA[0] = 1;
    denA[1] = -1;
  
    numB[0] = 0.7635;
    numB[1] = -0.7006;
    denB[0] = 1;
    denB[1] = -1;   
  }

  return true;
}

void Robot::control() {

  // Compute update of motor voltages if controller is enabled (triggered by
  // pushing 'Button 0' in QRoboticsCenter)
  if(controlEnabled()) {
    LED1(ON); // Set both leds to the ON state clearly indicating that we are in the control loop
    LED2(ON);
    r = readValue(0);// Fill your control law here to conditionally update the motor voltage...
    e_A[1]=e_A[0];
    e_B[1]=e_B[0];
    u_A[1]=u_A[0];
    u_B[1]=u_B[0];
    omegaA = getSpeedMotorA();
    omegaB = getSpeedMotorB();
    e_A[0] = r-omegaA;
    e_B[0] = r-omegaB;

    u_A[0]=(numA[0]/denA[0])*e_A[0]+(numA[1]/denA[0])*e_A[1]-(denA[1]/denB[0])*u_A[1];

    u_B[0]=(numB[0]/denB[0])*e_B[0]+(numB[1]/denB[0])*e_B[1]-(denB[1]/denB[0])*u_B[1];

  } else {
    
    LED1(OFF); // Set LEDS to off indicating that the control is not enabled
    LED2(OFF);
    omegaA=getSpeedMotorA();
    omegaB=getSpeedMotorB();
    for(int k=0; k<2; k++){
    e_A[k] = 0.0;   // Set all components of the vector (float array) x to 0 as initialization
    e_B[k]= 0.0;
    u_A[k]=0.0;
    u_B[k]=0.0;
    r=0;
  }
  }
  setVoltageMotorA(u_A[0]);
  setVoltageMotorB(u_B[0]);

  writeValue(0,u_A[0]);
  writeValue(1,u_B[0]);
  writeValue(2,e_A[0]);
  writeValue(3,e_B[0]);
  writeValue(4,omegaA);
  writeValue(5,omegaB);
  writeValue(6,r);
}

bool Robot::controlEnabled() {
  return _button_states[0];       // The control is enabled if the state of button 0 is true
}

void Robot::button0callback() {
  if(toggleButton(0)) {           // Switches the state of button 0 and checks if the new state is true
    message("Robot enabled.");    // Display a message in the status bar of QRoboticsCenter
  }
  else {
    message("Robot disabled.");
  }
}

void Robot::button1callback() {
  toggleButton(1);
  init();                         // Reset the MECOtron and reinitialize the Robot object
  message("Reset.");
}
