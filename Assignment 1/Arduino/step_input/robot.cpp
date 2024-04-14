/*
 * MECOTRON TUTORIAL
 *
 * This is a template to get you started in the course of the tutorial on the
 * control theory platforms, a.k.a. the MECOtrons.s
 * The tasks of the tutorial session will guide you through this template and
 * ask you to make use of the platform's capabilities in a step-by-step fashion.
 *
 * Every function in this template comes with an opening comment that describes
 * its purpose and functionality. Please also pay attention to the remarks that
 * are made in comment blocks.
 *
 */

#include "robot.h"

bool Robot::init() {
  MECOtron::init(); // Initialize the MECOtron

  // Initializing the robot's specific variables
  for(int k=0; k<2; k++){
    x[k] = 0.0;   // Set all components of the vector (float array) x to 0 as initialization
    y[k] = 0.0;
  }

  return true;
}


int index = 0;

void Robot::control() {
  
  // Compute update of motor voltages if controller is enabled (triggered by
  // pushing 'Button 0' in QRoboticsCenter)
  if(controlEnabled()) {
    // Fill your control law here to conditionally update the motor voltage...
    LED1(ON);
    LED2(OFF);
    // apply step input
    uA = 4;
    uB = 4; 
    setVoltageMotorA(uA); // Apply 2.0 volts to motor B if the control is enabled
    setVoltageMotorB(uB);
    writeValue(1, uA);
    writeValue(2, getSpeedMotorA());
    writeValue(3, getPositionMotorA());
    writeValue(4, uB);
    writeValue(5, getSpeedMotorB());
    writeValue(6, getPositionMotorB());
   
  } 
  else {
    // If the controller is disabled, you might want to do something else...
    LED1(OFF);
    LED2(ON);
    uA = 0.0;
    uB = 0.0;
    setVoltageMotorA(uA); // Apply 0.0 volts to motor A if the control is disabled
    setVoltageMotorB(uB); // Apply 0.0 volts to motor B if the control is disabled
    writeValue(1, uA);
    writeValue(2, getSpeedMotorA());
    writeValue(3, getPositionMotorA());
    writeValue(4, uB);
    writeValue(5, getSpeedMotorB());
    writeValue(6, getPositionMotorB());
   
  }

  float va = getSpeedMotorA();    // Get the wheel speed of motor A (in radians/second)
  x[1] = x[0]; x[0] = va;         // Memorize the last two samples of the speed of motor A (in fact, a shift register)
  float vb = getSpeedMotorB();
  y[1] = y[0]; y[0] = vb;
  
  float k = readValue(0); // Read the value you set on QRoboticsCenter's channel 0
  writeValue(0, k);       // Send the value of variable k to QRoboticsCenter's channel 0
  
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
    period = 0;
    state = 0;
    index = 0;
  }
}

void Robot::button1callback() {
  toggleButton(1);
  init();                         // Reset the MECOtron and reinitialize the Robot object
  message("Reset.");
}

void Robot::button2callback(){
  if(toggleButton(2)){
    message("button 2 pressed.");
  }
  else{
    message("button 2 unpressed.");
  }
}
