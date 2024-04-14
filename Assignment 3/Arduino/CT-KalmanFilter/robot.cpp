/*
 * KALMAN FILTER TEMPLATE
 *
 * This is a template to get you started with the implementation of the Kalman filter
 * on your own cart.
 *
 */

#include "robot.h"

bool Robot::init() {
  MECOtron::init(); // Initialize the MECOtron
  
  // Initializing the robot's specific variables
  for(int k=0; k<2; k++){
    u_A[k] = 0.0; 
    u_B[k] = 0.0;   
    e_A[k] = 0.0;   
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
  desired_velocity(0) = 0;
  return true;
}

void Robot::control() {
  float volt_A = 0.0;
  float volt_B = 0.0;
  
  Matrix<1> desired_velocity; //control signal
  desired_velocity.Fill(0); //Initialize matrix with zeros
  
  // Assignment 3: 3.C
  //writeValue(0,xref(0));
  //writeValue(1,getFrontDistance());
  //writeValue(2, _xhat(0));
  //writeValue(3, _Phat(0));
  //writeValue(4, _L(0));
  //writeValue(5, _nu(0));
  //writeValue(6, _S(0)); 
  
  // Assignment 3: 3.e & 3.f
  _L(0) = -1*readValue(2); // Assignment 3: 3.f (10 times slower than state feedback)
  writeValue(0,xref(0));
  writeValue(1,getFrontDistance());
  writeValue(4, _Phat(0));  
  writeValue(5, _L(0));
  writeValue(6, _nu(0));
  writeValue(7, _S(0));
  writeValue(10, _xhat(0));
  
  // Kalman filtering
  if(KalmanFilterEnabled()) {   // only do this if controller is enabled (triggered by pushing 'Button 1' in QRoboticsCenter)
    // Correction step
    Matrix<1> distance_measurement;                                     // define a vector of length 1
    distance_measurement(0) = getFrontDistance();                       
    CorrectionUpdate(distance_measurement, _xhat, _Phat, _nu, _S, _L);     // do the correction step -> update _xhat, _Phat, _nu, _S
  }

  if(controlEnabled()) {   // only do this if controller is enabled (triggered by pushing 'Button 0' in QRoboticsCenter)

    // // UNCOMMENT AND COMPLETE LINES BELOW TO IMPLEMENT POSITION CONTROLLER
    float desired_position = readValue(0);      // use channel 0 to provide the constant position reference
    //desired_position multiplied by * - 1 =  sensor measures -x in our referenceframe ==> we want to input a postive value
    xref(0) = (-1*desired_position)/100;              // transform desired_position to the state reference (make sure units are consistent)
    K(0) = 2;                                    // state feedback gain K, to design
    desired_velocity = K * (xref - _xhat);      // calculate the state feedback signal, (i.e. the input for the velocity controller)

    //// UNCOMMENT AND COMPLETE LINES BELOW TO IMPLEMENT VELOCITY CONTROLLER
    r = desired_velocity(0)/wheel_radius;                    // speed to angular velocity (0.0315m = radius of wheel)
    e_A[1]=e_A[0];
    e_B[1]=e_B[0];
    u_A[1]=u_A[0];
    u_B[1]=u_B[0];
    omegaA = getSpeedMotorA();
    omegaB = getSpeedMotorB();
    e_A[0] = r-omegaA;
    e_B[0] = r-omegaB;

    u_A[0]=1*((numA[0]/denA[0])*e_A[0]+(numA[1]/denA[0])*e_A[1]-(denA[1]/denB[0])*u_A[1]);

    u_B[0]=1*((numB[0]/denB[0])*e_B[0]+(numB[1]/denB[0])*e_B[1]-(denB[1]/denB[0])*u_B[1]);

    if(u_A[0]>12){
      u_A[0]=12;
    }
    if(u_A[0]<-12){
      u_A[0]=-12;
    }
    
    if(u_B[0]>12){
      u_B[0]=12;
    }
    if(u_B[0]<-12){
      u_B[0]=-12;
    }
    
    volt_A = u_A[0];
    volt_B = u_B[0];
    
    // Send wheel speed command
    setVoltageMotorA(volt_A);
    setVoltageMotorB(volt_B);
  }
  else                      // do nothing since control is disables
  {
    xref(0)=0.0;
    desired_velocity(0) = 0.0;
    setVoltageMotorA(0.0);
    setVoltageMotorB(0.0);
    
    omegaA=getSpeedMotorA();
    omegaB=getSpeedMotorB();
    for(int k=0; k<2; k++){
    e_A[k] = 0.0;   // Set all components of the vector (float array) x to 0 as initialization
    e_B[k]= 0.0;
    u_A[k]=0.0;
    u_B[k]=0.0;
    }
  }
  
  // Kalman filtering
  if(KalmanFilterEnabled()) {   // only do this if controller is enabled (triggered by pushing 'Button 1' in QRoboticsCenter)
    // Prediction step
    PredictionUpdate(desired_velocity, _xhat, _Phat);                    // do the prediction step -> update _xhat and _Phat
  }
  // writeValue(8, _xhat(0)); // a priori state estimate
  // writeValue(9, _Phat(0)); // a priori state covariance
  // Send useful outputs to QRC
  
  // Assignment 3: 3.b
  //writeValue(0,xref(0));
  //writeValue(1,getFrontDistance());
  //writeValue(2, volt_A);
  //writeValue(3, volt_B);
  //writeValue(4, desired_velocity(0));
  //writeValue(5, getSpeedMotorA());
  //writeValue(6, getSpeedMotorB());
}

void Robot::resetController(){
  // Set all errors and control signals in the memory back to 0
  for(int k=0; k<2; k++){
    e_A[k] = 0.0;   // Set all components of the vector (float array) x to 0 as initialization
    e_B[k]= 0.0;
    u_A[k]=0.0;
    u_B[k]=0.0;
   }
}

void Robot::resetKalmanFilter() {
  // // UNCOMMENT AND MODIFIES LINES BELOW TO IMPLEMENT THE RESET OF THE KALMAN FILTER
  // // Initialize state covariance matrix
  _Phat.Fill(0);      // Initialize the covariance matrix
  // Assignment 3: 3.b & e & f
  _Phat(0,0) = 8*pow(10,-6);     // Fill the initial covariance matrix, you can change this according to your experiments
  // Assignment 3: 3.c
  //_Phat(0,0) = 8*pow(10,-7);        // uncertainty on the initial state equivalent to the uncertainty in the measurement = R-value
  //
  // // Initialize state estimate
  // Assignment 3: 3.b & c
  //_xhat(0) = -1*getFrontDistance();     // Change this according to your experiments
  // Assignment 3: 3.e & f
  _xhat(0) = -0.25;     
}

bool Robot::controlEnabled() {
  return _button_states[0];       // The control is enabled if the state of button 0 is true
}

bool Robot::KalmanFilterEnabled() {
  return _button_states[0];       // Assignment 3 3C&E: Change 1 to 0 so that controller and kalmanfilter start at same time
}

void Robot::button0callback() {
  if(toggleButton(0)) {           // Switches the state of button 0 and checks if the new state is true
    resetController();
    resetKalmanFilter();          // Assignment 3 3C&E: controller and kalmanfilter reset at same time
    message("Controller resed and enabled.");    // Display a message in the status bar of QRoboticsCenter
  }
  else {
    message("Control disabled.");
  }
}

void Robot::button1callback() {
  if(toggleButton(0)){                // Assignment 3 3C&E: Change 1 to 0 so that controller and kalmanfilter start at same time
      resetKalmanFilter();            // Reset the Kalman filter
      message("Kalman filter reset and enabled.");
  }
  else
  {
    message("Kalman filter disabled.");
  }
}
