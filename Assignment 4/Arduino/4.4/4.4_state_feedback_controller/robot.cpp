/*
 * EXTENDED KALMAN FILTER TEMPLATE
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
    e_B[k] = 0.0;   
    e_B[k] = 0.0;
    theta[k]=0.0;

    numA[0] = 0.7386;
    numA[1] = -0.6741;
    denA[0] = 1;
    denA[1] = -1;
  
    numB[0] = 0.7635;
    numB[1] = -0.7006;
    denB[0] = 1;
    denB[1] = -1;
  }
  desiredVelocityCart(0) = 0.0;
  return true;
}

void Robot::control() {

  float volt_A = 0.0;
  float volt_B = 0.0;
  float desiredVelocityMotorA = 0.0;
  float desiredVelocityMotorB = 0.0;
  Matrix<1> desiredVelocityCart;  // control signal
  desiredVelocityCart.Fill(0); //Initialize matrix with zeros
  Matrix<2> measurements;

// Kalman filtering
  if(KalmanFilterEnabled()) { // only do this if Kalman Filter is enabled (triggered by pushing 'Button 1' in QRoboticsCenter)
    // // UNCOMMENT AND MODIFY LINES BELOW TO IMPLEMENT THE KALMAN FILTER
    // // Correction step
    measurements(0) = (getPositionMotorA()+getPositionMotorB())/2*wheel_radius; //transform encoders measurement (getPositionMotorA() and getPositionMotorB()) to cart position measurement
    measurements(1) = getPendulumAngle();
    CorrectionUpdate(measurements, _xhat, _Phat, _nu, _S);     // do the correction step -> update _xhat, _Phat, _nu, _S
  
    // // Useful outputs to QRC for assignment questions
    // writeValue(1, _xhat(0));
    // writeValue(2, _xhat(1));
    // writeValue(3, _xhat(2));
    // writeValue(3, _Phat(0,0));
    // writeValue(4, _Phat(1,0));
    // writeValue(4, _Phat(1,0));
    // writeValue(6, _Phat(1,1));
    // writeValue(7, _Phat(2,0));
    // writeValue(7, _Phat(2,1));
    // writeValue(9, _Phat(2,2));
    // writeValue(10, measurements(0));
    // writeValue(11, measurements(1));
  }

  if(controlEnabled()) {   // only do this if controller is enabled (triggered by pushing 'Button 0' in QRoboticsCenter)

    // UNCOMMENT AND COMPLETE LINES BELOW TO IMPLEMENT THE FEEDFORWARD INPUTS (ASSIGNMENT 4.2, no state feedback here)
    // COMMENT OR REMOVE LINES BELOW ONCE YOU IMPLEMENT THE POSITION STATE FEEDBACK CONTROLLER
    // Compute the feedforward input of the cart
    // The feedforward is here returned by the built-in trajectory: trajectory.v()
    //desiredVelocityCart(0) = trajectory.v();  //desired forward velocity of the cart (in m/s)
    // The trajectory must be started by pushing 'Button 2' in QRoboticsCenter, otherwise will return zeros
    // after any experiment the trajectory must be reset pushing 'Button 3' in QRoboticsCenter
    //
    // // apply the static transformation between velocity of the cart and velocity of the motors
    //desiredVelocityMotorA = desiredVelocityCart(0)/wheel_radius;        // calculate the angular velocity of the motor A using desiredVelocityCart
    // desiredVelocityMotorB = ?;        // calculate the angular velocity of the motor B using desiredVelocityCart


    // // UNCOMMENT AND COMPLETE LINES BELOW TO IMPLEMENT POSITION CONTROLLER (ASSIGNMENT 4.4)
    // // step reference in the desired pendulum mass position
    ref_position(0)= readValue(0);  
    //
    // // State feedback controller
    float arrayKfb[1][3]{{ 3.2100  , 0.5139 ,   0.2805}};  // state feedback gain Kfb, to design
    Matrix<1, 3> Kfb = arrayKfb;
    //
    // // Compute feedback signal ufb = -Kfb*x
    Matrix<1> ufb = -Kfb * _xhat;
    //
    // // Feedforward controller
    float arrayKff[1][1] = {{3.2100}};  // feedforward gain Kff, to design
    Matrix<1> Kff = arrayKff;
    //
    
    Matrix<1> uff = Kff*ref_position(0);

    desiredVelocityCart(0) = uff(0) + ufb(0);
    // Calculation of State 3
    theta[1]=theta[0];
    theta[0]=getPendulumAngle();
    thetadot=(theta[0]-theta[1])/Ts;
    x3=thetadot*pendulum_length+getSpeedMotorA()*wheel_radius*cos(theta[0]);

   //// UNCOMMENT AND COMPLETE LINES BELOW TO IMPLEMENT VELOCITY CONTROLLER
    r = desiredVelocityCart(0)/wheel_radius;                    // speed to angular velocity (0.0315m = radius of wheel)
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
    desiredVelocityCart(0) = 0.0;
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
  if(KalmanFilterEnabled()) { // only do this if Kalman Filter is enabled (triggered by pushing 'Button 1' in QRoboticsCenter)
    // Prediction step
    PredictionUpdate(desiredVelocityCart, _xhat, _Phat);                    // do the prediction step -> update _xhat and _Phat
  }

  // Send useful outputs to QRC
  // to check functioning of trajectory and feedforward
  //writeValue(0, trajectory.v());
  //writeValue(1, trajectory.X());
  //writeValue(3, getSpeedMotorA());
  //writeValue(4, getSpeedMotorB());
  //writeValue(5, volt_A);
 // writeValue(6, volt_B);

  // // Useful outputs to QRC for assignment questions
  writeValue(0, desiredVelocityCart(0));
  writeValue(1, getPendulumAngle());
  writeValue(2, _xhat(0));
  writeValue(3, _xhat(1));
  writeValue(4, _xhat(2));
  writeValue(4, _xhat(2));
  writeValue(5,measurements(0)+pendulum_length*sin(measurements(1)));
  writeValue(6,measurements(0));
  writeValue(7,measurements(1));
  writeValue(8,volt_A);
  writeValue(9,volt_B);


  //triggers the trajectory to return the next values during the next cycle
  trajectory.update();
}

void Robot::resetController(){
  for(int k=0; k<2; k++){
    u_A[k] = 0.0; 
    u_B[k] = 0.0;   
    e_B[k] = 0.0;   
    e_B[k] = 0.0;
    theta[k]=0.0;
  }
}

void Robot::resetKalmanFilter() {
  // // UNCOMMENT AND MODIFY LINES BELOW TO IMPLEMENT THE RESET OF THE KALMAN FILTER
  // // Initialize state covariance matrix
   _Phat.Fill(0);       // Initialize the covariance matrix
   _Phat(0,0) = pow(10,-6);      // Fill the initial covariance matrix, you can change this according to your experiments
   _Phat(1,1) = pow(10,-6);
   _Phat(2,2) = pow(10,-2);
  //
  // // Initialize state estimate
   _xhat(0) = (getPositionMotorA()+getPositionMotorA())/2*wheel_radius;       // Change this according to your experiments
   _xhat(1) = -getPendulumAngle();
   _xhat(2) = 0;
  //
  // // Reset innovation and its covariance matrix
   _S.Fill(0);
  _nu.Fill(0);
}

bool Robot::controlEnabled() {
  return _button_states[0];       // The control is enabled if the state of button 0 is true
}

bool Robot::KalmanFilterEnabled() {
  return _button_states[0];       // Assignment 4.4: Change 1 to 0 so that controller and kalmanfilter start at same time
}

void Robot::button0callback() {
  if(toggleButton(0)) {           // Switches the state of button 0 and checks if the new state is true
    resetController();
    message("Controller reset and enabled.");    // Display a message in the status bar of QRoboticsCenter
  }
  else {
    message("Control disabled.");
  }
}

void Robot::button1callback() {
  if(toggleButton(0)){                // Assignment 4.4: Change 1 to 0 so that controller and kalmanfilter start at same time
      resetKalmanFilter();            // Reset the Kalman filter
      message("Kalman filter reset and enabled.");
  }
  else
  {
    message("Kalman filter disabled.");
  }
}

void Robot::button2callback() {
  if(toggleButton(2)) {
    trajectory.start();
    message("Trajectory started/resumed.");
  } else {
    trajectory.stop();
    message("Trajectory stopped.");
  }
}

void Robot::button3callback() {
     _button_states[2] = 0;
    trajectory.reset();
    message("Trajectory reset.");
}
