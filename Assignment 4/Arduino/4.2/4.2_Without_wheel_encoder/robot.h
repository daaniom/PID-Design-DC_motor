#ifndef ROBOT_H
#define ROBOT_H

/*
 * ROBOT Class
 *
 * Class incorporating the robot. This class is used to define state machines,
 * control algorithms, sensor readings,...
 * It should be interfaced with the communicator to send data to the world.
 *
 */

#include "mecotron.h" // Include MECOTRON header
#include <BasicLinearAlgebra.h> // Include BasicLinearAlgebra to make matrix manipulations easier
#include "extended_kalman_filter.h" // Include template to make extended Kalman filter implementation easier
#define PENDULUM
#include <trajectory.h> // Include trajectory, for assignment 4

class Robot : public MECOtron {
  private:

    // Class variables
    Trajectory trajectory; // define the reference trajectory object

    // Kalman filter
    Matrix<3> _xhat;       // state estimate vector
    Matrix<3,3> _Phat;     // state estimate covariance
    Matrix<1> _nu;         // innovation vector
    Matrix<1,1> _S;        // innovation covariance

    // Position controller
    Matrix<3> xref;        // reference state
    Matrix<1> desiredVelocityCart;  // control signal
    float Ts = 0.01;
    float wheel_radius = 0.0315;
    float pendulum_length = 0.155;
    float e_A[2];
    float e_B[2];
    float u_A[2];
    float u_B[2];
    float r=0;
    float omegaA=0;
    float omegaB=0;
    float numA[2]; 
    float denA[2]; 
    float numB[2];
    float denB[2];
    float theta[2];
    float thetadot=0.0;
    float x3=0.0;
  public:
    // Constructor
    Robot() { }

    void control();

    // General functions
    bool init();  // Set up the robot

    bool controlEnabled();
    bool KalmanFilterEnabled();

    void resetController();
    void resetKalmanFilter();

    void button0callback();
    void button1callback();
    void button2callback();
    void button3callback();

};

#endif // ROBOT_H
