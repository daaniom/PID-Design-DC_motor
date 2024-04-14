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

class Robot : public MECOtron {
  private:
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


  public:
    // Constructor
    Robot() { }

    void control();

    // General functions
    bool init();  // Set up the robot

    bool controlEnabled();

    void button0callback();
    void button1callback();

};

#endif // ROBOT_H
