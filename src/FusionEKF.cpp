#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;
    ekf_.P_ = MatrixXd(4,4);
    ekf_.P_  << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 10000, 0, // very unsure of the velocities
                0, 0, 0, 10000;

    //the initial transition matrix Q_
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ <<  0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double angle, radius;
      radius = measurement_pack.raw_measurements_(0);
      angle  = measurement_pack.raw_measurements_(1);
      ekf_.x_ << cos(angle) * radius, sin(angle) * radius, 0, 0;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      double x, y;
      x = measurement_pack.raw_measurements_(0);
      y  = measurement_pack.raw_measurements_(1);
      ekf_.x_ << x, y, 0, 0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    std::cout << "Initialized " << std::endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //std::cout << "Updating pose" << std::endl;

   //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.;	//dt - expressed in seconds
  float noise_ax = 9;
  previous_timestamp_ = measurement_pack.timestamp_;
  //1. Modify the F matrix so that the time is integrated
  ekf_.F_(0,2) = ekf_.F_(1,3) = dt;
  //2. Set the process covariance matrix Q
	ekf_.Q_(0,0) = ekf_.Q_(1,1) = (dt*dt*dt*dt) / 4 * noise_ax;
	ekf_.Q_(2,0) = ekf_.Q_(3,1) = (dt*dt*dt) / 2 * noise_ax;
	ekf_.Q_(0,2) = ekf_.Q_(1,3) = (dt*dt*dt) / 2 * noise_ax;
	ekf_.Q_(2,2) = ekf_.Q_(3,3) = (dt*dt) * noise_ax;

  ekf_.Predict();

  //std::cout << "Prediction step done" << std::endl;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    //std::cout << "UpdateEKF step" << std::endl;
    //std::cout << "Radar" << std::endl;
    // std::cout << "Hj: "
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    //std::cout << "UpdateEKF step done" << std::endl;
  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    //std::cout << "Update step" << std::endl;
    //std::cout << "Laser" << std::endl;
    ekf_.Update(measurement_pack.raw_measurements_);
    //std::cout << "Update step done" << std::endl;
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
