#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "kalman_filter.h"
#include "tools.h"

class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  /**
  * Calculate the linear measurement matrix given the current predicted state
  */
  MatrixXd CalculateJacobian(const VectorXd& x_state) {

  	MatrixXd Hj(3,4);
  	//recover state parameters
  	float px = x_state(0);
  	float py = x_state(1);
  	float vx = x_state(2);
  	float vy = x_state(3);
      float p_norm2 = (px * px + py * py);
      float p_norm  = sqrt(p_norm2);
      float p_norm3 = p_norm2 * p_norm;

  	//check division by zero
  	if (abs(p_norm) < 0.000001) {
  	    std::cout << "division by zero" << std::endl;
  	    return Hj;
  	}

  	//compute the Jacobian matrix
      Hj << px / p_norm, py / p_norm, 0, 0,
           -py / p_norm2, px / p_norm2, 0, 0,
           py * (vx * py - vy * px) / p_norm3, px * (vy * px - vx * py) / p_norm3,
           px / p_norm, py / p_norm;


  	return Hj;
  }

  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;

private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd Hj_;
};

#endif /* FusionEKF_H_ */
