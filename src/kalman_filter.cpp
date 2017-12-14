#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  //std::cout << "State before prediction: " << x_ << std::endl;
  //std::cout << "F_: " << F_ << std::endl;
  x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
  //std::cout << "State after prediction: " << x_ << std::endl;

}

void KalmanFilter::Update(const VectorXd &z) {

  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  double x, y_, x_dot, y_dot, rho, phi, rho_dot;

  x = x_(0);
  y_ = x_(1);
  x_dot = x_(2);
  y_dot = x_(3);

  rho = sqrt(x * x + y_ * y_);
  phi = atan2(y_, x);
  rho_dot = (x * x_dot + y_dot * y_) / rho;

  //std::cout << "Correcting for measurement" << std::endl;

  VectorXd z_pred(3);
  z_pred << rho, phi, rho_dot;

  //std::cout << "zpred: " << z_pred << std::endl << "z: " << z << std::endl;

	VectorXd y = z - z_pred;
  y(1) = fmod(y(1), (2*3.14159265));
	MatrixXd Ht = H_.transpose();
  //std::cout << "H_\n" << H_ << "\nP_\n" << P_ << "\nR_\n" << R_ << std::endl;
	MatrixXd S = H_ * P_ * Ht + R_;
  //std::cout << "Calculated measurement space variance" << std::endl;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}
