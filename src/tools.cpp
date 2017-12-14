#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  // ... your code here

  if(estimations.size() != ground_truth.size()
  		|| estimations.size() == 0){
  	cout << "Invalid estimation or ground_truth data" << endl;
  	return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
        // ... your code here
        if (estimations[i].cols() == ground_truth[i].cols() && estimations[i].rows() == ground_truth[i].rows()) {
            VectorXd residual = estimations[i] - ground_truth[i];
            //std::cout << "residual " << residual << std::endl;
            residual = residual.array().square();
            //std::cout << "residual2 " << residual << std::endl;
            rmse += residual;
        }
  }

  //std::cout << rmse << std::endl;

  //calculate the mean
  // ... your code here
    rmse /= (estimations.size());

  //calculate the squared root
  // ... your code here
    rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);
    float p_norm2 = (px * px + py * py);
    float p_norm  = sqrt(p_norm2);
    float p_norm3 = p_norm2 * p_norm;

	//TODO: YOUR CODE HERE

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
