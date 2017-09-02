#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check validity of inputs
	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		cout << "Invalid estimation or ground truth data." << endl;
		return rmse;
	}

	// Accumulate squared residuals
	for (unsigned int i=0; i < estimations.size(); i++) {

		VectorXd residual = estimations[i] - ground_truth[i];

		// coefficient-wise multiplication
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	// Calculate the mean
	rmse = rmse/estimations.size();

	// Sqrt
	rmse = rmse.array().sqrt();

	// Done
	if (rmse(0) > 0.11 || rmse(1) > 0.11) {
		cout << "x/y RMSE > 0.11: " << rmse(0) << ", " << rmse(1) << endl;
	}
	if (rmse(2) > 0.52 || rmse(3) > 0.52) {
		cout << "vx/vy RMSE > 0.52: " << rmse(2) << ", " << rmse(3) << endl;
	}

	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float c1 = px*px + py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if (fabs(c1) < 0.0001) {
		cout << "Error: CalculateJacobian() - Division by Zero or Close to Zero" << endl;
		return Hj;
	}

	/*
	//compute the Jacobian MatrixXd
	float term_0_0 = px/(sqrt(px*px + py*py));
	float term_0_1 = py/(sqrt(px*px + py*py));
	float term_0_2 = 0.0;
	float term_0_3 = 0.0;
	
	float term_1_0 = -py/(px*px + py*py);
	float term_1_1 = px/(px*px + py*py);
	float term_1_2 = 0.0;
	float term_1_3 = 0.0;
	
	float term_2_0 = (py*(vx*py-vy*px))/(pow((px*px+py*py),3/2));
	float term_2_1 = (px*(vy*px-vx*py))/(pow((px*px+py*py), 3/2));
	float term_2_2 = px/(sqrt(px*px + py*py));
	float term_2_3 = py/(sqrt(px*px + py*py));
	
	Hj << term_0_0, term_0_1, term_0_2, term_0_3,
	        term_1_0, term_1_1, term_1_2, term_1_3,
	        term_2_0, term_2_1, term_2_2, term_2_3;
	*/

	Hj << (px/c2), (py/c2), 0, 0,
		  -(py/c1), (px/c1), 0, 0,
		  py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
	
	return Hj;

}
