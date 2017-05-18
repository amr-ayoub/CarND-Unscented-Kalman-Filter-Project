#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here. RMSE should be less than or equal to the values [.09, .10, .40, .30]. 
  */
	
	
	VectorXd rmse = VectorXd(4);
	rmse << 0, 0, 0, 0;
	
	// checking the estimation vector size is not zero
	if (estimations.size() == 0)
	{
		std::cout << "estimation vector size is zero" << std::endl;
		return rmse;
	}
	
	// checking estimation vector size equals the ground truth vector size
	if (estimations.size() != ground_truth.size())
	{
		std::cout << "estimation vector size doesn't equal the ground truth vector size" << std::endl;
		return rmse;
	}
	
	//squared residuals
	for (int i = 0; i < estimations.size(); ++i)
	{
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array()*residual.array();
		rmse += residual;
	}
	
	//mean
	rmse = rmse / estimations.size();
	
	
	//the squared root
	rmse = rmse.array().sqrt();
	
	return rmse;
	
}
