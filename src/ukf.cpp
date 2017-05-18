#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.25;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  
  is_initialized_ = false;
  
    
  // time in us
  time_us_ = 0.0;
  
  
  // state dimension
  n_x_ = 5;
  
  // Augmented state dimension
  n_aug_ = 7;
  
  lambda_ = 3 - n_x_;
  
  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
 

  
  //weights vector
  weights_ = VectorXd(2 * n_aug_ + 1);

  
  // laser NIS
  NIS_laser_ = 0.0;
  
  // radar NIS
  NIS_radar_ = 0.0;

}


UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
  
  /////////////// INITIALIZATION ///////////////////
  
  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) ||
      (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)) 
  {
    
    if (!is_initialized_) 
    {
      //Initialize the state x_ with the first measurement
      x_ << 1, 1, 1, 1, 0.1;
      
      //covariance matrix
      P_ << 0.15,    0, 0, 0, 0,
               0, 0.15, 0, 0, 0,
               0,    0, 1, 0, 0,
               0,    0, 0, 1, 0,
               0,    0, 0, 0, 1;

      //timestamp
      time_us_ = meas_package.timestamp_;
      
      
      if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) 
      {
        x_(0) = meas_package.raw_measurements_(0);
        x_(1) = meas_package.raw_measurements_(1);
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
      {
        float roh = meas_package.raw_measurements_(0);
        float phi = meas_package.raw_measurements_(1);
        float roh_dot = meas_package.raw_measurements_(2);
        
        x_(0) = roh * cos(phi);
        x_(1) = roh * sin(phi);
        
      }
      
      // Finished initialization and return
      is_initialized_ = true;
      
      return;
      
    }
    
    
    ///////////////// PREDICTION ////////////////////
    
    //calculate elapsed time
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; 
    time_us_ = meas_package.timestamp_;
    
    // call prediction function
    Prediction(dt);
    
    
    
    ///////////////// UPDATE ////////////////
    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
      UpdateLidar(meas_package);
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
      UpdateRadar(meas_package);
    
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  //////////////// SIGMA POINTS GENERATION //////////////////
  
  //sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //square root of P
  MatrixXd A = P_.llt().matrixL();

  //lambda
  lambda_ = 3 - n_x_;

  //first column
  Xsig.col(0) = x_;
  
  //sigma points
  
  for (int point = 0; point < n_x_; point++)
  {
    Xsig.col(point + 1) = x_ + sqrt(lambda_ + n_x_) * A.col(point);
    Xsig.col(point + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(point);
  }
  
  
  //////////////// SIGMA POINTS Augmuntation //////////////////

  //augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  
  //augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  
  //augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //lambda
  lambda_ = 3 - n_aug_;
  
  //augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  //augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_*std_a_;
  P_aug(6, 6) = std_yawdd_*std_yawdd_;
  
  //square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //augmented sigma points
  Xsig_aug.col(0) = x_aug;
  
  for (int point = 0; point< n_aug_; point++)
  {
    Xsig_aug.col(point + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(point);
    Xsig_aug.col(point + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(point);
  }
  
  

  //////////////// SIGMA POINTS Prediction //////////////////

  for (int point = 0; point < 2 * n_aug_ + 1; point++)
  {
    
    double p_x      = Xsig_aug(0, point);
    double p_y      = Xsig_aug(1, point);
    double v        = Xsig_aug(2, point);
    double yaw      = Xsig_aug(3, point);
    double yawd     = Xsig_aug(4, point);
    double nu_a     = Xsig_aug(5, point);
    double nu_yawdd = Xsig_aug(6, point);
    
    //predicted state values
    double px_p, py_p;
    
    //check division by zero
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else 
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }
    
    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;
    
    //noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;
    
    //assign predicted sigma point
    Xsig_pred_(0, point) = px_p;
    Xsig_pred_(1, point) = py_p;
    Xsig_pred_(2, point) = v_p;
    Xsig_pred_(3, point) = yaw_p;
    Xsig_pred_(4, point) = yawd_p;
    
    
  }
  
  
  //////////////// SIGMA POINTS conversion to mean and covariance //////////////////
  
  //weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  
  for (int i = 1; i < 2 * n_aug_ + 1; i++)
  {
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }
  
  //mean state
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
  { 
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  
  //covariance matrix
  P_.fill(0.0); 
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
    
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  
  //measurement vector
  VectorXd z = meas_package.raw_measurements_;
  
  //measurement dimensions
  int n_z = 2;
  
  //sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  //sigma points transformation
  for (int point = 0; point < 2 * n_aug_ + 1; point++)
  {
    double p_x = Xsig_pred_(0, point);
    double p_y = Xsig_pred_(1, point);
    Zsig(0, point) = p_x;
    Zsig(1, point) = p_y;
  }
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  for (int point = 0; point < 2 * n_aug_ + 1; point++)
    z_pred = z_pred + weights_(point) * Zsig.col(point);
    
  
  //measurement covariance matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  
  for (int point = 0; point < 2 * n_aug_ + 1; point++)
  {
    VectorXd z_diff = Zsig.col(point) - z_pred;
    S = S + weights_(point) * z_diff * z_diff.transpose();
  }
  
  //measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_*std_laspx_, 0,
       0, std_laspy_*std_laspy_;
  
  S = S + R;
  
  
  //cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  for (int point = 0; point < 2 * n_aug_ + 1; point++)
  {
    VectorXd z_diff = Zsig.col(point) - z_pred;
    VectorXd x_diff = Xsig_pred_.col(point) - x_;
    Tc = Tc + weights_(point) * x_diff * z_diff.transpose();
  }
  
  //Kalman gain
  MatrixXd K = Tc * S.inverse();
  
  //NIS
  VectorXd z_diff = z - z_pred;
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  
  //state mean
  x_ = x_ + K * z_diff;
  
  //covariance matrix
  P_ = P_ - K*S*K.transpose();
  
 
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  
  //measurement vector
  VectorXd z = meas_package.raw_measurements_;
  
  //measurement dimensions
  int n_z = 3;
  
  //sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  //sigma points into measurement space
  for (int point = 0; point < 2 * n_aug_ + 1; point++)
  {
    
    double p_x = Xsig_pred_(0, point);
    double p_y = Xsig_pred_(1, point);
    double v   = Xsig_pred_(2, point);
    double yaw = Xsig_pred_(3, point);
    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;
    
    //measurements
    //roh
    Zsig(0, point) = sqrt(p_x*p_x + p_y*p_y);                        
    //phi
    Zsig(1, point) = atan2(p_y, p_x);
    //r_dot
    Zsig(2, point) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);
  
  }
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) 
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  
  
  //measurement covariance matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  
  for (int point = 0; point < 2 * n_aug_ + 1; point++)
  {
    VectorXd z_diff = Zsig.col(point) - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
    
    S = S + weights_(point) * z_diff * z_diff.transpose();
    
    
  }
  
  //measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_*std_radr_,                       0,                     0,
                         0, std_radphi_*std_radphi_,                     0,
                         0,                       0, std_radrd_*std_radrd_;
  
  S = S + R;
  
  
  //cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  
  for (int point = 0; point < 2 * n_aug_ + 1; point++)
  {
    VectorXd z_diff = Zsig.col(point) - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
    
    VectorXd x_diff = Xsig_pred_.col(point) - x_;
    
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;
    
    Tc = Tc + weights_(point) * x_diff * z_diff.transpose();
  }
  
  //Kalman gain
  MatrixXd K = Tc * S.inverse();
  
  VectorXd z_diff = z - z_pred;
  
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
  
  //NIS
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  
  //state mean
  x_ = x_ + K * z_diff;
  
  //covariance matrix
  P_ = P_ - K*S*K.transpose();
  
  

  
}
