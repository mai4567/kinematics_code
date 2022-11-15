#ifndef UR3_CONTROLLER_H
#define UR3_CONTROLLER_H

#include <Eigen/Dense>
#include <QWidget>
#include <iostream>


#define PI acos(-1)


class UR3_CONTROLLER{
public:
  std::vector<double> thetax_;
  std::vector<double> theta_record_;  //用了记录上一组角度，做动画插值
  std::vector<double> dx_;
  std::vector<double> ax_;
  std::vector<double> alphax_;
  Eigen::Matrix4d T_;
  std::vector<std::vector<double>> result_;
  double nx;double ox;double ax;double px;
  double ny;double oy;double ay;double py;
  double nz;double oz;double az;double pz;
  void initDH();
  void initDH(std::vector<float> a,
              std::vector<float> d,
              std::vector<float> alpha);
  void printfDH();
  void setTheta(double theta1,double theta2,
                double theta3,double theta4,
                double theta5,double theta6);
  void setTheta(std::vector<double> thetas);
  void setAngle(std::vector<float> angles);
  void setT(Eigen::Matrix4d T);
  Eigen::Matrix4d getTx(double a,double alpha,double d,double theta);
  void positiveKinematics();  //正运动学
  Eigen::Matrix4d positiveKinematicsMat(std::vector<double> angle);
  std::vector<double> getTheta3(double theta1,double theta5,double theta6);
  double getTheta2(double theta1,double theta3,double theta6);
  double getTheta4(double theta1,double theta2,double theta3,double theta6);
  void negitiveKinematics();  //逆运动学
  void filterSlove(bool positive);  //挑选正确的解
};
#endif // UR3_CONTROLLER_H
