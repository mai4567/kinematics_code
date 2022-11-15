#include "./ur3_controller.h"
#include <QDebug>

void UR3_CONTROLLER::initDH(){
    //初始化DH参数列表
    dx_.push_back(151.9);
    dx_.push_back(0);
    dx_.push_back(0);
    dx_.push_back(110.4);
    dx_.push_back(83.4);
    dx_.push_back(81.4);
    ax_.push_back(0);
    ax_.push_back(0);
    ax_.push_back(243.65);
    ax_.push_back(213);
    ax_.push_back(0);
    ax_.push_back(0);
    alphax_.push_back(0);
    alphax_.push_back(PI/2);
    alphax_.push_back(0);
    alphax_.push_back(0);
    alphax_.push_back(PI/2);
    alphax_.push_back(-PI/2);
    printfDH();
}

void UR3_CONTROLLER::initDH(std::vector<float> a,
                            std::vector<float> d,
                            std::vector<float> alpha){
    //初始化DH参数列表
    for (int i=0;i<a.size();i++){
        ax_.push_back(a[i]*1000);
        dx_.push_back(d[i]*1000);
        alphax_.push_back(alpha[i]*PI/180);
    }
    printfDH();
}



void UR3_CONTROLLER::printfDH(){
    std::cout<<std::endl;
    std::cout<<"-------DH参数列表加载成功-------"<<std::endl;
    for(int i=0;i<ax_.size();i++){
        std::cout<<"d"<<i<<":"<<dx_[i]<<"  a"<<i<<":"<<ax_[i]
                 <<"  alpha"<<i<<":"<<alphax_[i];
        if (i != ax_.size()-1)
            std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"------------------------------"<<std::endl;
    std::cout<<std::endl;
}


void UR3_CONTROLLER::setTheta(double theta1,double theta2,
                              double theta3,double theta4,
                              double theta5,double theta6){
    if (thetax_.size()!=6){
        thetax_.clear();
        thetax_.push_back(theta1);
        thetax_.push_back(theta2);
        thetax_.push_back(theta3);
        thetax_.push_back(theta4);
        thetax_.push_back(theta5);
        thetax_.push_back(theta6);
    }
    else{
        thetax_[0] = theta1;
        thetax_[1] = theta2;
        thetax_[2] = theta3;
        thetax_[3] = theta4;
        thetax_[4] = theta5;
        thetax_[5] = theta6;
    }
}


void UR3_CONTROLLER::setTheta(std::vector<double> thetas){
    thetax_.clear();
    for (int i=0;i<thetas.size();i++)
        thetax_.push_back(thetas[i]);
}


void UR3_CONTROLLER::setAngle(std::vector<float> angles){
    thetax_.clear();
    for (int i=0;i<angles.size();i++)
        thetax_.push_back(angles[i]*PI/180);
}


Eigen::Matrix4d UR3_CONTROLLER::getTx(double a,double alpha,double d,double theta){
    Eigen::Matrix4d Tx;
    double T11 = cos(theta);
    double T12 = -1*sin(theta);
    double T14 = a;
    double T21 = sin(theta)*cos(alpha);
    double T22 = cos(theta)*cos(alpha);
    double T23 = -1*sin(alpha);
    double T24 = -1*d*sin(alpha);
    double T31 = sin(theta)*sin(alpha);
    double T32 = cos(theta)*sin(alpha);
    double T33 = cos(alpha);
    double T34 = d*cos(alpha);
    Tx << T11  ,T12  ,0    ,T14  ,
          T21  ,T22  ,T23  ,T24  ,
          T31  ,T32  ,T33  ,T34  ,
            0  ,  0  ,  0  ,  1;
    Eigen::Matrix4d mat1;
    mat1<<1,0,0,a,
          0,cos(alpha),-1*sin(alpha),0,
          0,sin(alpha),cos(alpha),0,
          0,0,0,1;
    Eigen::Matrix4d mat2;
    mat2<<cos(theta),-1*sin(theta),0,0,
              sin(theta),cos(theta),0,0,
              0,0,1,d,
              0,0,0,1;
    return mat1*mat2;
}

void UR3_CONTROLLER::setT(Eigen::Matrix4d T){
    T_ = T;
    nx = T_(0,0);ox = T_(0,1);ax = T_(0,2);px = T_(0,3);
    ny = T_(1,0);oy = T_(1,1);ay = T_(1,2);py = T_(1,3);
    nz = T_(2,0);oz = T_(2,1);az = T_(2,2);pz = T_(2,3);
}

//正运动学
void UR3_CONTROLLER::positiveKinematics(){
    Eigen::Matrix4d T10 = getTx(ax_[0],alphax_[0],dx_[0],thetax_[0]);
    Eigen::Matrix4d T21 = getTx(ax_[1],alphax_[1],dx_[1],thetax_[1]);
    Eigen::Matrix4d T32 = getTx(ax_[2],alphax_[2],dx_[2],thetax_[2]);
    Eigen::Matrix4d T43 = getTx(ax_[3],alphax_[3],dx_[3],thetax_[3]);
    Eigen::Matrix4d T54 = getTx(ax_[4],alphax_[4],dx_[4],thetax_[4]);
    Eigen::Matrix4d T65 = getTx(ax_[5],alphax_[5],dx_[5],thetax_[5]);
    T_ = T10*T21*T32*T43*T54*T65;
    nx = T_(0,0);ox = T_(0,1);ax = T_(0,2);px = T_(0,3);
    ny = T_(1,0);oy = T_(1,1);ay = T_(1,2);py = T_(1,3);
    nz = T_(2,0);oz = T_(2,1);az = T_(2,2);pz = T_(2,3);
    //std::cout<<"正运动学解:"<<T_<<std::endl;
}

//逆运动学
void UR3_CONTROLLER::negitiveKinematics(){
    result_.clear();
    //求解theta1
    double m = py-dx_[5]*ay;
    double n = dx_[5]*ax-px;
    double phi = atan2(m,n);
    double theta11 = std::atan2(-1*dx_[3],sqrt(m*m+n*n-dx_[3]*dx_[3]))-phi;
    double theta12 = std::atan2(-1*dx_[3],-1*sqrt(m*m+n*n-dx_[3]*dx_[3]))-phi;
    if (theta11<-1*PI){
        theta11 = theta11+2*PI;
    }
    else if(theta11>PI){
        theta11 = theta11-2*PI;
    }
    if (theta12<-1*PI){
        theta12 = theta12+2*PI;
    }
    else if (theta12>PI){
        theta12 = theta12-2*PI;
    }
    double angle11 = theta11*180/PI;
    double angle12 = theta12*180/PI;
//    std::cout<<"angle11:"<<angle11<<std::endl;
//    std::cout<<"angle12:"<<angle12<<std::endl;

    //求解theta5,一个theta1对应两个theta5
    double theta51 = acos(sin(theta11)*ax-cos(theta11)*ay);
    double theta52 = -1*acos(sin(theta11)*ax-cos(theta11)*ay);
    double theta53 = acos(sin(theta12)*ax-cos(theta12)*ay);
    double theta54 = -1*acos(sin(theta12)*ax-cos(theta12)*ay);
    double angle51 = theta51*180/PI;
    double angle52 = theta52*180/PI;
    double angle53 = theta53*180/PI;
    double angle54 = theta54*180/PI;
//    std::cout<<"angle51:"<<angle51<<std::endl;
//    std::cout<<"angle52:"<<angle52<<std::endl;
//    std::cout<<"angle53:"<<angle53<<std::endl;
//    std::cout<<"angle54:"<<angle54<<std::endl;

    //求解theta6,一个theta1对应一个theta6,由atan()计算所得的角都应该注意值域的问题
    double theta61,theta62,theta63,theta64;

    double fenzi1 = (-1*sin(theta11)*ox+cos(theta11)*oy)/sin(theta51);
    double fenmu1 = (sin(theta11)*nx-cos(theta11)*ny)/(sin(theta51));
//    std::cout<<"fenzi1:"<<fenzi1<<std::endl;
//    std::cout<<"fenmu1:"<<fenmu1<<std::endl;
    theta61 = std::atan2(fenzi1,fenmu1);

    double fenzi2 = (-1*sin(theta11)*ox+cos(theta11)*oy)/sin(theta52);
    double fenmu2 = (sin(theta11)*nx-cos(theta11)*ny)/(sin(theta52));
//    std::cout<<"fenzi2:"<<fenzi2<<std::endl;
//    std::cout<<"fenmu2:"<<fenmu2<<std::endl;
    theta62 = std::atan2(fenzi2,fenmu2);

    double fenzi3 = (-1*sin(theta12)*ox+cos(theta12)*oy)/sin(theta53);
    double fenmu3 = (sin(theta12)*nx-cos(theta12)*ny)/(sin(theta53));
    theta63 = std::atan2(fenzi3,fenmu3);

    double fenzi4 = (-1*sin(theta12)*ox+cos(theta12)*oy)/sin(theta54);
    double fenmu4 = (sin(theta12)*nx-cos(theta12)*ny)/(sin(theta54));
    theta64 = std::atan2(fenzi4,fenmu4);

    double angle61 = theta61*180/PI;
    double angle62 = theta62*180/PI;
    double angle63 = theta63*180/PI;
    double angle64 = theta64*180/PI;
//    std::cout<<"angle61:"<<angle61<<std::endl;
//    std::cout<<"angle62:"<<angle62<<std::endl;
//    std::cout<<"angle63:"<<angle63<<std::endl;
//    std::cout<<"angle64:"<<angle64<<std::endl;

    //求解theta3
    std::vector<double> theta3;  //中间变量暂存返回结果
    theta3 = getTheta3(theta11,theta51,theta61);
    double theta31 = theta3[0];
    double theta32 = theta3[1];
    double angle31 = theta31*180/PI;
    double angle32 = theta32*180/PI;
    theta3 = getTheta3(theta11,theta52,theta62);
    double theta33 = theta3[0];
    double theta34 = theta3[1];
    double angle33 = theta33*180/PI;
    double angle34 = theta34*180/PI;
    theta3 = getTheta3(theta12,theta53,theta63);
    double theta35 = theta3[0];
    double theta36 = theta3[1];
    double angle35 = theta35*180/PI;
    double angle36 = theta36*180/PI;
    theta3 = getTheta3(theta12,theta54,theta64);
    double theta37 = theta3[0];
    double theta38 = theta3[1];
    double angle37 = theta37*180/PI;
    double angle38 = theta38*180/PI;
//    std::cout<<"angle31:"<<angle31<<std::endl;
//    std::cout<<"angle32:"<<angle32<<std::endl;
//    std::cout<<"angle33:"<<angle33<<std::endl;
//    std::cout<<"angle34:"<<angle34<<std::endl;
//    std::cout<<"angle35:"<<angle35<<std::endl;
//    std::cout<<"angle36:"<<angle36<<std::endl;
//    std::cout<<"angle37:"<<angle37<<std::endl;
//    std::cout<<"angle38:"<<angle38<<std::endl;

    //求解theta2
    double theta21 = getTheta2(theta11,theta31,theta61);
    double theta22 = getTheta2(theta11,theta32,theta61);
    double theta23 = getTheta2(theta11,theta33,theta62);
    double theta24 = getTheta2(theta11,theta34,theta62);
    double theta25 = getTheta2(theta12,theta35,theta63);
    double theta26 = getTheta2(theta12,theta36,theta63);
    double theta27 = getTheta2(theta12,theta37,theta64);
    double theta28 = getTheta2(theta12,theta38,theta64);
    double angle21 = theta21*180/PI;
    double angle22 = theta22*180/PI;
    double angle23 = theta23*180/PI;
    double angle24 = theta24*180/PI;
    double angle25 = theta25*180/PI;
    double angle26 = theta26*180/PI;
    double angle27 = theta27*180/PI;
    double angle28 = theta28*180/PI;
//    std::cout<<"angle21:"<<angle21<<std::endl;
//    std::cout<<"angle22:"<<angle22<<std::endl;
//    std::cout<<"angle23:"<<angle23<<std::endl;
//    std::cout<<"angle24:"<<angle24<<std::endl;
//    std::cout<<"angle25:"<<angle25<<std::endl;
//    std::cout<<"angle26:"<<angle26<<std::endl;
//    std::cout<<"angle27:"<<angle27<<std::endl;
//    std::cout<<"angle28:"<<angle28<<std::endl;

    //求解theta4
    double theta41 = getTheta4(theta11,theta21,theta31,theta61);
    double theta42 = getTheta4(theta11,theta22,theta32,theta61);
    double theta43 = getTheta4(theta11,theta23,theta33,theta62);
    double theta44 = getTheta4(theta11,theta24,theta34,theta62);
    double theta45 = getTheta4(theta12,theta25,theta35,theta63);
    double theta46 = getTheta4(theta12,theta26,theta36,theta63);
    double theta47 = getTheta4(theta12,theta27,theta37,theta64);
    double theta48 = getTheta4(theta12,theta28,theta38,theta64);
    double angle41 = theta41*180/PI;
    double angle42 = theta42*180/PI;
    double angle43 = theta43*180/PI;
    double angle44 = theta44*180/PI;
    double angle45 = theta45*180/PI;
    double angle46 = theta46*180/PI;
    double angle47 = theta47*180/PI;
    double angle48 = theta48*180/PI;
//    std::cout<<"angle41:"<<angle41<<std::endl;
//    std::cout<<"angle42:"<<angle42<<std::endl;
//    std::cout<<"angle43:"<<angle43<<std::endl;
//    std::cout<<"angle44:"<<angle44<<std::endl;
//    std::cout<<"angle45:"<<angle45<<std::endl;
//    std::cout<<"angle46:"<<angle46<<std::endl;
//    std::cout<<"angle47:"<<angle47<<std::endl;
//    std::cout<<"angle48:"<<angle48<<std::endl;
    //压入结果容器
    result_.clear();
    std::vector<double> middle_result;
    middle_result.push_back(angle11);middle_result.push_back(angle21);middle_result.push_back(angle31);
    middle_result.push_back(angle41);middle_result.push_back(angle51);middle_result.push_back(angle61);
    result_.push_back(middle_result);
    middle_result.clear();

    middle_result.push_back(angle11);middle_result.push_back(angle22);middle_result.push_back(angle32);
    middle_result.push_back(angle42);middle_result.push_back(angle51);middle_result.push_back(angle61);
    result_.push_back(middle_result);
    middle_result.clear();

    middle_result.push_back(angle11);middle_result.push_back(angle23);middle_result.push_back(angle33);
    middle_result.push_back(angle43);middle_result.push_back(angle52);middle_result.push_back(angle62);
    result_.push_back(middle_result);
    middle_result.clear();

    middle_result.push_back(angle11);middle_result.push_back(angle24);middle_result.push_back(angle34);
    middle_result.push_back(angle44);middle_result.push_back(angle52);middle_result.push_back(angle62);
    result_.push_back(middle_result);
    middle_result.clear();

    middle_result.push_back(angle12);middle_result.push_back(angle25);middle_result.push_back(angle35);
    middle_result.push_back(angle45);middle_result.push_back(angle53);middle_result.push_back(angle63);
    result_.push_back(middle_result);
    middle_result.clear();

    middle_result.push_back(angle12);middle_result.push_back(angle26);middle_result.push_back(angle36);
    middle_result.push_back(angle46);middle_result.push_back(angle53);middle_result.push_back(angle63);
    result_.push_back(middle_result);
    middle_result.clear();

    middle_result.push_back(angle12);middle_result.push_back(angle27);middle_result.push_back(angle37);
    middle_result.push_back(angle47);middle_result.push_back(angle54);middle_result.push_back(angle64);
    result_.push_back(middle_result);
    middle_result.clear();

    middle_result.push_back(angle12);middle_result.push_back(angle28);middle_result.push_back(angle38);
    middle_result.push_back(angle48);middle_result.push_back(angle54);middle_result.push_back(angle64);
    result_.push_back(middle_result);
    middle_result.clear();

    //输出结果(可能含有不正确的解)
    for (int i=0;i<result_.size();i++){
        std::cout<<"第"<<i<<"组解:";
        for (int j=0;j<result_[i].size();j++){
            std::cout<<"   "<<result_[i][j];
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}


//求解theta3
std::vector<double> UR3_CONTROLLER::getTheta3(double theta1,double theta5,double theta6){
    std::vector<double> theta3;  // 存放运算结果
    double a = dx_[4]*(sin(theta6)*(cos(theta1)*nx+sin(theta1)*ny)+cos(theta6)*(cos(theta1)*ox+sin(theta1)*oy))-dx_[5]*(cos(theta1)*ax+sin(theta1)*ay)+cos(theta1)*px+sin(theta1)*py;
    double b = dx_[4]*(sin(theta6)*nz+cos(theta6)*oz)-dx_[5]*az+pz-dx_[0];

    double fenzi = (a*a+b*b-ax_[3]*ax_[3]-ax_[2]*ax_[2]);
    double fenmu = (2*ax_[2]*ax_[3]);

    //std::cout<<"fenzi2:"<<fenzi2/fenmu2<<std::endl;

    double theta31 = acos(fenzi/fenmu);
    double theta32 = -1*acos(fenzi/fenmu);
    theta3.clear();
    theta3.push_back(theta31);
    theta3.push_back(theta32);
    return theta3;
}


//求解theta2
double UR3_CONTROLLER::getTheta2(double theta1,double theta3,double theta6){
    double a = dx_[4]*(sin(theta6)*(cos(theta1)*nx+sin(theta1)*ny)+cos(theta6)*(cos(theta1)*ox+sin(theta1)*oy))-dx_[5]*(cos(theta1)*ax+sin(theta1)*ay)+cos(theta1)*px+sin(theta1)*py;
    double b = dx_[4]*(sin(theta6)*nz+cos(theta6)*oz)-dx_[5]*az+pz-dx_[0];
    double fenzi = (b*(ax_[3]*cos(theta3)+ax_[2])-ax_[3]*sin(theta3)*a)/(ax_[2]*ax_[2]+ax_[3]*ax_[3]+2*ax_[2]*ax_[3]*cos(theta3));
    double fenmu = (a*(ax_[3]*cos(theta3)+ax_[2])+ax_[3]*sin(theta3)*b)/(ax_[2]*ax_[2]+ax_[3]*ax_[3]+2*ax_[2]*ax_[3]*cos(theta3));
    double theta2 = std::atan2(fenzi,fenmu);
    return theta2;
}


//求解theta4
double UR3_CONTROLLER::getTheta4(double theta1,double theta2,double theta3,double theta6){
    double nx = T_(0,0);double ny = T_(1,0);double nz = T_(2,0);
    double ox = T_(0,1);double oy = T_(1,1);double oz = T_(2,1);
    double fenzi = -1*sin(theta6)*(cos(theta1)*nx+sin(theta1)*ny)-cos(theta6)*(cos(theta1)*ox+sin(theta1)*oy);
    double fenmu = nz*sin(theta6)+oz*cos(theta6);
    //因为atan(x)所求的结果落在-PI/2到PI/2上,theta2+theta3+theta4显然不一定落在此区域上
    double theta234;
    if ((fenzi>0 && fenmu<0)){
        theta234 = PI+atan(fenzi/fenmu);
    }
    else if (fenzi<0 && fenmu<0){
        theta234 = atan(fenzi/fenmu)-PI;
    }
    else{
        theta234 = atan(fenzi/fenmu);
    }
    //因为theta2+theta3+theta4的取值范围可以到[-540,540] (该范围可以在后期传入范围参数确定)
    double theta4;
    //尝试覆盖[-540,540]所有theta2+theta3+theta4的可能性
    if (2*PI+theta234-theta2-theta3<PI && 2*PI+theta234>-1*PI){
        theta4 = theta234+2*PI-theta2-theta3;
    }
    else if (-2*PI+theta234-theta2-theta3<PI && -2*PI+theta234>-1*PI){
        theta4 = theta234-2*PI-theta2-theta3;
    }
    else
        theta4 = theta234-theta2-theta3;
    return theta4;
}

Eigen::Matrix4d UR3_CONTROLLER::positiveKinematicsMat(std::vector<double> angle){
    double theta1 = angle[0]*PI/180;
    double theta2 = angle[1]*PI/180;
    double theta3 = angle[2]*PI/180;
    double theta4 = angle[3]*PI/180;
    double theta5 = angle[4]*PI/180;
    double theta6 = angle[5]*PI/180;
    Eigen::Matrix4d T10 = getTx(ax_[0],alphax_[0],dx_[0],theta1);
    Eigen::Matrix4d T21 = getTx(ax_[1],alphax_[1],dx_[1],theta2);
    Eigen::Matrix4d T32 = getTx(ax_[2],alphax_[2],dx_[2],theta3);
    Eigen::Matrix4d T43 = getTx(ax_[3],alphax_[3],dx_[3],theta4);
    Eigen::Matrix4d T54 = getTx(ax_[4],alphax_[4],dx_[4],theta5);
    Eigen::Matrix4d T65 = getTx(ax_[5],alphax_[5],dx_[5],theta6);
    return T10*T21*T32*T43*T54*T65;
}


void UR3_CONTROLLER::filterSlove(bool positive){
    if (positive == false){//如果是逆运动学
        std::vector<std::vector<double>> right_result;
        for (int i=0;i<result_.size();i++){
            bool status = true;
            for (int j=0;j<result_[i].size();j++){  //检查每个元素的取值
                if (status == false)
                    break;
                if (result_[i][j]<180&&result_[i][j]>-1*180){
                }
                else{
                    status = false;
                }
            }
            if (status == true){
                right_result.push_back(result_[i]);
            }
        }
        result_.clear();
        result_ = right_result;
    }
    else{  //如果是正运动学
        std::vector<std::vector<double>> right_result;
        std::cout<<"输出矩阵查看问题:"<<T_<<std::endl;
        for (int i=0;i<result_.size();i++){
            Eigen::Matrix4d T = positiveKinematicsMat(result_[i]);
            std::cout<<"对比矩阵:"<<T<<std::endl;
            //计算两矩阵偏差
            double dv = pow(T(0,0)-T_(0,0),2)+pow(T(0,1)-T_(0,1),2)+pow(T(0,2)-T_(0,2),2)+pow(T(0,3)-T_(0,3),2)
                       +pow(T(1,0)-T_(1,0),2)+pow(T(1,1)-T_(1,1),2)+pow(T(1,2)-T_(1,2),2)+pow(T(1,3)-T_(1,3),2)
                       +pow(T(2,0)-T_(2,0),2)+pow(T(2,1)-T_(2,1),2)+pow(T(2,2)-T_(2,2),2)+pow(T(2,3)-T_(2,3),2);
            if (dv<0.1){
//                std::cout<<"第"<<i<<"组解正确"<<std::endl;
                right_result.push_back(result_[i]);
            }
        }
        result_.clear();
        result_ = right_result;
    }
    //输出正确结果
    for (int i=0;i<result_.size();i++){
        std::cout<<"第"<<i<<"组解:";
        for (int j=0;j<result_[i].size();j++){
            std::cout<<"   "<<result_[i][j];
        }
        std::cout<<std::endl;
    }
}
