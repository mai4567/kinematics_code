function [T] = positiveKinematics(angle)
theta1 = angle(1)*3.14159/180;theta3 = angle(3)*3.14159/180;theta5 = angle(5)*3.14159/180;
theta2 = angle(2)*3.14159/180;theta4 = angle(4)*3.14159/180;theta6 = angle(6)*3.14159/180;
d1 = 151.9; a1 = 0; alpha1 = 0;
d2 = 0; a2 = 0; alpha2 = pi/2;
d3 = 0; a3 = 243.65; alpha3 = 0;
d4 = 110.4; a4 = 213; alpha4 = 0;
d5 = 83.4; a5 = 0; alpha5 = pi/2;
d6 = 81.4; a6 = 0; alpha6 = -pi/2;
T10 = vpa(getT(a1,alpha1,d1,theta1));
T21 = vpa(getT(a2,alpha2,d2,theta2));
T32 = vpa(getT(a3,alpha3,d3,theta3));
T43 = vpa(getT(a4,alpha4,d4,theta4));
T54 = vpa(getT(a5,alpha5,d5,theta5));
T65 = vpa(getT(a6,alpha6,d6,theta6));
T = vpa(T10*T21*T32*T43*T54*T65);
end

