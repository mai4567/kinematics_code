function [theta4] = getTheta4(theta1,theta2,theta3,theta6,T60)
nx = T60(1,1);ox = T60(1,2);ax = T60(1,3);px = T60(1,4);
ny = T60(2,1);oy = T60(2,2);ay = T60(2,3);py = T60(2,4);
nz = T60(3,1);oz = T60(3,2);az = T60(3,3);pz = T60(3,4);
d1 = 151.9; a1 = 0; alpha1 = 0;
d2 = 0; a2 = 0; alpha2 = pi/2;
d3 = 0; a3 = 243.65; alpha3 = 0;
d4 = 110.4; a4 = 213; alpha4 = 0;
d5 = 83.4; a5 = 0; alpha5 = pi/2;
d6 = 81.4; a6 = 0; alpha6 = -pi/2;
fenzi = -sin(theta6)*(cos(theta1)*nx+sin(theta1)*ny)-cos(theta6)*(cos(theta1)*ox+sin(theta1)*oy);
fenmu = nz*sin(theta6)+oz*cos(theta6);
theta234 = atan(fenzi/fenmu);
if (fenzi>0 && fenmu<0)
    theta234 = atan(fenzi/fenmu)+3.14159;
end
if(fenzi<0 && fenmu<0)
    theta234 = atan(fenzi/fenmu)-3.14159;
end
if (theta234-2*3.14159-theta2-theta3>-3.14159 && theta234-2*3.14159-theta2-theta3<3.14159)
    theta4 = theta234-2*3.14159-theta2-theta3;
elseif(theta234+2*3.14159-theta2-theta3>-3.14159 && theta234+2*3.14159-theta2-theta3<3.14159)
    theta4 = theta234+2*3.14159-theta2-theta3;
else
    theta4 = theta234-theta2-theta3;
end
end

