function [theta2] = getTheta2(theta1,theta3,theta6,T60)
theta2 = 0;
if (isreal(theta3))
    nx = T60(1,1);ox = T60(1,2);ax = T60(1,3);px = T60(1,4);
    ny = T60(2,1);oy = T60(2,2);ay = T60(2,3);py = T60(2,4);
    nz = T60(3,1);oz = T60(3,2);az = T60(3,3);pz = T60(3,4);
    d1 = 151.9; a1 = 0; alpha1 = 0;
    d2 = 0; a2 = 0; alpha2 = pi/2;
    d3 = 0; a3 = 243.65; alpha3 = 0;
    d4 = 110.4; a4 = 213; alpha4 = 0;
    d5 = 83.4; a5 = 0; alpha5 = pi/2;
    d6 = 81.4; a6 = 0; alpha6 = -pi/2;
    a = d5*sin(theta6)*(cos(theta1)*nx+sin(theta1)*ny)+d5*cos(theta6)*(cos(theta1)*ox+sin(theta1)*oy)-d6*(cos(theta1)*ax+sin(theta1)*ay)+cos(theta1)*px+sin(theta1)*py;
    b = d5*(nz*sin(theta6)+oz*cos(theta6))-az*d6+pz-d1;
    fenzi = (b*(a4*cos(theta3)+a3)-a4*sin(theta3)*a)/(a3*a3+a4*a4+2*a3*a4*cos(theta3));
    fenmu = (a*(a4*cos(theta3)+a3)+a4*sin(theta3)*b)/(a3*a3+a4*a4+2*a3*a4*cos(theta3));
    theta2 = atan2(fenzi,fenmu);
end
end

