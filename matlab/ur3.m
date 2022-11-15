%%
syms pi;
d1 = 151.9; a1 = 0; alpha1 = 0;
d2 = 0; a2 = 0; alpha2 = pi/2;
d3 = 0; a3 = 243.65; alpha3 = 0;
d4 = 110.4; a4 = 213; alpha4 = 0;
d5 = 83.4; a5 = 0; alpha5 = pi/2;
d6 = 81.4; a6 = 0; alpha6 = -pi/2;
%求解运动学反解
%推导一般式
syms theta1_s theta2_s theta3_s theta4_s theta5_s theta6_s;
syms d1_s d2_s d3_s d4_s d5_s d6_s;
syms a1_s a2_s a3_s a4_s a5_s a6_s;
T10_s = vpa(getT(0,alpha1,d1_s,theta1_s));
T21_s = vpa(getT(0,alpha2,0,theta2_s));
T32_s = vpa(getT(a3_s,alpha3,0,theta3_s));
T43_s = vpa(getT(a4_s,alpha4,d4_s,theta4_s));
T54_s = vpa(getT(0,alpha5,d5_s,theta5_s));
T65_s = vpa(getT(0,alpha6,d6_s,theta6_s));
%获取一些中间结果
T64_s = vpa(T54_s*T65_s);
T63_s = vpa(T43_s*T64_s);
T31_s = vpa(T21_s*T32_s);
T61_s = vpa(T31_s*T63_s);
T10_s_ = vpa(inv(T10_s));
T65_s_ = vpa(inv(T65_s));
T54_s_ = vpa(inv(T54_s));
T51_s = vpa(T21_s*T32_s*T43_s*T54_s);
%求解具体值
%---------------------------------------------------------------------
angles = [-150,150,-20,30,-30,-120];
%---------------------------------------------------------------------
T10 = vpa(getT(a1,alpha1,d1,angles(1)*3.14159/180));
T21 = vpa(getT(a2,alpha2,d2,angles(2)*3.14159/180));
T32 = vpa(getT(a3,alpha3,d3,angles(3)*3.14159/180));
T43 = vpa(getT(a4,alpha4,d4,angles(4)*3.14159/180));
T54 = vpa(getT(a5,alpha5,d5,angles(5)*3.14159/180));
T65 = vpa(getT(a6,alpha6,d6,angles(6)*3.14159/180));
% T60 = vpa(T10*T21*T32*T43*T54*T65);

T60 =   [0.84988,-0.4609,-0.255464,-265,
        -0.0894,0.351527,-0.931892,-200,
        0.51937,0.814819,0.257537,590.49,
        0,0,0,1];

nx = T60(1,1);ox = T60(1,2);ax = T60(1,3);px = T60(1,4);
ny = T60(2,1);oy = T60(2,2);ay = T60(2,3);py = T60(2,4);
nz = T60(3,1);oz = T60(3,2);az = T60(3,3);pz = T60(3,4);

%求解theta1
m = py - d6*ay; n = d6*ax-px;
phi = atan2(m,n);
theta11 = atan2(-d4,sqrt(m^2+n^2-d4^2))-phi;
theta12 = atan2(-d4,-sqrt(m^2+n^2-d4^2))-phi;
if (theta11<-pi)
    theta11 = theta11+2*pi;
elseif(theta11>pi)
    theta11 = theta11-2*pi;
end
if (theta12<-pi)
    theta12 = theta12+2*pi;
elseif(theta12>pi)
    theta12 = theta12-2*pi;
end
angle11 = vpa(theta11*180/pi,3);
angle12 = vpa(theta12*180/pi,3);

T31 = T21_s*T32_s;
T41 = T31*T43_s;

%求解theta5
theta51 = vpa(acos(sin(theta11)*ax-cos(theta11)*ay),5);
theta52 = vpa(-acos(sin(theta11)*ax-cos(theta11)*ay),5);
theta53 = vpa(acos(sin(theta12)*ax-cos(theta12)*ay),5);
theta54 = vpa(-acos(sin(theta12)*ax-cos(theta12)*ay),5);
angle51 = theta51*180/3.14159;
angle52 = theta52*180/3.14159;
angle53 = theta53*180/3.14159;
angle54 = theta54*180/3.14159;

% %求解theta6
s1 = (sin(theta11)*nx-cos(theta11)*ny)/(sin(theta51));
t1 = (-sin(theta11)*ox+cos(theta11)*oy)/(sin(theta51));
theta61 = vpa(atan2(t1,s1),5);
angle61 = theta61*180/3.14159;

s2 = (sin(theta11)*nx-cos(theta11)*ny)/(sin(theta52));
t2 = (-sin(theta11)*ox+cos(theta11)*oy)/(sin(theta52));
theta62 = vpa(atan2(t2,s2),5);
angle62 = theta62*180/3.14159;

s3 = (sin(theta12)*nx-cos(theta12)*ny)/(sin(theta53));
t3 = (-sin(theta12)*ox+cos(theta12)*oy)/(sin(theta53));
theta63 = vpa(atan2(t3,s3),5);
angle63 = theta63*180/3.14159;

s4 = (sin(theta12)*nx-cos(theta12)*ny)/(sin(theta54));
t4 = (-sin(theta12)*ox+cos(theta12)*oy)/(sin(theta54));
theta64 = vpa(atan2(t4,s4),5);
angle64 = theta64*180/3.14159;

% %求解theta3
[theta31,theta32] = getTheta3(theta11,theta61,T60);
[theta33,theta34] = getTheta3(theta11,theta62,T60);
[theta35,theta36] = getTheta3(theta12,theta63,T60);
[theta37,theta38] = getTheta3(theta12,theta64,T60);

angle31 = theta31*180/3.14159;
angle32 = theta32*180/3.14159;
angle33 = theta33*180/3.14159;
angle34 = theta34*180/3.14159;
angle35 = theta35*180/3.14159;
angle36 = theta36*180/3.14159;
angle37 = theta37*180/3.14159;
angle38 = theta38*180/3.14159;

%求theta2
theta21 = getTheta2(theta11,theta31,theta61,T60);
theta22 = getTheta2(theta11,theta32,theta61,T60);
theta23 = getTheta2(theta11,theta33,theta62,T60);
theta24 = getTheta2(theta11,theta34,theta62,T60);
theta25 = getTheta2(theta12,theta35,theta63,T60);
theta26 = getTheta2(theta12,theta36,theta63,T60);
theta27 = getTheta2(theta12,theta37,theta64,T60);
theta28 = getTheta2(theta12,theta38,theta64,T60);

angle21 = theta21*180/3.14159;
angle22 = theta22*180/3.14159;
angle23 = theta23*180/3.14159;
angle24 = theta24*180/3.14159;
angle25 = theta25*180/3.14159;
angle26 = theta26*180/3.14159;
angle27 = theta27*180/3.14159;
angle28 = theta28*180/3.14159;


%求解theta4
theta41 = getTheta4(theta11,theta21,theta31,theta61,T60);
theta42 = getTheta4(theta11,theta22,theta32,theta61,T60);
theta43 = getTheta4(theta11,theta23,theta33,theta62,T60);
theta44 = getTheta4(theta11,theta24,theta34,theta62,T60);
theta45 = getTheta4(theta12,theta25,theta35,theta63,T60);
theta46 = getTheta4(theta12,theta26,theta36,theta63,T60);
theta47 = getTheta4(theta12,theta27,theta37,theta64,T60);
theta48 = getTheta4(theta12,theta28,theta38,theta64,T60);


angle41 = theta41*180/3.14159;
angle42 = theta42*180/3.14159;
angle43 = theta43*180/3.14159;
angle44 = theta44*180/3.14159;
angle45 = theta45*180/3.14159;
angle46 = theta46*180/3.14159;
angle47 = theta47*180/3.14159;
angle48 = theta48*180/3.14159;

digits(3) ;
result1 = [angle11 angle21 angle31 angle41 angle51 angle61]
result2 = [angle11 angle22 angle32 angle42 angle51 angle61]
result3 = [angle11 angle23 angle33 angle43 angle52 angle62]
result4 = [angle11 angle24 angle34 angle44 angle52 angle62]
result5 = [angle12 angle25 angle35 angle45 angle53 angle63]
result6 = [angle12 angle26 angle36 angle46 angle53 angle63]
result7 = [angle12 angle27 angle37 angle47 angle54 angle64]
result8 = [angle12 angle28 angle38 angle48 angle54 angle64]


T_60 = positiveKinematics(result5);