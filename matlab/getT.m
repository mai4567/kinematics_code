function [T] = getT(a,alpha,d,theta)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
T=[ cos(theta)                -sin(theta)               0             a;
    sin(theta)*cos(alpha)  cos(theta)*cos(alpha)  -sin(alpha)  -d*sin(alpha);
    sin(theta)*sin(alpha)  cos(theta)*sin(alpha)   cos(alpha)   d*cos(alpha);
    0                              0                    0              1];
end

