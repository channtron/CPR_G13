%% Estimación de parámetros
clear all;
% syms T1 T2 T3 q1 qd1 qdd1 q2 qd2 qdd2 q3 qd3 qdd3 real
% PI = sym('pi');
syms m1 m2 m3 lc1 lc2 lc3 Ixx1 Ixx2 Ixx3 Iyy1 Iyy2 Iyy3 Izz1 Izz2 Izz3 Jm1 Jm2 Jm3 Bm1 Bm2 Bm3 real
Tm=0.001;

A11=0.5; A12=.5; F11=5; F12=3; Ph11=0.5; Ph12=0;

sim('sk_R3GDL_2017');
% Factores de reducción
  %R1=50; R2=30; R3=15;
  R1=1; R2=1; R3=1; %Accionamiento directo
 
% K
  K=diag([0.5,0.4,0.35]); %(N*m/A)


Sigma_red=[Iyy1+(R1^2)*Jm1+Iyy2-Izz2-(R2^2)*Jm2+Iyy3-Izz3 ;Bm1 ;m2-m2*lc2^2-Izz2-(R2^2)*Jm2+(16/9)*m3*lc3^2+(16/9)*Izz3;...
    m2*lc2-m2*lc2^2-Izz2-(R2^2)*Jm2;Ixx2-Iyy2+Izz2+(R2^2)*Jm2;Bm2;m3-(16/9)*m3*lc3^2-(16/9)*Izz3;...
    m3*lc3-(4/3)*m3*lc3^2-(4/3)*Izz3;Ixx3-Iyy3+Izz3;Jm3;Bm3];

% Gamma_red= [  qdd1, 2500*qd1, 0.5*qdd1 + 0.5*qdd1*cos(2.0*q2) - 1.0*qd1*qd2*sin(2.0*q2), 2.0*qd1*qd2*sin(2.0*q2) - 1.0*qdd1*cos(2.0*q2) - 1.0*qdd1, qdd1/2 - 0.5*qdd1*cos(2.0*q2) + qd1*qd2*sin(2.0*q2),       0, 0.78125*qdd1 + 0.75*qdd1*cos(2.0*q2 + q3) + 0.28125*qdd1*cos(2*q2 + 2*q3) + 0.5*qdd1*cos(2.0*q2) + 0.75*qdd1*cos(q3) - 0.75*qd1*qd3*sin(q3) - 1.5*qd1*qd2*sin(2.0*q2 + q3) - 0.75*qd1*qd3*sin(2.0*q2 + q3) - 0.5625*qd1*qd2*sin(2*q2 + 2*q3) - 0.5625*qd1*qd3*sin(2*q2 + 2*q3) - 1.0*qd1*qd2*sin(2.0*q2), qd1*qd3*sin(q3) - 1.0*qdd1*cos(2.0*q2 + q3) - 0.75*qdd1*cos(2*q2 + 2*q3) - 1.0*qdd1*cos(q3) - 0.75*qdd1 + 2.0*qd1*qd2*sin(2.0*q2 + q3) + qd1*qd3*sin(2.0*q2 + q3) + 1.5*qd1*qd2*sin(2*q2 + 2*q3) + 1.5*qd1*qd3*sin(2*q2 + 2*q3),  qdd1/2 - 0.5*qdd1*cos(2*q2 + 2*q3) + qd1*qd2*sin(2*q2 + 2*q3) + qd1*qd3*sin(2*q2 + 2*q3),        0,       0;...
%              0    ,        0,                  0.5*sin(2.0*q2)*qd1^2 + qdd2 + g*cos(q2),            - 1.0*sin(2.0*q2)*qd1^2 - 2.0*qdd2 - g*cos(q2),                              -0.5*qd1^2*sin(2.0*q2), 900*qd2,                                                                     1.5625*qdd2 + 0.5625*qdd3 - 0.75*qd3^2*sin(q3) + 0.75*g*cos(q2 + q3) + 0.75*qd1^2*sin(2.0*q2 + q3) + g*cos(q2) + 0.28125*qd1^2*sin(2*q2 + 2*q3) + 1.5*qdd2*cos(q3) + 0.75*qdd3*cos(q3) + 0.5*qd1^2*sin(2.0*q2) - 1.5*qd2*qd3*sin(q3),                                                 qd3^2*sin(q3) - 1.5*qdd3 - 1.5*qdd2 - 1.0*g*cos(q2 + q3) - 1.0*qd1^2*sin(2.0*q2 + q3) - 0.75*qd1^2*sin(2*q2 + 2*q3) - 2.0*qdd2*cos(q3) - 1.0*qdd3*cos(q3) + 2.0*qd2*qd3*sin(q3),                                                               -0.5*qd1^2*sin(2*q2 + 2*q3),        0,       0;...
%              0    ,        0,                                                         0,                                                         0,                                                   0,       0,                                                                                                                           0.5625*qdd2 + 0.5625*qdd3 + 0.375*qd1^2*sin(q3) + 0.75*qd2^2*sin(q3) + 0.75*g*cos(q2 + q3) + 0.375*qd1^2*sin(2.0*q2 + q3) + 0.28125*qd1^2*sin(2*q2 + 2*q3) + 0.75*qdd2*cos(q3),                                                                - 1.5*qdd2 - 1.5*qdd3 - 0.5*qd1^2*sin(q3) - 1.0*qd2^2*sin(q3) - 1.0*g*cos(q2 + q3) - 0.5*qd1^2*sin(2.0*q2 + q3) - 0.75*qd1^2*sin(2*q2 + 2*q3) - 1.0*qdd2*cos(q3),                                                               -0.5*qd1^2*sin(2*q2 + 2*q3), 225*qdd3, 225*qd3];
%  

<<<<<<< HEAD
Y=Im*diag([R1 R2 R3])*K;
g=9.81;
j=1;
for i=2000:10:size(t)
    
    Gamma(j:j+2,:)= [  qdd(i,1), (R1^2)*qd(i,1), 0.5*qdd(i,1) + 0.5*qdd(i,1)*cos(2.0*q(i,2)) - 1.0*qd(i,1)*qd(i,2)*sin(2.0*q(i,2)), 2.0*qd(i,1)*qd(i,2)*sin(2.0*q(i,2)) - 1.0*qdd(i,1)*cos(2.0*q(i,2)) - 1.0*qdd(i,1), qdd(i,1)/2 - 0.5*qdd(i,1)*cos(2.0*q(i,2)) + qd(i,1)*qd(i,2)*sin(2.0*q(i,2)),       0, 0.78125*qdd(i,1) + 0.75*qdd(i,1)*cos(2.0*q(i,2) + q(i,3)) + 0.28125*qdd(i,1)*cos(2*q(i,2) + 2*q(i,3)) + 0.5*qdd(i,1)*cos(2.0*q(i,2)) + 0.75*qdd(i,1)*cos(q(i,3)) - 0.75*qd(i,1)*qd(i,3)*sin(q(i,3)) - 1.5*qd(i,1)*qd(i,2)*sin(2.0*q(i,2) + q(i,3)) - 0.75*qd(i,1)*qd(i,3)*sin(2.0*q(i,2) + q(i,3)) - 0.5625*qd(i,1)*qd(i,2)*sin(2*q(i,2) + 2*q(i,3)) - 0.5625*qd(i,1)*qd(i,3)*sin(2*q(i,2) + 2*q(i,3)) - 1.0*qd(i,1)*qd(i,2)*sin(2.0*q(i,2)), qd(i,1)*qd(i,3)*sin(q(i,3)) - 1.0*qdd(i,1)*cos(2.0*q(i,2) + q(i,3)) - 0.75*qdd(i,1)*cos(2*q(i,2) + 2*q(i,3)) - 1.0*qdd(i,1)*cos(q(i,3)) - 0.75*qdd(i,1) + 2.0*qd(i,1)*qd(i,2)*sin(2.0*q(i,2) + q(i,3)) + qd(i,1)*qd(i,3)*sin(2.0*q(i,2) + q(i,3)) + 1.5*qd(i,1)*qd(i,2)*sin(2*q(i,2) + 2*q(i,3)) + 1.5*qd(i,1)*qd(i,3)*sin(2*q(i,2) + 2*q(i,3)),  qdd(i,1)/2 - 0.5*qdd(i,1)*cos(2*q(i,2) + 2*q(i,3)) + qd(i,1)*qd(i,2)*sin(2*q(i,2) + 2*q(i,3)) + qd(i,1)*qd(i,3)*sin(2*q(i,2) + 2*q(i,3)),        0,       0;...
             0    ,        0,                  0.5*sin(2.0*q(i,2))*qd(i,1)^2 + qdd(i,2) + g*cos(q(i,2)),            - 1.0*sin(2.0*q(i,2))*qd(i,1)^2 - 2.0*qdd(i,2) - g*cos(q(i,2)),                              -0.5*qd(i,1)^2*sin(2.0*q(i,2)), (R2^2)*qd(i,2),                                                                     1.5625*qdd(i,2) + 0.5625*qdd(i,3) - 0.75*qd(i,3)^2*sin(q(i,3)) + 0.75*g*cos(q(i,2) + q(i,3)) + 0.75*qd(i,1)^2*sin(2.0*q(i,2) + q(i,3)) + g*cos(q(i,2)) + 0.28125*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)) + 1.5*qdd(i,2)*cos(q(i,3)) + 0.75*qdd(i,3)*cos(q(i,3)) + 0.5*qd(i,1)^2*sin(2.0*q(i,2)) - 1.5*qd(i,2)*qd(i,3)*sin(q(i,3)),                                                 qd(i,3)^2*sin(q(i,3)) - 1.5*qdd(i,3) - 1.5*qdd(i,2) - 1.0*g*cos(q(i,2) + q(i,3)) - 1.0*qd(i,1)^2*sin(2.0*q(i,2) + q(i,3)) - 0.75*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)) - 2.0*qdd(i,2)*cos(q(i,3)) - 1.0*qdd(i,3)*cos(q(i,3)) + 2.0*qd(i,2)*qd(i,3)*sin(q(i,3)),                                                               -0.5*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)),        0,       0;...
             0    ,        0,                                                         0,                                                         0,                                                   0,       0,                                                                                                                           0.5625*qdd(i,2) + 0.5625*qdd(i,3) + 0.375*qd(i,1)^2*sin(q(i,3)) + 0.75*qd(i,2)^2*sin(q(i,3)) + 0.75*g*cos(q(i,2) + q(i,3)) + 0.375*qd(i,1)^2*sin(2.0*q(i,2) + q(i,3)) + 0.28125*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)) + 0.75*qdd(i,2)*cos(q(i,3)),                                                                - 1.5*qdd(i,2) - 1.5*qdd(i,3) - 0.5*qd(i,1)^2*sin(q(i,3)) - 1.0*qd(i,2)^2*sin(q(i,3)) - 1.0*g*cos(q(i,2) + q(i,3)) - 0.5*qd(i,1)^2*sin(2.0*q(i,2) + q(i,3)) - 0.75*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)) - 1.0*qdd(i,2)*cos(q(i,3)),                                                               -0.5*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)), (R3^2)*qdd(i,3), (R3^2)*qd(i,3)];
=======
Y=Im*K*R;
g=9.81;
j=1;
for i=2000:10:size(t)
   
    Gamma(j:j+2,:)= [  qdd(i,1), 2500*qd(i,1), 0.5*qdd(i,1) + 0.5*qdd(i,1)*cos(2.0*q(i,2)) - 1.0*qd(i,1)*qd(i,2)*sin(2.0*q(i,2)), 2.0*qd(i,1)*qd(i,2)*sin(2.0*q(i,2)) - 1.0*qdd(i,1)*cos(2.0*q(i,2)) - 1.0*qdd(i,1), qdd(i,1)/2 - 0.5*qdd(i,1)*cos(2.0*q(i,2)) + qd(i,1)*qd(i,2)*sin(2.0*q(i,2)),       0, 0.78125*qdd(i,1) + 0.75*qdd(i,1)*cos(2.0*q(i,2) + q(i,3)) + 0.28125*qdd(i,1)*cos(2*q(i,2) + 2*q(i,3)) + 0.5*qdd(i,1)*cos(2.0*q(i,2)) + 0.75*qdd(i,1)*cos(q(i,3)) - 0.75*qd(i,1)*qd(i,3)*sin(q(i,3)) - 1.5*qd(i,1)*qd(i,2)*sin(2.0*q(i,2) + q(i,3)) - 0.75*qd(i,1)*qd(i,3)*sin(2.0*q(i,2) + q(i,3)) - 0.5625*qd(i,1)*qd(i,2)*sin(2*q(i,2) + 2*q(i,3)) - 0.5625*qd(i,1)*qd(i,3)*sin(2*q(i,2) + 2*q(i,3)) - 1.0*qd(i,1)*qd(i,2)*sin(2.0*q(i,2)), qd(i,1)*qd(i,3)*sin(q(i,3)) - 1.0*qdd(i,1)*cos(2.0*q(i,2) + q(i,3)) - 0.75*qdd(i,1)*cos(2*q(i,2) + 2*q(i,3)) - 1.0*qdd(i,1)*cos(q(i,3)) - 0.75*qdd(i,1) + 2.0*qd(i,1)*qd(i,2)*sin(2.0*q(i,2) + q(i,3)) + qd(i,1)*qd(i,3)*sin(2.0*q(i,2) + q(i,3)) + 1.5*qd(i,1)*qd(i,2)*sin(2*q(i,2) + 2*q(i,3)) + 1.5*qd(i,1)*qd(i,3)*sin(2*q(i,2) + 2*q(i,3)),  qdd(i,1)/2 - 0.5*qdd(i,1)*cos(2*q(i,2) + 2*q(i,3)) + qd(i,1)*qd(i,2)*sin(2*q(i,2) + 2*q(i,3)) + qd(i,1)*qd(i,3)*sin(2*q(i,2) + 2*q(i,3)),        0,       0;...
             0    ,        0,                  0.5*sin(2.0*q(i,2))*qd(i,1)^2 + qdd(i,2) + g*cos(q(i,2)),            - 1.0*sin(2.0*q(i,2))*qd(i,1)^2 - 2.0*qdd(i,2) - g*cos(q(i,2)),                              -0.5*qd(i,1)^2*sin(2.0*q(i,2)), 900*qd(i,2),                                                                     1.5625*qdd(i,2) + 0.5625*qdd(i,3) - 0.75*qd(i,3)^2*sin(q(i,3)) + 0.75*g*cos(q(i,2) + q(i,3)) + 0.75*qd(i,1)^2*sin(2.0*q(i,2) + q(i,3)) + g*cos(q(i,2)) + 0.28125*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)) + 1.5*qdd(i,2)*cos(q(i,3)) + 0.75*qdd(i,3)*cos(q(i,3)) + 0.5*qd(i,1)^2*sin(2.0*q(i,2)) - 1.5*qd(i,2)*qd(i,3)*sin(q(i,3)),                                                 qd(i,3)^2*sin(q(i,3)) - 1.5*qdd(i,3) - 1.5*qdd(i,2) - 1.0*g*cos(q(i,2) + q(i,3)) - 1.0*qd(i,1)^2*sin(2.0*q(i,2) + q(i,3)) - 0.75*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)) - 2.0*qdd(i,2)*cos(q(i,3)) - 1.0*qdd(i,3)*cos(q(i,3)) + 2.0*qd(i,2)*qd(i,3)*sin(q(i,3)),                                                               -0.5*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)),        0,       0;...
             0    ,        0,                                                         0,                                                         0,                                                   0,       0,                                                                                                                           0.5625*qdd(i,2) + 0.5625*qdd(i,3) + 0.375*qd(i,1)^2*sin(q(i,3)) + 0.75*qd(i,2)^2*sin(q(i,3)) + 0.75*g*cos(q(i,2) + q(i,3)) + 0.375*qd(i,1)^2*sin(2.0*q(i,2) + q(i,3)) + 0.28125*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)) + 0.75*qdd(i,2)*cos(q(i,3)),                                                                - 1.5*qdd(i,2) - 1.5*qdd(i,3) - 0.5*qd(i,1)^2*sin(q(i,3)) - 1.0*qd(i,2)^2*sin(q(i,3)) - 1.0*g*cos(q(i,2) + q(i,3)) - 0.5*qd(i,1)^2*sin(2.0*q(i,2) + q(i,3)) - 0.75*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)) - 1.0*qdd(i,2)*cos(q(i,3)),                                                               -0.5*qd(i,1)^2*sin(2*q(i,2) + 2*q(i,3)), 225*qdd(i,3), 225*qd(i,3)];
>>>>>>> 1c8453d2d46543f2d5cf3f1d08537800c7900125
    Taux(j:j+2,1)=Y(i,:)';
    j= j+3;
end

theta=lscov(Gamma,Taux); %pseudo inversa del primero por el segundo
%theta=inv(Gamma'*Gamma)*Gamma'*Taux;
[rn,m]=size(Gamma);
sigma_cuad= (norm(Taux-Gamma*theta)^2)/(rn-m);
csig=sigma_cuad*inv(Gamma'*Gamma);
size(csig);
var=sqrt(diag(csig));
fprintf('El error al estimar los parametros en porcentaje es: \n');
est_err=(100*(var./theta))'




