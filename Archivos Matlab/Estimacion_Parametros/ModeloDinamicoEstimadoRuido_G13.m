function [qdd] = ModeloDinamicoEstimadoRuido_G13(in)

% Variables de entrada en la funcion: [q(2)  qp(2)  Imotor(2)]
q1        = in(1);
q2        = in(2);
q3        = in(3);
qd1       = in(4);
qd2       = in(5);
qd3       = in(6);
Im1      = in(7);
Im2      = in(8);
Im3      = in(9);

R1=1; R2=1; R3=1; %Accionamiento directo
R=diag([R1,R2,R3]);
% K
  K=diag([0.5,0.4,0.35]); %(N*m/A)
  g=9.81;

Ma= [ 1.9283*cos(2.0*q2 + q3) + 2.9292*cos(2.0*q2) + 1.9283*cos(q3) + 0.54422*cos(2.0*q2 + 2.0*q3) + 3.4915,  0,  0;
    0, 3.8566*cos(q3) + 6.9901,   1.9283*cos(q3) + 1.0898;
    0, 1.9283*cos(q3) + 1.0898, 0.000093518*R3^2 + 1.0898];
 

Va= [ -8.4703e-22*qd1*(- 6.4183e15*R1^2 + 4.553e21*qd2*sin(2.0*q2 + q3) + 2.2765e21*qd3*sin(2.0*q2 + q3) + 6.9164e21*qd2*sin(2.0*q2) + 2.2765e21*qd3*sin(q3) + 1.285e21*qd2*sin(2.0*q2 + 2.0*q3) + 1.285e21*qd3*sin(2.0*q2 + 2.0*q3));
    0.54422*qd1^2*sin(2.0*q2 + 2.0*q3) - 1.9283*qd3^2*sin(q3) + 8.693e-6*R2^2*qd2 + 1.9283*qd1^2*sin(2.0*q2 + q3) + 2.9292*qd1^2*sin(2.0*q2) - 3.8566*qd2*qd3*sin(q3);
	0.96414*qd1^2*sin(q3) + 1.9283*qd2^2*sin(q3) + 0.54422*qd1^2*sin(2.0*q2 + 2.0*q3) + 0.000060046*R3^2*qd3 + 0.96414*qd1^2*sin(2.0*q2 + q3)];

Ga = [ 0;
    g*(1.9283*cos(q2 + q3) + 6.742*cos(q2));
	1.9283*g*cos(q2 + q3)];
 
   Im=[Im1;Im2;Im3];
  Tau=[Im1*R1*K(1,1);Im2*R2*K(2,2);Im3*R3*K(3,3)];
% Por lo que:  
% Aceleraciones
  qdd = inv(Ma)*(Tau-Va-Ga);