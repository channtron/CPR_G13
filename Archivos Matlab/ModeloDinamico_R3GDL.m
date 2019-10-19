function [qdd] = ModeloDinamico_R3GDL(in)

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

% Factores de reducción
  R=diag([50,30,15]) ;
% K
  K=diag([0.5,0.4,0.35]); %(N*m/A)
  g=9.81;

% Matriz de Inercias
  M =[ 2.0808*cos(2.0*q2 + q3) + 2.7549*cos(2.0*q2) + 2.0808*cos(q3) + 0.61697*cos(2.0*q2 + 2.0*q3) + 29.839,  0,   0;
    0, 4.1616*cos(q3) + 48.364, 2.0808*cos(q3) + 1.2355;
  0, 2.0808*cos(q3) + 1.2355,                   6.722];
 
  
% Matriz de aceleraciones centrípetas y de Coriolis
  V = [ -2.1684e-19*qd1*(1.9192e19*qd2*sin(2.0*q2 + q3) + 9.5961e18*qd3*sin(2.0*q2 + q3) + 2.541e19*qd2*sin(2.0*q2) + 9.5961e18*qd3*sin(q3) + 5.6906e18*qd2*sin(2.0*q2 + 2.0*q3) + 5.6906e18*qd3*sin(2.0*q2 + 2.0*q3) - 1.3835e20);
  7.65*qd2 - 2.0808*qd3^2*sin(q3) + 0.61697*qd1^2*sin(2.0*q2 + 2.0*q3) + 2.0808*qd1^2*sin(2.0*q2 + q3) + 2.7549*qd1^2*sin(2.0*q2) - 4.1616*qd2*qd3*sin(q3);
   3.375*qd3 + 1.0404*qd1^2*sin(q3) + 2.0808*qd2^2*sin(q3) + 0.61697*qd1^2*sin(2.0*q2 + 2.0*q3) + 1.0404*qd1^2*sin(2.0*q2 + q3)];
 
% Par gravitatorio                
  G = [  0;
 g*(2.0808*cos(q2 + q3) + 6.3812*cos(q2));
                    2.0808*g*cos(q2 + q3)];

% Ecuación del robot
%    Tau = M*qpp + V + G
  Im=[Im1;Im2;Im3];
  Tau=[Im1*R(1)*K(1);Im2*R(2)*K(2);Im3*R(3)*K(3)];
% Por lo que:  
% Aceleraciones
  qdd = inv(M)*(Tau-V-G);
  