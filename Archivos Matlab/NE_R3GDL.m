%Estimacion_Parametros;

clear g;
syms T1 T2 T3 q1 qd1 qdd1 q2 qd2 qdd2 q3 qd3 qdd3 g real  
PI = sym('pi');


Gamma_red= [  qdd1, 2500*qd1, 0.5*qdd1 + 0.5*qdd1*cos(2.0*q2) - 1.0*qd1*qd2*sin(2.0*q2), 2.0*qd1*qd2*sin(2.0*q2) - 1.0*qdd1*cos(2.0*q2) - 1.0*qdd1, qdd1/2 - 0.5*qdd1*cos(2.0*q2) + qd1*qd2*sin(2.0*q2),       0, 0.78125*qdd1 + 0.75*qdd1*cos(2.0*q2 + q3) + 0.28125*qdd1*cos(2*q2 + 2*q3) + 0.5*qdd1*cos(2.0*q2) + 0.75*qdd1*cos(q3) - 0.75*qd1*qd3*sin(q3) - 1.5*qd1*qd2*sin(2.0*q2 + q3) - 0.75*qd1*qd3*sin(2.0*q2 + q3) - 0.5625*qd1*qd2*sin(2*q2 + 2*q3) - 0.5625*qd1*qd3*sin(2*q2 + 2*q3) - 1.0*qd1*qd2*sin(2.0*q2), qd1*qd3*sin(q3) - 1.0*qdd1*cos(2.0*q2 + q3) - 0.75*qdd1*cos(2*q2 + 2*q3) - 1.0*qdd1*cos(q3) - 0.75*qdd1 + 2.0*qd1*qd2*sin(2.0*q2 + q3) + qd1*qd3*sin(2.0*q2 + q3) + 1.5*qd1*qd2*sin(2*q2 + 2*q3) + 1.5*qd1*qd3*sin(2*q2 + 2*q3),  qdd1/2 - 0.5*qdd1*cos(2*q2 + 2*q3) + qd1*qd2*sin(2*q2 + 2*q3) + qd1*qd3*sin(2*q2 + 2*q3),        0,       0;...
             0    ,        0,                  0.5*sin(2.0*q2)*qd1^2 + qdd2 + g*cos(q2),            - 1.0*sin(2.0*q2)*qd1^2 - 2.0*qdd2 - g*cos(q2),                              -0.5*qd1^2*sin(2.0*q2), 900*qd2,                                                                     1.5625*qdd2 + 0.5625*qdd3 - 0.75*qd3^2*sin(q3) + 0.75*g*cos(q2 + q3) + 0.75*qd1^2*sin(2.0*q2 + q3) + g*cos(q2) + 0.28125*qd1^2*sin(2*q2 + 2*q3) + 1.5*qdd2*cos(q3) + 0.75*qdd3*cos(q3) + 0.5*qd1^2*sin(2.0*q2) - 1.5*qd2*qd3*sin(q3),                                                 qd3^2*sin(q3) - 1.5*qdd3 - 1.5*qdd2 - 1.0*g*cos(q2 + q3) - 1.0*qd1^2*sin(2.0*q2 + q3) - 0.75*qd1^2*sin(2*q2 + 2*q3) - 2.0*qdd2*cos(q3) - 1.0*qdd3*cos(q3) + 2.0*qd2*qd3*sin(q3),                                                               -0.5*qd1^2*sin(2*q2 + 2*q3),        0,       0;...
             0    ,        0,                                                         0,                                                         0,                                                   0,       0,                                                                                                                           0.5625*qdd2 + 0.5625*qdd3 + 0.375*qd1^2*sin(q3) + 0.75*qd2^2*sin(q3) + 0.75*g*cos(q2 + q3) + 0.375*qd1^2*sin(2.0*q2 + q3) + 0.28125*qd1^2*sin(2*q2 + 2*q3) + 0.75*qdd2*cos(q3),                                                                - 1.5*qdd2 - 1.5*qdd3 - 0.5*qd1^2*sin(q3) - 1.0*qd2^2*sin(q3) - 1.0*g*cos(q2 + q3) - 0.5*qd1^2*sin(2.0*q2 + q3) - 0.75*qd1^2*sin(2*q2 + 2*q3) - 1.0*qdd2*cos(q3),                                                               -0.5*qd1^2*sin(2*q2 + 2*q3), 225*qdd3, 225*qd3];
  Tau=Gamma_red*theta;
  T1=Tau(1);
  T2=Tau(2);
  T3=Tau(3);
  
 
  T1=simplify(T1);
  T2=simplify(T2);
  T3=simplify(T3);
  
  %T=M(q)*qdd + V(q,qd) + G(q)
  
  %Primera ecuacion
  M11=diff(T1,qdd1);
  T1aux=simplify(T1-M11*qdd1);
  M12=diff(T1aux,qdd2);
  T1aux=simplify(T1aux-M12*qdd2);
  M13=diff(T1aux,qdd3);
  T1aux=simplify(T1aux-M13*qdd3);  
  G1=diff(T1aux,g)*g;
  T1aux=simplify(T1aux-G1);  
  V1=T1aux;
  
  %Segunda ecuacion
  M21=diff(T2,qdd1);    
  T2aux=simplify(T2-M21*qdd1);
  M22=diff(T2aux,qdd2);
  T2aux=simplify(T2aux-M22*qdd2);
  M23=diff(T2aux,qdd3);
  T2aux=simplify(T2aux-M23*qdd3);  
  G2=diff(T2aux,g)*g;
  T2aux=simplify(T2aux-G2);  
  V2=T2aux;
  
  %Segunda ecuacion
  M31=diff(T3,qdd1);
  T3aux=simplify(T3-M31*qdd1);
  M32=diff(T3aux,qdd2);
  T3aux=simplify(T3aux-M32*qdd2);
  M33=diff(T3aux,qdd3);
  T3aux=simplify(T3aux-M33*qdd3);
  G3=diff(T3aux,g)*g;
  T3aux=simplify(T3aux-G3);  
  V3=T3aux;

  %Simplificamos todos los términos
  M11=simplify(M11);    M12=simplify(M12);  M13=simplify(M13);
  M21=simplify(M21);    M22=simplify(M22);  M23=simplify(M23);
  M31=simplify(M31);    M32=simplify(M32);  M33=simplify(M33);
  V1=simplify(V1);  V2=simplify(V2);    V3=simplify(V3);
  G1=simplify(G1);  G2=simplify(G2);    G3=simplify(G3);
  
  %Creamos matrices y vectores
  M=[M11 M12 M13; M21 M22 M23; M31 M32 M33];
  V=[V1 V2 V3]';
  G=[G1 G2 G3]';
  
  
  %Redondeamos
  M=vpa(M,5)
  V=vpa(V,5)
  G=vpa(G,5)
  