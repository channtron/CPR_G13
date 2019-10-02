%Necesitamos conocer:
    %vcdm,i (velocidad lineal de cada cdm respecto a U)
    %wi respecto a i
%La formulacion se basa en sacar las fuerzas y pares aplicados a partir de
%la lagrangiana del conjunto L=K-U

syms TauL1 TauL2 TauL3 q1 qd1 qdd1 q2 qd2 qdd2 q3 qd3 qdd3 g real  
PI = sym('pi');

% DATOS CINEMÁTICOS DEL BRAZO DEL ROBOT
% Dimensiones (m)
  L0=0;  % Base
  L1=1;  % Eslabón 1 
  L2=0.6557;  % Eslabón 2
  L3=0.5;  % Eslabón 3
  
% Parámetros de Denavit-Hartenberg (utilizado en primera regla de Newton-Euler)
% Eslabón base (no utilizado)
  theta0=0; d0=L0; a0=0; alpha0=0;
% Eslabón 1:
  theta1=q1; d1=1; a1=0 ; alpha1=0;
% Eslabón 2:
  theta2=PI/2 ; d2=q2+0.3 ; a2=0.3; alpha2=PI/2;
% Eslabón 3:
  theta3=0; d3=1+q3 ; a3=0; alpha3=0 ;
% Entre eslabón 3 y marco donde se ejerce la fuerza (a definir según
% experimento)
  theta4=0; d4=0 ; a4=0 ; alpha4=0 ;

% DATOS DINÁMICOS DEL BRAZO DEL ROBOT
% Eslabón 1
  m1=5.5 ; % kg
  s11 = [0 , 0 ,-0.5 ]'; % m
  I11=[ 0.4583,0  ,0  ;0 ,0.4583  ,0 ;0 ,0 ,0.0113 ]; % kg.m2

% Eslabón 2
  m2=4 ; % kg
  s22 = [-0.1226, -0.0409 , 0.1138]'; % m
  I22=[0.1281 -0.0290 -0.0558 ; -0.0290 0.1728 -0.0186 ; -0.0558 -0.0186 0.1004 ]; % kg.m2
 
% Eslabón 3
  m3=3.5 ; % kg
  s33 = [0, 0 ,-0.25 ]'; % m
  I33=[ 0.0729,  0, 0 ;0 , 0.0729 ,0 ;0 ,0 , 0.0072]; % kg.m2

% DATOS DE LOS MOTORES
% Inercias
  Jm1=0.25; Jm2=0.25; Jm3=0.08; % kg.m2% Coeficientes de fricción viscosa
% Coeficientes de fricción viscosa
  Bm1=3.6e-5 ; Bm2=3.6e-5 ; Bm3=3.6e-5 ; % N.m / (rad/s)
% Factores de reducción
  R1=25 ; R2=20 ; R3=25 ;
  
 %L1: Obtención de las matrices de rotación (i)R(i-1) y de sus inversas
  R01=[cos(theta1) -cos(alpha1)*sin(theta1) sin(alpha1)*sin(theta1);
      sin(theta1)  cos(alpha1)*cos(theta1)  -sin(alpha1)*cos(theta1);
      0            sin(alpha1)                cos(alpha1)           ];
  R10= R01';

  R12=[cos(theta2) -cos(alpha2)*sin(theta2) sin(alpha2)*sin(theta2);
      sin(theta2)  cos(alpha2)*cos(theta2)  -sin(alpha2)*cos(theta2);
      0            sin(alpha2)              cos(alpha2)           ];
  R21= R12';

  R23=[cos(theta3) -cos(alpha3)*sin(theta3) sin(alpha3)*sin(theta3);
      sin(theta3)  cos(alpha3)*cos(theta3)  -sin(alpha3)*cos(theta3);
      0            sin(alpha3)              cos(alpha3)           ];
  R32= R23';

  R34=[cos(theta4) -cos(alpha4)*sin(theta4) sin(alpha4)*sin(theta4);
      sin(theta4)  cos(alpha4)*cos(theta4)  -sin(alpha4)*cos(theta4);
      0            sin(alpha4)              cos(alpha4)           ];
  R43= R34';
  
  %Matrices de transformación a partir de los paarametros de Denavit
  %Hartenberg
  A01=MDH(q1,L1,0,0);
  A12=MDH(PI/2,L2+q2,L3,PI/2);
  A23=MDH(0,1+q3,0,0);
  
 % L2: Asignación a cada eslabón de sistema de referencia de acuerdo con las normas de D-H.
  % Eslabón 1:
    p11 = [a1, d1*sin(alpha1), d1*cos(alpha1)]';   
  % Eslabón 2:
    p22 = [a2, d2*sin(alpha2), d2*cos(alpha2)]'; 
  % Eslabón 3:
    p33 = [a3, d3*sin(alpha3), d3*cos(alpha3)]'; 
  % Entre eslabón 2 y marco donde se ejerce la fuerza (supongo que el mismo
  % que el Z0
    p44 = [a4, d4*sin(alpha4), d4*cos(alpha4)]'; 
    
    
  %Calculamos las velocidades lineales absolutas (cdm i) y velocidades angulares
  %absolutas del eslabón (i) en función del (i-1) --> RECURSIVIDAD
  
  %L3: Condiciones iniciales de la base
  w00=[0 0 0]';
  wd00 = [0 0 0]';
  v00 = [0 0 0]';
  vd00 = [0 0 g]'; % Aceleración de la gravedad en el eje Z0 negativo

    % Condiciones iniciales para el extremo del robot
         f44= [0 0 0]';
         n44= [0 0 0]';

    % Definición de vector local Z
          Z=[0 0 1]';
  
 %L4: Obtención velocidades angulares absolutas
    w11= R10*(w00+Z*qd1);  %Rotación
    w22 = R21*w11;      %Translación
    w33 = R32*w22;      %Translación
    
 %L5: Obtencion velocidades lineares absolutas
    v11=R10*v00+cross(w11,p11);     %Rotación
    v22=R21*(v11+Z*qd2)+cross(w22,p22);     %Translación
    v33=R32*(v22+Z*qd3)+cross(w33,p33);     %Translación
    
 %L6: Obtencióm v1elocidades lineares absolutas (cdm i)
    vcdm11=v11+cross(w11,s11);
    vcdm22=v22+cross(w22,s22);
    vcdm33=v33+cross(w33,s33);
    
 %L7:Calculo Posiciones cdm respecto a U
    AuxPcdmU1=A01*[s11;1];  PcdmU1=[AuxPcdmU1(1);AuxPcdmU1(2);AuxPcdmU1(3)];
    AuxPcdmU2=A01*A12*[s22;1];  PcdmU2=[AuxPcdmU2(1);AuxPcdmU2(2);AuxPcdmU2(3)];
    AuxPcdmU3=A01*A12*A23*[s33;1];  PcdmU3=[AuxPcdmU3(1);AuxPcdmU3(2);AuxPcdmU3(3)];
            
 %L8: Cálculo Energía Cinética y Potencial
    K=(0.5*m1*vcdm11'*vcdm11+0.5*w11'*I11*w11)+(0.5*m2*vcdm22'*vcdm22+0.5*w22'*I22*w22)+(0.5*m3*vcdm33'*vcdm33+0.5*w33'*I33*w33);
    U=(m1*g*Z'*PcdmU1)+(m2*g*Z'*PcdmU2)+(m3*g*Z'*PcdmU3);
       
 %L9: Cálculo Lagrangiana
    L=K-U;
    L=simplify(L);
  
 %L10: Calculamos las fuerzas aplicadas
    %Derivada de L respecto a qd
    dL_dqd1=diff(L,qd1);
    dL_dqd2=diff(L,qd2);
    dL_dqd3=diff(L,qd3);
    
    %Derivada de la derivada de L respecto a qd respecto a t->
    %diff(diff(L,qd),t)
    dL_dqd1_dt= diff(dL_dqd1,q1)*qd1 + diff(dL_dqd1,q2)*qd2 + diff(dL_dqd1,q3)*qd3 + diff(dL_dqd1,qd1)*qdd1 + diff(dL_dqd1,qd2)*qdd2 + diff(dL_dqd1,qd3)*qdd3;
    dL_dqd2_dt= diff(dL_dqd2,q1)*qd1 + diff(dL_dqd2,q2)*qd2 + diff(dL_dqd2,q3)*qd3 + diff(dL_dqd2,qd1)*qdd1 + diff(dL_dqd2,qd2)*qdd2 + diff(dL_dqd2,qd3)*qdd3;
    dL_dqd3_dt= diff(dL_dqd3,q1)*qd1 + diff(dL_dqd3,q2)*qd2 + diff(dL_dqd3,q3)*qd3 + diff(dL_dqd3,qd1)*qdd1 + diff(dL_dqd3,qd2)*qdd2 + diff(dL_dqd3,qd3)*qdd3;
    
    %Calculamos las Tau   
    TauL1=dL_dqd1_dt-diff(L,q1);
    TauL2=dL_dqd2_dt-diff(L,q2);
    TauL3=dL_dqd3_dt-diff(L,q3);
    
    %Simplificamos y redondeamos
    TauL1=vpa(simplify(TauL1),5)
    TauL2=vpa(simplify(TauL2),5)
    TauL3=vpa(simplify(TauL3),5)
    
    %Tal y como se puede observar, las fuerzas y pares aplicados
    %calculados por el Algoritmo de Lagrange son exactamente iguales que
    %los calculados por el metodo de Newton-Euler
    
    %Sin embargo, Lagrange no nos da los esfuerzos sobre las articulaciones
    %mientras que Newton-Euler si
    
  T1=simplify(TauL1)
  T2=simplify(TauL2)
  T3=simplify(TauL3)
  
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
  
  %Metemos la reductora
  R=diag([R1 R2 R3]);
  Jm=diag([Jm1 Jm2 Jm3]);
  Bm=diag([Bm1 Bm2 Bm3]);
  
  %Construimos la matrices finales
  Ma=M+R*R*Jm;      Ma=simplify(Ma);
  Va=V+R*R*Bm*[qd1;qd2;qd3];        Va=simplify(Va);
  Ga=G;
  
  %Redondeamos
  Ma=vpa(Ma,5)
  Va=vpa(Va,5)
  Ga=vpa(Ga,5)
  