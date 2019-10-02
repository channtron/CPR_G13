% Ejemplo de la utilización del algoritmo de Newton Euler para la dinámica
% de un robot de 3 DGL
% M.G. Ortega (2017)

syms T1 T2 T3 q1 qd1 qdd1 q2 qd2 qdd2 q3 qd3 qdd3 g real  
PI = sym('pi');



% DATOS CINEMÁTICOS DEL BRAZO DEL ROBOT
% Dimensiones (m)
    L1=0.65;
    L2=1;
    L0=0.5;
    L3=0.75;
  
% Parámetros de Denavit-Hartenberg (utilizado en primera regla de Newton-Euler)
% Eslabón base (no utilizado)
  theta0=0; d0=L0; a0=0; alpha0=0;
% Eslabón 1:
  theta1=q1; d1=L1; a1=0 ; alpha1=PI/2;
% Eslabón 2:
  theta2=q2 ; d2=0 ; a2=L2; alpha2=0;
% Eslabón 3:
  theta3=q3; d3=0 ; a3=L3; alpha3=0 ;
% Entre eslabón 3 y marco donde se ejerce la fuerza (a definir según
% experimento)
  theta4=0; d4=0 ; a4=0 ; alpha4=0 ;

% DATOS DINÁMICOS DEL BRAZO DEL ROBOT
    %TODOS SIMBOLICOS
    syms m1 m2 m3 lc1 lc2 lc3 Ixx1 Ixx2 Ixx3 Iyy1 Iyy2 Iyy3 Izz1 Izz2 Izz3 Jm1 Jm2 Jm3 Bm1 Bm2 Bm3 real
% Eslabón 1
%   m1=5.5 ; % kg
   s11 = [0 , -lc1,0]'; % m
   I11=[ Ixx1,0  ,0  ;0 ,Iyy1  ,0 ;0 ,0 ,Izz1]; % kg.m2

% Eslabón 2
%   m2=4 ; % kg
   s22 = [-lc2, 0 , 0]'; % m
   I22=[Ixx2 ,0  ,0  ;0 ,Iyy2  ,0 ;0 ,0 ,Izz2 ]; % kg.m2

% Eslabón 3
%   m3=3.5 ; % kg
   s33 = [-lc3, 0 ,0]'; % m
   I33=[Ixx3 ,0  ,0  ;0 ,Iyy3  ,0 ;0 ,0 ,Izz3 ]; % kg.m2

% DATOS DE LOS MOTORES
% Inercias
  %Jm1=0.25; Jm2=0.25; Jm3=0.08; % kg.m2
% Coeficientes de fricción viscosa
  %Bm1=3.6e-5 ; Bm2=3.6e-5 ; Bm3=3.6e-5 ; % N.m / (rad/s)
% Factores de reducción
  R1=50 ; R2=30 ; R3=15 ;
% K
  K1=0.5; K2=0.4; K3=0.35; %(N*m/A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORÍTMO DE NEWTON-EULER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% wij : velocidad angular absoluta de eje j expresada en i
% wdij : aceleración angular absoluta de eje j expresada en i
% vij : velocidad lineal absoluta del origen del marco j expresada en i
% vdij : aceleración lineal absoluta del origen del marco j expresada en i
% aii : aceleración del centro de gravedad del eslabón i, expresado en i?

% fij : fuerza ejercida sobre la articulación j-1 (unión barra j-1 con j),
% expresada en i-1
%
% nij : par ejercido sobre la articulación j-1 (unión barra j-1 con j),
% expresada en i-1

% pii : vector (libre) que une el origen de coordenadas de i-1 con el de i,
% expresadas en i : [ai, di*sin(alphai), di*cos(alphai)] (a,d,aplha: parámetros de DH)
%
% sii : coordenadas del centro de masas del eslabón i, expresada en el sistema
% i

% Iii : matriz de inercia del eslabón i expresado en un sistema paralelo al
% i y con el origen en el centro de masas del eslabón
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% N-E 1: Asignación a cada eslabón de sistema de referencia de acuerdo con las normas de D-H.
  % Eslabón 1:
    p11 = [a1, d1*sin(alpha1), d1*cos(alpha1)]';   
  % Eslabón 2:
    p22 = [a2, d2*sin(alpha2), d2*cos(alpha2)]'; 
  % Eslabón 3:
    p33 = [a3, d3*sin(alpha3), d3*cos(alpha3)]'; 
  % Entre eslabón 2 y marco donde se ejerce la fuerza (supongo que el mismo
  % que el Z0
    p44 = [a4, d4*sin(alpha4), d4*cos(alpha4)]'; 

% N-E 2: Condiciones iniciales de la base
  w00=[0 0 0]';
  wd00 = [0 0 0]';
  v00 = [0 0 0]';
  vd00 = [0 0 g]'; % Aceleración de la gravedad en el eje Z0 negativo

% Condiciones iniciales para el extremo del robot
  f44= [0 0 0]';
  n44= [0 0 0]';

% Definición de vector local Z
  Z=[0 0 1]';


% N-E 3: Obtención de las matrices de rotación (i)R(i-1) y de sus inversas
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


%%%%%%% ITERACIÓN HACIA EL EXTERIOR (CINEMÁTICA)

% N-E 4: Obtención de las velocidades angulares absolutas
 % Articulación 1
    w11= R10*(w00+Z*qd1);  % Si es de rotación
   % w11 = R10*w00;      % Si es de translación
 % Articulación 2
    w22= R21*(w11+Z*qd2);  % Si es de rotación
    %w22 = R21*w11;      % Si es de translación
 % Articulación 3
   w33= R32*(w22+Z*qd3);  % Si es de rotación
    %w33 = R32*w22;      % Si es de translación

% N-E 5: Obtención de las aceleraciones angulares absolutas
 % Articulación 1
    wd11 = R10*(wd00+Z*qdd1+cross(w00,Z*qd1));  % si es de rotación
   % wd11 = R10*wd00;                                % si es de translación
 % Articulación 2
   wd22 = R21*(wd11+Z*qdd2+cross(w11,Z*qd2));  % si es de rotación
%     wd22 = R21*wd11;                                % si es de translación
 % Articulación 3
   wd33 = R32*(wd22+Z*qdd3+cross(w22,Z*qd3));  % si es de rotación
%     wd33 = R32*wd22;                                % si es de translación

% N-E 6: Obtención de las aceleraciones lineales de los orígenes de los
% sistemas
 % Articulación 1
    vd11 = cross(wd11,p11)+cross(w11,cross(w11,p11))+R10*vd00;  % si es de rotación
   % vd11 = R10*(Z*qdd1+vd00)+cross(wd11,p11)+2*cross(w11,R10*Z*qd1) + cross(w11,cross(w11,p11));    % si es de translación
 % Articulación 2
   vd22 = cross(wd22,p22)+cross(w22,cross(w22,p22))+R21*vd11;  % si es de rotación
%     vd22 = R21*(Z*qdd2+vd11)+cross(wd22,p22)+2*cross(w22,R21*Z*qd2) + cross(w22,cross(w22,p22));    % si es de translación
 % Articulación 3
   vd33 = cross(wd33,p33)+cross(w33,cross(w33,p33))+R32*vd22;  % si es de rotación
%     vd33 = R32*(Z*qdd3+vd22)+cross(wd33,p33)+2*cross(w33,R32*Z*qd3) + cross(w33,cross(w33,p33));    % si es de translación

% N-E 7: Obtención de las aceleraciones lineales de los centros de gravedad
    a11 = cross(wd11,s11)+cross(w11,cross(w11,s11))+vd11;
    a22 = cross(wd22,s22)+cross(w22,cross(w22,s22))+vd22;
    a33 = cross(wd33,s33)+cross(w33,cross(w33,s33))+vd33;

%%%%%%% ITERACIÓN HACIA EL INTERIOR (DINÁMICA)

% N-E 8: Obtención de las fuerzas ejercidas sobre los eslabones
  f33=R34*f44+m3*a33;
  f22=R23*f33+m2*a22;
  f11=R12*f22+m1*a11;

% N-E 9: Obtención de los pares ejercidas sobre los eslabones
  n33 = R34*(n44+cross(R43*p33,f44))+cross(p33+s33,m3*a33)+I33*wd33+cross(w33,I33*w33);
  n22 = R23*(n33+cross(R32*p22,f33))+cross(p22+s22,m2*a22)+I22*wd22+cross(w22,I22*w22);
  n11 = R12*(n22+cross(R21*p11,f22))+cross(p11+s11,m1*a11)+I11*wd11+cross(w11,I11*w11);

% N-E 10: Obtener la fuerza o par aplicado sobre la articulación
  N3z = n33'*R32*Z;  % Si es de rotación
  N3  = n33'*R32;    % Para ver todos los pares, no solo el del eje Z
%   F3z = f33'*R32*Z;  % Si es de translacion;
%   F3  = f33'*R32;    % Para ver todas las fuerzas, no solo la del eje Z
  N2z = n22'*R21*Z;  % Si es de rotación
  N2  = n22'*R21;    % Para ver todos los pares, no solo el del eje Z
%   F2z = f22'*R21*Z;  % Si es de translacion;
%   F2  = f22'*R21;    % Para ver todas las fuerzas, no solo la del eje Z
  N1z = n11'*R10*Z;  % Si es de rotación
  N1  = n11'*R10;    % Para ver todos los pares, no solo el del eje Z
  %F1z = f11'*R10*Z;  % Si es de translacion;
  %F1  = f11'*R10;    % Para ver todas las fuerzas, no solo la del eje Z

% Robot ??? (descomentar los que procedan)
  % T1=F1z;
    T1=N1z;
   %T2=F2z;
   T2=N2z;
  % T3=F3z;
   T3=N3z;
  
  %A PARTIR DE AQUÍ ES CÓDIGO NUESTRO
  %---------------------------------------------------------------------
  T1=simplify(T1);
  T2=simplify(T2);
  T3=simplify(T3);
  
  %Sacamos los parametros derivando
  p_aux1=diff(T1,m1);
  p_aux2=diff(p_aux1,lc1);
  p_m1aa1=diff(p_aux2,lc1)/2; %T1 m1*lc1^2 
  T1_aux=simplify(T1-p_m1aa1*m1*lc1^2);
  
  p_aux1=diff(T1_aux,m2);
  p_aux2=diff(p_aux1,lc2);
  p_m2aa1=diff(p_aux2,lc2)/2; %T1 m2*lc2^2 
  T1_aux=simplify(T1_aux-p_m2aa1*m2*lc2^2);
  
  p_aux1=diff(T1_aux,m3);
  p_aux2=diff(p_aux1,lc3);
  p_m3aa1=diff(p_aux2,lc3)/2; %T1 m3*lc3^2 
  T1_aux=simplify(T1_aux-p_m3aa1*m3*lc3^2);
  
  p_aux1=diff(T1_aux,m1);
  p_m1a1=diff(p_aux1,lc1);  %T1 m1*lc1
  T1_aux=simplify(T1_aux-p_m1a1*m1*lc1);
  
  p_aux1=diff(T1_aux,m2);
  p_m2a1=diff(p_aux1,lc2);  %T1 m2*lc2
  T1_aux=simplify(T1_aux-p_m2a1*m2*lc2);
  
  p_aux1=diff(T1_aux,m3);
  p_m3a1=diff(p_aux1,lc3);  %T1 m3*lc3
  T1_aux=simplify(T1_aux-p_m3a1*m3*lc3);
  
  p_m11=diff(T1_aux,m1);     %T1 m1
  T1_aux=simplify(T1_aux-p_m11*m1);
  
  p_m21=diff(T1_aux,m2);     %T1 m2
  T1_aux=simplify(T1_aux-p_m21*m2);
  
  p_m31=diff(T1_aux,m3);     %T1 m3
  T1_aux=simplify(T1_aux-p_m31*m3);

  p_Ixx11=diff(T1_aux,Ixx1); %T1 Ixx1
  T1_aux=simplify(T1_aux-p_Ixx11*Ixx1);
  
  p_Ixx21=diff(T1_aux,Ixx2); %T1 Ixx2
  T1_aux=simplify(T1_aux-p_Ixx21*Ixx2);
  
  p_Ixx31=diff(T1_aux,Ixx3); %T1 Ixx3
  T1_aux=simplify(T1_aux-p_Ixx31*Ixx3);
  
  p_Iyy11=diff(T1_aux,Iyy1); %T1 Iyy1
  T1_aux=simplify(T1_aux-p_Iyy11*Iyy1);
  
  p_Iyy21=diff(T1_aux,Iyy2); %T1 Iyy2
  T1_aux=simplify(T1_aux-p_Iyy21*Iyy2);
  
  p_Iyy31=diff(T1_aux,Iyy3); %T1 Iyy3
  T1_aux=simplify(T1_aux-p_Iyy31*Iyy3);
  
  p_Izz11=diff(T1_aux,Izz1); %T1 Izz1
  T1_aux=simplify(T1_aux-p_Izz11*Izz1);
  
  p_Izz21=diff(T1_aux,Izz2); %T1 Izz2
  T1_aux=simplify(T1_aux-p_Izz21*Izz2);
  
  p_Izz31=diff(T1_aux,Izz3); %T1 Izz3
  T1_aux=simplify(T1_aux-p_Izz31*Izz3);
  
  