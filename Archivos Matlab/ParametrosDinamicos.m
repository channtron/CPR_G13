
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
  % T2=F2z;
   T2=N2z;
  % T3=F3z;
   T3=N3z;
  
  %A PARTIR DE AQUÍ ES CÓDIGO NUESTRO
  %---------------------------------------------------------------------
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
  
  %%
  
  %Ecuación con las matrices ampliadas
  T=simplify(Ma*[qdd1;qdd2;qdd3]+Va+Ga); 
  T1=T(1);
  T2=T(2);
  T3=T(3);
 
  %Sacamos los parametros derivando
  
  %Eslabón 1
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
  
  p_Jm11=diff(T1_aux,Jm1);   %T1 Jm1
  T1_aux=simplify(T1_aux-p_Jm11*Jm1);
  
  p_Jm21=diff(T1_aux,Jm2);   %T1 Jm2
  T1_aux=simplify(T1_aux-p_Jm21*Jm2);
  
  p_Jm31=diff(T1_aux,Jm3);   %T1 Jm3
  T1_aux=simplify(T1_aux-p_Jm31*Jm3);
  
  p_Bm11=diff(T1_aux,Bm1);   %T1 Bm1
  T1_aux=simplify(T1_aux-p_Bm11*Bm1);
  
  p_Bm21=diff(T1_aux,Bm2);   %T1 Bm2
  T1_aux=simplify(T1_aux-p_Bm21*Bm2);
  
  p_Bm31=diff(T1_aux,Bm3);   %T1 Bm3
  T1_aux=simplify(T1_aux-p_Bm31*Bm3);
  
  %Vector de los parametros para el eslabón 1
  p_1=[p_m11, p_m21, p_m31, p_m1a1, p_m2a1, p_m3a1, p_m1aa1, p_m2aa1, p_m3aa1,...
      p_Ixx11, p_Ixx21, p_Ixx31, p_Iyy11, p_Iyy21, p_Iyy31, p_Izz11, p_Izz21,p_Izz31,...
      p_Jm11, p_Jm21, p_Jm31, p_Bm11, p_Bm21, p_Bm31];
  
 
  %Eslabón 2
  p_aux1=diff(T2,m1);
  p_aux2=diff(p_aux1,lc1);
  p_m1aa2=diff(p_aux2,lc1)/2; %T2 m1*lc1^2 
  T2_aux=simplify(T2-p_m1aa2*m1*lc1^2);
  
  p_aux1=diff(T2_aux,m2);
  p_aux2=diff(p_aux1,lc2);
  p_m2aa2=diff(p_aux2,lc2)/2; %T2 m2*lc2^2 
  T2_aux=simplify(T2_aux-p_m2aa2*m2*lc2^2);
  
  p_aux1=diff(T2_aux,m3);
  p_aux2=diff(p_aux1,lc3);
  p_m3aa2=diff(p_aux2,lc3)/2; %T2 m3*lc3^2 
  T2_aux=simplify(T2_aux-p_m3aa2*m3*lc3^2);
  
  p_aux1=diff(T2_aux,m1);
  p_m1a2=diff(p_aux1,lc1);  %T2 m1*lc1
  T2_aux=simplify(T2_aux-p_m1a2*m1*lc1);
  
  p_aux1=diff(T2_aux,m2);
  p_m2a2=diff(p_aux1,lc2);  %T2 m2*lc2
  T2_aux=simplify(T2_aux-p_m2a2*m2*lc2);
  
  p_aux1=diff(T2_aux,m3);
  p_m3a2=diff(p_aux1,lc3);  %T2 m3*lc3
  T2_aux=simplify(T2_aux-p_m3a2*m3*lc3);
  
  p_m12=diff(T2_aux,m1);     %T2 m1
  T2_aux=simplify(T2_aux-p_m12*m1);
  
  p_m22=diff(T2_aux,m2);     %T2 m2
  T2_aux=simplify(T2_aux-p_m22*m2);
  
  p_m32=diff(T2_aux,m3);     %T2 m3
  T2_aux=simplify(T2_aux-p_m32*m3);
  
  p_Ixx12=diff(T2_aux,Ixx1); %T2 Ixx1
  T2_aux=simplify(T2_aux-p_Ixx12*Ixx1);
  
  p_Ixx22=diff(T2_aux,Ixx2); %T2 Ixx2
  T2_aux=simplify(T2_aux-p_Ixx22*Ixx2);
  
  p_Ixx32=diff(T2_aux,Ixx3); %T2 Ixx3
  T2_aux=simplify(T2_aux-p_Ixx32*Ixx3);
  
  p_Iyy12=diff(T2_aux,Iyy1); %T2 Iyy1
  T2_aux=simplify(T2_aux-p_Iyy12*Iyy1);
  
  p_Iyy22=diff(T2_aux,Iyy2); %T2 Iyy2
  T2_aux=simplify(T2_aux-p_Iyy22*Iyy2);
  
  p_Iyy32=diff(T2_aux,Iyy3); %T2 Iyy3
  T2_aux=simplify(T2_aux-p_Iyy32*Iyy3);
  
  p_Izz12=diff(T2_aux,Izz1); %T2 Izz1
  T2_aux=simplify(T2_aux-p_Izz12*Izz1);
  
  p_Izz22=diff(T2_aux,Izz2); %T2 Izz2
  T2_aux=simplify(T2_aux-p_Izz22*Izz2);
  
  p_Izz32=diff(T2_aux,Izz3); %T2 Izz3
  T2_aux=simplify(T2_aux-p_Izz32*Izz3);
  
  p_Jm12=diff(T2_aux,Jm1);   %T2 Jm1
  T2_aux=simplify(T2_aux-p_Jm12*Jm1);
  
  p_Jm22=diff(T2_aux,Jm2);   %T2 Jm2
  T2_aux=simplify(T2_aux-p_Jm22*Jm2);
  
  p_Jm32=diff(T2_aux,Jm3);   %T2 Jm3
  T2_aux=simplify(T2_aux-p_Jm32*Jm3);
  
  p_Bm12=diff(T2_aux,Bm1);   %T2 Bm1
  T2_aux=simplify(T2_aux-p_Bm12*Bm1);
  
  p_Bm22=diff(T2_aux,Bm2);   %T2 Bm2
  T2_aux=simplify(T2_aux-p_Bm22*Bm2);
  
  p_Bm32=diff(T2_aux,Bm3);   %T2 Bm3
  T2_aux=simplify(T2_aux-p_Bm32*Bm3);
  
  %Vector de los parametros para el eslabón 2
  p_2=[p_m12, p_m22, p_m32, p_m1a2, p_m2a2, p_m3a2, p_m1aa2, p_m2aa2, p_m3aa2,...
      p_Ixx12, p_Ixx22, p_Ixx32, p_Iyy12, p_Iyy22, p_Iyy32, p_Izz12, p_Izz22, p_Izz32,...
      p_Jm12, p_Jm22, p_Jm32, p_Bm12, p_Bm22, p_Bm32];
  
  
  %Eslabón 3
  p_aux1=diff(T3,m1);
  p_aux2=diff(p_aux1,lc1);
  p_m1aa3=diff(p_aux2,lc1)/2; %T3 m1*lc1^2 
  T3_aux=simplify(T3-p_m1aa3*m1*lc1^2);
  
  p_aux1=diff(T3_aux,m2);
  p_aux2=diff(p_aux1,lc2);
  p_m2aa3=diff(p_aux2,lc2)/2; %T3 m2*lc2^2 
  T3_aux=simplify(T3_aux-p_m2aa3*m2*lc2^2);
  
  p_aux1=diff(T3_aux,m3);
  p_aux2=diff(p_aux1,lc3);
  p_m3aa3=diff(p_aux2,lc3)/2; %T3 m3*lc3^2 
  T3_aux=simplify(T3_aux-p_m3aa3*m3*lc3^2);
  
  p_aux1=diff(T3_aux,m1);
  p_m1a3=diff(p_aux1,lc1);  %T3 m1*lc1
  T3_aux=simplify(T3_aux-p_m1a3*m1*lc1);
  
  p_aux1=diff(T3_aux,m2);
  p_m2a3=diff(p_aux1,lc2);  %T3 m2*lc2
  T3_aux=simplify(T3_aux-p_m2a3*m2*lc2);
  
  p_aux1=diff(T3_aux,m3);
  p_m3a3=diff(p_aux1,lc3);  %T3 m3*lc3
  T3_aux=simplify(T3_aux-p_m3a3*m3*lc3);
  
  p_m13=diff(T3_aux,m1);     %T3 m1
  T3_aux=simplify(T3_aux-p_m13*m1);
  
  p_m23=diff(T3_aux,m2);     %T3 m2
  T3_aux=simplify(T3_aux-p_m23*m2);
  
  p_m33=diff(T3_aux,m3);     %T3 m3
  T3_aux=simplify(T3_aux-p_m33*m3);
  
  p_Ixx13=diff(T3_aux,Ixx1); %T3 Ixx1
  T3_aux=simplify(T3_aux-p_Ixx13*Ixx1);
  
  p_Ixx23=diff(T3_aux,Ixx2); %T3 Ixx2
  T3_aux=simplify(T3_aux-p_Ixx23*Ixx2);
  
  p_Ixx33=diff(T3_aux,Ixx3); %T3 Ixx3
  T3_aux=simplify(T3_aux-p_Ixx33*Ixx3);
  
  p_Iyy13=diff(T3_aux,Iyy1); %T3 Iyy1
  T3_aux=simplify(T3_aux-p_Iyy13*Iyy1);
  
  p_Iyy23=diff(T3_aux,Iyy2); %T3 Iyy2
  T3_aux=simplify(T3_aux-p_Iyy23*Iyy2);
  
  p_Iyy33=diff(T3_aux,Iyy3); %T3 Iyy3
  T3_aux=simplify(T3_aux-p_Iyy33*Iyy3);
  
  p_Izz13=diff(T3_aux,Izz1); %T3 Izz1
  T3_aux=simplify(T3_aux-p_Izz13*Izz1);
  
  p_Izz23=diff(T3_aux,Izz2); %T3 Izz2
  T3_aux=simplify(T3_aux-p_Izz23*Izz2);
  
  p_Izz33=diff(T3_aux,Izz3); %T3 Izz3
  T3_aux=simplify(T3_aux-p_Izz33*Izz3);
  
  p_Jm13=diff(T3_aux,Jm1);   %T3 Jm1
  T3_aux=simplify(T3_aux-p_Jm13*Jm1);
  
  p_Jm23=diff(T3_aux,Jm2);   %T3 Jm2
  T3_aux=simplify(T3_aux-p_Jm23*Jm2);
  
  p_Jm33=diff(T3_aux,Jm3);   %T3 Jm3
  T3_aux=simplify(T3_aux-p_Jm33*Jm3);
  
  p_Bm13=diff(T3_aux,Bm1);   %T3 Bm1
  T3_aux=simplify(T3_aux-p_Bm13*Bm1);
  
  p_Bm23=diff(T3_aux,Bm2);   %T3 Bm2
  T3_aux=simplify(T3_aux-p_Bm23*Bm2);
  
  p_Bm33=diff(T3_aux,Bm3);   %T3 Bm3
  T3_aux=simplify(T3_aux-p_Bm33*Bm3);
  
  %Vector de los parametros para el eslabón 3
  p_3=[p_m13, p_m23, p_m33, p_m1a3, p_m2a3, p_m3a3, p_m1aa3, p_m2aa3, p_m3aa3,...
      p_Ixx13, p_Ixx23, p_Ixx33, p_Iyy13, p_Iyy23, p_Iyy33, p_Izz13, p_Izz23,p_Izz33,...
      p_Jm13, p_Jm23, p_Jm33, p_Bm13, p_Bm23, p_Bm33];
  
  
  %K*R*I=[.....]3x24 * [..]24x1
  
  %Formamos la matriz 3x24 con los valores que multiplican a los parámetros a identificar
  Mult=simplify([p_1;p_2;p_3]);
  
  %Formamos el vector con los 24 parámetros
  Param=[m1; m2; m3; m1*lc1; m2*lc2; m3*lc3; m1*lc1^2; m2*lc2^2; m3*lc3^2; Ixx1; Ixx2; Ixx3;...
      Iyy1; Iyy2; Iyy3; Izz1; Izz2; Izz3; Jm1; Jm2; Jm3; Bm1; Bm2; Bm3];
  
  
  %Comprobamos que es correcto
    %Mult*Param - Ma*qdd - Va - Ga=
    simplify(Mult*Param-T)
  
  
  %Matrices
  Mult=[0, 0.5*qdd1 + 0.5*qdd1*cos(2.0*q2) - 1.0*qd1*qd2*sin(2.0*q2), 0.78125*qdd1 + 0.75*qdd1*cos(2.0*q2 + q3) + 0.28125*qdd1*cos(2*q2 + 2*q3) + 0.5*qdd1*cos(2.0*q2) + 0.75*qdd1*cos(q3) - 0.75*qd1*qd3*sin(q3) - 1.5*qd1*qd2*sin(2.0*q2 + q3) - 0.75*qd1*qd3*sin(2.0*q2 + q3) - 0.5625*qd1*qd2*sin(2*q2 + 2*q3) - 0.5625*qd1*qd3*sin(2*q2 + 2*q3) - 1.0*qd1*qd2*sin(2.0*q2), 0, 2.0*qd1*qd2*sin(2.0*q2) - 1.0*qdd1*cos(2.0*q2) - 1.0*qdd1, qd1*qd3*sin(q3) - 1.0*qdd1*cos(2.0*q2 + q3) - 0.75*qdd1*cos(2*q2 + 2*q3) - 1.0*qdd1*cos(q3) - 0.75*qdd1 + 2.0*qd1*qd2*sin(2.0*q2 + q3) + qd1*qd3*sin(2.0*q2 + q3) + 1.5*qd1*qd2*sin(2*q2 + 2*q3) + 1.5*qd1*qd3*sin(2*q2 + 2*q3), 0, (qdd1*(1.0*cos(2.0*q2) + 1.0))/2 - 1.0*qd1*qd2*sin(2.0*q2), 0.5*qdd1 + 0.5*qdd1*cos(2*q2 + 2*q3) - 1.0*qd1*qd2*sin(2*q2 + 2*q3) - 1.0*qd1*qd3*sin(2*q2 + 2*q3), 0, qdd1/2 - 0.5*qdd1*cos(2.0*q2) + qd1*qd2*sin(2.0*q2), qdd1/2 - 0.5*qdd1*cos(2*q2 + 2*q3) + qd1*qd2*sin(2*q2 + 2*q3) + qd1*qd3*sin(2*q2 + 2*q3), qdd1, 0.5*qdd1 + 0.5*qdd1*cos(2.0*q2) - 1.0*qd1*qd2*sin(2.0*q2), 0.5*qdd1 + 0.5*qdd1*cos(2*q2 + 2*q3) - 1.0*qd1*qd2*sin(2*q2 + 2*q3) - 1.0*qd1*qd3*sin(2*q2 + 2*q3), 0,0,0, 2500*qdd1,0,0, 2500*qd1,0,0;...
        0,0.5*sin(2.0*q2)*qd1^2 + qdd2 + g*cos(q2),1.5625*qdd2 + 0.5625*qdd3 - 0.75*qd3^2*sin(q3) + 0.75*g*cos(q2 + q3) + 0.75*qd1^2*sin(2.0*q2 + q3) + g*cos(q2) + 0.28125*qd1^2*sin(2*q2 + 2*q3) + 1.5*qdd2*cos(q3) + 0.75*qdd3*cos(q3) + 0.5*qd1^2*sin(2.0*q2) - 1.5*qd2*qd3*sin(q3), 0, - 1.0*sin(2.0*q2)*qd1^2 - 2.0*qdd2 - g*cos(q2),qd3^2*sin(q3) - 1.5*qdd3 - 1.5*qdd2 - 1.0*g*cos(q2 + q3) - 1.0*qd1^2*sin(2.0*q2 + q3) - 0.75*qd1^2*sin(2*q2 + 2*q3) - 2.0*qdd2*cos(q3) - 1.0*qdd3*cos(q3) + 2.0*qd2*qd3*sin(q3), 0,0.5*sin(2.0*q2)*qd1^2 + qdd2, 0.5*sin(2*q2 + 2*q3)*qd1^2 + qdd2 + qdd3, 0,-0.5*qd1^2*sin(2.0*q2),-0.5*qd1^2*sin(2*q2 + 2*q3),0, 0.5*qd1^2*sin(2.0*q2),                                                                         0.5*qd1^2*sin(2*q2 + 2*q3), 0, qdd2, qdd2 + qdd3,0, 900*qdd2,0,0, 900*qd2,0;...
        0,0,0.5625*qdd2 + 0.5625*qdd3 + 0.375*qd1^2*sin(q3) + 0.75*qd2^2*sin(q3) + 0.75*g*cos(q2 + q3) + 0.375*qd1^2*sin(2.0*q2 + q3) + 0.28125*qd1^2*sin(2*q2 + 2*q3) + 0.75*qdd2*cos(q3), 0,0,- 1.5*qdd2 - 1.5*qdd3 - 0.5*qd1^2*sin(q3) - 1.0*qd2^2*sin(q3) - 1.0*g*cos(q2 + q3) - 0.5*qd1^2*sin(2.0*q2 + q3) - 0.75*qd1^2*sin(2*q2 + 2*q3) - 1.0*qdd2*cos(q3), 0, 0, 0.5*sin(2*q2 + 2*q3)*qd1^2 + qdd2 + qdd3, 0,0,-0.5*qd1^2*sin(2*q2 + 2*q3),    0,0,0.5*qd1^2*sin(2*q2 + 2*q3), 0,0, qdd2 + qdd3, 0,0, 225*qdd3,0,0, 225*qd3]
  
  
  %%
  %Son nulas m1, m1*lc1, m1*lc1^2, Ixx1 e Izz1
  Gamma=[  qdd1/2 + (qdd1*cos(2*q2))/2 - qd1*qd2*sin(2*q2), (25*qdd1)/32 + (qdd1*cos(2*q2))/2 + (9*qdd1*cos(2*q2 + 2*q3))/32 + (3*qdd1*cos(q3))/4 + (3*qdd1*cos(2*q2 + q3))/4 - (3*qd1*qd3*sin(q3))/4 - (3*qd1*qd2*sin(2*q2 + q3))/2 - (3*qd1*qd3*sin(2*q2 + q3))/4 - qd1*qd2*sin(2*q2) - (9*qd1*qd2*sin(2*q2 + 2*q3))/16 - (9*qd1*qd3*sin(2*q2 + 2*q3))/16,  2*qd1*qd2*sin(2*q2) - qdd1*cos(2*q2) - qdd1, qd1*qd3*sin(q3) - (3*qdd1*cos(2*q2 + 2*q3))/4 - qdd1*cos(q3) - qdd1*cos(2*q2 + q3) - (3*qdd1)/4 + 2*qd1*qd2*sin(2*q2 + q3) + qd1*qd3*sin(2*q2 + q3) + (3*qd1*qd2*sin(2*q2 + 2*q3))/2 + (3*qd1*qd3*sin(2*q2 + 2*q3))/2,  (qdd1*(cos(2*q2) + 1))/2 - qd1*qd2*sin(2*q2), qdd1/2 + (qdd1*cos(2*q2 + 2*q3))/2 - qd1*qd2*sin(2*q2 + 2*q3) - qd1*qd3*sin(2*q2 + 2*q3),  qdd1/2 - (qdd1*cos(2*q2))/2 + qd1*qd2*sin(2*q2), qdd1/2 - (qdd1*cos(2*q2 + 2*q3))/2 + qd1*qd2*sin(2*q2 + 2*q3) + qd1*qd3*sin(2*q2 + 2*q3), qdd1, qdd1/2 + (qdd1*cos(2*q2))/2 - qd1*qd2*sin(2*q2), qdd1/2 + (qdd1*cos(2*q2 + 2*q3))/2 - qd1*qd2*sin(2*q2 + 2*q3) - qd1*qd3*sin(2*q2 + 2*q3),     0,           0, 2500*qdd1,        0,        0, 2500*qd1,       0,       0;
            (sin(2*q2)*qd1^2)/2 + qdd2 + g*cos(q2),                                                        (25*qdd2)/16 + (9*qdd3)/16 - (3*qd3^2*sin(q3))/4 + (3*qd1^2*sin(2*q2 + q3))/4 + (3*g*cos(q2 + q3))/4 + (qd1^2*sin(2*q2))/2 + g*cos(q2) + (9*qd1^2*sin(2*q2 + 2*q3))/32 + (3*qdd2*cos(q3))/2 + (3*qdd3*cos(q3))/4 - (3*qd2*qd3*sin(q3))/2,       - sin(2*q2)*qd1^2 - 2*qdd2 - g*cos(q2),                                                    qd3^2*sin(q3) - (3*qdd3)/2 - (3*qdd2)/2 - qd1^2*sin(2*q2 + q3) - g*cos(q2 + q3) - (3*qd1^2*sin(2*q2 + 2*q3))/4 - 2*qdd2*cos(q3) - qdd3*cos(q3) + 2*qd2*qd3*sin(q3),                    (sin(2*q2)*qd1^2)/2 + qdd2,                                                 (sin(2*q2 + 2*q3)*qd1^2)/2 + qdd2 + qdd3,                             -(qd1^2*sin(2*q2))/2,                                                              -(qd1^2*sin(2*q2 + 2*q3))/2,    0,                             (qd1^2*sin(2*q2))/2,                                                               (qd1^2*sin(2*q2 + 2*q3))/2,  qdd2, qdd2 + qdd3,         0, 900*qdd2,        0,        0, 900*qd2,       0;
                                                0,                                                                                                                  (9*qdd2)/16 + (9*qdd3)/16 + (3*qd1^2*sin(q3))/8 + (3*qd2^2*sin(q3))/4 + (3*qd1^2*sin(2*q2 + q3))/8 + (3*g*cos(q2 + q3))/4 + (9*qd1^2*sin(2*q2 + 2*q3))/32 + (3*qdd2*cos(q3))/4,                                            0,                                                               - sin(q3)*qd2^2 - (3*qdd2)/2 - (3*qdd3)/2 - (qd1^2*sin(q3))/2 - (qd1^2*sin(2*q2 + q3))/2 - g*cos(q2 + q3) - (3*qd1^2*sin(2*q2 + 2*q3))/4 - qdd2*cos(q3),                                             0,                                                 (sin(2*q2 + 2*q3)*qd1^2)/2 + qdd2 + qdd3,                                                0,                                                              -(qd1^2*sin(2*q2 + 2*q3))/2,    0,                                               0,                                                               (qd1^2*sin(2*q2 + 2*q3))/2,     0, qdd2 + qdd3,         0,        0, 225*qdd3,        0,       0, 225*qd3];

  clear Gamma_numerica;
  for ind=1:3:19
    qdd1=rand; qdd2=rand; qdd3=rand; qd1=rand; qd2=rand; qd3=rand; q1=rand; q2=rand; q3=rand; g=rand;
    Gamma_numerica(ind:ind+2,:)=eval(Gamma);
  end
  % Gamma_numerica=Gamma_numerica_aux(1:19,:);
    size(Gamma_numerica);

  rank(Gamma_numerica);
  [AA,ColumnasIndep]=rref(Gamma_numerica)
  % Linealmente Independientes 1 2 3 4 7 8 9 16 17 18 19
  
  %%
    clear all
    syms T1 T2 T3 q1 qd1 qdd1 q2 qd2 qdd2 q3 qd3 qdd3 g real
    PI = sym('pi');
    syms m1 m2 m3 lc1 lc2 lc3 Ixx1 Ixx2 Ixx3 Iyy1 Iyy2 Iyy3 Izz1 Izz2 Izz3 Jm1 Jm2 Jm3 Bm1 Bm2 Bm3 real
    Gamma=[  qdd1/2 + (qdd1*cos(2*q2))/2 - qd1*qd2*sin(2*q2), (25*qdd1)/32 + (qdd1*cos(2*q2))/2 + (9*qdd1*cos(2*q2 + 2*q3))/32 + (3*qdd1*cos(q3))/4 + (3*qdd1*cos(2*q2 + q3))/4 - (3*qd1*qd3*sin(q3))/4 - (3*qd1*qd2*sin(2*q2 + q3))/2 - (3*qd1*qd3*sin(2*q2 + q3))/4 - qd1*qd2*sin(2*q2) - (9*qd1*qd2*sin(2*q2 + 2*q3))/16 - (9*qd1*qd3*sin(2*q2 + 2*q3))/16,  2*qd1*qd2*sin(2*q2) - qdd1*cos(2*q2) - qdd1, qd1*qd3*sin(q3) - (3*qdd1*cos(2*q2 + 2*q3))/4 - qdd1*cos(q3) - qdd1*cos(2*q2 + q3) - (3*qdd1)/4 + 2*qd1*qd2*sin(2*q2 + q3) + qd1*qd3*sin(2*q2 + q3) + (3*qd1*qd2*sin(2*q2 + 2*q3))/2 + (3*qd1*qd3*sin(2*q2 + 2*q3))/2,  (qdd1*(cos(2*q2) + 1))/2 - qd1*qd2*sin(2*q2), qdd1/2 + (qdd1*cos(2*q2 + 2*q3))/2 - qd1*qd2*sin(2*q2 + 2*q3) - qd1*qd3*sin(2*q2 + 2*q3),  qdd1/2 - (qdd1*cos(2*q2))/2 + qd1*qd2*sin(2*q2), qdd1/2 - (qdd1*cos(2*q2 + 2*q3))/2 + qd1*qd2*sin(2*q2 + 2*q3) + qd1*qd3*sin(2*q2 + 2*q3), qdd1, qdd1/2 + (qdd1*cos(2*q2))/2 - qd1*qd2*sin(2*q2), qdd1/2 + (qdd1*cos(2*q2 + 2*q3))/2 - qd1*qd2*sin(2*q2 + 2*q3) - qd1*qd3*sin(2*q2 + 2*q3),     0,           0, 2500*qdd1,        0,        0, 2500*qd1,       0,       0;
    (sin(2*q2)*qd1^2)/2 + qdd2 + g*cos(q2),                                                        (25*qdd2)/16 + (9*qdd3)/16 - (3*qd3^2*sin(q3))/4 + (3*qd1^2*sin(2*q2 + q3))/4 + (3*g*cos(q2 + q3))/4 + (qd1^2*sin(2*q2))/2 + g*cos(q2) + (9*qd1^2*sin(2*q2 + 2*q3))/32 + (3*qdd2*cos(q3))/2 + (3*qdd3*cos(q3))/4 - (3*qd2*qd3*sin(q3))/2,       - sin(2*q2)*qd1^2 - 2*qdd2 - g*cos(q2),                                                    qd3^2*sin(q3) - (3*qdd3)/2 - (3*qdd2)/2 - qd1^2*sin(2*q2 + q3) - g*cos(q2 + q3) - (3*qd1^2*sin(2*q2 + 2*q3))/4 - 2*qdd2*cos(q3) - qdd3*cos(q3) + 2*qd2*qd3*sin(q3),                    (sin(2*q2)*qd1^2)/2 + qdd2,                                                 (sin(2*q2 + 2*q3)*qd1^2)/2 + qdd2 + qdd3,                             -(qd1^2*sin(2*q2))/2,                                                              -(qd1^2*sin(2*q2 + 2*q3))/2,    0,                             (qd1^2*sin(2*q2))/2,                                                               (qd1^2*sin(2*q2 + 2*q3))/2,  qdd2, qdd2 + qdd3,         0, 900*qdd2,        0,        0, 900*qd2,       0;
    0,                                                                                                                  (9*qdd2)/16 + (9*qdd3)/16 + (3*qd1^2*sin(q3))/8 + (3*qd2^2*sin(q3))/4 + (3*qd1^2*sin(2*q2 + q3))/8 + (3*g*cos(q2 + q3))/4 + (9*qd1^2*sin(2*q2 + 2*q3))/32 + (3*qdd2*cos(q3))/4,                                            0,                                                               - sin(q3)*qd2^2 - (3*qdd2)/2 - (3*qdd3)/2 - (qd1^2*sin(q3))/2 - (qd1^2*sin(2*q2 + q3))/2 - g*cos(q2 + q3) - (3*qd1^2*sin(2*q2 + 2*q3))/4 - qdd2*cos(q3),                                             0,                                                 (sin(2*q2 + 2*q3)*qd1^2)/2 + qdd2 + qdd3,                                                0,                                                              -(qd1^2*sin(2*q2 + 2*q3))/2,    0,                                               0,                                                               (qd1^2*sin(2*q2 + 2*q3))/2,     0, qdd2 + qdd3,         0,        0, 225*qdd3,        0,       0, 225*qd3];
    AA=rref(Gamma);
    %AA=simplify(AA)
%     AA =
%  
% [ 1, 0, 0, ((9*qd1^4*qdd1*cos(q3))/2 - (9*qd1^4*qdd1*cos(4*q2 + 3*q3))/2 + 36*qdd1*qdd3^2*cos(q3) - 72*qdd1*qdd2^2*cos(2*q2 + q3) + 50*g^2*qdd1*cos(q2 + q3)*cos(q2) + 144*qd1*qd2*qdd2^2*sin(2*q2 + q3) + 72*qd1*qd3*qdd2^2*sin(2*q2 + q3) + 36*qd1^2*qdd1*qdd2*sin(2*q2 + q3) + 18*qd1^2*qdd1*qdd3*sin(2*q2 + q3) + 36*g*qdd1*qdd2*cos(q2 + q3) - 18*qd1^3*qd2*qdd2*cos(4*q2 + 3*q3) - 9*qd1^2*qdd1*qdd2*sin(4*q2 + 3*q3) + 36*qdd1*qdd2^2*cos(2*q2)*cos(q3) + 36*qdd1*qdd3^2*cos(2*q2)*cos(q3) - 72*qdd1*qdd2*qdd3*cos(2*q2 + q3) + 9*qd1^4*qdd1*sin(2*q2)*sin(q3) + 36*qdd1*qdd2^2*cos(2*q2 + 2*q3)*cos(q3) + 96*g*qdd1*qdd2*cos(q2)^3 + 96*g*qdd1*qdd3*cos(q2)^3 + 18*qd1^3*qd2*qdd2*cos(q3) + 72*qd1*qd3*qdd2^2*sin(q3) - 27*qd1^2*qdd1*qdd2*sin(q3) - 18*qd1^2*qdd1*qdd3*sin(q3) - 36*qd2^2*qdd1*qdd3*sin(q3) - 36*qd3^2*qdd1*qdd2*sin(q3) - 36*qd3^2*qdd1*qdd3*sin(q3) + 9*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2) - 9*qd1^4*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 72*qd2*qd3*qdd1*qdd2*sin(q3) + 72*qd1*qd3*qdd2*qdd3*sin(q3) - 72*qd2*qd3*qdd1*qdd3*sin(q3) + 32*g*qd1^2*qdd1*cos(q2)^3*sin(q3) + 64*g*qd2^2*qdd1*cos(q2)^3*sin(q3) - 72*qd1*qd2*qdd2^2*sin(2*q2)*cos(q3) - 72*qd1*qd2*qdd3^2*sin(2*q2)*cos(q3) - 18*qd1^2*qdd1*qdd2*cos(2*q2)*sin(q3) - 18*qd1^2*qdd1*qdd2*sin(2*q2)*cos(q3) - 18*qd1^2*qdd1*qdd3*cos(2*q2)*sin(q3) - 36*qd1^2*qdd1*qdd3*sin(2*q2)*cos(q3) - 36*qd2^2*qdd1*qdd2*cos(2*q2)*sin(q3) - 36*qd2^2*qdd1*qdd3*cos(2*q2)*sin(q3) - 36*qd3^2*qdd1*qdd2*cos(2*q2)*sin(q3) - 36*qd3^2*qdd1*qdd3*cos(2*q2)*sin(q3) - 9*qd1^4*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 9*qd1^4*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) - 18*qd1^4*qdd1*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 128*g*qd1*qd2^3*sin(q2)^3*sin(q3) + 64*g*qd1^3*qd2*sin(q2)^3*sin(q3) + 144*qd1*qd2*qdd2*qdd3*sin(2*q2 + q3) + 72*qd1*qd3*qdd2*qdd3*sin(2*q2 + q3) + 72*qd1*qd2^3*qdd2*sin(2*q2)*sin(q3) + 36*qd1^3*qd2*qdd2*sin(2*q2)*sin(q3) + 72*qd1*qd2^3*qdd3*sin(2*q2)*sin(q3) + 36*qd1^3*qd2*qdd3*sin(2*q2)*sin(q3) + 36*qd1^3*qd3*qdd2*sin(2*q2)*sin(q3) + 36*qd1^3*qd3*qdd3*sin(2*q2)*sin(q3) + 32*g^2*qdd1*cos(2*q2)*cos(q2 + q3)*cos(q2) + 24*g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2) - 18*g*qdd1*qdd2*cos(q2)*cos(q3) - 36*g*qdd1*qdd3*cos(q2)*cos(q3) - 72*qd1*qd2*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) - 72*qd1*qd3*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) - 36*qd1^2*qdd1*qdd2*cos(2*q2 + q3)*sin(2*q2) + 18*qd1^2*qdd1*qdd2*sin(2*q2 + q3)*cos(2*q2) + 18*qd1^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) - 18*qd1^2*qdd1*qdd2*sin(2*q2 + 2*q3)*cos(q3) - 36*qd1^2*qdd1*qdd3*cos(2*q2 + q3)*sin(2*q2) + 18*qd1^2*qdd1*qdd3*sin(2*q2 + q3)*cos(2*q2) + 18*qd1^2*qdd1*qdd3*sin(2*q2 + 2*q3)*cos(q3) + 36*qd2^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) - 18*qd1^4*qdd1*cos(2*q2 + q3)*sin(2*q2)*sin(2*q2 + 2*q3) + 9*qd1^4*qdd1*sin(2*q2 + q3)*cos(2*q2)*sin(2*q2 + 2*q3) + 9*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2)*cos(2*q2 + 2*q3) + 192*g*qd1*qd2*qdd2*sin(q2)^3 + 192*g*qd1*qd2*qdd3*sin(q2)^3 + 36*qd1^3*qd2*qdd2*sin(2*q2 + q3)*sin(2*q2) - 72*qd1*qd2^3*qdd2*sin(2*q2 + 2*q3)*sin(q3) - 36*qd1^3*qd2*qdd2*sin(2*q2 + 2*q3)*sin(q3) + 36*qd1^3*qd2*qdd3*sin(2*q2 + q3)*sin(2*q2) + 36*qd1^3*qd3*qdd2*sin(2*q2 + q3)*sin(2*q2) + 36*qd1^3*qd3*qdd3*sin(2*q2 + q3)*sin(2*q2) + 18*g^2*qdd1*cos(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 36*g*qdd1*qdd2*cos(2*q2 + q3)*cos(q2) - 36*g*qdd1*qdd3*cos(2*q2 + q3)*cos(q2) + 18*qd1^2*qd2^2*qdd1*sin(2*q2)*sin(q3) + 64*g*qdd1*qdd2*cos(q2)^3*cos(q3) + 36*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2 + q3) + 72*qdd1*qdd2*qdd3*cos(2*q2)*cos(q3) + 9*g*qd1^2*qdd1*cos(q2)*sin(q3) + 18*g*qd2^2*qdd1*cos(q2)*sin(q3) - 18*qd1^2*qd2^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 18*qd1^2*qd3^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 128*g*qd1*qd2^3*sin(q2)*sin(q3) - 64*g*qd1^3*qd2*sin(q2)*sin(q3) + 25*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(q2) + 18*g*qd1^2*qdd1*sin(2*q2)*cos(q2 + q3) - 192*g*qd1*qd2*qdd2*sin(q2) - 192*g*qd1*qd2*qdd3*sin(q2) + 72*g*qd1*qd2*qdd2*sin(2*q2 + q3)*cos(q2) + 72*g*qd1*qd2*qdd3*sin(2*q2 + q3)*cos(q2) + 36*g*qd1*qd3*qdd2*sin(2*q2 + q3)*cos(q2) + 36*g*qd1*qd3*qdd3*sin(2*q2 + q3)*cos(q2) + 36*qd1^3*qd2*qd3^2*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 36*qd1^3*qd2^2*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 16*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2)*cos(q2) + 9*g*qd1^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) - 18*g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) + 18*g*qd2^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) - 72*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) - 72*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) - 32*g*qd1^3*qd2*sin(2*q2 + q3)*sin(2*q2)*cos(q2) - 36*g*qd1*qd2^3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) - 18*g*qd1^3*qd2*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) - 144*qd1*qd2*qdd2*qdd3*sin(2*q2)*cos(q3) - 72*qd2*qd3*qdd1*qdd2*cos(2*q2)*sin(q3) - 72*qd2*qd3*qdd1*qdd3*cos(2*q2)*sin(q3) - 18*g*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) + 9*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2 + 2*q3)*cos(q2) + 18*g*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q2 + q3) + 18*g*qd1^3*qd2*sin(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) - 36*g*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) - 36*g*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) + 24*g*qd1^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q2) + 18*qd1^2*qdd1*qdd2*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 18*qd1^2*qdd1*qdd2*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q3) + 18*qd1^2*qdd1*qdd3*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 48*g*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2) - 72*qd1^3*qd2*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 36*qd1^3*qd2*qdd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 36*qd1^3*qd3*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 72*qd1*qd2*qd3^2*qdd2*sin(2*q2)*sin(q3) + 144*qd1*qd2^2*qd3*qdd2*sin(2*q2)*sin(q3) + 72*qd1*qd2*qd3^2*qdd3*sin(2*q2)*sin(q3) + 144*qd1*qd2^2*qd3*qdd3*sin(2*q2)*sin(q3) + 18*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2)*cos(q3) - 64*g^2*qd1*qd2*sin(2*q2)*cos(q2 + q3)*cos(q2) + 36*g*qd1*qd3*qdd2*cos(q2)*sin(q3) + 36*g*qd1*qd3*qdd3*cos(q2)*sin(q3) - 36*qd1^2*qd2*qd3*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 72*qd1*qd2^2*qd3*qdd2*sin(2*q2 + 2*q3)*sin(q3) - 36*g^2*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 36*g^2*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 18*qd1^2*qd2^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 18*qd1^2*qd2^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) - 18*qd1^2*qd3^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 128*g*qd1*qd2*qdd2*cos(q2)^2*cos(q3)*sin(q2) - 36*qd1^2*qd2*qd3*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 36*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) - 36*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) - 36*g*qd1*qd2^2*qd3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3))/(6*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))), -1, -(2*((3*qd1^4*qdd1*cos(q3))/2 - (3*qd1^4*qdd1*cos(4*q2 + 3*q3))/2 + 12*qdd1*qdd3^2*cos(q3) - 24*qdd1*qdd2^2*cos(2*q2 + q3) + 6*g^2*qdd1*cos(q2 + q3)*cos(q2) + 48*qd1*qd2*qdd2^2*sin(2*q2 + q3) + 24*qd1*qd3*qdd2^2*sin(2*q2 + q3) + 12*qd1^2*qdd1*qdd2*sin(2*q2 + q3) + 6*qd1^2*qdd1*qdd3*sin(2*q2 + q3) + 12*g*qdd1*qdd2*cos(q2 + q3) - 6*qd1^3*qd2*qdd2*cos(4*q2 + 3*q3) - 3*qd1^2*qdd1*qdd2*sin(4*q2 + 3*q3) + 12*qdd1*qdd2^2*cos(2*q2)*cos(q3) + 12*qdd1*qdd3^2*cos(2*q2)*cos(q3) - 24*qdd1*qdd2*qdd3*cos(2*q2 + q3) + 3*qd1^4*qdd1*sin(2*q2)*sin(q3) + 12*qdd1*qdd2^2*cos(2*q2 + 2*q3)*cos(q3) + 16*g*qdd1*qdd2*cos(q2)^3 + 16*g*qdd1*qdd3*cos(q2)^3 + 6*qd1^3*qd2*qdd2*cos(q3) + 24*qd1*qd3*qdd2^2*sin(q3) - 9*qd1^2*qdd1*qdd2*sin(q3) - 6*qd1^2*qdd1*qdd3*sin(q3) - 12*qd2^2*qdd1*qdd3*sin(q3) - 12*qd3^2*qdd1*qdd2*sin(q3) - 12*qd3^2*qdd1*qdd3*sin(q3) + 3*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2) - 3*qd1^4*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 24*qd2*qd3*qdd1*qdd2*sin(q3) + 24*qd1*qd3*qdd2*qdd3*sin(q3) - 24*qd2*qd3*qdd1*qdd3*sin(q3) - 24*qd1*qd2*qdd2^2*sin(2*q2)*cos(q3) - 24*qd1*qd2*qdd3^2*sin(2*q2)*cos(q3) - 6*qd1^2*qdd1*qdd2*cos(2*q2)*sin(q3) - 6*qd1^2*qdd1*qdd2*sin(2*q2)*cos(q3) - 6*qd1^2*qdd1*qdd3*cos(2*q2)*sin(q3) - 12*qd1^2*qdd1*qdd3*sin(2*q2)*cos(q3) - 12*qd2^2*qdd1*qdd2*cos(2*q2)*sin(q3) - 12*qd2^2*qdd1*qdd3*cos(2*q2)*sin(q3) - 12*qd3^2*qdd1*qdd2*cos(2*q2)*sin(q3) - 12*qd3^2*qdd1*qdd3*cos(2*q2)*sin(q3) - 3*qd1^4*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 3*qd1^4*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) - 6*qd1^4*qdd1*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 48*qd1*qd2*qdd2*qdd3*sin(2*q2 + q3) + 24*qd1*qd3*qdd2*qdd3*sin(2*q2 + q3) + 24*qd1*qd2^3*qdd2*sin(2*q2)*sin(q3) + 12*qd1^3*qd2*qdd2*sin(2*q2)*sin(q3) + 24*qd1*qd2^3*qdd3*sin(2*q2)*sin(q3) + 12*qd1^3*qd2*qdd3*sin(2*q2)*sin(q3) + 12*qd1^3*qd3*qdd2*sin(2*q2)*sin(q3) + 12*qd1^3*qd3*qdd3*sin(2*q2)*sin(q3) + 4*g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2) - 6*g*qdd1*qdd2*cos(q2)*cos(q3) - 12*g*qdd1*qdd3*cos(q2)*cos(q3) - 24*qd1*qd2*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) - 24*qd1*qd3*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) - 12*qd1^2*qdd1*qdd2*cos(2*q2 + q3)*sin(2*q2) + 6*qd1^2*qdd1*qdd2*sin(2*q2 + q3)*cos(2*q2) + 6*qd1^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) - 6*qd1^2*qdd1*qdd2*sin(2*q2 + 2*q3)*cos(q3) - 12*qd1^2*qdd1*qdd3*cos(2*q2 + q3)*sin(2*q2) + 6*qd1^2*qdd1*qdd3*sin(2*q2 + q3)*cos(2*q2) + 6*qd1^2*qdd1*qdd3*sin(2*q2 + 2*q3)*cos(q3) + 12*qd2^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) - 6*qd1^4*qdd1*cos(2*q2 + q3)*sin(2*q2)*sin(2*q2 + 2*q3) + 3*qd1^4*qdd1*sin(2*q2 + q3)*cos(2*q2)*sin(2*q2 + 2*q3) + 3*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2)*cos(2*q2 + 2*q3) + 32*g*qd1*qd2*qdd2*sin(q2)^3 + 32*g*qd1*qd2*qdd3*sin(q2)^3 + 12*qd1^3*qd2*qdd2*sin(2*q2 + q3)*sin(2*q2) - 24*qd1*qd2^3*qdd2*sin(2*q2 + 2*q3)*sin(q3) - 12*qd1^3*qd2*qdd2*sin(2*q2 + 2*q3)*sin(q3) + 12*qd1^3*qd2*qdd3*sin(2*q2 + q3)*sin(2*q2) + 12*qd1^3*qd3*qdd2*sin(2*q2 + q3)*sin(2*q2) + 12*qd1^3*qd3*qdd3*sin(2*q2 + q3)*sin(2*q2) + 6*g^2*qdd1*cos(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 12*g*qdd1*qdd2*cos(2*q2 + q3)*cos(q2) - 12*g*qdd1*qdd3*cos(2*q2 + q3)*cos(q2) + 6*qd1^2*qd2^2*qdd1*sin(2*q2)*sin(q3) + 12*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2 + q3) + 24*qdd1*qdd2*qdd3*cos(2*q2)*cos(q3) + 3*g*qd1^2*qdd1*cos(q2)*sin(q3) + 6*g*qd2^2*qdd1*cos(q2)*sin(q3) - 6*qd1^2*qd2^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 6*qd1^2*qd3^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 3*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(q2) + 6*g*qd1^2*qdd1*sin(2*q2)*cos(q2 + q3) - 32*g*qd1*qd2*qdd2*sin(q2) - 32*g*qd1*qd2*qdd3*sin(q2) + 24*g*qd1*qd2*qdd2*sin(2*q2 + q3)*cos(q2) + 24*g*qd1*qd2*qdd3*sin(2*q2 + q3)*cos(q2) + 12*g*qd1*qd3*qdd2*sin(2*q2 + q3)*cos(q2) + 12*g*qd1*qd3*qdd3*sin(2*q2 + q3)*cos(q2) + 12*qd1^3*qd2*qd3^2*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 12*qd1^3*qd2^2*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 3*g*qd1^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) - 6*g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) + 6*g*qd2^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) - 24*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) - 24*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) - 12*g*qd1*qd2^3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) - 6*g*qd1^3*qd2*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) - 48*qd1*qd2*qdd2*qdd3*sin(2*q2)*cos(q3) - 24*qd2*qd3*qdd1*qdd2*cos(2*q2)*sin(q3) - 24*qd2*qd3*qdd1*qdd3*cos(2*q2)*sin(q3) - 6*g*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) + 3*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2 + 2*q3)*cos(q2) + 6*g*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q2 + q3) + 6*g*qd1^3*qd2*sin(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) - 12*g*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) - 12*g*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) + 4*g*qd1^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q2) + 6*qd1^2*qdd1*qdd2*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 6*qd1^2*qdd1*qdd2*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q3) + 6*qd1^2*qdd1*qdd3*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 8*g*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2) - 24*qd1^3*qd2*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 12*qd1^3*qd2*qdd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 12*qd1^3*qd3*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 24*qd1*qd2*qd3^2*qdd2*sin(2*q2)*sin(q3) + 48*qd1*qd2^2*qd3*qdd2*sin(2*q2)*sin(q3) + 24*qd1*qd2*qd3^2*qdd3*sin(2*q2)*sin(q3) + 48*qd1*qd2^2*qd3*qdd3*sin(2*q2)*sin(q3) + 6*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2)*cos(q3) + 12*g*qd1*qd3*qdd2*cos(q2)*sin(q3) + 12*g*qd1*qd3*qdd3*cos(q2)*sin(q3) - 12*qd1^2*qd2*qd3*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 24*qd1*qd2^2*qd3*qdd2*sin(2*q2 + 2*q3)*sin(q3) - 12*g^2*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 12*g^2*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 6*qd1^2*qd2^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 6*qd1^2*qd2^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) - 6*qd1^2*qd3^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 12*qd1^2*qd2*qd3*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 12*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) - 12*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) - 12*g*qd1*qd2^2*qd3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3)))/(3*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))), -(2*qdd1*qdd2 + g*qdd1*cos(q2) - g*qdd1*cos(q2)^3 - 2*qdd1*qdd2*cos(q2)^2 + qd1^2*qdd1*sin(2*q2) + 2*g*qd1*qd2*(sin(q2) - sin(q2)^3) + 2*qd1*qd2*qdd2*sin(2*q2))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))), -(18*qdd1*qdd2^2 + 3*qd1^4*qdd1*cos(4*q2 + 3*q3) - 18*qdd1*qdd2^2*cos(2*q2 + 2*q3) - 3*qd1^4*qdd1*cos(q3) + 24*qdd1*qdd2^2*cos(q3) + 18*qdd1*qdd2*qdd3 + 12*g^2*qdd1*cos(q2 + q3)*cos(q2) - 18*qdd1*qdd2*qdd3*cos(2*q2 + 2*q3) + 12*qd1^2*qdd1*qdd2*sin(2*q2 + q3) + 24*g*qdd1*qdd2*cos(q2 + q3) + 9*qd1^2*qdd1*qdd2*sin(2*q2) + 9*qd1^2*qdd1*qdd3*sin(2*q2) + 9*qd1^4*qdd1*sin(2*q2)*sin(2*q2 + 2*q3) + 12*qd1^3*qd2*qdd2*cos(4*q2 + 3*q3) + 9*g*qdd1*qdd2*cos(q2) + 9*g*qdd1*qdd3*cos(q2) + 36*qd1*qd2*qdd2^2*sin(2*q2 + 2*q3) + 36*qd1*qd3*qdd2^2*sin(2*q2 + 2*q3) + 18*qd1^2*qdd1*qdd2*sin(2*q2 + 2*q3) + 6*qd1^2*qdd1*qdd2*sin(4*q2 + 3*q3) + 6*qd1^4*qdd1*sin(2*q2)*sin(q3) - 24*qdd1*qdd2^2*cos(2*q2 + 2*q3)*cos(q3) - 12*qd1^3*qd2*qdd2*cos(q3) + 30*qd1^2*qdd1*qdd2*sin(q3) + 24*qd2^2*qdd1*qdd2*sin(q3) + 6*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2) + 6*qd1^4*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 12*qd1^2*qdd1*qdd2*sin(2*q2)*cos(q3) + 6*qd1^4*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 6*qd1^4*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) + 12*qd1^4*qdd1*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2) + 12*g*qdd1*qdd2*cos(q2)*cos(q3) + 48*qd1*qd2*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) + 48*qd1*qd3*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) - 12*qd1^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) + 12*qd1^2*qdd1*qdd2*sin(2*q2 + 2*q3)*cos(q3) - 12*qd1^2*qdd1*qdd3*sin(2*q2 + 2*q3)*cos(q3) - 24*qd2^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) + 12*qd1^4*qdd1*cos(2*q2 + q3)*sin(2*q2)*sin(2*q2 + 2*q3) - 6*qd1^4*qdd1*sin(2*q2 + q3)*cos(2*q2)*sin(2*q2 + 2*q3) - 6*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2)*cos(2*q2 + 2*q3) + 48*qd1*qd2^3*qdd2*sin(2*q2 + 2*q3)*sin(q3) + 24*qd1^3*qd2*qdd2*sin(2*q2 + 2*q3)*sin(q3) - 12*g^2*qdd1*cos(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) + 36*qd1*qd2*qdd2*qdd3*sin(2*q2 + 2*q3) + 36*qd1*qd3*qdd2*qdd3*sin(2*q2 + 2*q3) + 12*qd1^2*qd2^2*qdd1*sin(2*q2)*sin(q3) - 24*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2 + q3) - 9*qd1^2*qdd1*qdd2*sin(2*q2)*cos(2*q2 + 2*q3) - 9*qd1^2*qdd1*qdd3*sin(2*q2)*cos(2*q2 + 2*q3) + 6*g*qd1^2*qdd1*cos(q2)*sin(q3) + 12*g*qd2^2*qdd1*cos(q2)*sin(q3) + 18*qd1^3*qd2*qdd2*sin(2*q2)*sin(2*q2 + 2*q3) + 18*qd1^3*qd2*qdd3*sin(2*q2)*sin(2*q2 + 2*q3) + 18*qd1^3*qd3*qdd2*sin(2*q2)*sin(2*q2 + 2*q3) + 18*qd1^3*qd3*qdd3*sin(2*q2)*sin(2*q2 + 2*q3) + 12*qd1^2*qd2^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 12*qd1^2*qd3^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 9*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2) - 9*g*qdd1*qdd3*cos(2*q2 + 2*q3)*cos(q2) + 6*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(q2) + 12*g*qd1^2*qdd1*sin(2*q2)*cos(q2 + q3) - 24*qd1^3*qd2*qd3^2*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 24*qd1^3*qd2^2*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 6*g*qd1^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) + 12*g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) - 12*g*qd2^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) + 48*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) + 48*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) + 24*g*qd1*qd2^3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) + 12*g*qd1^3*qd2*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) + 12*g*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) - 6*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2 + 2*q3)*cos(q2) - 12*g*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q2 + q3) + 18*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2) + 18*g*qd1*qd2*qdd3*sin(2*q2 + 2*q3)*cos(q2) + 18*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2) + 18*g*qd1*qd3*qdd3*sin(2*q2 + 2*q3)*cos(q2) - 12*g*qd1^3*qd2*sin(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) + 24*g*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) + 24*g*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) - 8*g*qd1^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q2) - 12*qd1^2*qdd1*qdd2*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 12*qd1^2*qdd1*qdd2*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q3) - 12*qd1^2*qdd1*qdd3*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 16*g*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2) + 48*qd1^3*qd2*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 24*qd1^3*qd2*qdd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 24*qd1^3*qd3*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 12*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2)*cos(q3) + 24*qd1^2*qd2*qd3*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 48*qd1*qd2^2*qd3*qdd2*sin(2*q2 + 2*q3)*sin(q3) + 24*g^2*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) + 24*g^2*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) + 12*qd1^2*qd2^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 12*qd1^2*qd2^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) + 12*qd1^2*qd3^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 24*qd1^2*qd2*qd3*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 24*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) + 24*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) + 24*g*qd1*qd2^2*qd3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3))/(3*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))), -(qdd1*(2*cos(q2)*sin(q2)*qd1^2 + 2*qdd2 + g*cos(q2)))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))), -(2*qdd2 + g*cos(q2))/(g*cos(q2)), -(18*qdd1*qdd2^2 - 3*qd1^4*qdd1*cos(4*q2 + 3*q3) + 18*qdd1*qdd2^2*cos(2*q2 + 2*q3) + 3*qd1^4*qdd1*cos(q3) + 24*qdd1*qdd2^2*cos(q3) + 18*qdd1*qdd2*qdd3 + 12*g^2*qdd1*cos(q2 + q3)*cos(q2) + 18*qdd1*qdd2*qdd3*cos(2*q2 + 2*q3) + 12*qd1^2*qdd1*qdd2*sin(2*q2 + q3) + 24*g*qdd1*qdd2*cos(q2 + q3) + 9*qd1^2*qdd1*qdd2*sin(2*q2) + 9*qd1^2*qdd1*qdd3*sin(2*q2) - 12*qd1^3*qd2*qdd2*cos(4*q2 + 3*q3) + 9*g*qdd1*qdd2*cos(q2) + 9*g*qdd1*qdd3*cos(q2) - 36*qd1*qd2*qdd2^2*sin(2*q2 + 2*q3) - 36*qd1*qd3*qdd2^2*sin(2*q2 + 2*q3) - 6*qd1^2*qdd1*qdd2*sin(4*q2 + 3*q3) + 6*qd1^4*qdd1*sin(2*q2)*sin(q3) + 24*qdd1*qdd2^2*cos(2*q2 + 2*q3)*cos(q3) + 12*qd1^3*qd2*qdd2*cos(q3) - 6*qd1^2*qdd1*qdd2*sin(q3) + 24*qd2^2*qdd1*qdd2*sin(q3) + 6*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2) - 6*qd1^4*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 12*qd1^2*qdd1*qdd2*sin(2*q2)*cos(q3) - 6*qd1^4*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 6*qd1^4*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) - 12*qd1^4*qdd1*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 8*g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2) + 12*g*qdd1*qdd2*cos(q2)*cos(q3) - 48*qd1*qd2*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) - 48*qd1*qd3*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) + 12*qd1^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) - 12*qd1^2*qdd1*qdd2*sin(2*q2 + 2*q3)*cos(q3) + 12*qd1^2*qdd1*qdd3*sin(2*q2 + 2*q3)*cos(q3) + 24*qd2^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) - 12*qd1^4*qdd1*cos(2*q2 + q3)*sin(2*q2)*sin(2*q2 + 2*q3) + 6*qd1^4*qdd1*sin(2*q2 + q3)*cos(2*q2)*sin(2*q2 + 2*q3) + 6*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2)*cos(2*q2 + 2*q3) - 48*qd1*qd2^3*qdd2*sin(2*q2 + 2*q3)*sin(q3) - 24*qd1^3*qd2*qdd2*sin(2*q2 + 2*q3)*sin(q3) + 12*g^2*qdd1*cos(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 36*qd1*qd2*qdd2*qdd3*sin(2*q2 + 2*q3) - 36*qd1*qd3*qdd2*qdd3*sin(2*q2 + 2*q3) + 12*qd1^2*qd2^2*qdd1*sin(2*q2)*sin(q3) + 24*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2 + q3) + 9*qd1^2*qdd1*qdd2*sin(2*q2)*cos(2*q2 + 2*q3) + 9*qd1^2*qdd1*qdd3*sin(2*q2)*cos(2*q2 + 2*q3) + 6*g*qd1^2*qdd1*cos(q2)*sin(q3) + 12*g*qd2^2*qdd1*cos(q2)*sin(q3) - 18*qd1^3*qd2*qdd2*sin(2*q2)*sin(2*q2 + 2*q3) - 18*qd1^3*qd2*qdd3*sin(2*q2)*sin(2*q2 + 2*q3) - 18*qd1^3*qd3*qdd2*sin(2*q2)*sin(2*q2 + 2*q3) - 18*qd1^3*qd3*qdd3*sin(2*q2)*sin(2*q2 + 2*q3) - 12*qd1^2*qd2^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 12*qd1^2*qd3^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 9*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2) + 9*g*qdd1*qdd3*cos(2*q2 + 2*q3)*cos(q2) + 6*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(q2) + 12*g*qd1^2*qdd1*sin(2*q2)*cos(q2 + q3) + 24*qd1^3*qd2*qd3^2*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 24*qd1^3*qd2^2*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 6*g*qd1^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) - 12*g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) + 12*g*qd2^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) - 48*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) - 48*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) - 24*g*qd1*qd2^3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) - 12*g*qd1^3*qd2*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) - 12*g*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) + 6*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2 + 2*q3)*cos(q2) + 12*g*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q2 + q3) - 18*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2) - 18*g*qd1*qd2*qdd3*sin(2*q2 + 2*q3)*cos(q2) - 18*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2) - 18*g*qd1*qd3*qdd3*sin(2*q2 + 2*q3)*cos(q2) + 12*g*qd1^3*qd2*sin(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) - 24*g*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) - 24*g*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) + 8*g*qd1^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q2) + 12*qd1^2*qdd1*qdd2*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 12*qd1^2*qdd1*qdd2*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q3) + 12*qd1^2*qdd1*qdd3*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 16*g*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2) - 48*qd1^3*qd2*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 24*qd1^3*qd2*qdd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 24*qd1^3*qd3*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 12*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2)*cos(q3) - 24*qd1^2*qd2*qd3*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 48*qd1*qd2^2*qd3*qdd2*sin(2*q2 + 2*q3)*sin(q3) - 24*g^2*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 24*g^2*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 12*qd1^2*qd2^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 12*qd1^2*qd2^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) - 12*qd1^2*qd3^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 24*qd1^2*qd2*qd3*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 24*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) - 24*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) - 24*g*qd1*qd2^2*qd3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3))/(3*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))), (2*qdd2)/(g*cos(q2)), -((qdd2 + qdd3)*(24*qdd1*qdd3*cos(q3) - 9*g*qdd1*cos(q2) - 24*qdd1*qdd2*cos(q3) - 18*qdd1*qdd2 - 48*qdd1*qdd2*cos(2*q2 + q3) + 32*g*qdd1*cos(q2)^3 - 12*qd1^2*qdd1*sin(q3) - 24*qd2^2*qdd1*sin(q3) - 24*qd3^2*qdd1*sin(q3) - 18*qdd1*qdd2*cos(2*q2 + 2*q3) + 12*qd1^2*qdd1*sin(2*q2 + q3) - 9*qd1^2*qdd1*sin(2*q2) - 24*g*qdd1*cos(2*q2 + q3)*cos(q2) + 36*qd1*qd2*qdd2*sin(2*q2 + 2*q3) + 36*qd1*qd3*qdd2*sin(2*q2 + 2*q3) - 9*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3) + 24*qdd1*qdd2*cos(2*q2)*cos(q3) + 24*qdd1*qdd3*cos(2*q2)*cos(q3) + 18*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3) + 18*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3) - 9*g*qdd1*cos(2*q2 + 2*q3)*cos(q2) - 64*g*qd1*qd2*sin(q2) + 48*qd1*qd3*qdd2*sin(q3) - 48*qd2*qd3*qdd1*sin(q3) - 12*qd1^2*qdd1*cos(2*q2)*sin(q3) - 24*qd1^2*qdd1*sin(2*q2)*cos(q3) - 24*qd2^2*qdd1*cos(2*q2)*sin(q3) - 24*qd3^2*qdd1*cos(2*q2)*sin(q3) + 96*qd1*qd2*qdd2*sin(2*q2 + q3) + 48*qd1*qd3*qdd2*sin(2*q2 + q3) + 48*qd1*qd2^3*sin(2*q2)*sin(q3) + 24*qd1^3*qd2*sin(2*q2)*sin(q3) + 24*qd1^3*qd3*sin(2*q2)*sin(q3) - 24*g*qdd1*cos(q2)*cos(q3) - 24*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2) + 12*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2) + 64*g*qd1*qd2*sin(q2)^3 + 24*qd1^3*qd2*sin(2*q2 + q3)*sin(2*q2) + 24*qd1^3*qd3*sin(2*q2 + q3)*sin(2*q2) + 48*qd1*qd2*qd3^2*sin(2*q2)*sin(q3) + 96*qd1*qd2^2*qd3*sin(2*q2)*sin(q3) + 24*g*qd1*qd3*cos(q2)*sin(q3) + 48*g*qd1*qd2*sin(2*q2 + q3)*cos(q2) + 24*g*qd1*qd3*sin(2*q2 + q3)*cos(q2) - 48*qd1*qd2*qdd2*sin(2*q2)*cos(q3) - 48*qd1*qd2*qdd3*sin(2*q2)*cos(q3) - 48*qd2*qd3*qdd1*cos(2*q2)*sin(q3) + 18*g*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2) + 18*g*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2)))/(3*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))), -(2500*qdd1*(2*cos(q2)*sin(q2)*qd1^2 + 2*qdd2 + g*cos(q2)))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))), (1800*qdd2)/(g*cos(q2)), -(75*qdd3*(18*qdd1*qdd3 - 9*g*qdd1*cos(q2) + 9*qd1^2*qdd1*sin(2*q2 + 2*q3) + 24*qdd1*qdd3*cos(q3) - 48*qdd1*qdd2*cos(2*q2 + q3) + 32*g*qdd1*cos(q2)^3 + 18*qdd1*qdd2*cos(2*q2) + 18*qdd1*qdd3*cos(2*q2) - 24*qd3^2*qdd1*sin(q3) - 18*qdd1*qdd2*cos(2*q2 + 2*q3) + 24*qd1^2*qdd1*sin(2*q2 + q3) + 24*g*qdd1*cos(q2 + q3) - 9*qd1^2*qdd1*sin(2*q2) - 24*g*qdd1*cos(2*q2 + q3)*cos(q2) + 24*g*qdd1*cos(2*q2)*cos(q2 + q3) + 36*qd1*qd2*qdd2*sin(2*q2 + 2*q3) + 36*qd1*qd3*qdd2*sin(2*q2 + 2*q3) + 9*qd1^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3) - 9*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3) + 48*qdd1*qdd2*cos(2*q2)*cos(q3) + 24*qdd1*qdd3*cos(2*q2)*cos(q3) + 18*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3) - 9*g*qdd1*cos(2*q2 + 2*q3)*cos(q2) - 64*g*qd1*qd2*sin(q2) + 48*qd1*qd3*qdd2*sin(q3) - 48*qd2*qd3*qdd1*sin(q3) - 24*qd1^2*qdd1*sin(2*q2)*cos(q3) - 24*qd3^2*qdd1*cos(2*q2)*sin(q3) + 96*qd1*qd2*qdd2*sin(2*q2 + q3) + 48*qd1*qd3*qdd2*sin(2*q2 + q3) + 24*qd1^3*qd3*sin(2*q2)*sin(q3) - 24*g*qdd1*cos(q2)*cos(q3) - 24*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2) + 24*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2) + 64*g*qd1*qd2*sin(q2)^3 - 36*qd1*qd2*qdd2*sin(2*q2) - 36*qd1*qd2*qdd3*sin(2*q2) + 24*qd1^3*qd3*sin(2*q2 + q3)*sin(2*q2) + 48*qd1*qd2*qd3^2*sin(2*q2)*sin(q3) + 96*qd1*qd2^2*qd3*sin(2*q2)*sin(q3) + 24*g*qd1*qd3*cos(q2)*sin(q3) + 48*g*qd1*qd2*sin(2*q2 + q3)*cos(q2) - 48*g*qd1*qd2*sin(2*q2)*cos(q2 + q3) + 24*g*qd1*qd3*sin(2*q2 + q3)*cos(q2) - 96*qd1*qd2*qdd2*sin(2*q2)*cos(q3) - 48*qd1*qd2*qdd3*sin(2*q2)*cos(q3) - 48*qd2*qd3*qdd1*cos(2*q2)*sin(q3) + 18*g*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2) + 18*g*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2)))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))), -(2500*qd1*(2*cos(q2)*sin(q2)*qd1^2 + 2*qdd2 + g*cos(q2)))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))), (1800*qd2)/(g*cos(q2)), -(75*qd3*(18*qdd1*qdd3 - 9*g*qdd1*cos(q2) + 9*qd1^2*qdd1*sin(2*q2 + 2*q3) + 24*qdd1*qdd3*cos(q3) - 48*qdd1*qdd2*cos(2*q2 + q3) + 32*g*qdd1*cos(q2)^3 + 18*qdd1*qdd2*cos(2*q2) + 18*qdd1*qdd3*cos(2*q2) - 24*qd3^2*qdd1*sin(q3) - 18*qdd1*qdd2*cos(2*q2 + 2*q3) + 24*qd1^2*qdd1*sin(2*q2 + q3) + 24*g*qdd1*cos(q2 + q3) - 9*qd1^2*qdd1*sin(2*q2) - 24*g*qdd1*cos(2*q2 + q3)*cos(q2) + 24*g*qdd1*cos(2*q2)*cos(q2 + q3) + 36*qd1*qd2*qdd2*sin(2*q2 + 2*q3) + 36*qd1*qd3*qdd2*sin(2*q2 + 2*q3) + 9*qd1^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3) - 9*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3) + 48*qdd1*qdd2*cos(2*q2)*cos(q3) + 24*qdd1*qdd3*cos(2*q2)*cos(q3) + 18*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3) - 9*g*qdd1*cos(2*q2 + 2*q3)*cos(q2) - 64*g*qd1*qd2*sin(q2) + 48*qd1*qd3*qdd2*sin(q3) - 48*qd2*qd3*qdd1*sin(q3) - 24*qd1^2*qdd1*sin(2*q2)*cos(q3) - 24*qd3^2*qdd1*cos(2*q2)*sin(q3) + 96*qd1*qd2*qdd2*sin(2*q2 + q3) + 48*qd1*qd3*qdd2*sin(2*q2 + q3) + 24*qd1^3*qd3*sin(2*q2)*sin(q3) - 24*g*qdd1*cos(q2)*cos(q3) - 24*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2) + 24*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2) + 64*g*qd1*qd2*sin(q2)^3 - 36*qd1*qd2*qdd2*sin(2*q2) - 36*qd1*qd2*qdd3*sin(2*q2) + 24*qd1^3*qd3*sin(2*q2 + q3)*sin(2*q2) + 48*qd1*qd2*qd3^2*sin(2*q2)*sin(q3) + 96*qd1*qd2^2*qd3*sin(2*q2)*sin(q3) + 24*g*qd1*qd3*cos(q2)*sin(q3) + 48*g*qd1*qd2*sin(2*q2 + q3)*cos(q2) - 48*g*qd1*qd2*sin(2*q2)*cos(q2 + q3) + 24*g*qd1*qd3*sin(2*q2 + q3)*cos(q2) - 96*qd1*qd2*qdd2*sin(2*q2)*cos(q3) - 48*qd1*qd2*qdd3*sin(2*q2)*cos(q3) - 48*qd2*qd3*qdd1*cos(2*q2)*sin(q3) + 18*g*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2) + 18*g*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2)))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3)))]
% [ 0, 1, 0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              -(8*(6*qdd2 + 6*qdd3 + 2*qd1^2*sin(q3) + 4*qd2^2*sin(q3) + 2*qd1^2*sin(2*q2 + q3) + 4*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 4*qdd2*cos(q3)))/(3*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))),  0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               (32*((sin(2*q2 + 2*q3)*qd1^2)/2 + qdd2 + qdd3))/(3*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))),                                                                                                                                                                                                                 0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       -(16*qd1^2*sin(2*q2 + 2*q3))/(3*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))),                                                                                                       0,                                 0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         (16*qd1^2*sin(2*q2 + 2*q3))/(3*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))),                    0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               (32*(qdd2 + qdd3))/(3*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))),                                                                                                            0,                       0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         (2400*qdd3)/(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3)),                                                                                                           0,                      0,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         (2400*qd3)/(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))]
% [ 0, 0, 1,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             (3*(2*qdd1*qdd2^2*cos(2*q2 - q3) + 2*qdd1*qdd3^2*cos(2*q2 - q3) + 2*qdd1*qdd2^2*cos(2*q2 + 3*q3) + 3*g^2*qdd1*cos(q3) + 4*qdd1*qdd3^2*cos(q3) + 3*g^2*qdd1*cos(2*q2 + q3) - 4*qdd1*qdd2^2*cos(2*q2 + q3) + 2*qdd1*qdd3^2*cos(2*q2 + q3) + g^2*qdd1*cos(2*q2 + 3*q3) + g^2*qdd1*cos(4*q2 + 3*q3) + 3*g*qdd1*qdd2*cos(3*q2 + 3*q3) - g*qd1^3*qd3*cos(q2 - q3) + 2*g*qd1*qd2^3*cos(q2 + 3*q3) - 2*g*qd1*qd2^3*cos(3*q2 + q3) - g*qd1^3*qd3*cos(q2 + 3*q3) + g*qd1^3*qd3*cos(3*q2 + q3) + 4*qdd1*qdd2*qdd3*cos(2*q2 - q3) - 8*qd1*qd2^3*qdd2*cos(2*q2 + q3) - 4*qd1^3*qd2*qdd2*cos(2*q2 + q3) - 4*qd1*qd2^3*qdd3*cos(2*q2 + q3) - 2*qd1^3*qd2*qdd3*cos(2*q2 + q3) - 2*qd1^3*qd3*qdd2*cos(2*q2 + q3) - 2*qd1^3*qd3*qdd3*cos(2*q2 + q3) - qd1^3*qd2*qdd3*cos(4*q2 + q3) - qd1^3*qd3*qdd2*cos(4*q2 + q3) - 2*qd1^3*qd3*qdd3*cos(4*q2 + q3) - 2*g^2*qd1*qd2*sin(2*q2 + q3) - 2*g^2*qd1*qd3*sin(2*q2 + q3) + 2*g*qd1^2*qdd1*sin(q2 - q3) - 2*g*qd2^2*qdd1*sin(q2 - q3) - g*qd1^2*qdd1*sin(q2 + 3*q3) + g*qd1^2*qdd1*sin(3*q2 + q3) + g*qd2^2*qdd1*sin(q2 + 3*q3) - g*qd2^2*qdd1*sin(3*q2 + q3) + 8*qd1*qd2*qdd2^2*sin(2*q2 + q3) - 4*qd1*qd2*qdd3^2*sin(2*q2 + q3) + 4*qd1*qd3*qdd2^2*sin(2*q2 + q3) - 4*qd2^2*qdd1*qdd2*sin(2*q2 + q3) - 2*qd2^2*qdd1*qdd3*sin(2*q2 + q3) - 2*qd3^2*qdd1*qdd2*sin(2*q2 + q3) - (qd1^2*qdd1*qdd3*sin(4*q2 + q3))/2 - 2*qd3^2*qdd1*qdd3*sin(2*q2 + q3) - 2*qd1^3*qd2*qdd2*cos(3*q3) - qd1^3*qd2*qdd3*cos(3*q3) - qd1^3*qd3*qdd2*cos(3*q3) - qd1^2*qd2^2*qdd1*cos(q3) - (qd1^2*qd3^2*qdd1*cos(q3))/2 + g*qdd1*qdd2*cos(q2 + q3) - 8*g*qdd1*qdd3*cos(q2 + q3) + (qd1^2*qdd1*qdd3*sin(3*q3))/2 - qd1^3*qd2*qd3^2*sin(q3) - qd1^3*qd2^2*qd3*sin(q3) + 2*g*qd1*qd2^3*cos(3*q2 + 3*q3) + g*qd1^3*qd3*cos(5*q2 + 3*q3) + 4*qd1*qd2^3*qdd2*cos(2*q2 - q3) + 2*qd1^3*qd2*qdd2*cos(2*q2 - q3) + 4*qd1*qd2^3*qdd3*cos(2*q2 - q3) + 2*qd1^3*qd2*qdd3*cos(2*q2 - q3) + 2*qd1^3*qd3*qdd2*cos(2*q2 - q3) + 4*qd1*qd2^3*qdd2*cos(2*q2 + 3*q3) + 2*qd1^3*qd2*qdd2*cos(2*q2 + 3*q3) + 2*qd1^3*qd3*qdd3*cos(2*q2 - q3) + qd1^3*qd2*qdd3*cos(4*q2 + 3*q3) + qd1^3*qd3*qdd2*cos(4*q2 + 3*q3) - 2*qd1^2*qd2^2*qdd1*cos(2*q2 + q3) - qd1^2*qd3^2*qdd1*cos(2*q2 + q3) - (qd1^2*qd3^2*qdd1*cos(4*q2 + q3))/2 - 2*g^2*qd1*qd2*sin(2*q2 + 3*q3) - 2*g^2*qd1*qd3*sin(2*q2 + 3*q3) - 2*g^2*qd1*qd2*sin(4*q2 + 3*q3) - 2*g^2*qd1*qd3*sin(4*q2 + 3*q3) - g*qd1^2*qdd1*sin(3*q2 + 3*q3) + g*qd2^2*qdd1*sin(3*q2 + 3*q3) - 4*qd1*qd2*qdd2^2*sin(2*q2 - q3) - 4*qd1*qd2*qdd3^2*sin(2*q2 - q3) - 4*qd1*qd2*qdd2^2*sin(2*q2 + 3*q3) - 4*qd1*qd3*qdd2^2*sin(2*q2 + 3*q3) - qd1^2*qdd1*qdd3*sin(2*q2 - q3) + 2*qd2^2*qdd1*qdd2*sin(2*q2 - q3) + 2*qd2^2*qdd1*qdd3*sin(2*q2 - q3) + 2*qd3^2*qdd1*qdd2*sin(2*q2 - q3) + qd1^2*qdd1*qdd3*sin(2*q2 + 3*q3) + 2*qd2^2*qdd1*qdd2*sin(2*q2 + 3*q3) + 2*qd3^2*qdd1*qdd3*sin(2*q2 - q3) + (qd1^2*qdd1*qdd3*sin(4*q2 + 3*q3))/2 + qd1^3*qd2*qd3^2*sin(4*q2 + q3) + qd1^3*qd2^2*qd3*sin(4*q2 + q3) + qd1^2*qd2^2*qdd1*cos(3*q3) + (qd1^2*qd3^2*qdd1*cos(3*q3))/2 - 2*g*qdd1*qdd2*cos(q2 - q3) - 4*g*qdd1*qdd3*cos(q2 - q3) + g*qdd1*qdd2*cos(q2 + 3*q3) - 3*g*qdd1*qdd2*cos(3*q2 + q3) - 4*g*qdd1*qdd3*cos(3*q2 + q3) - 2*g*qd1*qd2^3*cos(q2 + q3) + qd1^3*qd2*qd3^2*sin(3*q3) + qd1^3*qd2^2*qd3*sin(3*q3) - 4*qdd1*qdd2*qdd3*cos(2*q2 + q3) - g*qd1^2*qdd1*sin(q2 + q3) + g*qd2^2*qdd1*sin(q2 + q3) + qd1^2*qd2^2*qdd1*cos(2*q2 - q3) + qd1^2*qd2^2*qdd1*cos(2*q2 + 3*q3) + qd1^2*qd3^2*qdd1*cos(2*q2 + 3*q3) + (qd1^2*qd3^2*qdd1*cos(4*q2 + 3*q3))/2 - qd1^3*qd2*qd3^2*sin(4*q2 + 3*q3) - qd1^3*qd2^2*qd3*sin(4*q2 + 3*q3) + 2*qd1^3*qd2*qdd2*cos(q3) + qd1^3*qd2*qdd3*cos(q3) + qd1^3*qd3*qdd2*cos(q3) + 2*qd1^3*qd3*qdd3*cos(q3) - 2*g^2*qd1*qd2*sin(q3) - 2*g^2*qd1*qd3*sin(q3) + 8*qd1*qd3*qdd2^2*sin(q3) + (3*qd1^2*qdd1*qdd3*sin(q3))/2 - 4*qd2^2*qdd1*qdd3*sin(q3) - 4*qd3^2*qdd1*qdd2*sin(q3) - 4*qd3^2*qdd1*qdd3*sin(q3) - 8*qd2*qd3*qdd1*qdd2*sin(q3) + 8*qd1*qd3*qdd2*qdd3*sin(q3) - 8*qd2*qd3*qdd1*qdd3*sin(q3) - 2*g*qd1*qd2^2*qd3*cos(q2 + q3) - 4*g*qd1*qd3*qdd2*sin(q2 - q3) - 2*g*qd1*qd2*qdd2*sin(q2 + 3*q3) + 6*g*qd1*qd2*qdd2*sin(3*q2 + q3) - 4*g*qd1*qd3*qdd3*sin(q2 - q3) + 8*g*qd1*qd2*qdd3*sin(3*q2 + q3) - 2*g*qd1*qd3*qdd2*sin(q2 + 3*q3) + 2*g*qd1*qd3*qdd2*sin(3*q2 + q3) + 4*g*qd1*qd3*qdd3*sin(3*q2 + q3) + 8*qd1*qd2*qdd2*qdd3*sin(2*q2 + q3) - 4*qd2*qd3*qdd1*qdd2*sin(2*q2 + q3) + 8*qd1*qd3*qdd2*qdd3*sin(2*q2 + q3) - 4*qd2*qd3*qdd1*qdd3*sin(2*q2 + q3) - qd1^2*qd2*qd3*qdd1*cos(q3) + 2*g*qd1*qd2^2*qd3*cos(q2 + 3*q3) - 2*g*qd1*qd2^2*qd3*cos(3*q2 + q3) - 2*qd1^2*qd2*qd3*qdd1*cos(2*q2 + q3) - 4*qd1*qd2*qd3^2*qdd2*cos(2*q2 + q3) - 12*qd1*qd2^2*qd3*qdd2*cos(2*q2 + q3) - 4*qd1*qd2*qd3^2*qdd3*cos(2*q2 + q3) - 8*qd1*qd2^2*qd3*qdd3*cos(2*q2 + q3) - qd1^2*qd2*qd3*qdd1*cos(4*q2 + q3) - 6*g*qd1*qd2*qdd2*sin(3*q2 + 3*q3) - 6*g*qd1*qd3*qdd2*sin(3*q2 + 3*q3) - 8*qd1*qd2*qdd2*qdd3*sin(2*q2 - q3) + 4*qd2*qd3*qdd1*qdd2*sin(2*q2 - q3) + 4*qd2*qd3*qdd1*qdd3*sin(2*q2 - q3) + qd1^2*qd2*qd3*qdd1*cos(3*q3) + 2*g*qd1*qd2*qdd2*sin(q2 + q3) + 8*g*qd1*qd2*qdd3*sin(q2 + q3) + 2*g*qd1*qd3*qdd2*sin(q2 + q3) + 8*g*qd1*qd3*qdd3*sin(q2 + q3) + 2*g*qd1*qd2^2*qd3*cos(3*q2 + 3*q3) + 4*qd1*qd2*qd3^2*qdd2*cos(2*q2 - q3) + 8*qd1*qd2^2*qd3*qdd2*cos(2*q2 - q3) + 4*qd1*qd2*qd3^2*qdd3*cos(2*q2 - q3) + 8*qd1*qd2^2*qd3*qdd3*cos(2*q2 - q3) + 2*qd1^2*qd2*qd3*qdd1*cos(2*q2 + 3*q3) + 4*qd1*qd2^2*qd3*qdd2*cos(2*q2 + 3*q3) + qd1^2*qd2*qd3*qdd1*cos(4*q2 + 3*q3)))/(4*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))), -1,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                -(2*qdd1*qdd2^2*cos(2*q2 - q3) + 2*qdd1*qdd3^2*cos(2*q2 - q3) + 2*qdd1*qdd2^2*cos(2*q2 + 3*q3) + 3*g^2*qdd1*cos(q3) + 4*qdd1*qdd3^2*cos(q3) + 3*g^2*qdd1*cos(2*q2 + q3) - 4*qdd1*qdd2^2*cos(2*q2 + q3) + 2*qdd1*qdd3^2*cos(2*q2 + q3) + g^2*qdd1*cos(2*q2 + 3*q3) + g^2*qdd1*cos(4*q2 + 3*q3) + 3*g*qdd1*qdd2*cos(3*q2 + 3*q3) - g*qd1^3*qd3*cos(q2 - q3) + 2*g*qd1*qd2^3*cos(q2 + 3*q3) - 2*g*qd1*qd2^3*cos(3*q2 + q3) - g*qd1^3*qd3*cos(q2 + 3*q3) + g*qd1^3*qd3*cos(3*q2 + q3) + 4*qdd1*qdd2*qdd3*cos(2*q2 - q3) - 8*qd1*qd2^3*qdd2*cos(2*q2 + q3) - 4*qd1^3*qd2*qdd2*cos(2*q2 + q3) - 4*qd1*qd2^3*qdd3*cos(2*q2 + q3) - 2*qd1^3*qd2*qdd3*cos(2*q2 + q3) - 2*qd1^3*qd3*qdd2*cos(2*q2 + q3) - 2*qd1^3*qd3*qdd3*cos(2*q2 + q3) - qd1^3*qd2*qdd3*cos(4*q2 + q3) - qd1^3*qd3*qdd2*cos(4*q2 + q3) - 2*qd1^3*qd3*qdd3*cos(4*q2 + q3) - 2*g^2*qd1*qd2*sin(2*q2 + q3) - 2*g^2*qd1*qd3*sin(2*q2 + q3) + 2*g*qd1^2*qdd1*sin(q2 - q3) - 2*g*qd2^2*qdd1*sin(q2 - q3) - g*qd1^2*qdd1*sin(q2 + 3*q3) + g*qd1^2*qdd1*sin(3*q2 + q3) + g*qd2^2*qdd1*sin(q2 + 3*q3) - g*qd2^2*qdd1*sin(3*q2 + q3) + 8*qd1*qd2*qdd2^2*sin(2*q2 + q3) - 4*qd1*qd2*qdd3^2*sin(2*q2 + q3) + 4*qd1*qd3*qdd2^2*sin(2*q2 + q3) - 4*qd2^2*qdd1*qdd2*sin(2*q2 + q3) - 2*qd2^2*qdd1*qdd3*sin(2*q2 + q3) - 2*qd3^2*qdd1*qdd2*sin(2*q2 + q3) - (qd1^2*qdd1*qdd3*sin(4*q2 + q3))/2 - 2*qd3^2*qdd1*qdd3*sin(2*q2 + q3) - 2*qd1^3*qd2*qdd2*cos(3*q3) - qd1^3*qd2*qdd3*cos(3*q3) - qd1^3*qd3*qdd2*cos(3*q3) - qd1^2*qd2^2*qdd1*cos(q3) - (qd1^2*qd3^2*qdd1*cos(q3))/2 + g*qdd1*qdd2*cos(q2 + q3) - 8*g*qdd1*qdd3*cos(q2 + q3) + (qd1^2*qdd1*qdd3*sin(3*q3))/2 - qd1^3*qd2*qd3^2*sin(q3) - qd1^3*qd2^2*qd3*sin(q3) + 2*g*qd1*qd2^3*cos(3*q2 + 3*q3) + g*qd1^3*qd3*cos(5*q2 + 3*q3) + 4*qd1*qd2^3*qdd2*cos(2*q2 - q3) + 2*qd1^3*qd2*qdd2*cos(2*q2 - q3) + 4*qd1*qd2^3*qdd3*cos(2*q2 - q3) + 2*qd1^3*qd2*qdd3*cos(2*q2 - q3) + 2*qd1^3*qd3*qdd2*cos(2*q2 - q3) + 4*qd1*qd2^3*qdd2*cos(2*q2 + 3*q3) + 2*qd1^3*qd2*qdd2*cos(2*q2 + 3*q3) + 2*qd1^3*qd3*qdd3*cos(2*q2 - q3) + qd1^3*qd2*qdd3*cos(4*q2 + 3*q3) + qd1^3*qd3*qdd2*cos(4*q2 + 3*q3) - 2*qd1^2*qd2^2*qdd1*cos(2*q2 + q3) - qd1^2*qd3^2*qdd1*cos(2*q2 + q3) - (qd1^2*qd3^2*qdd1*cos(4*q2 + q3))/2 - 2*g^2*qd1*qd2*sin(2*q2 + 3*q3) - 2*g^2*qd1*qd3*sin(2*q2 + 3*q3) - 2*g^2*qd1*qd2*sin(4*q2 + 3*q3) - 2*g^2*qd1*qd3*sin(4*q2 + 3*q3) - g*qd1^2*qdd1*sin(3*q2 + 3*q3) + g*qd2^2*qdd1*sin(3*q2 + 3*q3) - 4*qd1*qd2*qdd2^2*sin(2*q2 - q3) - 4*qd1*qd2*qdd3^2*sin(2*q2 - q3) - 4*qd1*qd2*qdd2^2*sin(2*q2 + 3*q3) - 4*qd1*qd3*qdd2^2*sin(2*q2 + 3*q3) - qd1^2*qdd1*qdd3*sin(2*q2 - q3) + 2*qd2^2*qdd1*qdd2*sin(2*q2 - q3) + 2*qd2^2*qdd1*qdd3*sin(2*q2 - q3) + 2*qd3^2*qdd1*qdd2*sin(2*q2 - q3) + qd1^2*qdd1*qdd3*sin(2*q2 + 3*q3) + 2*qd2^2*qdd1*qdd2*sin(2*q2 + 3*q3) + 2*qd3^2*qdd1*qdd3*sin(2*q2 - q3) + (qd1^2*qdd1*qdd3*sin(4*q2 + 3*q3))/2 + qd1^3*qd2*qd3^2*sin(4*q2 + q3) + qd1^3*qd2^2*qd3*sin(4*q2 + q3) + qd1^2*qd2^2*qdd1*cos(3*q3) + (qd1^2*qd3^2*qdd1*cos(3*q3))/2 - 2*g*qdd1*qdd2*cos(q2 - q3) - 4*g*qdd1*qdd3*cos(q2 - q3) + g*qdd1*qdd2*cos(q2 + 3*q3) - 3*g*qdd1*qdd2*cos(3*q2 + q3) - 4*g*qdd1*qdd3*cos(3*q2 + q3) - 2*g*qd1*qd2^3*cos(q2 + q3) + qd1^3*qd2*qd3^2*sin(3*q3) + qd1^3*qd2^2*qd3*sin(3*q3) - 4*qdd1*qdd2*qdd3*cos(2*q2 + q3) - g*qd1^2*qdd1*sin(q2 + q3) + g*qd2^2*qdd1*sin(q2 + q3) + qd1^2*qd2^2*qdd1*cos(2*q2 - q3) + qd1^2*qd2^2*qdd1*cos(2*q2 + 3*q3) + qd1^2*qd3^2*qdd1*cos(2*q2 + 3*q3) + (qd1^2*qd3^2*qdd1*cos(4*q2 + 3*q3))/2 - qd1^3*qd2*qd3^2*sin(4*q2 + 3*q3) - qd1^3*qd2^2*qd3*sin(4*q2 + 3*q3) + 2*qd1^3*qd2*qdd2*cos(q3) + qd1^3*qd2*qdd3*cos(q3) + qd1^3*qd3*qdd2*cos(q3) + 2*qd1^3*qd3*qdd3*cos(q3) - 2*g^2*qd1*qd2*sin(q3) - 2*g^2*qd1*qd3*sin(q3) + 8*qd1*qd3*qdd2^2*sin(q3) + (3*qd1^2*qdd1*qdd3*sin(q3))/2 - 4*qd2^2*qdd1*qdd3*sin(q3) - 4*qd3^2*qdd1*qdd2*sin(q3) - 4*qd3^2*qdd1*qdd3*sin(q3) - 8*qd2*qd3*qdd1*qdd2*sin(q3) + 8*qd1*qd3*qdd2*qdd3*sin(q3) - 8*qd2*qd3*qdd1*qdd3*sin(q3) - 2*g*qd1*qd2^2*qd3*cos(q2 + q3) - 4*g*qd1*qd3*qdd2*sin(q2 - q3) - 2*g*qd1*qd2*qdd2*sin(q2 + 3*q3) + 6*g*qd1*qd2*qdd2*sin(3*q2 + q3) - 4*g*qd1*qd3*qdd3*sin(q2 - q3) + 8*g*qd1*qd2*qdd3*sin(3*q2 + q3) - 2*g*qd1*qd3*qdd2*sin(q2 + 3*q3) + 2*g*qd1*qd3*qdd2*sin(3*q2 + q3) + 4*g*qd1*qd3*qdd3*sin(3*q2 + q3) + 8*qd1*qd2*qdd2*qdd3*sin(2*q2 + q3) - 4*qd2*qd3*qdd1*qdd2*sin(2*q2 + q3) + 8*qd1*qd3*qdd2*qdd3*sin(2*q2 + q3) - 4*qd2*qd3*qdd1*qdd3*sin(2*q2 + q3) - qd1^2*qd2*qd3*qdd1*cos(q3) + 2*g*qd1*qd2^2*qd3*cos(q2 + 3*q3) - 2*g*qd1*qd2^2*qd3*cos(3*q2 + q3) - 2*qd1^2*qd2*qd3*qdd1*cos(2*q2 + q3) - 4*qd1*qd2*qd3^2*qdd2*cos(2*q2 + q3) - 12*qd1*qd2^2*qd3*qdd2*cos(2*q2 + q3) - 4*qd1*qd2*qd3^2*qdd3*cos(2*q2 + q3) - 8*qd1*qd2^2*qd3*qdd3*cos(2*q2 + q3) - qd1^2*qd2*qd3*qdd1*cos(4*q2 + q3) - 6*g*qd1*qd2*qdd2*sin(3*q2 + 3*q3) - 6*g*qd1*qd3*qdd2*sin(3*q2 + 3*q3) - 8*qd1*qd2*qdd2*qdd3*sin(2*q2 - q3) + 4*qd2*qd3*qdd1*qdd2*sin(2*q2 - q3) + 4*qd2*qd3*qdd1*qdd3*sin(2*q2 - q3) + qd1^2*qd2*qd3*qdd1*cos(3*q3) + 2*g*qd1*qd2*qdd2*sin(q2 + q3) + 8*g*qd1*qd2*qdd3*sin(q2 + q3) + 2*g*qd1*qd3*qdd2*sin(q2 + q3) + 8*g*qd1*qd3*qdd3*sin(q2 + q3) + 2*g*qd1*qd2^2*qd3*cos(3*q2 + 3*q3) + 4*qd1*qd2*qd3^2*qdd2*cos(2*q2 - q3) + 8*qd1*qd2^2*qd3*qdd2*cos(2*q2 - q3) + 4*qd1*qd2*qd3^2*qdd3*cos(2*q2 - q3) + 8*qd1*qd2^2*qd3*qdd3*cos(2*q2 - q3) + 2*qd1^2*qd2*qd3*qdd1*cos(2*q2 + 3*q3) + 4*qd1*qd2^2*qd3*qdd2*cos(2*q2 + 3*q3) + qd1^2*qd2*qd3*qdd1*cos(4*q2 + 3*q3))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))),   -(qdd1*qdd2 + g*qdd1*cos(q2) - g*qdd1*cos(q2)^3 - qdd1*qdd2*cos(q2)^2 + (qd1^2*qdd1*sin(2*q2))/2 + 2*g*qd1*qd2*(sin(q2) - sin(q2)^3) + qd1*qd2*qdd2*sin(2*q2))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))),                                                                                                                                                                 -(6*qdd1*qdd2^2 + qd1^4*qdd1*cos(4*q2 + 3*q3) - 6*qdd1*qdd2^2*cos(2*q2 + 2*q3) - qd1^4*qdd1*cos(q3) + 8*qdd1*qdd2^2*cos(q3) + 6*qdd1*qdd2*qdd3 + 8*g^2*qdd1*cos(q2 + q3)*cos(q2) - 6*qdd1*qdd2*qdd3*cos(2*q2 + 2*q3) + 4*qd1^2*qdd1*qdd2*sin(2*q2 + q3) + 8*g*qdd1*qdd2*cos(q2 + q3) + 3*qd1^2*qdd1*qdd2*sin(2*q2) + 3*qd1^2*qdd1*qdd3*sin(2*q2) + 3*qd1^4*qdd1*sin(2*q2)*sin(2*q2 + 2*q3) + 4*qd1^3*qd2*qdd2*cos(4*q2 + 3*q3) + 6*g*qdd1*qdd2*cos(q2) + 6*g*qdd1*qdd3*cos(q2) + 12*qd1*qd2*qdd2^2*sin(2*q2 + 2*q3) + 12*qd1*qd3*qdd2^2*sin(2*q2 + 2*q3) + 6*qd1^2*qdd1*qdd2*sin(2*q2 + 2*q3) + 2*qd1^2*qdd1*qdd2*sin(4*q2 + 3*q3) + 2*qd1^4*qdd1*sin(2*q2)*sin(q3) - 8*qdd1*qdd2^2*cos(2*q2 + 2*q3)*cos(q3) - 4*qd1^3*qd2*qdd2*cos(q3) + 10*qd1^2*qdd1*qdd2*sin(q3) + 8*qd2^2*qdd1*qdd2*sin(q3) + 2*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2) + 2*qd1^4*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 4*qd1^2*qdd1*qdd2*sin(2*q2)*cos(q3) + 2*qd1^4*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 2*qd1^4*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) + 4*qd1^4*qdd1*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 6*g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2) + 8*g*qdd1*qdd2*cos(q2)*cos(q3) + 16*qd1*qd2*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) + 16*qd1*qd3*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) - 4*qd1^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) + 4*qd1^2*qdd1*qdd2*sin(2*q2 + 2*q3)*cos(q3) - 4*qd1^2*qdd1*qdd3*sin(2*q2 + 2*q3)*cos(q3) - 8*qd2^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) + 4*qd1^4*qdd1*cos(2*q2 + q3)*sin(2*q2)*sin(2*q2 + 2*q3) - 2*qd1^4*qdd1*sin(2*q2 + q3)*cos(2*q2)*sin(2*q2 + 2*q3) - 2*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2)*cos(2*q2 + 2*q3) + 16*qd1*qd2^3*qdd2*sin(2*q2 + 2*q3)*sin(q3) + 8*qd1^3*qd2*qdd2*sin(2*q2 + 2*q3)*sin(q3) - 8*g^2*qdd1*cos(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) + 12*qd1*qd2*qdd2*qdd3*sin(2*q2 + 2*q3) + 12*qd1*qd3*qdd2*qdd3*sin(2*q2 + 2*q3) + 4*qd1^2*qd2^2*qdd1*sin(2*q2)*sin(q3) - 8*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2 + q3) - 3*qd1^2*qdd1*qdd2*sin(2*q2)*cos(2*q2 + 2*q3) - 3*qd1^2*qdd1*qdd3*sin(2*q2)*cos(2*q2 + 2*q3) + 4*g*qd1^2*qdd1*cos(q2)*sin(q3) + 8*g*qd2^2*qdd1*cos(q2)*sin(q3) + 6*qd1^3*qd2*qdd2*sin(2*q2)*sin(2*q2 + 2*q3) + 6*qd1^3*qd2*qdd3*sin(2*q2)*sin(2*q2 + 2*q3) + 6*qd1^3*qd3*qdd2*sin(2*q2)*sin(2*q2 + 2*q3) + 6*qd1^3*qd3*qdd3*sin(2*q2)*sin(2*q2 + 2*q3) + 4*qd1^2*qd2^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 4*qd1^2*qd3^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 6*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2) - 6*g*qdd1*qdd3*cos(2*q2 + 2*q3)*cos(q2) + 4*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(q2) + 4*g*qd1^2*qdd1*sin(2*q2)*cos(q2 + q3) - 8*qd1^3*qd2*qd3^2*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 8*qd1^3*qd2^2*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 4*g*qd1^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) + 8*g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) - 8*g*qd2^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) + 16*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) + 16*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) + 16*g*qd1*qd2^3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) + 8*g*qd1^3*qd2*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) + 8*g*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) - 4*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2 + 2*q3)*cos(q2) - 4*g*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q2 + q3) + 12*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2) + 12*g*qd1*qd2*qdd3*sin(2*q2 + 2*q3)*cos(q2) + 12*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2) + 12*g*qd1*qd3*qdd3*sin(2*q2 + 2*q3)*cos(q2) - 8*g*qd1^3*qd2*sin(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) + 8*g*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) + 8*g*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) - 4*qd1^2*qdd1*qdd2*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 4*qd1^2*qdd1*qdd2*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q3) - 4*qd1^2*qdd1*qdd3*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 16*qd1^3*qd2*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 8*qd1^3*qd2*qdd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 8*qd1^3*qd3*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 8*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2)*cos(q3) + 8*qd1^2*qd2*qd3*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 16*qd1*qd2^2*qd3*qdd2*sin(2*q2 + 2*q3)*sin(q3) + 16*g^2*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) + 16*g^2*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) + 4*qd1^2*qd2^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 4*qd1^2*qd2^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) + 4*qd1^2*qd3^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 8*qd1^2*qd2*qd3*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 16*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) + 16*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) + 16*g*qd1*qd2^2*qd3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3))/(2*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))),     -(qdd1*(cos(q2)*sin(q2)*qd1^2 + qdd2 + g*cos(q2)))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))),   -(qdd2 + g*cos(q2))/(g*cos(q2)),                                                                                                                                                                                                            -(6*qdd1*qdd2^2 - qd1^4*qdd1*cos(4*q2 + 3*q3) + 6*qdd1*qdd2^2*cos(2*q2 + 2*q3) + qd1^4*qdd1*cos(q3) + 8*qdd1*qdd2^2*cos(q3) + 6*qdd1*qdd2*qdd3 + 8*g^2*qdd1*cos(q2 + q3)*cos(q2) + 6*qdd1*qdd2*qdd3*cos(2*q2 + 2*q3) + 4*qd1^2*qdd1*qdd2*sin(2*q2 + q3) + 8*g*qdd1*qdd2*cos(q2 + q3) + 3*qd1^2*qdd1*qdd2*sin(2*q2) + 3*qd1^2*qdd1*qdd3*sin(2*q2) - 4*qd1^3*qd2*qdd2*cos(4*q2 + 3*q3) + 6*g*qdd1*qdd2*cos(q2) + 6*g*qdd1*qdd3*cos(q2) - 12*qd1*qd2*qdd2^2*sin(2*q2 + 2*q3) - 12*qd1*qd3*qdd2^2*sin(2*q2 + 2*q3) - 2*qd1^2*qdd1*qdd2*sin(4*q2 + 3*q3) + 2*qd1^4*qdd1*sin(2*q2)*sin(q3) + 8*qdd1*qdd2^2*cos(2*q2 + 2*q3)*cos(q3) + 4*qd1^3*qd2*qdd2*cos(q3) - 2*qd1^2*qdd1*qdd2*sin(q3) + 8*qd2^2*qdd1*qdd2*sin(q3) + 2*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2) - 2*qd1^4*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 4*qd1^2*qdd1*qdd2*sin(2*q2)*cos(q3) - 2*qd1^4*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 2*qd1^4*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) - 4*qd1^4*qdd1*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 8*g*qdd1*qdd2*cos(q2)*cos(q3) - 16*qd1*qd2*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) - 16*qd1*qd3*qdd2^2*sin(2*q2 + 2*q3)*cos(q3) + 4*qd1^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) - 4*qd1^2*qdd1*qdd2*sin(2*q2 + 2*q3)*cos(q3) + 4*qd1^2*qdd1*qdd3*sin(2*q2 + 2*q3)*cos(q3) + 8*qd2^2*qdd1*qdd2*cos(2*q2 + 2*q3)*sin(q3) - 4*qd1^4*qdd1*cos(2*q2 + q3)*sin(2*q2)*sin(2*q2 + 2*q3) + 2*qd1^4*qdd1*sin(2*q2 + q3)*cos(2*q2)*sin(2*q2 + 2*q3) + 2*qd1^4*qdd1*sin(2*q2 + q3)*sin(2*q2)*cos(2*q2 + 2*q3) - 16*qd1*qd2^3*qdd2*sin(2*q2 + 2*q3)*sin(q3) - 8*qd1^3*qd2*qdd2*sin(2*q2 + 2*q3)*sin(q3) + 8*g^2*qdd1*cos(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 12*qd1*qd2*qdd2*qdd3*sin(2*q2 + 2*q3) - 12*qd1*qd3*qdd2*qdd3*sin(2*q2 + 2*q3) + 4*qd1^2*qd2^2*qdd1*sin(2*q2)*sin(q3) + 8*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2 + q3) + 3*qd1^2*qdd1*qdd2*sin(2*q2)*cos(2*q2 + 2*q3) + 3*qd1^2*qdd1*qdd3*sin(2*q2)*cos(2*q2 + 2*q3) + 4*g*qd1^2*qdd1*cos(q2)*sin(q3) + 8*g*qd2^2*qdd1*cos(q2)*sin(q3) - 6*qd1^3*qd2*qdd2*sin(2*q2)*sin(2*q2 + 2*q3) - 6*qd1^3*qd2*qdd3*sin(2*q2)*sin(2*q2 + 2*q3) - 6*qd1^3*qd3*qdd2*sin(2*q2)*sin(2*q2 + 2*q3) - 6*qd1^3*qd3*qdd3*sin(2*q2)*sin(2*q2 + 2*q3) - 4*qd1^2*qd2^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 4*qd1^2*qd3^2*qdd1*sin(2*q2 + 2*q3)*sin(q3) + 6*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2) + 6*g*qdd1*qdd3*cos(2*q2 + 2*q3)*cos(q2) + 4*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(q2) + 4*g*qd1^2*qdd1*sin(2*q2)*cos(q2 + q3) + 8*qd1^3*qd2*qd3^2*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 8*qd1^3*qd2^2*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 4*g*qd1^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) - 8*g*qd1^2*qdd1*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) + 8*g*qd2^2*qdd1*cos(2*q2 + 2*q3)*cos(q2)*sin(q3) - 16*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) - 16*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2 + q3) - 16*g*qd1*qd2^3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) - 8*g*qd1^3*qd2*sin(2*q2 + 2*q3)*cos(q2)*sin(q3) - 8*g*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) + 4*g*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2 + 2*q3)*cos(q2) + 4*g*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q2 + q3) - 12*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2) - 12*g*qd1*qd2*qdd3*sin(2*q2 + 2*q3)*cos(q2) - 12*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2) - 12*g*qd1*qd3*qdd3*sin(2*q2 + 2*q3)*cos(q2) + 8*g*qd1^3*qd2*sin(2*q2 + q3)*sin(2*q2 + 2*q3)*cos(q2) - 8*g*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) - 8*g*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q2 + q3) + 4*qd1^2*qdd1*qdd2*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 4*qd1^2*qdd1*qdd2*sin(2*q2)*cos(2*q2 + 2*q3)*cos(q3) + 4*qd1^2*qdd1*qdd3*cos(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 16*qd1^3*qd2*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 8*qd1^3*qd2*qdd3*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) - 8*qd1^3*qd3*qdd2*sin(2*q2)*sin(2*q2 + 2*q3)*cos(q3) + 8*g*qdd1*qdd2*cos(2*q2 + 2*q3)*cos(q2)*cos(q3) - 8*qd1^2*qd2*qd3*qdd1*sin(2*q2 + 2*q3)*sin(q3) - 16*qd1*qd2^2*qd3*qdd2*sin(2*q2 + 2*q3)*sin(q3) - 16*g^2*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 16*g^2*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2 + q3)*cos(q2) - 4*qd1^2*qd2^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) + 4*qd1^2*qd2^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3)*sin(q3) - 4*qd1^2*qd3^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 8*qd1^2*qd2*qd3*qdd1*cos(2*q2)*sin(2*q2 + 2*q3)*sin(q3) - 16*g*qd1*qd2*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) - 16*g*qd1*qd3*qdd2*sin(2*q2 + 2*q3)*cos(q2)*cos(q3) - 16*g*qd1*qd2^2*qd3*sin(2*q2 + 2*q3)*cos(q2)*sin(q3))/(2*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))),     qdd2/(g*cos(q2)),                                                                                             -((qdd2 + qdd3)*(8*qdd1*qdd3*cos(q3) - 6*g*qdd1*cos(q2) - 8*qdd1*qdd2*cos(q3) - 6*qdd1*qdd2 - 16*qdd1*qdd2*cos(2*q2 + q3) - 4*qd1^2*qdd1*sin(q3) - 8*qd2^2*qdd1*sin(q3) - 8*qd3^2*qdd1*sin(q3) - 6*qdd1*qdd2*cos(2*q2 + 2*q3) + 4*qd1^2*qdd1*sin(2*q2 + q3) - 3*qd1^2*qdd1*sin(2*q2) - 16*g*qdd1*cos(2*q2 + q3)*cos(q2) + 12*qd1*qd2*qdd2*sin(2*q2 + 2*q3) + 12*qd1*qd3*qdd2*sin(2*q2 + 2*q3) - 3*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3) + 8*qdd1*qdd2*cos(2*q2)*cos(q3) + 8*qdd1*qdd3*cos(2*q2)*cos(q3) + 6*qd1^3*qd2*sin(2*q2)*sin(2*q2 + 2*q3) + 6*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3) - 6*g*qdd1*cos(2*q2 + 2*q3)*cos(q2) + 16*qd1*qd3*qdd2*sin(q3) - 16*qd2*qd3*qdd1*sin(q3) - 4*qd1^2*qdd1*cos(2*q2)*sin(q3) - 8*qd1^2*qdd1*sin(2*q2)*cos(q3) - 8*qd2^2*qdd1*cos(2*q2)*sin(q3) - 8*qd3^2*qdd1*cos(2*q2)*sin(q3) + 32*qd1*qd2*qdd2*sin(2*q2 + q3) + 16*qd1*qd3*qdd2*sin(2*q2 + q3) + 16*qd1*qd2^3*sin(2*q2)*sin(q3) + 8*qd1^3*qd2*sin(2*q2)*sin(q3) + 8*qd1^3*qd3*sin(2*q2)*sin(q3) - 16*g*qdd1*cos(q2)*cos(q3) - 8*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2) + 4*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2) + 8*qd1^3*qd2*sin(2*q2 + q3)*sin(2*q2) + 8*qd1^3*qd3*sin(2*q2 + q3)*sin(2*q2) + 16*qd1*qd2*qd3^2*sin(2*q2)*sin(q3) + 32*qd1*qd2^2*qd3*sin(2*q2)*sin(q3) + 16*g*qd1*qd3*cos(q2)*sin(q3) + 32*g*qd1*qd2*sin(2*q2 + q3)*cos(q2) + 16*g*qd1*qd3*sin(2*q2 + q3)*cos(q2) - 16*qd1*qd2*qdd2*sin(2*q2)*cos(q3) - 16*qd1*qd2*qdd3*sin(2*q2)*cos(q3) - 16*qd2*qd3*qdd1*cos(2*q2)*sin(q3) + 12*g*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2) + 12*g*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2)))/(2*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))),     -(2500*qdd1*(cos(q2)*sin(q2)*qd1^2 + qdd2 + g*cos(q2)))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))),  (900*qdd2)/(g*cos(q2)),                                                                                     -(225*qdd3*(6*qdd1*qdd3 - 6*g*qdd1*cos(q2) + 3*qd1^2*qdd1*sin(2*q2 + 2*q3) + 8*qdd1*qdd3*cos(q3) - 16*qdd1*qdd2*cos(2*q2 + q3) + 6*qdd1*qdd2*cos(2*q2) + 6*qdd1*qdd3*cos(2*q2) - 8*qd3^2*qdd1*sin(q3) - 6*qdd1*qdd2*cos(2*q2 + 2*q3) + 8*qd1^2*qdd1*sin(2*q2 + q3) + 8*g*qdd1*cos(q2 + q3) - 3*qd1^2*qdd1*sin(2*q2) - 16*g*qdd1*cos(2*q2 + q3)*cos(q2) + 8*g*qdd1*cos(2*q2)*cos(q2 + q3) + 12*qd1*qd2*qdd2*sin(2*q2 + 2*q3) + 12*qd1*qd3*qdd2*sin(2*q2 + 2*q3) + 3*qd1^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3) - 3*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3) + 16*qdd1*qdd2*cos(2*q2)*cos(q3) + 8*qdd1*qdd3*cos(2*q2)*cos(q3) + 6*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3) - 6*g*qdd1*cos(2*q2 + 2*q3)*cos(q2) + 16*qd1*qd3*qdd2*sin(q3) - 16*qd2*qd3*qdd1*sin(q3) - 8*qd1^2*qdd1*sin(2*q2)*cos(q3) - 8*qd3^2*qdd1*cos(2*q2)*sin(q3) + 32*qd1*qd2*qdd2*sin(2*q2 + q3) + 16*qd1*qd3*qdd2*sin(2*q2 + q3) + 8*qd1^3*qd3*sin(2*q2)*sin(q3) - 16*g*qdd1*cos(q2)*cos(q3) - 8*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2) + 8*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2) - 12*qd1*qd2*qdd2*sin(2*q2) - 12*qd1*qd2*qdd3*sin(2*q2) + 8*qd1^3*qd3*sin(2*q2 + q3)*sin(2*q2) + 16*qd1*qd2*qd3^2*sin(2*q2)*sin(q3) + 32*qd1*qd2^2*qd3*sin(2*q2)*sin(q3) + 16*g*qd1*qd3*cos(q2)*sin(q3) + 32*g*qd1*qd2*sin(2*q2 + q3)*cos(q2) - 16*g*qd1*qd2*sin(2*q2)*cos(q2 + q3) + 16*g*qd1*qd3*sin(2*q2 + q3)*cos(q2) - 32*qd1*qd2*qdd2*sin(2*q2)*cos(q3) - 16*qd1*qd2*qdd3*sin(2*q2)*cos(q3) - 16*qd2*qd3*qdd1*cos(2*q2)*sin(q3) + 12*g*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2) + 12*g*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2)))/(2*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3))),     -(2500*qd1*(cos(q2)*sin(q2)*qd1^2 + qdd2 + g*cos(q2)))/(g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))),  (900*qd2)/(g*cos(q2)),                                                                                     -(225*qd3*(6*qdd1*qdd3 - 6*g*qdd1*cos(q2) + 3*qd1^2*qdd1*sin(2*q2 + 2*q3) + 8*qdd1*qdd3*cos(q3) - 16*qdd1*qdd2*cos(2*q2 + q3) + 6*qdd1*qdd2*cos(2*q2) + 6*qdd1*qdd3*cos(2*q2) - 8*qd3^2*qdd1*sin(q3) - 6*qdd1*qdd2*cos(2*q2 + 2*q3) + 8*qd1^2*qdd1*sin(2*q2 + q3) + 8*g*qdd1*cos(q2 + q3) - 3*qd1^2*qdd1*sin(2*q2) - 16*g*qdd1*cos(2*q2 + q3)*cos(q2) + 8*g*qdd1*cos(2*q2)*cos(q2 + q3) + 12*qd1*qd2*qdd2*sin(2*q2 + 2*q3) + 12*qd1*qd3*qdd2*sin(2*q2 + 2*q3) + 3*qd1^2*qdd1*cos(2*q2)*sin(2*q2 + 2*q3) - 3*qd1^2*qdd1*sin(2*q2)*cos(2*q2 + 2*q3) + 16*qdd1*qdd2*cos(2*q2)*cos(q3) + 8*qdd1*qdd3*cos(2*q2)*cos(q3) + 6*qd1^3*qd3*sin(2*q2)*sin(2*q2 + 2*q3) - 6*g*qdd1*cos(2*q2 + 2*q3)*cos(q2) + 16*qd1*qd3*qdd2*sin(q3) - 16*qd2*qd3*qdd1*sin(q3) - 8*qd1^2*qdd1*sin(2*q2)*cos(q3) - 8*qd3^2*qdd1*cos(2*q2)*sin(q3) + 32*qd1*qd2*qdd2*sin(2*q2 + q3) + 16*qd1*qd3*qdd2*sin(2*q2 + q3) + 8*qd1^3*qd3*sin(2*q2)*sin(q3) - 16*g*qdd1*cos(q2)*cos(q3) - 8*qd1^2*qdd1*cos(2*q2 + q3)*sin(2*q2) + 8*qd1^2*qdd1*sin(2*q2 + q3)*cos(2*q2) - 12*qd1*qd2*qdd2*sin(2*q2) - 12*qd1*qd2*qdd3*sin(2*q2) + 8*qd1^3*qd3*sin(2*q2 + q3)*sin(2*q2) + 16*qd1*qd2*qd3^2*sin(2*q2)*sin(q3) + 32*qd1*qd2^2*qd3*sin(2*q2)*sin(q3) + 16*g*qd1*qd3*cos(q2)*sin(q3) + 32*g*qd1*qd2*sin(2*q2 + q3)*cos(q2) - 16*g*qd1*qd2*sin(2*q2)*cos(q2 + q3) + 16*g*qd1*qd3*sin(2*q2 + q3)*cos(q2) - 32*qd1*qd2*qdd2*sin(2*q2)*cos(q3) - 16*qd1*qd2*qdd3*sin(2*q2)*cos(q3) - 16*qd2*qd3*qdd1*cos(2*q2)*sin(q3) + 12*g*qd1*qd2*sin(2*q2 + 2*q3)*cos(q2) + 12*g*qd1*qd3*sin(2*q2 + 2*q3)*cos(q2)))/(2*g*cos(q2)^2*(qdd1*cos(q2) - 2*qd1*qd2*sin(q2))*(6*qdd2 + 6*qdd3 + 4*qd1^2*sin(q3) + 8*qd2^2*sin(q3) + 4*qd1^2*sin(2*q2 + q3) + 8*g*cos(q2 + q3) + 3*qd1^2*sin(2*q2 + 2*q3) + 8*qdd2*cos(q3)))]
%  
