clear all;
syms T1 T2 T3 q1 qd1 qdd1 q2 qd2 qdd2 q3 qd3 qdd3 R1 R2 R3 real  
PI = sym('pi');
g= 9.81;

% theta_err=[1, -0.008449050110542,    -5.9862;
%            2,  0.000034811978541,     2.1984;
%            3,  0.081299020277041,     4.0684;
%            4,  0.040884638254871,     3.8348;
%            5,  0.007078196497825,     7.7617;
%            6,  0.000032714576887,    15.7319;
%            7,  0.006693026242915,    15.5411;
%            8,  0.004295166132318,    14.4064;
%            9,  0.001665474017116,    20.1662;
%           10,  0.000010549988836,    13.8915;
%           11,  0.000116833679087,     3.2593];
      
theta_err=[0.025143665463728
   0.000005436552851
   4.378912858948338
   0.841650362316433
   0.041912895110567
   0.000008692987532
   3.204748093330908
   0.475275653184315
   0.001313979835617
   0.000093518420424
   0.000060046467256];
      
% theta_est=theta_err(:,2);

R1=1; R2=1; R3=1; %Accionamiento directo
R=diag([R1,R2,R3]);
% K
  K=diag([0.5,0.4,0.35]); %(N*m/A)
  g=9.81;
  
Tau=[R1*K(1,1);R2*K(2,2);R3*K(3,3)];

Ma= [ 1.9283*cos(2.0*q2 + q3) + 2.9292*cos(2.0*q2) + 1.9283*cos(q3) + 0.54422*cos(2.0*q2 + 2.0*q3) + 3.4915,  0,  0;
    0, 3.8566*cos(q3) + 6.9901,   1.9283*cos(q3) + 1.0898;
    0, 1.9283*cos(q3) + 1.0898, 0.000093518*R3^2 + 1.0898];
 

Va= [ -8.4703e-22*qd1*(- 6.4183e15*R1^2 + 4.553e21*qd2*sin(2.0*q2 + q3) + 2.2765e21*qd3*sin(2.0*q2 + q3) + 6.9164e21*qd2*sin(2.0*q2) + 2.2765e21*qd3*sin(q3) + 1.285e21*qd2*sin(2.0*q2 + 2.0*q3) + 1.285e21*qd3*sin(2.0*q2 + 2.0*q3));
    0.54422*qd1^2*sin(2.0*q2 + 2.0*q3) - 1.9283*qd3^2*sin(q3) + 8.693e-6*R2^2*qd2 + 1.9283*qd1^2*sin(2.0*q2 + q3) + 2.9292*qd1^2*sin(2.0*q2) - 3.8566*qd2*qd3*sin(q3);
	0.96414*qd1^2*sin(q3) + 1.9283*qd2^2*sin(q3) + 0.54422*qd1^2*sin(2.0*q2 + 2.0*q3) + 0.000060046*R3^2*qd3 + 0.96414*qd1^2*sin(2.0*q2 + q3)];

Ga = [ 0;
    g*(1.9283*cos(q2 + q3) + 6.742*cos(q2));
	1.9283*g*cos(q2 + q3)];
 
 %linealizando la M y la V
%  Tau=[Im1*R1*K(1,1);Im2*R2*K(2,2);Im3*R3*K(3,3)];
  G1s=tf(Tau(1),[10.8215  5.4365e-6 0]);
  G2s=tf(Tau(2),[10.8467  8.6930e-6 0]);
  G3s=tf(Tau(3),[1.0899  6.0046e-5 0]);


  G1s=tf(Tau(1)/5.4365e-6,[10.8215 5.4365e-6 0]/5.4365e-6);
  G2s=tf(Tau(2)/8.6930e-6,[10.8467  8.6930e-6 0]/8.6930e-6);
  G3s=tf(Tau(3)/6.0046e-5,[1.0899  6.0046e-5 0]/6.0046e-5);
  
  %%
  %MODELO LINEALIZANDO M Y V
  %PD ajustado con rltool para que tsbc de 10ms polos de GBc en 30rad/s
  Tauc=1/15;
  K1=19479;
  K2=24405;
  K3=2802.6;
  C1s=tf(K1*[Tauc 1],1);
  C2s=tf(K2*[Tauc 1],1);
  C3s=tf(K3*[Tauc 1],1);
%   Tauc=1/25;
%   K1=54107;
%   K2=67792;
%   K3=7785;
%   C1s=tf(K1*[Tauc 1],1);
%   C2s=tf(K2*[Tauc 1],1);
%   C3s=tf(K3*[Tauc 1],1);
  
  Td1=Tauc;  
  Td2=Tauc;
  Td3=Tauc;
  Kp1=K1*Tauc;
  Kp2=K2*Tauc;
  Kp3=K3*Tauc;
  
  %%
  %MODELO LINEALIZANDO M Y V
  %PID ajustado con rltool para que tsbc de 10ms polos de GBc en 30rad/s
  Tauc=1/10;
  K1=146090;
  K2=183040;
  K3=21019;
  C1s=tf(K1*conv([Tauc 1],[Tauc 1]),[1 0]);
  C2s=tf(K2*conv([Tauc 1],[Tauc 1]),[1 0]);
  C3s=tf(K3*conv([Tauc 1],[Tauc 1]),[1 0]);
  
  Ti1=2*Tauc;
  Ti2=2*Tauc;
  Ti3=2*Tauc;
  Td1=(Tauc^2)/Ti1;  
  Td2=(Tauc^2)/Ti2;
  Td3=(Tauc^2)/Ti3;
  Kp1=K1*Ti1;
  Kp2=K2*Ti2;
  Kp3=K3*Ti3;
  
  %%
  %PI linealizando solo la M
  G1s=tf(Tau(1),[10.8215 0 0]);
  G2s=tf(Tau(2),[10.8467 0 0]);
  G3s=tf(Tau(3),[1.0899  0 0]);
  
  G1s=tf(Tau(1)/10.8215,[10.8215 0 0]/10.8215);
  G2s=tf(Tau(2)/10.8467,[10.8467 0 0]/10.8467);
  G3s=tf(Tau(3)/1.0899,[1.0899  0 0]/1.0899);
  
  Tauc=1/15;
  K1=19479;
  K2=24405;
  K3=2802.6;
  C1s=tf(K1*[Tauc 1],1);
  C2s=tf(K2*[Tauc 1],1);
  C3s=tf(K3*[Tauc 1],1);
  
  Td1=Tauc;  
  Td2=Tauc;
  Td3=Tauc;
  Kp1=K1*Tauc;
  Kp2=K2*Tauc;
  Kp3=K3*Tauc;

  
  %%
  % PID linealizando la M 
% Tau=[Im1*R1*K(1,1);Im2*R2*K(2,2);Im3*R3*K(3,3)];

  G1s=tf(Tau(1)/10.8215,[10.8215 0 0]/10.8215);
  G2s=tf(Tau(2)/10.8467,[10.8467 0 0]/10.8467);
  G3s=tf(Tau(3)/1.0899,[1.0899  0 0]/1.0899);
  
  Tauc=1/10;
  K1=146090;
  K2=183040;
  K3=21019;
  C1s=tf(K1*conv([Tauc 1],[Tauc 1]),[1 0]);
  C2s=tf(K2*conv([Tauc 1],[Tauc 1]),[1 0]);
  C3s=tf(K3*conv([Tauc 1],[Tauc 1]),[1 0]);
  
  Ti1=2*Tauc;
  Ti2=2*Tauc;
  Ti3=2*Tauc;
  Td1=(Tauc^2)/Ti1;  
  Td2=(Tauc^2)/Ti2;
  Td3=(Tauc^2)/Ti3;
  Kp1=K1*Ti1;
  Kp2=K2*Ti2;
  Kp3=K3*Ti3;
  
  %%
  %PI PAR CALCULADO (LINEALIZAMOS TODO)
  G1s=tf(Tau(1),[1 0 0]);
  G2s=tf(Tau(2),[1 0 0]);
  G3s=tf(Tau(3),[1 0 0]);
  
  Tauc=1/10;
  K1=13500;
  K2=16875;
  K3=19286;
  C1s=tf(K1*[Tauc 1],1);
  C2s=tf(K2*[Tauc 1],1);
  C3s=tf(K3*[Tauc 1],1);
  
  Td1=Tauc;  
  Td2=Tauc;
  Td3=Tauc;
  Kp1=K1*Tauc;
  Kp2=K2*Tauc;
  Kp3=K3*Tauc;
  
  %%
  %PID PAR CALCULADO (LINEALIZAMOS TODO)
  
  G1s=tf(Tau(1),[1 0 0]);
  G2s=tf(Tau(2),[1 0 0]);
  G3s=tf(Tau(3),[1 0 0]);

  
  Tauc=1/10;
  K1=13500;
  K2=16875;
  K3=19286;
  C1s=tf(K1*conv([Tauc 1],[Tauc 1]),[1 0]);
  C2s=tf(K2*conv([Tauc 1],[Tauc 1]),[1 0]);
  C3s=tf(K3*conv([Tauc 1],[Tauc 1]),[1 0]);
  
  Ti1=2*Tauc;
  Ti2=2*Tauc;
  Ti3=2*Tauc;
  Td1=(Tauc^2)/Ti1;  
  Td2=(Tauc^2)/Ti2;
  Td3=(Tauc^2)/Ti3;
  Kp1=K1*Ti1;
  Kp2=K2*Ti2;
  Kp3=K3*Ti3;