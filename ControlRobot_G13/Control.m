function [I] = Control(in)
%Esta funcion implementa un controlador, tanto continuo como discreto
%   Detailed explanation goes here
qr1        = in(1);
qr2        = in(2);
qr3        = in(3);
qdr1       = in(4);
qdr2       = in(5);
qdr3       = in(6);
qddr1      = in(7);
qddr2      = in(8);
qddr3      = in(9);
q1         = in(10);
q2         = in(11);
q3         = in(12);
qd1        = in(13);
qd2        = in(14);
qd3        = in(15);
tiempo     = in(16);

qr=[qr1;qr2;qr3];
qdr=[qdr1;qdr2;qdr3];
qddr=[qddr1;qddr2;qddr3];
q=[q1;q2;q3];
qd=[qd1;qd2;qd3];

%Definimos Variables Estáticas y tiempo de muestreo
persistent Int_err;
persistent numf;
persistent denf;
Tm=0.001;
IMax=1000;
IMin=-1000;

%Nos aseguramos de que el error se inicializa cada vez que simulo
if(tiempo<1e-5)
    Int_err=[0;0;0];
%     [numf,denf]=butter(2,(2*pi*60)/(pi/Tm),'low');
end

% qd=filter(numf,denf,qd_aux);



%CONTROLADORES
%------------------------------------------------------------------------
R1=1; R2=1; R3=1; 
g=9.81;
%Definimos las constantes Kp Ki Kd
%PD Sin Precompensar (hecho pero sale meh)
% Kp=1.0e3*[1.2986;1.627;0.18684];
% Ti=[0;0;0];
% Td=[ 0.0667; 0.0667; 0.0667];

%PD Precompensando V y G (igual que sin precompensar)
% Kp=1.0e3*[1.298;1.627;0.18684];
% Ti=[0;0;0];
% Td=[0.0667;0.0667;0.0667];

%PD por par calculado (este es el unico PD que funciona bien de verdad)
% Kp=[1350;1687.5;1928.6];
% Ti=[0;0;0];
% Td=[0.1;0.1;0.1];

%PID Sin Precompensar (hecho y sale guay)
% Kp=[29218;36608;4203.8];
% Ti=[0.2;0.2;0.2];
% Td=[0.05;0.05;0.05];

%PID compensando Va y Ga (hecho y sale)
% Kp=[29218;36608;4203.8];
% Ti=[0.2;0.2;0.2];
% Td=[0.05;0.05;0.05];

%PID por par calculado (listo c:)
Kp=[2700;3375;3857.2];
Ti=[0.2;0.2;0.2];
Td=[0.05;0.05;0.05];

%PD por diseño en frecuencia
% Kp=[3.2278e+03;4665.600;1868.400];
% Ti=[0;0;0];
% Td=[0.0667;0.0667;0.0667];
%-----------------------------------------------------------------------

%Calculamos los errores
Err_q=qr-q;
Err_qd=qdr-qd;

%Calculamos el control
Int_err=Int_err+Tm*Err_q;

%PID
U=Kp.*Err_q + Kp.*Td.*Err_qd + Kp./Ti.*Int_err;

%PD
% U=Kp.*Err_q + Kp.*Td.*Err_qd;



%Precompensar
Ma= [ 1.9283*cos(2.0*q2 + q3) + 2.9292*cos(2.0*q2) + 1.9283*cos(q3) + 0.54422*cos(2.0*q2 + 2.0*q3) + 3.4915,  0,  0;
    0, 3.8566*cos(q3) + 6.9901,   1.9283*cos(q3) + 1.0898;
    0, 1.9283*cos(q3) + 1.0898, 0.000093518*R3^2 + 1.0898];
 

Va= [ -8.4703e-22*qdr1*(- 6.4183e15*R1^2 + 4.553e21*qdr2*sin(2.0*q2 + q3) + 2.2765e21*qdr3*sin(2.0*q2 + q3) + 6.9164e21*qdr2*sin(2.0*q2) + 2.2765e21*qdr3*sin(q3) + 1.285e21*qdr2*sin(2.0*q2 + 2.0*q3) + 1.285e21*qdr3*sin(2.0*q2 + 2.0*q3));
    0.54422*qdr1^2*sin(2.0*q2 + 2.0*q3) - 1.9283*qdr3^2*sin(q3) + 8.693e-6*R2^2*qdr2 + 1.9283*qdr1^2*sin(2.0*q2 + q3) + 2.9292*qdr1^2*sin(2.0*q2) - 3.8566*qdr2*qdr3*sin(q3);
	0.96414*qdr1^2*sin(q3) + 1.9283*qdr2^2*sin(q3) + 0.54422*qdr1^2*sin(2.0*q2 + 2.0*q3) + 0.000060046*R3^2*qdr3 + 0.96414*qdr1^2*sin(2.0*q2 + q3)];

Ga = [ 0;
    g*(1.9283*cos(q2 + q3) + 6.742*cos(q2));
	1.9283*g*cos(q2 + q3)];

% I=U;    %Sin precompensar nada
% I=U+Ga; %Precompensando la gravedad
% I=U+Va+Ga; %Precompensando Va y Ga
I=Va+Ga+Ma*(U+qddr); %Precompensando todo

% Saturacion de las intensidades y antiwindup 
if I(1)>IMax
    I(1)=IMax;
    Int_err(1)=Int_err(1)-Tm*Err_q(1);
elseif I(1)<IMin
    I(1)=IMin;
    Int_err(1)=Int_err(1)-Tm*Err_q(1);
end

if I(2)>IMax
    I(2)=IMax;
    Int_err(2)=Int_err(2)-Tm*Err_q(2);
elseif I(2)<IMin
    I(2)=IMin;
    Int_err(2)=Int_err(2)-Tm*Err_q(2);
end

if I(3)>IMax
    I(3)=IMax;
    Int_err(3)=Int_err(3)-Tm*Err_q(3);
elseif I(3)<IMin
    I(3)=IMin;
    Int_err(3)=Int_err(3)-Tm*Err_q(3);
end
end