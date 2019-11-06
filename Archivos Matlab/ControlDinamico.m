function [I] = ControlDinamico(in)
%Esta funcion implementa un controlador, tanto continuo como discreto
%   Detailed explanation goes here
qr1        = in(1);
qr2        = in(2);
qr3        = in(3);
qdr1       = in(4);
qdr2       = in(5);
qdr3       = in(6);
q1         = in(7);
q2         = in(8);
q3         = in(9);
qd1        = in(10);
qd2        = in(11);
qd3        = in(12);
qddr1      = in(13);
qddr2      = in(14);
qddr3      = in(15);
tiempo     = in(16);

qr=[qr1;qr2;qr3];
qdr=[qdr1;qdr2;qdr3];
qddr=[qddr1;qddr2;qddr3];
q=[q1;q2;q3];
qd=[qd1;qd2;qd3];

%Definimos Variables Estáticas y tiempo de muestreo
persistent Int_err;
Tm=0.001;

%Nos aseguramos de que el error se inicializa cada vez que simulo
if(tiempo<1e-5)
    Int_err=[0;0;0];
end

%CONTROLADORES
%------------------------------------------------------------------------
R1=1; R2=1; R3=1; 
g=9.81;
%Definimos las constantes Kp Ki Kd
%PI Sin Precompensar
Kp=1.0e3*[1.2986;1.627;0.18684];
Ki=[0;0;0];
Kd=[0.0667;0.0667;0.0667];

%PI Precompensando V y G
% Kp=1.0e6*[3.5481;2.5119;1.2589];
% Ki=[0;0;0];
% Kd=1.0e4*[6.4989;4.6009;2.3059];

%PI por par calculado
% Kp=1.0e3*[7.6736 7.6736 7.6736];
% Ki=[0;0;0];
% Kd=[140.5539 140.5539 140.5539];

%PID Sin Precompensar
% Kp=1e6*[1.2264;0.8484;0.4252];
% Ki=1e6*[8.1091;5.6101;2.8117];
% Kd=1e4*[2.3184;1.6039;0.8039];

%PID compensando Va y Ga
% Kp=1e6*[1.1985;0.8484;0.4252];
% Ki=1e6*[7.9245;5.6101;2.8117];
% Kd=1e4*[2.2656;1.6039;0.8039];

%PID por par calculado
% Kp=1e3*[7.5617;7.5617;7.5617];
% Ki=[50000;50000;50000];
% Kd=[142.9486;142.9486;142.9486];
%-----------------------------------------------------------------------

%Calculamos los errores
Err_q=qr-q;
Err_qd=qdr-qd;

%Calculamos el control
Int_err=Int_err+Tm*Err_q;
U=Kp.*Err_q + Kd.*Err_qd + Ki.*Int_err;

%Precompensar
Ma= [ 1.9283*cos(2.0*q2 + q3) + 2.9292*cos(2.0*q2) + 1.9283*cos(q3) + 0.54422*cos(2.0*q2 + 2.0*q3) + 3.4915,  0,  0;
    0, 3.8566*cos(q3) + 6.9901,   1.9283*cos(q3) + 1.0898;
    0, 1.9283*cos(q3) + 1.0898, 0.000093518*R3^2 + 1.0898];
 

Va= [ -8.4703e-22*qd1*(- 6.4183e15*R1^2 + 4.553e21*qd2*sin(2.0*q2 + q3) + 2.2765e21*qd3*sin(2.0*q2 + q3) + 6.9164e21*qd2*sin(2.0*q2) + 2.2765e21*qd3*sin(q3) + 1.285e21*qd2*sin(2.0*q2 + 2.0*q3) + 1.285e21*qd3*sin(2.0*q2 + 2.0*q3));
    0.54422*qd1^2*sin(2.0*q2 + 2.0*q3) - 1.9283*qd3^2*sin(q3) + 8.693e-6*R2^2*qd2 + 1.9283*qd1^2*sin(2.0*q2 + q3) + 2.9292*qd1^2*sin(2.0*q2) - 3.8566*qd2*qd3*sin(q3);
	0.96414*qd1^2*sin(q3) + 1.9283*qd2^2*sin(q3) + 0.54422*qd1^2*sin(2.0*q2 + 2.0*q3) + 0.000060046*R3^2*qd3 + 0.96414*qd1^2*sin(2.0*q2 + q3)];

Ga = [ 0;
    g*(1.9283*cos(q2 + q3) + 6.742*cos(q2));
	1.9283*g*cos(q2 + q3)];

%I=U;    %Sin precompensar nada
I=U+Ga; %Precompensando la gravedad
%I=U+Va+Ga; %Precompensando Va y Ga
%I=Va+Ga+Ma*(U+qddr); %Precompensando todo

