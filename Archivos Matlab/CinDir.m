function [p] = CinDir(q)
L1=0.65;
L2=1;
L0=0.5;
L3=0.75;

q1=q(1);
q2=q(2);
q3=q(3);

%Modelo cinematico directo
px=cos(q1)*(L3*cos(q2 + q3) + L2*cos(q2));
py=sin(q1)*(L3*cos(q2 + q3) + L2*cos(q2));
pz=L0 + L1 + L3*sin(q2 + q3) + L2*sin(q2);

p=[px;py;pz]
end

