%Matrices TH
syms q1 q2 q3
L1=0.65;
L2=1;
L0=0.5;
L3=0.75;
A01=MTH([0 L0 0 0]);
A12=MTH([q1 L1 0 pi/2]);
A23=MTH([q2 0 L2 0]);
A34=MTH([q3 0 L3 0]);
A40=A01*A12*A23*A34;
simplify(A40)
%Cinematica Directa
%x=cos(q1)*(L3*cos(q2 + q3) + L2*cos(q2))
%y=sin(q1)*(L3*cos(q2 + q3) + L2*cos(q2))
%z=L0 + L1 + L3*sin(q2 + q3) + L2*sin(q2)