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

%Cinematica Inversa
A01^(-1)
simplify(A12*A23*A34)

simplify(A12^(-1))*(A01^(-1))
simplify(A23*A34)

simplify(A23^(-1)*(A12^(-1))*(A01^(-1)))
A34

syms x y z

%Radio=1 centro=(0,0)
%altura z=1.3
for angulo=0:0.1:2*pi
 x=1*cos(angulo);
 y=1*sin(angulo);
 z=1.3;
 [q1 q2 q3]=CinematicaInversa(x,y,z)
 
end


eq=[q1==atan(y/x), x*cos(q1)+y*sin(q1)==3/4*cos(q3+q2) +cos(q2), -x*cos(q1)*sin(q2)-y*sin(q1)-sin(q2) + z*cos(q2) - (23/20)*cos(q2)==3/4 *sin(q3)]
S=solve(eq,[q1 q2 q3])
