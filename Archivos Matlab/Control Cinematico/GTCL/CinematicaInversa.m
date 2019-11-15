function q = CinematicaInversa(in)

x       = in(1);           % Posición cartesianas
y       = in(2);           % 
z       = in(3);           % 

% A rellenar por el alumno
L1=0.65;
L2=1;
L0=0.5;
L3=0.75;
q3cos=(x^2 + y^2 + (z-L1-L0)^2 - (L2)^2 - (L3)^2)/(2*L2*L3);

q1=atan2(y,x);
q3=atan2(sqrt(1-(q3cos)^2),q3cos);
q2=atan2(z-L1-L0,sqrt(x^2 + y^2)) - atan2((L3*sin(q3)),(L2+L3*cos(q3)));

q=[q1;q2;q3];
