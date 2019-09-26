function [A] = MTH(den)

theta=den(1);
d=den(2);
a=den(3);
alpha=den(4);

A=[ cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha),  a*cos(theta);
    sin(theta), cos(theta)*cos(alpha),   -cos(theta)*sin(alpha), a*sin(theta);
    0,          sin(alpha),             cos(alpha),             d;
    0,          0,                      0,                      1];

end

