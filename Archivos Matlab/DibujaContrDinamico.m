% Se dibujan 7 graficas:
% -Referencias de artuculacion vs articulacion
% -Errores articulacion(referencia - real)
% -Referencias de velocidad vs velocidad
% -Errores velocidad(referencia - real)
% -Fuerzas y pares
% -Movimiento en XYZ
% -Error en XYZ

%q vs q ref
figure();subplot(3,1,1);plot(ts,q(:,1),ts,qr(:,1));grid;xlabel('Tiempo(s)');ylabel('q1 vs q1 ref(rad)');legend('q1','q1ref','location','southwest');...
    subplot(3,1,2);plot(ts,q(:,2),ts,qr(:,2));grid;xlabel('Tiempo(s)');ylabel('q2 vs q2 ref(m)');legend('q2','q2ref','location','northwest');...
    subplot(3,1,3);plot(ts,q(:,3),ts,qr(:,3));grid;xlabel('Tiempo(s)');ylabel('q3 vs q3 ref(m)');legend('q3','q3ref','location','southwest');

%Error q (qref-q)
figure();subplot(3,1,1);plot(ts,(qr(:,1)-q(:,1)));grid;xlabel('Tiempo(s)');ylabel('Err q1(rad)');...
    subplot(3,1,2);plot(ts,(qr(:,2)-q(:,2)));grid;xlabel('Tiempo(s)');ylabel('Err q2(m)');...
    subplot(3,1,3);plot(ts,(qr(:,3)-q(:,3)));grid;xlabel('Tiempo(s)');ylabel('Err q3(m)');

%qd vs qd ref
figure();subplot(3,1,1);plot(ts,qd(:,1),ts,qdr(:,1));grid;xlabel('Tiempo(s)');ylabel('qd1 vs qd1 ref(rad/s)');legend('qd1','qd1ref','location','southwest');...
    subplot(3,1,2);plot(ts,qd(:,2),ts,qdr(:,2));grid;xlabel('Tiempo(s)');ylabel('qd2 vs qd2 ref(m/s)');legend('qd2','qd2ref','location','northwest');...
    subplot(3,1,3);plot(ts,qd(:,3),ts,qdr(:,3));grid;xlabel('Tiempo(s)');ylabel('qd3 vs qd3 ref(m/s)');legend('qd3','qd3ref','location','northwest');

%Err qd (qd ref-q
figure();subplot(3,1,1);plot(ts,(qdr(:,1)-qd(:,1)));grid;xlabel('Tiempo(s)');ylabel('Err qd1(rad/s)');...
    subplot(3,1,2);plot(ts,(qdr(:,2)-qd(:,2)));grid;xlabel('Tiempo(s)');ylabel('Err qd2(rad/s)');...
    subplot(3,1,3);plot(ts,(qdr(:,3)-qd(:,3)));grid;xlabel('Tiempo(s)');ylabel('Err qd3(rad/s)');

%Intensidades
figure();subplot(3,1,1);plot(ts,Ims(:,1));grid;xlabel('Tiempo(s)');ylabel('N1(N*m)');...
    subplot(3,1,2);plot(ts,Ims(:,2));grid;xlabel('Tiempo(s)');ylabel('F2(N)');...
    subplot(3,1,3);plot(ts,Ims(:,3));grid;xlabel('Tiempo(s)');ylabel('F3(N)');

%Movimiento XYZ
figure();plot3(x,y,z,xr,yr,zr,'.-');grid;xlabel('EJE X');ylabel('EJE Y');zlabel('EJE Z');

%Error valor absoluto XYZ
DifMod=sqrt((xr-x).^2 + (yr-y).^2 + (zr-z).^2);
figure();plot(ts,DifMod);grid;title('Módulo del error');xlabel('Tiempo(s)');ylabel('Error(m)');