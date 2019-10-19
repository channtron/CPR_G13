% Se dibujan 7 graficas:
% -Referencias de artuculacion vs articulacion
% -Errores articulacion(referencia - real)
% -Referencias de velocidad vs velocidad
% -Errores velocidad(referencia - real)
% -Fuerzas y pares
% -Movimiento en XYZ
% -Error en XYZ

%q vs q ref
figure();subplot(3,1,1);plot(t,q(:,1),t,qr(:,1));grid;xlabel('Tiempo(s)');ylabel('q1 vs q1 ref(rad)');legend('q1','q1ref','location','southwest');...
    subplot(3,1,2);plot(t,q(:,2),t,qr(:,2));grid;xlabel('Tiempo(s)');ylabel('q2 vs q2 ref(m)');legend('q2','q2ref','location','northwest');...
    subplot(3,1,3);plot(t,q(:,3),t,qr(:,3));grid;xlabel('Tiempo(s)');ylabel('q3 vs q3 ref(m)');legend('q3','q3ref','location','southwest');

%qd vs qd ref
figure();subplot(3,1,1);plot(t,qd(:,1),t,qdr(:,1));grid;xlabel('Tiempo(s)');ylabel('qd1 vs qd1 ref(rad/s)');legend('qd1','qd1ref','location','southwest');...
    subplot(3,1,2);plot(t,qd(:,2),t,qdr(:,2));grid;xlabel('Tiempo(s)');ylabel('qd2 vs qd2 ref(m/s)');legend('qd2','qd2ref','location','northwest');...
    subplot(3,1,3);plot(t,qd(:,3),t,qdr(:,3));grid;xlabel('Tiempo(s)');ylabel('qd3 vs qd3 ref(m/s)');legend('qd3','qd3ref','location','northwest');

%qdd vs qdd ref
figure();subplot(3,1,1);plot(t,qdd(:,1),t,qddr(:,1));grid;xlabel('Tiempo(s)');ylabel('qdd1 vs qdd1 ref(rad/s)');legend('qdd1','qdd1ref','location','southwest');...
    subplot(3,1,2);plot(t,qdd(:,2),t,qddr(:,2));grid;xlabel('Tiempo(s)');ylabel('qdd2 vs qdd2 ref(m/s)');legend('qdd2','qdd2ref','location','northwest');...
    subplot(3,1,3);plot(t,qdd(:,3),t,qddr(:,3));grid;xlabel('Tiempo(s)');ylabel('qdd3 vs qdd3 ref(m/s)');legend('qdd3','qdd3ref','location','northwest');

%Fuerzas y pares
figure();subplot(3,1,1);plot(t,Im(:,1));grid;xlabel('Tiempo(s)');ylabel('N1(N*m)');...
    subplot(3,1,2);plot(t,Im(:,2));grid;xlabel('Tiempo(s)');ylabel('F2(N)');...
    subplot(3,1,3);plot(t,Im(:,3));grid;xlabel('Tiempo(s)');ylabel('F3(N)');
