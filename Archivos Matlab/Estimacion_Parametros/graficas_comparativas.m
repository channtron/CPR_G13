sim('sk_R3GDL.mdl');

% Graficas comparativas entre modelo estimados y resultados reales
figure(3);
subplot(3,1,1); plot(ts,qpps(:,1),ts,qpps_est(:,1)); xlabel('t (s)'); ylabel('qpp_1 (rad/s^2)'); grid; legend('Real', 'Estimado');
subplot(3,1,2); plot(ts,qpps(:,2),ts,qpps_est(:,2)); xlabel('t (s)'); ylabel('qpp_2 (rad/s^2)'); grid; legend('Real', 'Estimado');
subplot(3,1,3); plot(ts,qpps(:,3),ts,qpps_est(:,3)); xlabel('t (s)'); ylabel('qpp_3 (rad/s^2)'); grid; legend('Real', 'Estimado');

figure(2);
subplot(3,1,1); plot(ts,qps(:,1),ts,qps_est(:,1)); xlabel('t (s)'); ylabel('qp_1 (rad/s)'); grid; legend('Real', 'Estimado');
subplot(3,1,2); plot(ts,qps(:,2),ts,qps_est(:,2)); xlabel('t (s)'); ylabel('qp_2 (rad/s)'); grid; legend('Real', 'Estimado');
subplot(3,1,3); plot(ts,qps(:,3),ts,qps_est(:,3)); xlabel('t (s)'); ylabel('qp_3 (rad/s)'); grid; legend('Real', 'Estimado');

figure(1);
subplot(3,1,1); plot(ts,qs(:,1),ts,qs_est(:,1)); xlabel('t (s)'); ylabel('q_1 (rad)'); grid; legend('Real', 'Estimado');
subplot(3,1,2); plot(ts,qs(:,2),ts,qs_est(:,2)); xlabel('t (s)'); ylabel('q_2 (rad)'); grid; legend('Real', 'Estimado');
subplot(3,1,3); plot(ts,qs(:,3),ts,qs_est(:,3)); xlabel('t (s)'); ylabel('q_3 (rad)'); grid; legend('Real', 'Estimado');