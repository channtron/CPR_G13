% sim('sk_R3GDL.mdl');

% Graficas comparativas entre modelo estimados y resultados reales
figure(3);
subplot(3,1,1); plot(t_m,qdd_s(:,1),t_m,qddrs(:,1)); xlabel('t (s)'); ylabel('qpp_1 (rad/s^2)'); grid; legend('Real', 'Estimado');
subplot(3,1,2); plot(t_m,qdd_s(:,2),t_m,qddrs(:,2)); xlabel('t (s)'); ylabel('qpp_2 (rad/s^2)'); grid; legend('Real', 'Estimado');
subplot(3,1,3); plot(t_m,qdd_s(:,3),t_m,qddrs(:,3)); xlabel('t (s)'); ylabel('qpp_3 (rad/s^2)'); grid; legend('Real', 'Estimado');

figure(2);
subplot(3,1,1); plot(t_m,qd_s(:,1),t_m,qdrs(:,1)); xlabel('t (s)'); ylabel('qp_1 (rad/s)'); grid; legend('Real', 'Estimado');
subplot(3,1,2); plot(t_m,qd_s(:,2),t_m,qdrs(:,2)); xlabel('t (s)'); ylabel('qp_2 (rad/s)'); grid; legend('Real', 'Estimado');
subplot(3,1,3); plot(t_m,qd_s(:,3),t_m,qdrs(:,3)); xlabel('t (s)'); ylabel('qp_3 (rad/s)'); grid; legend('Real', 'Estimado');

figure(1);
subplot(3,1,1); plot(t_m,q_s(:,1),t_m,qrs(:,1)); xlabel('t (s)'); ylabel('q_1 (rad)'); grid; legend('Real', 'Estimado');
subplot(3,1,2); plot(t_m,q_s(:,2),t_m,qrs(:,2)); xlabel('t (s)'); ylabel('q_2 (rad)'); grid; legend('Real', 'Estimado');
subplot(3,1,3); plot(t_m,q_s(:,3),t_m,qrs(:,3)); xlabel('t (s)'); ylabel('q_3 (rad)'); grid; legend('Real', 'Estimado');
% 
% figure(4);
% subplot(3,1,1); plot(t_m,qd_ms(:,1),t_m,qd_s(:,1)); xlabel('t (s)'); ylabel('qp_1 (rad/s)'); grid; legend('Ruido','Ideal');
% subplot(3,1,2); plot(t_m,qd_ms(:,2),t_m,qd_s(:,2)); xlabel('t (s)'); ylabel('qp_2 (rad/s)'); grid; legend('Ruido','Ideal');
% subplot(3,1,3); plot(t_m,qd_ms(:,3),t_m,qd_s(:,3)); xlabel('t (s)'); ylabel('qp_3 (rad/s)'); grid; legend('Ruido','Ideal');