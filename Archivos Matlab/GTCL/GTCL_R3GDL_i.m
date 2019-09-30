    function trayectoria = GTCL_R3GDL(in)
                       
Xini       = in(1);           % Posición cartesiana inicial
Yini       = in(2);           % 
Zini       = in(3);           % 
Xfin       = in(4);           % Posición cartesiana final
Yfin       = in(5);           % 
Zfin       = in(6);           % 
N          = in(7);           % Número de puntos intermedios
t_ini      = in(8);           % Tiempo de inicio del movimiento
T          = in(9);           % Duración del movimiento
t          = in(10);          % Reloj

inc_t=(T/(N+1));
inc_x=(Xfin-Xini)/(N+1);
inc_y=(Yfin-Yini)/(N+1);
inc_z=(Zfin-Zini)/(N+1);

puntos=[Xini, Yini, Zini, t_ini];
qlist=[CinematicaInversa([Xini, Yini, Zini])'];
for i=2:N+2
    % puntos incrementales
    wv=puntos(i-1,:);
    xi=wv(1)+inc_x;
    yi=wv(2)+inc_y;
    zi=wv(3)+inc_z;
    ti=wv(4)+inc_t;
    puntos=[puntos;xi, yi, zi, ti]; 
    % saltos de q
    qlist=[qlist;CinematicaInversa([xi, yi, zi])'];

end

% Velocidades euristica

qdlist=[0,0,0];
for i=2:N+1 %recorremos los puntos interiores
    wvp=qlist(i-1,:);   %working vector past
    wv=qlist(i,:);      %working vector 
    wvf=qlist(i+1,:);   %working vector future
    for j=1:3 % recorremos las tres articulaciones
        wn1=wvf(j)-wv(j);
        wn2=wv(j)-wvp(j);
        if(sign(wn1)==sign(wn2))
            qdlist(i,j)=((wn1)/inc_t+(wn2)/inc_t)/2;
        else
            qdlist(i,j)=0;
        end
    end    
end
qdlist(N+2,:)=[0,0,0];
    
% i=1;
wv=puntos(:,4);
q=qlist(1,:);
qd=qdlist(1,:);
qdd=[0,0,0];
% while wv(i)<=t
%     if wv(i)<t %&& wv(i+1)>t
%         for j=1:3
%             ar(j)=qlist(i,j);
%             br(j)=qdlist(i,j);
%             cr(j)=(3/inc_t^2)*(qlist(i+1,j)-qlist(i,j))-(qdlist(i+1,j)+2*qdlist(i,j))/inc_t;
%             dr(j)=(-2/inc_t^3)*(qlist(i+1,j)-qlist(i,j))+(qdlist(i+1,j)+qdlist(i,j))/inc_t^2;
%         end
%         q=ar+br*(t-wv(i))+cr*(t-wv(i))^2+dr*(t-wv(i))^3;
%         qd=br+2*cr*(t-wv(i))+3*dr*(t-wv(i))^2;
%         qdd=2*cr+6*dr*(t-wv(i));
%         break
%     else
%         if i>N+1  % nos encontramos en reposo tras finalizar el movimiento
%             q=qlist(N+2,:);
%             qd=qdlist(N+2,:);
%             qdd=[0,0,0];
%             break
%         end
%     end 
%     i= i+1;
% end
if t<=t_ini
    q=qlist(1,:);
    qd=qdlist(1,:);
    qdd=[0,0,0];
elseif t_ini<t && t<t_ini+T
    for i=1:N+1
        if t<wv(i)
            i=i-1;
            break
        end
    end
    for j=1:3
        ar(j)=qlist(i,j);
        br(j)=qdlist(i,j);
        cr(j)=(3/inc_t^2)*(qlist(i+1,j)-qlist(i,j))-(qdlist(i+1,j)+2*qdlist(i,j))/inc_t;
        dr(j)=(-2/inc_t^3)*(qlist(i+1,j)-qlist(i,j))+(qdlist(i+1,j)+qdlist(i,j))/inc_t^2;
    end
    q=ar+br*(t-wv(i))+cr*(t-wv(i))^2+dr*(t-wv(i))^3;
    qd=br+2*cr*(t-wv(i))+3*dr*(t-wv(i))^2;
    qdd=2*cr+6*dr*(t-wv(i));
else
    q=qlist(N+2,:);
    qd=qdlist(N+2,:);
    qdd=[0,0,0];
end

q1=q(1);        q2=q(2);        q3=q(3);

qp1=qd(1);      qp2=qd(2);      qp3=qd(3);

qpp1=qdd(1);    qpp2=qdd(2);    qpp3=qdd(3);


    % FIN DE INTERPORLACIÓN MEDIANTE POLINOMIOS DE TERCER ORDEN
trayectoria =  [q1; q2; q3; qp1; qp2; qp3; qpp1; qpp2; qpp3];
return;
