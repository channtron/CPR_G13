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
for i=2:N+2
    % puntos incrementales
    wv=puntos(i-1,:);
    xi=wv(1)+inc_x;
    yi=wv(2)+inc_y;
    zi=wv(3)+inc_z;
    ti=wv(4)+inc_t;
    puntos=[puntos;xi, yi, zi, ti]; 
    % saltos de q

end

if t<=t_ini
    q=CinematicaInversa([puntos(1,1),puntos(1,2),puntos(1,3)]');
    qd=[0,0,0];
    qdd=[0,0,0];
elseif t_ini<t && t<t_ini+T
    wvt=puntos(:,4); %working vector time
    for i=1:N+1
        if t<wvt(i)
            i=i-1;
            break
        end
    end
    % Calculamos unicamente las q y qd relevantes para el segmento, qi-1 qi
    % qi+1 qi+2, qdi qdi+1
    switch i 
        case 1 % Nos encontramos en el primer segmento
            for j=i:i+2
            qlist(j,:)=CinematicaInversa([puntos(j,1),puntos(j,2),puntos(j,3)]');
            end
            qdlist(i,:)=[0,0,0];
            %velocidad en el punto j=i+1
            j=i+1;
            wvp=qlist(j-1,:);   %working vector past
            wv=qlist(j,:);      %working vector 
            wvf=qlist(j+1,:);   %working vector future
            for k=1:3 % recorremos las tres articulaciones
                wn1=wvf(k)-wv(k);
                wn2=wv(k)-wvp(k);
                if(sign(wn1)==sign(wn2))
                    qdlist(j,k)=((wn1)/inc_t+(wn2)/inc_t)/2;
                else
                    qdlist(j,k)=0;
                end
            end 
            
        case N+1 % Nos encontramos en el ultimo segmento
            for j=i-1:i+1
            qlist(j,:)=CinematicaInversa([puntos(j,1),puntos(j,2),puntos(j,3)]');
            end
            qdlist(i+1,:)=[0,0,0];
            %velocidades en el punto j=i
            j=i;
            wvp=qlist(j-1,:);   %working vector past
            wv=qlist(j,:);      %working vector 
            wvf=qlist(j+1,:);   %working vector future
            for k=1:3 % recorremos las tres articulaciones
                wn1=wvf(k)-wv(k);
                wn2=wv(k)-wvp(k);
                if(sign(wn1)==sign(wn2))
                    qdlist(j,k)=((wn1)/inc_t+(wn2)/inc_t)/2;
                else
                    qdlist(j,k)=0;
                end
            end
            
        otherwise % Nos encontramos en algun segmento intermedio
            for j=i-1:i+2
            qlist(j,:)=CinematicaInversa([puntos(j,1),puntos(j,2),puntos(j,3)]');
            end
            for j=i:i+1
                wvp=qlist(j-1,:);   %working vector past
                wv=qlist(j,:);      %working vector 
                wvf=qlist(j+1,:);   %working vector future
                for k=1:3 % recorremos las tres articulaciones
                    wn1=wvf(k)-wv(k);
                    wn2=wv(k)-wvp(k);
                    if(sign(wn1)==sign(wn2))
                        qdlist(j,k)=((wn1)/inc_t+(wn2)/inc_t)/2;
                    else
                        qdlist(j,k)=0;
                    end
                end
            end
    end
    
    for j=1:3
        ar(j)=qlist(i,j);
        br(j)=qdlist(i,j);
        cr(j)=(3/inc_t^2)*(qlist(i+1,j)-qlist(i,j))-(qdlist(i+1,j)+2*qdlist(i,j))/inc_t;
        dr(j)=(-2/inc_t^3)*(qlist(i+1,j)-qlist(i,j))+(qdlist(i+1,j)+qdlist(i,j))/inc_t^2;
    end
    q=ar+br*(t-wvt(i))+cr*(t-wvt(i))^2+dr*(t-wvt(i))^3;
    qd=br+2*cr*(t-wvt(i))+3*dr*(t-wvt(i))^2;
    qdd=2*cr+6*dr*(t-wvt(i));
else
    q=CinematicaInversa([puntos(N+2,1),puntos(N+2,2),puntos(N+2,3)]');
    qd=[0,0,0];
    qdd=[0,0,0];
end

% q1=q(1);        q2=q(2);        q3=q(3);
% 
% qp1=qd(1);      qp2=qd(2);      qp3=qd(3);
% 
% qpp1=qdd(1);    qpp2=qdd(2);    qpp3=qdd(3);


    % FIN DE INTERPORLACIÓN MEDIANTE POLINOMIOS DE TERCER ORDEN
trayectoria =  [q(1); q(2); q(3); qd(1); qd(2); qd(3); qdd(1); qdd(2); qdd(3)];
return;
