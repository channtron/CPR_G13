    function trayectoria = GTCL_R3GDL(in)
                       
Xini       = in(1);           % Posici�n cartesiana inicial
Yini       = in(2);           % 
Zini       = in(3);           % 
Xfin       = in(4);           % Posici�n cartesiana final
Yfin       = in(5);           % 
Zfin       = in(6);           % 
N          = in(7);           % N�mero de puntos intermedios
t_ini      = in(8);           % Tiempo de inicio del movimiento
T          = in(9);           % Duraci�n del movimiento
t          = in(10);          % Reloj

inc_t=(T/(N+1));
inc_x=(Xfin-Xini)/(N+1);
inc_y=(Yfin-Yini)/(N+1);
inc_z=(Zfin-Zini)/(N+1);

puntos=[Xini, Yini, Zini, t_ini];
qlist=[CinematicaInversa([Xini, Yini, Zini])'];
vtramlist=[0, 0, 0];
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
    
    
    qd2=((wvf(2)-wv(2))/inc_t+(wv(2)-wvp(2))/inc_t)/2;
    
    
    qd3=((wvf(3)-wv(3))/inc_t+(wv(3)-wvp(3))/inc_t)/2;

% q1=qlist(:,1);
% q2=qlist(:,2);
% q3=qlist(:,3);

%% Metodo Euristico

qdlist=

q1=0; q2=0; q3=0;

qp1=0; qp2=0; qp3=0;

qpp1=0; qpp2=0; qpp3=0;


    % FIN DE INTERPORLACI�N MEDIANTE POLINOMIOS DE TERCER ORDEN
trayectoria =  [q1; q2; q3; qp1; qp2; qp3; qpp1; qpp2; qpp3];
return;