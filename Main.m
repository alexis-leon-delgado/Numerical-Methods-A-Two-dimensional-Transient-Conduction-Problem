%% A Two-dimensional Transient Conduction Problem

% Autor: Alexis Leon

close all; clear; clc;
format long;

%% 1. Entrada de dades

% Dades físiques:
T_inicial=8; % [ºC]
x_p1=0.50; % [m]
y_p1=0.40; % [m]
x_p2=0.50; % [m]
y_p2=0.70; % [m]
x_p3=1.10; % [m]
y_p3=0.80; % [m]
rho_P=[1500.00, 1600.00, 1900.00, 2500.00]; % [kg/m^3]
C_p=[750.00, 770.00, 810.00, 930.00]; %[J/(kgK)]
k=[170.00, 140.00, 200.00, 140.00]; %[W/(mK)]
T_bottom=23.00; %[ºC]
Q_flow=60.00; %[W/m]
T_g=33.00; %[ºC]
alpha=9.00; %[W/(m^2 K)]
T_right=@(temps) 8.00+0.005*temps; %[ºC]

% Dades numèriques:
densitat_node=50;
N_x1=round(densitat_node*x_p1);
N_x2=round(densitat_node*(x_p3-x_p1));

N_y1=round(densitat_node*y_p1);
N_y2=round(densitat_node*(y_p2-y_p1));
N_y3=round(densitat_node*(y_p3-y_p2));

delta_t=0.1; %[s]
t_final=5000; %[s]
delta=1e-7;
beta=1;
fr=1;

%% 2. Càlculs previs

% Generació de la malla:
N_x=N_x1+N_x2;
N_y=N_y1+N_y2+N_y3;
delta_x1=x_p1/N_x1; % [m]
delta_x2=(x_p3-x_p1)/N_x2; % [m]
delta_y1=y_p1/N_y1; % [m]
delta_y2=(y_p2-y_p1)/N_y2; % [m]
delta_y3=(y_p3-y_p2)/N_y3; % [m]
x_CV=zeros(1,N_x+2);
for i=2:(N_x+2)
    if i<=N_x1+2
        x_CV(i)=(i-2)*delta_x1;
    else
        x_CV(i)=x_CV(N_x1+2)+(i-N_x1-2)*delta_x2;
    end
end
y_CV=zeros(1,N_y+2);
for j=2:(N_y+2)
    if j<=N_y1+2
        y_CV(j)=(j-2)*delta_y1;
    elseif (j>=N_y1+3) && (j<=N_y1+N_y2+2)
        y_CV(j)=y_CV(N_y1+2)+(j-N_y1-2)*delta_y2;
    else
        y_CV(j)=y_CV(N_y1+N_y2+2)+(j-N_y1-N_y2-2)*delta_y3;
    end
end
x_P=zeros(1,N_x+2);
for i=2:(N_x+2)
    if i<=N_x1+1
        x_P(i)=x_CV(i)+delta_x1/2;
    elseif i==N_x+2
        x_P(i)=x_p3;
    else
        x_P(i)=x_CV(i)+delta_x2/2;
    end 
end
y_P=zeros(1,N_y+2);
for j=2:(N_y+2)
    if j<=N_y1+1
        y_P(j)=y_CV(j)+delta_y1/2;
    elseif (j>=N_y1+2) && (j<=N_y1+N_y2+1)
        y_P(j)=y_CV(j)+delta_y2/2;
    elseif j==N_y+2
        y_P(j)=y_p3;
    else
        y_P(j)=y_CV(j)+delta_y3/2;
    end
end

% Geometria:
S_w=zeros(N_x+2,N_y+2);
S_w(:,2:N_y1+1)=delta_y1;
S_w(:,N_y1+2:N_y1+N_y2+1)=delta_y2;
S_w(:,N_y1+N_y2+2:N_y+1)=delta_y3;
S_e=S_w;

S_s=zeros(N_x+2,N_y+2);
S_s(2:N_x1+1,:)=delta_x1;
S_s(N_x1+2:N_x+1,:)=delta_x2;
S_n=S_s;

V_P=zeros(N_x+2,N_y+2);
V_P(2:N_x1+1,2:N_y1+1)=delta_x1*delta_y1*1; % Es considera una amplada de 1 [m]
V_P(2:N_x1+1,N_y1+2:N_y1+N_y2+1)=delta_x1*delta_y2*1;
V_P(2:N_x1+1,N_y1+N_y2+2:N_y+1)=delta_x1*delta_y3*1;
V_P(N_x1+2:N_x+1,2:N_y1+1)=delta_x2*delta_y1*1;
V_P(N_x1+2:N_x+1,N_y1+2:N_y1+N_y2+1)=delta_x2*delta_y2*1;
V_P(N_x1+2:N_x+1,N_y1+N_y2+2:N_y+1)=delta_x2*delta_y3*1;

mat=zeros(N_x+2,N_y+2);
mat(1:N_x1+1,1:N_y1+1)=1;
mat(1:N_x1+1,N_y1+2:end)=3;
mat(N_x1+2:end,1:N_y1+N_y2+1)=2;
mat(N_x1+2:end,N_y1+N_y2+2:end)=4;


% Definició de matrius
a_P=zeros(N_x+2,N_y+2);
a_W=zeros(N_x+2,N_y+2);
a_E=zeros(N_x+2,N_y+2);
a_S=zeros(N_x+2,N_y+2);
a_N=zeros(N_x+2,N_y+2);
b_P=zeros(N_x+2,N_y+2);
t=zeros(1,t_final/delta_t+1);
T=zeros(N_x+2,N_y+2,length(t));
Q_abs=zeros(N_x+2,N_y+2,length(t));
Q_ced=zeros(N_x+2,N_y+2,length(t));
Q_total_posterior=zeros(N_x+2,N_y+2);
terme_acumulatiu=zeros(N_x+2,N_y+2);
diferencia_conservacio_energia=zeros(N_x+2,N_y+2);

%% 3. Mapa inicial de temperatures

T(:,:,1)=T_inicial*ones(N_x+2,N_y+2);

%% 4. Càlcul de temperatures a t^(n+1)=t^n+delta_t

tic;
iteracio=zeros(1,length(t));  % Comptador d'iteracions
for n=1:(length(t)-1) % Inici del for temporal que supleix el punt 5 de l'algorisme
    t(n+1)=t(n)+delta_t;
    
% 4.1. Estimació del camp de temperatures:
T_temps_posterior_est=T(:,:,n);



desviacio_max=10;   % Inicialització de la desviació màxima
while desviacio_max>delta % Inici del while per avaluar la convergència
    iteracio(n)=iteracio(n)+1;
    T_temps_posterior_est_iteracio_anterior=T_temps_posterior_est;
% 4.2. Avaluació dels coeficients de discretització:
for i=1:(N_x+2)
    for j=1:(N_y+2)
        if (i==1) && (j>=2 && j<=N_y+1) % Nodes del contorn esquerre
            a_S(i,j)=0;
            a_N(i,j)=0;
            a_W(i,j)=0;
            d_PE=x_P(i+1)-x_P(i);
            a_E(i,j)=k(mat(i,j))*S_e(i,j)/d_PE;
            a_P(i,j)=a_E(i,j)+alpha*S_w(i,j)+rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))/delta_t;
            b_P(i,j)=alpha*S_w(i,j)*T_g+rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))*T(i,j,n)/delta_t;
        elseif (i==N_x+2) && (j>=2 && j<=N_y+1) % Nodes del contorn dret
            a_S(i,j)=0;
            a_N(i,j)=0;
            a_E(i,j)=0;
            %d_PW=x_P(i)-x_P(i-1);
            a_W(i,j)=0;%k(mat(i,j))*S_w(i,j)/d_PW;
            a_P(i,j)=1;%a_W(i,j)+rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))/delta_t;
            b_P(i,j)=T_right(t(n+1));%rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))*T_right(t(n))/delta_t; %T(i,j,n)
        elseif (i>=2 && i<=N_x+1) && (j==1) % Nodes del contorn inferior
            a_S(i,j)=0;
            a_N(i,j)=0;
            a_W(i,j)=0;
            a_E(i,j)=0;
            a_P(i,j)=1;
            b_P(i,j)=T_bottom;
        elseif (i>=2 && i<=N_x+1) && (j==N_y+2) % Nodes del contorn superior
            a_N(i,j)=0;
            a_W(i,j)=0;
            a_E(i,j)=0;
            d_PS=y_P(j)-y_P(j-1);
            a_S(i,j)=k(mat(i,j))*S_s(i,j)/d_PS;
            a_P(i,j)=a_S(i,j)+rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))/delta_t;
            b_P(i,j)=rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))*T(i,j,n)/delta_t+Q_flow*S_n(i,j);
        elseif (i>=2 && i<=N_x+1) && (j>=2 && j<=N_y+1) % Nodes interns
            % Cara NORTH:
            d_PN=y_P(j+1)-y_P(j);
            if ((i>=2 && i<=N_x1+1) && (j==N_y1+1)) || ((i>=N_x1+2 && i<=N_x+1) && (j==N_y1+N_y2+1)) % Nodes amb cara NORTH compartida entre els materials 1&3 i 2&4
                k_n=conductivitat_mitja_harmonica('n',d_PN,i,j,x_P,y_P,x_CV,y_CV,k,mat);
            else
                k_n=k(mat(i,j));
            end
            a_N(i,j)=k_n*S_n(i,j)/d_PN;
            
            % Cara SOUTH:
            d_PS=y_P(j)-y_P(j-1);
            if ((i>=2 && i<=N_x1+1) && (j==N_y1+2)) || ((i>=N_x1+2 && i<=N_x+1) && (j==N_y1+N_y2+2)) % Nodes amb cara SOUTH compartida entre els materials 1&3 i 2&4
                k_s=conductivitat_mitja_harmonica('s',d_PS,i,j,x_P,y_P,x_CV,y_CV,k,mat);
            else
                k_s=k(mat(i,j));
            end
            a_S(i,j)=k_s*S_s(i,j)/d_PS;
            
            % Cara EAST:
            d_PE=x_P(i+1)-x_P(i);
            if (i==N_x1+1) && (j>=2 && j<=N_y+1) % Nodes amb cara EAST compartida entre els materials 1&2, 2&3 i 3&4
                k_e=conductivitat_mitja_harmonica('e',d_PE,i,j,x_P,y_P,x_CV,y_CV,k,mat);
            else
                k_e=k(mat(i,j));
            end
            a_E(i,j)=k_e*S_e(i,j)/d_PE;
            
            % Cara WEST:
            d_PW=x_P(i)-x_P(i-1);
            if (i==N_x1+2) && (j>=2 && j<=N_y+1) % Nodes amb cara WEST compartida entre els materials 1&2, 2&3 i 3&4
                k_w=conductivitat_mitja_harmonica('w',d_PW,i,j,x_P,y_P,x_CV,y_CV,k,mat);
            else
                k_w=k(mat(i,j));
            end
            a_W(i,j)=k_w*S_w(i,j)/d_PW;
            
            % Resta de coeficients respecte P:
            a_P(i,j)=a_N(i,j)+a_W(i,j)+a_E(i,j)+a_S(i,j)+rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))/delta_t;
            b_P(i,j)=rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))*T(i,j,n)/delta_t;
                
        end % Fi del for dels tipus de nodes contorn o interns
    end % Fi del for de nodes j
end % Fi del for de nodes i
    

% 4.3. Resolució del sistema d'equacions: line-by-line

P=zeros(N_x+2,N_y+2);
Q=zeros(N_x+2,N_y+2);
b_P_linebyline=zeros(N_x+2,N_y+2);
for j=1:N_y+2 % For que canvia la línia o fila a resoldre
    % Càlcul de P i Q:
    for i=1:N_x+2
        if (i==1 && (j==1||j==N_y+2)) || (i==N_x+2 && (j==1||j==N_y+2))
            % Si es tracta d'un dels nodes de les cantonades, no es
            % realitza cap acció
        else
        if j==1
            b_P_linebyline(i,j)=b_P(i,j)+a_N(i,j)*T(i,j+1,n+1);
        elseif j==N_y+2
            b_P_linebyline(i,j)=b_P(i,j)+a_S(i,j)*T(i,j-1,n+1);
        else
            b_P_linebyline(i,j)=b_P(i,j)+a_N(i,j)*T(i,j+1,n+1)+a_S(i,j)*T(i,j-1,n+1);
        end
        if i==1% Node inicial
            P(i,j)=a_E(i,j)/a_P(i,j);    
            Q(i,j)=b_P_linebyline(i,j)/a_P(i,j);
        else % Resta de nodes
            P(i,j)=a_E(i,j)/(a_P(i,j)-a_W(i,j)*P(i-1,j));
            Q(i,j)=(b_P_linebyline(i,j)+a_W(i,j)*Q(i-1,j))/(a_P(i,j)-a_W(i,j)*P(i-1,j));  
        end
        end
    end
    
    % Càlcul de T:
    for i=N_x+2:-1:1 % Es recorren els nodes en ordre invers
        if i==N_x+2 % Cas de l'últim node
            T(i,j,n+1)=Q(i,j);
        else
            T(i,j,n+1)=P(i,j)*T(i+1,j,n+1)+Q(i,j);
        end
        T_temps_posterior_est(i,j)=T_temps_posterior_est_iteracio_anterior(i,j)+fr*(T(i,j,n+1)-T_temps_posterior_est_iteracio_anterior(i,j));
    end
    
end % Fi del for que canvia la línia o fila a resoldre


% 4.4. Condició de convergència
desviacio_max=max(max(abs(T_temps_posterior_est_iteracio_anterior-T(:,:,n+1))));
end % Fi del while per avaluar la convergència
end % Fi del for temporal
elapsedTime=toc;

%% 5. Càlculs finals i impressió de resultats

% instant=5000-delta_t;
% n=find(t==instant);
% for i=1:(N_x+2)
%     for j=1:(N_y+2)
%         if (i==1) && (j>=2 && j<=N_y+1) % Nodes del contorn esquerre
%             d_PE=x_P(i+1)-x_P(i);
%             Q_w=alpha*(T_g-T(i,j,n+1))*S_w(i,j);
%             Q_e=-k(mat(i,j))*(T(i+1,j,n+1)-T(i,j,n+1))*S_e(i,j)/d_PE;
%             Q_abs(i,j,n+1)=Q_w;
%             Q_ced(i,j,n+1)=Q_e;
%         elseif (i==N_x+2) && (j>=2 && j<=N_y+1) % Nodes del contorn dret
%             Q_abs(i,j,n+1)=0;
%             Q_ced(i,j,n+1)=0;
%         elseif (i>=2 && i<=N_x+1) && (j==1) % Nodes del contorn inferior
%             Q_abs(i,j,n+1)=0;
%             Q_ced(i,j,n+1)=0;
%         elseif (i>=2 && i<=N_x+1) && (j==N_y+2) % Nodes del contorn superior
%             d_PS=y_P(j)-y_P(j-1);
%             Q_n=Q_flow*S_n(i,j);
%             Q_s=-k(mat(i,j))*(T(i,j,n+1)-T(i,j-1,n))*S_s(i,j)/d_PS;
%             Q_abs(i,j,n+1)=Q_n+Q_s;
%             Q_ced(i,j,n+1)=0;
%         elseif (i>=2 && i<=N_x+1) && (j>=2 && j<=N_y+1) % Nodes interns
%             % Cara NORTH:
%             d_PN=y_P(j+1)-y_P(j);
%             if ((i>=2 && i<=N_x1+1) && (j==N_y1+1)) || ((i>=N_x1+2 && i<=N_x+1) && (j==N_y1+N_y2+1)) % Nodes amb cara NORTH compartida entre els materials 1&3 i 2&4
%                 k_n=conductivitat_mitja_harmonica('n',d_PN,i,j,x_P,y_P,x_CV,y_CV,k,mat);
%             else
%                 k_n=k(mat(i,j));
%             end
%             
%             % Cara SOUTH:
%             d_PS=y_P(j)-y_P(j-1);
%             if ((i>=2 && i<=N_x1+1) && (j==N_y1+2)) || ((i>=N_x1+2 && i<=N_x+1) && (j==N_y1+N_y2+2)) % Nodes amb cara SOUTH compartida entre els materials 1&3 i 2&4
%                 k_s=conductivitat_mitja_harmonica('s',d_PS,i,j,x_P,y_P,x_CV,y_CV,k,mat);
%             else
%                 k_s=k(mat(i,j));
%             end
%             
%             % Cara EAST:
%             d_PE=x_P(i+1)-x_P(i);
%             if (i==N_x1+1) && (j>=2 && j<=N_y+1) % Nodes amb cara EAST compartida entre els materials 1&2, 2&3 i 3&4
%                 k_e=conductivitat_mitja_harmonica('e',d_PE,i,j,x_P,y_P,x_CV,y_CV,k,mat);
%             else
%                 k_e=k(mat(i,j));
%             end
%             
%             % Cara WEST:
%             d_PW=x_P(i)-x_P(i-1);
%             if (i==N_x1+2) && (j>=2 && j<=N_y+1) % Nodes amb cara WEST compartida entre els materials 1&2, 2&3 i 3&4
%                 k_w=conductivitat_mitja_harmonica('w',d_PW,i,j,x_P,y_P,x_CV,y_CV,k,mat);
%             else
%                 k_w=k(mat(i,j));
%             end
%             
%             % Calors:
%             Q_w=-k(mat(i,j))*(T(i,j,n+1)-T(i-1,j,n+1))*S_w(i,j)/d_PW;
%             Q_e=-k(mat(i,j))*(T(i+1,j,n+1)-T(i,j,n+1))*S_e(i,j)/d_PE;
%             Q_s=-k(mat(i,j))*(T(i,j,n+1)-T(i,j-1,n))*S_s(i,j)/d_PS;
%             Q_n=-k(mat(i,j))*(T(i,j+1,n+1)-T(i,j,n+1))*S_n(i,j)/d_PN;
%             Q_abs(i,j,n+1)=Q_w+Q_s;
%             Q_ced(i,j,n+1)=Q_e+Q_n;
%                 
%         end % Fi del for dels tipus de nodes contorn o interns
%           
%         terme_acumulatiu(i,j)=rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))*(T(i,j,n+1)-T(i,j,n))/delta_t;
%         Q_total_posterior(i,j)=Q_abs(i,j,n+1)-Q_ced(i,j,n+1);
%         diferencia_conservacio_energia(i,j)=terme_acumulatiu(i,j)-Q_total_posterior(i,j);
%     end % Fi del for de nodes j
% end % Fi del for de nodes i
% Q_total=sum(sum(diferencia_conservacio_energia));
% diferencia_energia_total=sum(sum(diferencia_conservacio_energia));

% Obtenció de la temperatura als punts A i B a estudiar:
instant_a_estudiar=5000;
n=find(abs(t-instant_a_estudiar)==min(abs(t-instant_a_estudiar)));
% Punts A i B
x_punt_A=0.65; % [m]
y_punt_A=0.56; % [m]
x_punt_B=0.74; % [m]
y_punt_B=0.72; % [m]

x_punts=[x_punt_A x_punt_B];
y_punts=[y_punt_A y_punt_B];

T_punts=zeros(1,2);
for p=1:2
    
for i=1:N_x+1
    for j=1:N_y+1
       if (y_punts(p)<y_CV(j+1) && y_punts(p)>y_CV(j)) % El punt es troba entre la cara SUD i NORD del node [i,j]
           if (x_punts(p)<x_CV(i+1) && x_punts(p)>x_CV(i)) % El punt es troba entre la cara WEST i EAST del node [i,j]
               T_punts(p)=T(i,j,n);
           elseif x_punts(p)==x_CV(i) % El punt es troba a la cara WEST del node [i,j]
               T_punts(p)=1/2*(T(i-1,j,n)+T(i,j,n)); % Es fa la mitjana entre la T del node [i-1,j] i el [i,j]
           end
       elseif y_punts(p)==y_CV(j) % El punt es troba a la cara SUD del node [i,j]
           if (x_punts(p)<x_CV(i+1) && x_punts(p)>x_CV(i)) % El punt es troba entre la cara WEST i EAST del node [i,j]
               T_punts(p)=1/2*(T(i,j-1,n)+T(i,j,n)); % Es fa la mitjana entre la T del node [i,j-1] i el [i,j]
           elseif x_punts(p)==x_CV(i) % El punt es troba a la cara WEST del node [i,j]
               T_punts(p)=1/4*(T(i-1,j,n)+T(i,j,n)+T(i,j-1,n)+T(i-1,j-1,n)); % Es fa la mitjana entre la T dels nodes [i-1,j], [i,j], [i,j-1] i [i-1,j-1]
           end
       end
    end
end

end


%% Grafica

[X,Y]=meshgrid(x_CV,y_CV);

instants_temp=[5000];

%instants_temp=[1000,2000,3000,4000,5000,6000,7000,8000,9000,10000];
for index=1:length(instants_temp)
    
n=find(abs(t-instants_temp(index))==min(abs(t-instants_temp(index))));
T_actual=T(:,:,n);

figure(index);
set(figure(index),'Renderer', 'painters', 'Position',[50 50 1300 900])
colormap jet;
grafica=pcolor(X,Y,T_actual(:,:)');
ax=gca;
ax.TickDir = 'out';
grafica.EdgeAlpha=0.25;
xlabel('$X\;\mathrm{(m)}$','interpreter','latex','Fontsize',18);
ylabel('$Y\;\mathrm{(m)}$','interpreter','latex','Fontsize',18);
title(['Instant $t=$ ', num2str(instants_temp(index)),' s'],'interpreter','latex','Fontsize',18);
xticks([0:0.1:1.1]);
cbar=colorbar;
Tmin=min(T_actual(T_actual(:,:,end)~=0));
Tmax=max(max(T_actual(:,:,end)));
caxis([Tmin,Tmax]);
title(cbar,'$T\;\mathrm{(^{\circ}C)}$','interpreter','latex','Fontsize',12);
hold on;
isotermes=round(Tmin):Tmax;
[C,isot]=contour(X,Y,T_actual(:,:)',isotermes,'LineColor','k');
 clabel(C,isot,'Color','k','LabelSpacing',450,'Margin',0.1,'BackgroundColor',[0.7 0.7 0.7],'interpreter','latex');
end 



