function k_i=conductivitat_mitja_harmonica(cara,d_PI,i,j,x_P,y_P,x_CV,y_CV,k,mat)

if cara=='s'
d_Pi=y_P(j)-y_CV(j);
d_iI=y_CV(j)-y_P(j-1);
k_i=d_PI/(d_Pi/k(mat(i,j))+d_iI/k(mat(i,j-1)));
elseif cara=='n'
d_Pi=y_CV(j+1)-y_P(j);
d_iI=y_P(j+1)-y_CV(j+1);
k_i=d_PI/(d_Pi/k(mat(i,j))+d_iI/k(mat(i,j+1)));
elseif cara=='w'
d_Pi=x_P(i)-x_CV(i);
d_iI=x_CV(i)-x_P(i-1);
k_i=d_PI/(d_Pi/k(mat(i,j))+d_iI/k(mat(i-1,j)));
elseif cara=='e'
d_Pi=x_CV(i+1)-x_P(i);
d_iI=x_P(i+1)-x_CV(i+1);
k_i=d_PI/(d_Pi/k(mat(i,j))+d_iI/k(mat(i+1,j)));
end

end

            
%             if (i>=2 && i<=N_x1+1) && (j==N_y1+1) % Nodes amb cara NORTH compartida entre els materials 1 i 3
%                 d_PN=y_P(j+1)-y_P(j);
%                 k_n=conductivitat_mitja_harmonica('n',d_PN,i,j,x_P,y_P,x_CV,y_CV);
%                 a_N(i,j)=k_n*S_n(i,j)/d_PN;
%             elseif (i>=2 && i<=N_x1+1) && (j==N_y1+2) % Nodes amb cara SOUTH compartida entre els materials 1 i 3
%                 d_PS=y_P(j)-y_P(j-1);
%                 k_s=conductivitat_mitja_harmonica('s',d_PS,i,j,x_P,y_P,x_CV,y_CV);
%                 a_S(i,j)=k_s*S_s(i,j)/d_PS;
%                 
%             elseif (i==N_x1+1) && (j>=2 && j<=N_y1+1) % Nodes amb cara EAST compartida entre els materials 1 i 2
%                 d_PE=x_P(i+1)-x_P(i);
%                 k_e=conductivitat_mitja_harmonica('e',d_PE,i,j,x_P,y_P,x_CV,y_CV);
%                 a_E(i,j)=k_e*S_e(i,j)/d_PE;
%             elseif (i==N_x1+2) && (j>=2 && j<=N_y1+1) % Nodes amb cara WEST compartida entre els materials 1 i 2
%                 d_PW=x_P(i)-x_P(i-1);
%                 k_w=conductivitat_mitja_harmonica('w',d_PW,i,j,x_P,y_P,x_CV,y_CV);
%                 a_W(i,j)=k_w*S_w(i,j)/d_PW;
%                 
%             elseif (i==N_x1+1 || i==N_x1+2) && (j>=N_y1+2 && j<=N_y1+N_y2+1) % Nodes amb cara EAST compartida entre els materials 2 i 3
%                 d_PE=x_P(i+1)-x_P(i);
%                 k_e=conductivitat_mitja_harmonica('e',d_PE,i,j,x_P,y_P,x_CV,y_CV);
%                 a_E(i,j)=k_e*S_e(i,j)/d_PE;
%             elseif (i==N_x1+1 || i==N_x1+2) && (j>=N_y1+2 && j<=N_y1+N_y2+1) % Nodes amb cara WEST compartida entre els materials 2 i 3
%                 
%                 
%             elseif (i==N_x1+1) && (j>=N_y1+N_y2+2 && j<=N_y+1) % Nodes amb cara EAST compartida entre els materials 3 i 4
%                 
%             elseif (i==N_x1+2) && (j>=N_y1+N_y2+2 && j<=N_y+1) % Nodes amb cara WEST compartida entre els materials 3 i 4
%                 
%             
%             elseif (i>=N_x1+2 && i<=N_x+1) && (j==N_y1+N_y2+1) % Nodes amb cara NORTH compartida entre els materials 2 i 4
%             
%             elseif (i>=N_x1+2 && i<=N_x+1) && (j==N_y1+N_y2+2) % Nodes amb cara SOUTH compartida entre els materials 2 i 4



%             
%             
%             
%             
%             if ((i>=2 && i<=N_x1+1) && (j==N_y1+1)) || ((i>=N_x1+2 && i<=N_x+1) && (j==N_y1+N_y2+1)) % Nodes amb cara NORTH compartida entre els materials 1&3 i 2&4
%                 d_PN=y_P(j+1)-y_P(j);
%                 k_n=conductivitat_mitja_harmonica('n',d_PN,i,j,x_P,y_P,x_CV,y_CV);
%                 a_N(i,j)=k_n*S_n(i,j)/d_PN;
%             elseif ((i>=2 && i<=N_x1+1) && (j==N_y1+2)) || ((i>=N_x1+2 && i<=N_x+1) && (j==N_y1+N_y2+2)) % Nodes amb cara SOUTH compartida entre els materials 1&3 i 2&4
%                 d_PS=y_P(j)-y_P(j-1);
%                 k_s=conductivitat_mitja_harmonica('s',d_PS,i,j,x_P,y_P,x_CV,y_CV);
%                 a_S(i,j)=k_s*S_s(i,j)/d_PS;
%             end
%             if (i==N_x1+1) && (j>=2 && j<=N_y+1) % Nodes amb cara EAST compartida entre els materials 1&2, 2&3 i 3&4
%                 d_PE=x_P(i+1)-x_P(i);
%                 k_e=conductivitat_mitja_harmonica('e',d_PE,i,j,x_P,y_P,x_CV,y_CV);
%                 a_E(i,j)=k_e*S_e(i,j)/d_PE;
%             elseif (i==N_x1+2) && (j>=2 && j<=N_y+1) % Nodes amb cara WEST compartida entre els materials 1&2, 2&3 i 3&4
%                 d_PW=x_P(i)-x_P(i-1);
%                 k_w=conductivitat_mitja_harmonica('w',d_PW,i,j,x_P,y_P,x_CV,y_CV);
%                 a_W(i,j)=k_w*S_w(i,j)/d_PW;
%             end
%             % Nodes interns que no comparteixen cara/es amb altres materials
%                 d_PS=y_P(j)-y_P(j-1);
%                 d_PN=y_P(j+1)-y_P(j);
%                 d_PW=x_P(i)-x_P(j-1);
%                 d_PE=x_P(i+1)-x_P(i);
%                 a_N(i,j)=k(mat(i,j))*S_n(i,j)/d_PN;
%                 a_W(i,j)=k(mat(i,j))*S_w(i,j)/d_PW;
%                 a_E(i,j)=k(mat(i,j))*S_e(i,j)/d_PE;
%                 a_S(i,j)=k(mat(i,j))*S_s(i,j)/d_PS;
%                 a_P(i,j)=a_N(i,j)+a_W(i,j)+a_E(i,j)+a_S(i,j)+rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))/delta_t;
%                 b_P(i,j)=rho_P(mat(i,j))*V_P(i,j)*C_p(mat(i,j))*T(i,j,n)/delta_t;
