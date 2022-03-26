
function Q=MP_Hexagon(a,t,theta,chiral)
n=2^4;                                                                     % Element number for two beams (2^N)
% theta=30;                                                                   % Bending angle (degree)
if chiral==1
    type=1; % 1 for chiral, 2 for achiral
else
    type=2;
end

%% Material properties
E=71e9;                                                                    % Young's modulus (Pa)
G=27e9;                                                                    % Shear modulus
v=0.3;                                                                     % Poisson's ratio
ro=2700;                                                                   % Density (kg/m^3)
L0=a;                                                                  % Beam length (m)
% t=0.001;                                                                   % Beam thickness
h=1;                                                                       % Beam height
Iz=h*t^3/12;                                                               % Moment of inertia (m^4)
A=t*h;                                                                     % Cross section area (m^2)

%% Element stiffness matrix
syms t_var
L=L0/n;                                                                    % Element length
a0=L/2;
K_elm=[A*E/2/a0,             0,             0, -A*E/2/a0,             0,             0;
          0,  3*E*Iz/2/a0^3,  3*E*Iz/2/a0^2,        0, -3*E*Iz/2/a0^3,  3*E*Iz/2/a0^2;
          0,  3*E*Iz/2/a0^2,      2*E*Iz/a0,        0, -3*E*Iz/2/a0^2,        E*Iz/a0;
   -A*E/2/a0,             0,             0,  A*E/2/a0,             0,             0;
          0, -3*E*Iz/2/a0^3, -3*E*Iz/2/a0^2,        0,  3*E*Iz/2/a0^3, -3*E*Iz/2/a0^2;
          0,  3*E*Iz/2/a0^2,        E*Iz/a0,        0, -3*E*Iz/2/a0^2,     2*E*Iz/a0];           % Euler beam


%% K_OP, K_OG in 0 degree
ori=zeros(n,1);
if type==2
    for j=1:n     % achial
        ori(j)=theta*(1-1/n-(j-1)*2/n);
    end
elseif type==1
    for j=1:n/2  % chial
        ori(j)=theta*(1-2/n-(j-1)*4/n);
        ori(j+n/2)=-theta*(1-2/n-(j-1)*4/n);
    end                                                                        % Element orientation
else
    print('wrong type')
end
ori=ori/180*pi;
K_OP1=zeros(n*6);
for i=1:n
    T=[cos(ori(i)), sin(ori(i)), 0, 0, 0, 0;
        -sin(ori(i)), cos(ori(i)), 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0;
        0, 0, 0, cos(ori(i)), sin(ori(i)), 0;
        0, 0, 0, -sin(ori(i)), cos(ori(i)), 0;
        0, 0, 0, 0, 0, 1];

    Ke=T'*subs(K_elm,t_var,t)*T;
    K_OP1(6*(i-1)+1:6*(i-1)+6,6*(i-1)+1:6*(i-1)+6)=Ke;
end                                                                        % Diagonized Stiffness (not assembled)
I=[1,0,0;0,1,0;0,0,1];
TM_K=zeros(n*2*3,(n+1)*3);  % Combine the same nodes in K_OP1
TM_K(1:3,1:3)=I;
for j=2:n*2
    jj=floor(j/2);
    TM_K((j-1)*3+1:(j-1)*3+3,jj*3+1:jj*3+3)=I;
end                                                                       % Assembling matrix
K_OP2=TM_K'*K_OP1*TM_K;

K1(:,1:3)=K_OP2(:,1:3); K1(:,4:6)=K_OP2(:,3*n+1:3*n+3); K1(:,7:3*n+3)=K_OP2(:,4:3*n);
K2(1:3,:)=K1(1:3,:);K2(4:6,:)=K1(3*n+1:3*n+3,:);K2(7:3*n+3,:)=K1(4:3*n,:);
aa=K2(1:6,1:6);
ab=K2(1:6,7:end);
ba=K2(7:end,1:6);
bb=K2(7:end,7:end);
K_OP=aa-ab/bb*ba;
K_OP(abs(K_OP)<1e-3)=0;

%% K_OQ, half thickness
K_OQ1=zeros(n*6);
for i=1:n
    T=[cos(ori(i)), sin(ori(i)), 0, 0, 0, 0;
        -sin(ori(i)), cos(ori(i)), 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0;
        0, 0, 0, cos(ori(i)), sin(ori(i)), 0;
        0, 0, 0, -sin(ori(i)), cos(ori(i)), 0;
        0, 0, 0, 0, 0, 1];

    Ke=T'*subs(K_elm,t_var,t/2)*T;
    K_OQ1(6*(i-1)+1:6*(i-1)+6,6*(i-1)+1:6*(i-1)+6)=Ke;
end                                                                        % Diagonized Stiffness (not assembled)
I=[1,0,0;0,1,0;0,0,1];
TM_K=zeros(n*2*3,(n+1)*3);  % Combine the same nodes in K_OP1
TM_K(1:3,1:3)=I;
for j=2:n*2
    jj=floor(j/2);
    TM_K((j-1)*3+1:(j-1)*3+3,jj*3+1:jj*3+3)=I;
end                                                                       % Assembling matrix
K_OQ2=TM_K'*K_OQ1*TM_K;

K1(:,1:3)=K_OQ2(:,1:3); K1(:,4:6)=K_OQ2(:,3*n+1:3*n+3); K1(:,7:3*n+3)=K_OQ2(:,4:3*n);
K2(1:3,:)=K1(1:3,:);K2(4:6,:)=K1(3*n+1:3*n+3,:);K2(7:3*n+3,:)=K1(4:3*n,:);
aa=K2(1:6,1:6);
ab=K2(1:6,7:end);
ba=K2(7:end,1:6);
bb=K2(7:end,7:end);
K_OQ=aa-ab/bb*ba;
K_OQ(abs(K_OQ)<1e-3)=0;


%% K_OABCDEFGHI
% beam orientation, [CB,DC,ED,FE,AF,BA]
ori=[0, -2*pi/3, 2*pi/3, 0, -2*pi/3, 2*pi/3, 0, -2*pi/3, 2*pi/3, 0, -2*pi/3, 2*pi/3];
% connection list, ABCDEF=1,2,3,4,5,6,7
cnct=[2,1; 2,3; 4,3; 4,5; 6,5; 6,1; 6,7; 8,1; 2,9; 10,3; 4,11; 12,5;];
K_AtoL=zeros(12*3); % 12 nodes in total
for i=1:length(ori)
    T=[cos(ori(i)), sin(ori(i)), 0, 0, 0, 0;
        -sin(ori(i)), cos(ori(i)), 0, 0, 0, 0;
        0, 0, 1, 0, 0, 0;
        0, 0, 0, cos(ori(i)), sin(ori(i)), 0;
        0, 0, 0, -sin(ori(i)), cos(ori(i)), 0;
        0, 0, 0, 0, 0, 1];

    Ke=T'*K_OP*T;
    Ke(abs(Ke)<1e-3)=0;
    node1=cnct(i,1);
    node2=cnct(i,2);
    K_AtoL(3*(node1-1)+1:3*(node1-1)+3,3*(node1-1)+1:3*(node1-1)+3)=...
    K_AtoL(3*(node1-1)+1:3*(node1-1)+3,3*(node1-1)+1:3*(node1-1)+3)+Ke(1:3,1:3);
    K_AtoL(3*(node1-1)+1:3*(node1-1)+3,3*(node2-1)+1:3*(node2-1)+3)=...
    K_AtoL(3*(node1-1)+1:3*(node1-1)+3,3*(node2-1)+1:3*(node2-1)+3)+Ke(1:3,4:6);
    K_AtoL(3*(node2-1)+1:3*(node2-1)+3,3*(node1-1)+1:3*(node1-1)+3)=...
    K_AtoL(3*(node2-1)+1:3*(node2-1)+3,3*(node1-1)+1:3*(node1-1)+3)+Ke(4:6,1:3);
    K_AtoL(3*(node2-1)+1:3*(node2-1)+3,3*(node2-1)+1:3*(node2-1)+3)=...
    K_AtoL(3*(node2-1)+1:3*(node2-1)+3,3*(node2-1)+1:3*(node2-1)+3)+Ke(4:6,4:6);

end                                                                        % Diagonized Stiffness & Mass matrix (not assembled)

%% Global stiffness & mass matrix: K_OABCDEF
K_AtoL(abs(K_AtoL)<1e-3)=0;

aa=K_AtoL(1:15,1:15);
ab=K_AtoL(1:15,16:end);
ba=K_AtoL(16:end,1:15);
bb=K_AtoL(16:end,16:end);
K_FtoL=bb-ba*inv(aa)*ab;

Kuu=K_FtoL(1:3,1:3); Kuk=K_FtoL(1:3,4:end);
Kp=-Kuu\Kuk;
Kp(abs(Kp)<1e-3)=0;

aa=K_FtoL(1:3,1:3);
ab=K_FtoL(1:3,4:end);
ba=K_FtoL(4:end,1:3);
bb=K_FtoL(4:end,4:end);
K_GtoL=bb-ba*inv(aa)*ab;

%% Taylor expansion
V=3*sqrt(3)*L0^2*2;
ori=[0,pi/3,2/3*pi,pi,4/3*pi,5/3*pi];
Lxy=2*L0;
syms e11 e22 e12 e21 k13 k23 phi0 u0x u0y real
% A12=(e12-e21)/2;
% u12=e21-phi0;
% u21=e12+phi0;
u0x=0; u0y=0;
u12=(e21+e12)/2;
u21=(e21+e12)/2;
phi0=(e21-e12)/2; % Omega==0
psi0=phi0;
psi1=psi0+(k13*Lxy*cos(ori(1))+k23*Lxy*sin(ori(1)));
psi2=psi0+(k13*Lxy*cos(ori(2))+k23*Lxy*sin(ori(2)));
psi3=psi0+(k13*Lxy*cos(ori(3))+k23*Lxy*sin(ori(3)));
psi4=psi0+(k13*Lxy*cos(ori(4))+k23*Lxy*sin(ori(4)));
psi5=psi0+(k13*Lxy*cos(ori(5))+k23*Lxy*sin(ori(5)));
psi6=psi0+(k13*Lxy*cos(ori(6))+k23*Lxy*sin(ori(6)));
u1x=u0x+e11*Lxy*cos(ori(1))+u12*Lxy*sin(ori(1));
u1y=u0y+e22*Lxy*sin(ori(1))+u21*Lxy*cos(ori(1));
u2x=u0x+e11*Lxy*cos(ori(2))+u12*Lxy*sin(ori(2));
u2y=u0y+e22*Lxy*sin(ori(2))+u21*Lxy*cos(ori(2));
u3x=u0x+e11*Lxy*cos(ori(3))+u12*Lxy*sin(ori(3));
u3y=u0y+e22*Lxy*sin(ori(3))+u21*Lxy*cos(ori(3));
u4x=u0x+e11*Lxy*cos(ori(4))+u12*Lxy*sin(ori(4));
u4y=u0y+e22*Lxy*sin(ori(4))+u21*Lxy*cos(ori(4));
u5x=u0x+e11*Lxy*cos(ori(5))+u12*Lxy*sin(ori(5));
u5y=u0y+e22*Lxy*sin(ori(5))+u21*Lxy*cos(ori(5));
u6x=u0x+e11*Lxy*cos(ori(6))+u12*Lxy*sin(ori(6));
u6y=u0y+e22*Lxy*sin(ori(6))+u21*Lxy*cos(ori(6));


d=[u1x;u1y;psi1;u2x;u2y;psi2;u3x;u3y;psi3;u4x;u4y;psi4;u5x;u5y;psi5;u6x;u6y;psi6;];
w=d'*K_GtoL*d/(2*V);
c1=diff(w,e11,e11);c2=diff(w,e11,e22);c3=diff(w,e11,e12);
c4=diff(w,e11,e21);c5=diff(w,e11,k13);c6=diff(w,e11,k23);
c7=diff(w,e22,e22);c8=diff(w,e22,e12);c9=diff(w,e22,e21);
c10=diff(w,e22,k13);c11=diff(w,e22,k23);c12=diff(w,e12,e12);
c13=diff(w,e12,e21);c14=diff(w,e12,k13);c15=diff(w,e12,k23);
c16=diff(w,e21,e21);c17=diff(w,e21,k13);c18=diff(w,e21,k23);
c19=diff(w,k13,k13);c20=diff(w,k13,k23);c21=diff(w,k23,k23);
Q=[c1, c2, c3, c4, c5, c6;
   c2, c7, c8, c9,c10,c11;
   c3, c8,c12,c13,c14,c15;
   c4, c9,c13,c16,c17,c18;
   c5,c10,c14,c17,c19,c20;
   c6,c11,c15,c18,c20,c21];
Q=double(Q);
Q(abs(Q)<1e-3)=0;

q = [c1, c2, c3 + c4, c3 - c4, c5, c6;
   c2, c7, c8 + c9, c8 - c9, c10, c11;
   (c3 + c4)/2, (c8 + c9)/2, (c12 + c16)/2 + c13, (c12 - c16)/2, (c14 + c17)/2, (c15 + c18)/2;
   (c3 - c4)/2, (c8 - c9)/2, (c12 - c16)/2, (c12 + c16)/2 - c13, (c14 - c17)/2, (c15 - c18)/2;
   c5, c10, c14 + c17, c14 - c17, c19, c20;
   c6, c11, c15 + c18, c15 - c18, c20, c21];
q=double(q);
q(:,3:4)=q(:,3:4)/2;
q(abs(q)<1e-3)=0;
% q(4,:)=[];q(:,4)=[];
qi=inv(q);

C=Q(1:4,1:4);
% D=Q(5:6,5:6);
% B=Q(1:4,5:6);
Qi=inv(Q(1:4,1:4));
% Qi=inv(Q);
Em=1/Qi(1,1)
vm=-Qi(2,1)/Qi(1,1)
Gm=1/(Qi(3,3)+Qi(3,4))/2
Bm=1/(Qi(1,1)+Qi(1,2))/2
yita=(Qi(1,3)+Qi(1,4))/Qi(1,1);
end