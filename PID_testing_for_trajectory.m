clear; close; clc;

syms t
theta = load('theta.mat');
theta_new = theta.theta;
theta_dot = load('theta_dot.mat');
theta_dot_new = theta_dot.theta_dot;
theta_dot_new=[0 theta_dot_new']
starting_matrix = readmatrix("data_1.xlsx");
pt_matrix = starting_matrix(2:end,2:end);
X_ = pt_matrix(1,:);
Y_ = pt_matrix(2,:);
% plot(X',Y')
% hold on
points = [X_;Y_];
spline = (cscvn(points));
fnplt(spline)
way_points=cscvn(points(:,[1:end]));
num_of_points=100;
hold on
x1=[];
y1=[];
z1=[];
vx=[];
vy=[];
vz=[];
ax=[];
ay=[];
for iter=1:length(points)-1
    t1=linspace(way_points.breaks(iter), way_points.breaks(iter+1), num_of_points);
    eqn_x=way_points.coefs((iter-1)*2+1,4)+(way_points.coefs((iter-1)*2+1, 3)*(t-way_points.breaks(iter)))+(way_points.coefs((iter-1)*2+1, 2)*((t-way_points.breaks(iter))^2))+(way_points.coefs((iter-1)*2+1, 1)*(t-way_points.breaks(iter))^3);
    eqn_y=way_points.coefs((iter-1)*2+2,4)+(way_points.coefs((iter-1)*2+2, 3)*(t-way_points.breaks(iter)))+(way_points.coefs((iter-1)*2+2, 2)*((t-way_points.breaks(iter))^2))+(way_points.coefs((iter-1)*2+2, 1)*(t-way_points.breaks(iter))^3);   
    %eqn_z=way_points.coefs((iter-1)*3+3,4)+(way_points.coefs((iter-1)*3+3, 3)*(t-way_points.breaks(iter)))+(way_points.coefs((iter-1)*3+3, 2)*((t-way_points.breaks(iter))^2))+(way_points.coefs((iter-1)*3+3, 1)*(t-way_points.breaks(iter))^3);   
    vel_x=diff(eqn_x,t);
    vel_y=diff(eqn_y,t);
    %vel_z=diff(eqn_z,t);
    acc_x=diff(vel_x,t);
    acc_y=diff(vel_y,t);
    xxx=subs(eqn_x,t,t1);
    xxx=double(xxx);
    yyy=subs(eqn_y,t,t1);
    yyy=double(yyy);
    %zzz=subs(eqn_z,t,t1);
    x1=[x1 xxx];
    y1 =[y1 yyy];
    %z1= [z1 zzz];
     xxx=subs(vel_x,t,t1);
    yyy=subs(vel_y,t,t1);
    %zzz=subs(vel_z,t,t1);
    x1=double(x1);
    y1=double(y1);
    %z1=double(z1)
    vx=double([vx xxx]);
    vy =double([vy yyy]);
    xxx=subs(acc_x,t,t1);
    yyy=subs(acc_y,t,t1);
    ax=double([ax xxx]);
    ay=double([ay yyy]);
    %vz= double([vz zzz]);
end
theta_new = zeros(1,length(x1))
theta_dot_new=zeros(1,length(x1))
for i=1:length(x1)-1
    theta_new(1,i+1)= atan(y1(:,i+1)-y1(:,i)/x1(:,i+1)-x1(:,i))
    theta_dot_new(1,i+1)=theta_new(1,i+1)-theta_new(1,i)
end
%plot(x1,y1,'ok')
clear t
% dt = 0.1;
% ts = 5;
% t = 0:dt:ts;
 d = 0.2;
tot_points=num_of_points*(length(points)-1);
t=linspace(0, way_points.breaks(length(points)), tot_points);
dt=way_points.breaks(length(points))/tot_points;
l1=length(x1);
eta_d=[x1;y1;theta_new(1:l1)];
eta_dot_d=[vx;vy;theta_dot_new(1:l1)];
etad_dot_d=[ax;ay;zeros(l1,1)'];

% Simulation Parameters
%tuning parameters
Kp = diag([50,50,50]);%proportional gain
Kd = diag([30,30,30]);%differential gain
Ki= diag([3,3,3]);%integrator gain

% dt = 0.1;
% ts = 130;
% t = 0:dt:ts;
 d = 0.2;
% Initial Conditions
eta0 = zeros(3,1);
eta_dot0 = zeros(3,1);
zeta0 = zeros(3,1);
zetad0 = zeros(3,1);
zeta_dot0 = zeros(3,1);
zetad_dot0 = zeros(3,1);
error=zeros(3,1);

eta(:,1) = eta0;
eta_dot(:,1) = eta_dot0;
zeta(:,1) = zeta0;
zetad(:,1) = zetad0;
zeta_dot(:,1) = zeros(3,1);
zetad_dot0(:,1) = zeros(3,1);

% Robot Parameters
m = 10;
Iz = 0.05; % Inertia of robot
xbc=0; ybc= 0; % Coordinates of mass centre
Gamma = [1 1;0 0;-d d];


for i = 1:l1


    % Inertia Matrix D
    D = [   m       0          -ybc*m
            0       m           xbc*m
         -ybc*m   xbc*m   Iz+m*(xbc^2+ybc^2)];
     Gamma = [1 1;0 0;-d d];

      % Jacobian Matrix Based on Desired Values
    psid = eta_d(3,i);
    psid_dot = eta_dot_d(3,i);
    J_etad = [cos(psid),-sin(psid),0 ;sin(psid),cos(psid),0;0,0,1];
    J_dot_etad = [-sin(psid)*psid_dot,-cos(psid)*psid_dot,0;cos(psid)*psid_dot,-sin(psid)*psid_dot,0;0,0,0];
    ud = zetad(1);
    vd = zetad(2);
    rd = zetad(3);

    % Other Vectors
    n_vd = [-m*rd*(vd+xbc.*rd)
             m*rd*(ud-ybc.*rd)*0
             m*rd*((xbc.*ud)+(ybc.*vd))];



    % Input Vector
    if i == 1
        error_sum=zeros(3,1);
        tau(:,i) = zeros(3,1);
        error(:,i)=zeros(3,1);
        psi = 0;
        J_eta = [cos(psi),-sin(psi),0 ;sin(psi),cos(psi),0;0,0,1];
        u = zeta(1,i);
        v = zeta(2,i);
        r = zeta(3,i);

    else
        %%PD controller updated to PID
       
        tau(:,i) = D * ( inv(J_etad) * (etad_dot_d(:,i) - J_dot_etad*inv(J_etad)*eta_dot_d(:,i) )  ...
             + inv(J_eta)*( Kp*(eta_d(:,i) - eta(:,i-1)) + Kd*(eta_dot_d(:,i) - eta_dot(:,i-1)) +Ki*(error(:,i-1))) ) + n_vd;%Ki*(error_sum(:,i-1))
        zeta(:,i) = zeta(:,i-1)+dt*zeta_dot(:,i-1); %velocity update
        u(i) = zeta(1,i);
        v(i) = zeta(2,i);
        r(i) = zeta(3,i);
        n_v = [-m * r(i-1) * (v(i-1) + xbc.*r(i-1))
                m * r(i-1) * (u(i-1) - ybc.*r(i-1)) * 0
                m * r(i-1) * ((xbc.*u(i-1)) + (ybc.*v(i-1)))];
        zeta_dot(:,i) = inv(D)*(tau(:,i) - n_v - 0.5*zeta(:,i));
        wheels(:,i)= (pinv(Gamma))*zeta(:,i);
        zeta(:,i)=Gamma*wheels(:,i);

        % Jacobian Matrix
        psi = eta(3,i-1);
        J_eta = [cos(psi),-sin(psi),0 ; sin(psi),cos(psi),0 ; 0,0,1];

        eta_dot(:,i) = (J_eta*(zeta(:,i)+dt*zeta_dot(:,i)));
        eta(:,i) = eta(:,i-1)+ dt*(J_eta*(zeta(:,i)+dt*zeta_dot(:,i))); %state update
        error(:,i) = error(:,i-1)+(eta_d(:,i) - eta(:,i-1));%error update
         
    end
end

figure
tiledlayout(4,1)
nexttile
plot(eta_d(1,:))
hold on
plot(eta(1,:),'r--')
title('Desired vs Actual X Position')

nexttile
plot(eta_d(2,:))
hold on
plot(eta(2,:),'r--')
title('Desired vs Actual Y Position')

nexttile
plot(eta_d(3,:))
hold on
plot(eta(3,:),'r--')
title('Desired vs Actual \theta Position')

temp1 = eta_d(3,1:1301);
temp2 = eta(3,:);

figure
tiledlayout(1,4)
nexttile
plot(eta_d(1,:)-eta(1,:))
title('Error in X Position')

nexttile
plot(eta_d(2,:)-eta(2,:))
title('Error in Y Position')

nexttile
plot (eta_d(3,:)-eta(3,:))
title('Error in theta Position')
xlabel('Time(s)')
ylabel('Theta (rads)')
nexttile
%Animation
l = 0.4;
w = 0.2;
Lw=0
Rw=0
mr_co = [l/2,l/2,l/2,-l/2,-l/2;-w/2,-w/2,w/2,w/2,-w/2];


for i = 1:length(t)
    psi = eta_d(3,i);
    R = [cos(psi),-sin(psi) ; sin(psi),cos(psi)];
    v_pos = R*mr_co;
    fill(v_pos(1,:) + eta(1,i),v_pos(2,:) + eta(2,i),'g')
    hold on ; grid on;axis([0 150 0 150]),axis square
    plot(eta(1,1:i),eta(2,1:i),'b-')

    pause(0.01)
    hold off
end



