%clear; close; clc;
syms t
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
num_of_points=50;
hold off
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
    acc_x=diff(vel_x,t)
    acc_y=diff(vel_y,t)
    xxx=subs(eqn_x,t,t1);
    xxx=double(xxx)
    yyy=subs(eqn_y,t,t1);
    yyy=double(yyy)
    %zzz=subs(eqn_z,t,t1);
    x1=[x1 xxx];
    y1 =[y1 yyy];
    %z1= [z1 zzz];
     xxx=subs(vel_x,t,t1);
    yyy=subs(vel_y,t,t1);
    %zzz=subs(vel_z,t,t1);
    x1=double(x1)
    y1=double(y1)
    %z1=double(z1)
    vx=double([vx xxx]);
    vy =double([vy yyy]);
    xxx=subs(acc_x,t,t1);
    yyy=subs(acc_y,t,t1)
    ax=double([ax xxx])
    ay=double([ay yyy])
    %vz= double([vz zzz]);
end
theta_new = zeros(1,length(x1))
theta_dot_new=zeros(1,length(x1))
for i=1:length(x1)-1
    theta_new(1,i+1)= atan(y1(:,i+1)-y1(:,i)/x1(:,i+1)-x1(:,i))
    theta_dot_new(1,i+1)=theta_new(1,i+1)-theta_new(1,i)
end
tot_points=num_of_points*(length(points)-1);
t=linspace(0, way_points.breaks(length(points)), tot_points);
dt=way_points.breaks(length(points))/tot_points;
l1=length(x1);
eta_d=[x1;y1;theta_new(1:l1)];
eta_dot_d=[vx;vy;theta_dot_new(1:l1)];
etad_dot_d=[ax;ay;zeros(l1,1)'];

% normal sliding mode control




q_inertial = zeros(3,length(t)); % Xinertial,Yinertial, Theta inertial

q_inertial_dot = zeros(3,length(t));

q_inertial_ddot = zeros(3,length(t));

q_robotic = zeros(3,length(t));

wheel_matrix = zeros(2,length(t)); % angular velocity phi right wheel;angular velocity phi left wheel


% wheel_matrixdot = zeros(2,length(t));


q_robotic_dot = zeros(3,length(t));

Rotational_matrix = zeros(3,3,length(t));

control_velocity = zeros(2,length(t)); % v linear ; w angular

control_velocity_desired = zeros(2,length(t));

control_velocity_dot = zeros(2,length(t));

error = zeros(3,length(t));

error_dot = zeros(3,length(t));

error_ddot = zeros(3,length(t));

p_d = zeros(3,length(t));

p_d_dot_d = zeros(3,length(t));


% constant parameter

mc = 17;

mw = 0.5;

Ic = 0.537;

Iw = 0.0023;

d = 0.05;

L = 0.24;

Im = 0.0011;

R = 0.095;

I = Ic + mc*((d)^2) + 2*(mw*L*L) + 2*(Im);

m = mc + 2*(mw);
% Robot Parameters
m = 10;
Iz = 0.05; % Inertia of robot
xbc=0; ybc= 0; % Coordinates of mass centre


 D = [   m       0          -ybc*m
            0       m           xbc*m
         -ybc*m   xbc*m   Iz+m*(xbc^2+ybc^2)];;%H

B=[1 0 ; 0 1];

S = zeros(3,length(t));

S_dot = zeros(3,length(t));

lambda = 1;

tau = zeros(3,1);

K = 1;
ki=1

error_sum=zeros(3,1)
for i = 1 : length(x1)

    % Desired Trajectory

%     eta_d(:,i) = [0.5*cos(t(i)) ; 0.5*sin(t(i)) ; 0];
% 
%     eta_dot_d(:,i) = [-0.5*sin(t(i)) ; 0.5*cos(t(i)) ; rad(t(i))];
% 
%     etad_dot_d(:,i) = [-0.5*cos(t(i)) ; -0.5*sin(t(i)) ; 0];

   

    V = [-m*q_inertial(3,i)*(q_inertial(1,i)+xbc.*q_inertial(3,i))
             m*q_inertial(3,i)*(q_inertial(2,i)-ybc.*q_inertial(3,i))*0
             m*q_inertial(3,i)*((xbc.*q_inertial(2,i))+(ybc.*q_inertial(1,i)))];

       

    error(:,i) = eta_d(:,i) - q_inertial(:,i);

    error_dot(:,i) = eta_dot_d(:,i) - q_inertial_dot(:,i);

    error_ddot(:,i) = eta_dot_d(:,i) - q_inertial_ddot(:,i);

   

    S(:,i) = error_dot(:,i) + lambda*error(:,i);

    S_dot(:,i) = error_ddot(:,i) + lambda*error_dot(:,i);

   
    error_sum(:,1)=error_sum(:,1)+error(:,i)
    tau(:,i) = K.*D*sign(S(:,i)) + V +ki*(error_sum(:,1)) +D*(etad_dot_d(:,i) - lambda.*error_dot(:,i));

   

    R = [cos(q_inertial(3,i)) sin(q_inertial(3,i)) 0;-sin(q_inertial(3,i)) cos(q_inertial(3,i)) 0 ; 0 0 1];

    T = [cos(q_inertial(3,i)) 0;sin(q_inertial(3,i)) 0;0 1];

    W = [1,1;0 ,0;0.2,-0.2];

    J = R*W;

    psi = q_inertial(3,i);

    J_q = [cos(psi) -sin(psi) 0 ; sin(psi) cos(psi) 0 ; 0 0 1];

   

    q_inertial_ddot(:,i) = inv(D)*(tau(:,i) - V) ;

    % wheel_matrix(:,i) = inv(P)*J*q_inertial_dot(:,i);
    

    q_inertial_dot(:,i) = J_q*q_inertial_ddot(:,i) + dt*q_inertial_ddot(:,i);
     wheel_matrix(:,i)=pinv(J)* q_inertial_dot(:,i);
     q_inertial(:,i+1)= q_inertial(:,i)+(dt*q_inertial_dot(:,i));

    q_inertial(:,i+1) = q_inertial(:,i) + dt.*q_inertial_dot(:,i);

   

    %error(:,i+1) = eta_d(:,i) - q_inertial(:,i);

end
%% Animation

figure
tiledlayout(3,1)
nexttile
plot(eta_d(1,:),'color','bl')
hold on
plot(q_inertial(1,:),'color','r')
hold off

nexttile
plot(eta_d(2,:),'color','bl')
hold on
plot(q_inertial(2,:),'color','r')
hold off

nexttile
plot(temp1)
hold on
plot(temp2,'r--')
title('Desired vs Actual \theta Position')
% plot(eta_d(3,:),'color','bl')
% hold on
% plot(q_inertial(3,:),'color','r')
% hold off


% l = 0.20;
% w = 0.10;
% Lw=0
% Rw=0
% mr_co = [l/2,l/2,l/2,-l/2,-l/2;-w/2,-w/2,w/2,w/2,-w/2];
% for i = 1:length(t)
%     psi = q_inertial(3,i);
%     R = [cos(psi),-sin(psi) ; sin(psi),cos(psi)];
%     v_pos = R*mr_co;
%     fill(v_pos(1,:) + q_inertial(1,i),v_pos(2,:) + q_inertial(2,i),'g')
%     hold on ; grid on;axis([-1 3 -1 3]),axis square
%     plot(q_inertial(1,1:i),q_inertial(2,1:i),'b-')
% 
%     pause(0.1)
%     hold off
% end