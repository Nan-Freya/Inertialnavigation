clc
clear all

load('mustrain.mat', 'imudata')


%% Mean of rotation rate from the given IMU dataset in rad/s 
w_x= mean(imudata(:,2));
w_y= mean(imudata(:,3));
w_z= mean(imudata(:,4));
w= [w_x;w_y;w_z];
w_e= norm(w);  % Earth rotation rate in rad/s

%% Mean of specific force from the given IMU dataset in m/s^2
f_x= mean(imudata(:,5));
f_y= mean(imudata(:,6));
f_z= mean(imudata(:,7));
f= [f_x;f_y;f_z];
g= norm(f);


%% Intial Alignment - Euler Angels in degree
roll = atan2(-f_y,-f_z)*180/pi;
pitch = atand(-f_x/sqrt(f_y^2 + f_z^2));
yaw = atand(w_y/w_x);

%% Rotations Matrix/DCM
c1 = [1    0          0;
      0  cosd(roll) -sind(roll);
      0  sind(roll)  cosd(roll)];

c2 = [cosd(pitch)       0      sind(pitch);
        0               1         0;   
      -sind(pitch)      0      cosd(pitch)]; 
  
c3 = [cosd(yaw)  -sind(yaw)  0;
      sind(yaw)  cosd(yaw)   0; 
        0          0         1];
  
Cb_n = c1*c2*c3; % Geamte Rotationmatrix
Cn_b = Cb_n';


%% Trace of the matrix
tr = trace(Cb_n);

%% Rotation angel from the Trace of the matrix in degree
rotation_angel = acosd((tr-1)/2); % unit degree

%% Eigenvector for eigenvalue =1 of matrix Cb_n
h = [Cb_n(3,2)-Cb_n(2,3);
     Cb_n(1,3)-Cb_n(3,1);
     Cb_n(2,1)-Cb_n(1,2)];

u = h/norm(h);

%% Quaternion
a = 0.5*sqrt(1+tr);
b = 1/(4*a)*(Cb_n(3,2)-Cb_n(2,3));
c = 1/(4*a)*(Cb_n(1,3)-Cb_n(3,1));
d = 1/(4*a)*(Cb_n(2,1)-Cb_n(1,2));
q=[ a;
    b;
    c;
    d];




