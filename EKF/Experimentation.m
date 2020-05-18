close all;
clear all;
clc;

dt=0.2; %time step
T_f=5; %final time
T=0:dt:T_f;
test_case = 1; %1=IR, 2=thermocouple, 3=both
position=1; %1=northwest, 2=southeast

x_actual=[8.5;1.5];
x_pred=[5.5;2]; 
thermo=[1, 9;
               1, 1;
               5, 5;
               9, 9;
               9, 1]; %thermocouple position    

A=[1 0;0 1];
Q=[0.01 0;0 0.01]; %uncertainty related to system model

time=struct2cell(open('time.mat'));
time_data=time{1}';

%northwest
if position==1;
    nw=struct2cell(open('block_north_west.mat'));
    nw_data=nw{1};
    if test_case==1
        nw_voltage_lab=[nw_data(:,6), nw_data(:,8)]';
    elseif test_case==2
        nw_voltage_lab=[nw_data(:,1), nw_data(:,2), nw_data(:,3), nw_data(:,4), nw_data(:,5)]';
        coeff1=coeffvalues(fit(time_data', nw_voltage_lab(1,:)', 'exp1'));
        coeff2=coeffvalues(fit(time_data', nw_voltage_lab(2,:)', 'exp1'));
        coeff3=coeffvalues(fit(time_data', nw_voltage_lab(3,:)', 'exp1'));
        coeff4=coeffvalues(fit(time_data', nw_voltage_lab(4,:)', 'exp1'));
        coeff5=coeffvalues(fit(time_data', nw_voltage_lab(5,:)', 'exp1'));
    else
        nw_voltage_lab=[nw_data(:,1), nw_data(:,2), nw_data(:,3), nw_data(:,4), nw_data(:,5), nw_data(:,6), nw_data(:,8)]';
        coeff1=coeffvalues(fit(time_data', nw_voltage_lab(1,:)', 'exp1'));
        coeff2=coeffvalues(fit(time_data', nw_voltage_lab(2,:)', 'exp1'));
        coeff3=coeffvalues(fit(time_data', nw_voltage_lab(3,:)', 'exp1'));
        coeff4=coeffvalues(fit(time_data', nw_voltage_lab(4,:)', 'exp1'));
        coeff5=coeffvalues(fit(time_data', nw_voltage_lab(5,:)', 'exp1'));
    end
elseif position==2
    se=struct2cell(open('block_south_east.mat'));
    se_data=se{1};
    if test_case==1
        se_voltage_lab=[se_data(:,7), se_data(:,9)]';
    elseif test_case==2
        se_voltage_lab=[se_data(:,1), se_data(:,2), se_data(:,3), se_data(:,4), se_data(:,5)]';
        coeff1=coeffvalues(fit(time_data', se_voltage_lab(1,:)', 'exp1'));
        coeff2=coeffvalues(fit(time_data', se_voltage_lab(2,:)', 'exp1'));
        coeff3=coeffvalues(fit(time_data', se_voltage_lab(3,:)', 'exp1'));
        coeff4=coeffvalues(fit(time_data', se_voltage_lab(4,:)', 'exp1'));
        coeff5=coeffvalues(fit(time_data', se_voltage_lab(5,:)', 'exp1'));
    else
        se_voltage_lab=[se_data(:,1), se_data(:,2), se_data(:,3), se_data(:,4), se_data(:,5), se_data(:,7), se_data(:,9)]';
        coeff1=coeffvalues(fit(time_data', se_voltage_lab(1,:)', 'exp1'));
        coeff2=coeffvalues(fit(time_data', se_voltage_lab(2,:)', 'exp1'));
        coeff3=coeffvalues(fit(time_data', se_voltage_lab(3,:)', 'exp1'));
        coeff4=coeffvalues(fit(time_data', se_voltage_lab(4,:)', 'exp1'));
        coeff5=coeffvalues(fit(time_data', se_voltage_lab(5,:)', 'exp1'));
    end
end

if test_case==1
    R=[3.49e-6, 0, 0, 0; 0, 5.67e-5, 0, 0; 0, 0, 5.75e-6, 0; 0, 0, 0, 5.75e-6]; %uncertainty related to measurement model
elseif test_case==2
    R=[7.09e-7, 0, 0, 0, 0;
         0, 7.09e-7, 0, 0, 0;
         0, 0, 7.09e-7, 0, 0;
         0, 0, 0, 7.09e-7, 0;
         0, 0, 0, 0, 7.09e-7];
elseif test_case==3
    R=[3.49e-6, 0, 0, 0, 0, 0, 0, 0, 0;
          0, 5.67e-5, 0, 0, 0, 0, 0, 0, 0;
          0, 0, 5.75e-6, 0, 0, 0, 0, 0, 0;
          0, 0, 0, 5.75e-6, 0, 0, 0, 0, 0;
          0, 0, 0, 0, 7.09e-7, 0, 0, 0, 0;
          0, 0, 0, 0, 0, 7.09e-7, 0, 0, 0;
          0, 0, 0, 0, 0, 0, 7.09e-7, 0, 0;
          0, 0, 0, 0, 0, 0, 0, 7.09e-7, 0;
          0, 0, 0, 0, 0, 0, 0, 0, 7.09e-7];
end
    
P_0=1; %initial guess for covariance

x=zeros(2,length(T)+1); %initialize state vector
x(:,1)=x_actual;
x_hat=x_pred; %initialize prediction state vector

for k=1 : length(T)
    w_k=[Q(1,1)*randn(1);Q(2,2)*randn(1)]; %system noise
    if test_case==1
        v_k=[R(1,1)*randn(1); R(2,2)*randn(1); R(3,3)*randn(1); R(4,4)*randn(1)]; %measurement noise
    elseif test_case==2
        v_k=[R(1,1)*randn(1); R(2,2)*randn(1); R(3,3)*randn(1); R(4,4)*randn(1); R(5,5)*randn(1)];
    elseif test_case==3
        v_k=[R(1,1)*randn(1); R(2,2)*randn(1); R(3,3)*randn(1); R(4,4)*randn(1);
                 R(5,5)*randn(1); R(6,6)*randn(1); R(7,7)*randn(1); R(8,8)*randn(1); R(9,9)*randn(1)];
    end

    x(:,k+1)=A*x(:,k)+w_k; %system model equation

    %measurement models
    if test_case==1
        y(:,k)=[7.4E-4*x(1,k)^2 - 0.0932*x(1,k) + 2.96;
               2.18E-4*x(1,k)^2 - 0.0527*x(1,k) + 3.23;
               4.55E-3*x(2,k)^2 - 0.236*x(2,k) + 3.34;
               4.55E-3*x(2,k)^2 - 0.236*x(2,k) + 3.34]+v_k; 
    elseif test_case==2
        thermo_dist=[sqrt((abs(x(1,k)-thermo(1,1)))^2+(abs(x(2,k)-thermo(1,2)))^2);
                              sqrt((abs(x(1,k)-thermo(2,1)))^2+(abs(x(2,k)-thermo(2,2)))^2);
                              sqrt((abs(x(1,k)-thermo(3,1)))^2+(abs(x(2,k)-thermo(3,2)))^2);
                              sqrt((abs(x(1,k)-thermo(4,1)))^2+(abs(x(2,k)-thermo(4,2)))^2);
                              sqrt((abs(x(1,k)-thermo(5,1)))^2+(abs(x(2,k)-thermo(5,2)))^2)];       
                          
        y(:,k)=[0.0006*thermo_dist(1)^2 - 0.0094*thermo_dist(1) + 1.4499;
                    0.0006*thermo_dist(2)^2 - 0.0094*thermo_dist(2) + 1.4499;
                    0.0006*thermo_dist(3)^2 - 0.0094*thermo_dist(3) + 1.4499;
                    0.0006*thermo_dist(4)^2 - 0.0094*thermo_dist(4) + 1.4499;
                    0.0006*thermo_dist(5)^2 - 0.0094*thermo_dist(5) + 1.4499] + v_k;
    elseif test_case==3
        thermo_dist=[sqrt((abs(x(1,k)-thermo(1,1)))^2+(abs(x(2,k)-thermo(1,2)))^2);
                              sqrt((abs(x(1,k)-thermo(2,1)))^2+(abs(x(2,k)-thermo(2,2)))^2);
                              sqrt((abs(x(1,k)-thermo(3,1)))^2+(abs(x(2,k)-thermo(3,2)))^2);
                              sqrt((abs(x(1,k)-thermo(4,1)))^2+(abs(x(2,k)-thermo(4,2)))^2);
                              sqrt((abs(x(1,k)-thermo(5,1)))^2+(abs(x(2,k)-thermo(5,2)))^2)];               
        
        y(:,k)=[7.4E-4*x(1,k)^2 - 0.0932*x(1,k) + 2.96;
               2.18E-4*x(1,k)^2 - 0.0527*x(1,k) + 3.23;
               4.55E-3*x(2,k)^2 - 0.236*x(2,k) + 3.34;
               4.55E-3*x(2,k)^2 - 0.236*x(2,k) + 3.34;
               0.0006*thermo_dist(1)^2 - 0.0094*thermo_dist(1) + 1.4499;
               0.0006*thermo_dist(2)^2 - 0.0094*thermo_dist(2) + 1.4499;
               0.0006*thermo_dist(3)^2 - 0.0094*thermo_dist(3) + 1.4499;
               0.0006*thermo_dist(4)^2 - 0.0094*thermo_dist(4) + 1.4499;
               0.0006*thermo_dist(5)^2 - 0.0094*thermo_dist(5) + 1.4499] +v_k ;
    end
                
    %prediction
    x_hat_k=A*x_hat;
    P_k=A*P_0*A'+Q;
    
    if test_case==1
        h=[7.4E-4*x_hat_k(1)^2 - 0.0932*x_hat_k(1) + 2.96;
             2.18E-4*x_hat_k(1)^2 - 0.0527*x_hat_k(1) + 3.23;
             4.55E-3*x_hat_k(2)^2 - 0.236*x_hat_k(2) + 3.34;
             4.55E-3*x_hat_k(2)^2 - 0.236*x_hat_k(2) + 3.34];
         
             %update equation   
        H_k=[1.48E-3*x_hat_k(1)-0.0932, 0;
                 4.36E-4*x_hat_k(1)-0.0527, 0;
                 0, 9.1E-3*x_hat_k(2)-0.236;
                 0, 9.1E-3*x_hat_k(2)-0.236];
             
    elseif test_case==2
        thermo_dist=[sqrt((abs(x_hat_k(1)-thermo(1,1)))^2+(abs(x_hat_k(2)-thermo(1,2)))^2);
                              sqrt((abs(x_hat_k(1)-thermo(2,1)))^2+(abs(x_hat_k(2)-thermo(2,2)))^2);
                              sqrt((abs(x_hat_k(1)-thermo(3,1)))^2+(abs(x_hat_k(2)-thermo(3,2)))^2);
                              sqrt((abs(x_hat_k(1)-thermo(4,1)))^2+(abs(x_hat_k(2)-thermo(4,2)))^2);
                              sqrt((abs(x_hat_k(1)-thermo(5,1)))^2+(abs(x_hat_k(2)-thermo(5,2)))^2)];       
                          
        h=[0.0006*thermo_dist(1)^2 - 0.0094*thermo_dist(1) + 1.4499;
             0.0006*thermo_dist(2)^2 - 0.0094*thermo_dist(2) + 1.4499;
             0.0006*thermo_dist(3)^2 - 0.0094*thermo_dist(3) + 1.4499;
             0.0006*thermo_dist(4)^2 - 0.0094*thermo_dist(4) + 1.4499;
             0.0006*thermo_dist(5)^2 - 0.0094*thermo_dist(5) + 1.4499];
         
        %update equation
        H_k=[-0.0012*abs(x_hat_k(1)-thermo(1,1))+0.0094*abs(x_hat_k(1)-thermo(1,1))/thermo_dist(1), -0.0012*abs(x_hat_k(2)-thermo(1,2))+0.0094*abs(x_hat_k(2)-thermo(1,2))/thermo_dist(1);
                  -0.0012*abs(x_hat_k(1)-thermo(2,1))+0.0094*abs(x_hat_k(1)-thermo(2,1))/thermo_dist(2), -0.0012*abs(x_hat_k(2)-thermo(2,2))+0.0094*abs(x_hat_k(2)-thermo(2,2))/thermo_dist(2);
                  -0.0012*abs(x_hat_k(1)-thermo(3,1))+0.0094*abs(x_hat_k(1)-thermo(3,1))/thermo_dist(3), -0.0012*abs(x_hat_k(2)-thermo(3,2))+0.0094*abs(x_hat_k(2)-thermo(3,2))/thermo_dist(3);
                  -0.0012*abs(x_hat_k(1)-thermo(4,1))+0.0094*abs(x_hat_k(1)-thermo(4,1))/thermo_dist(4), -0.0012*abs(x_hat_k(2)-thermo(4,2))+0.0094*abs(x_hat_k(2)-thermo(4,2))/thermo_dist(4);
                  -0.0012*abs(x_hat_k(1)-thermo(5,1))+0.0094*abs(x_hat_k(1)-thermo(5,1))/thermo_dist(5), -0.0012*abs(x_hat_k(2)-thermo(5,2))+0.0094*abs(x_hat_k(2)-thermo(5,2))/thermo_dist(5)];
    elseif test_case==3
        thermo_dist=[sqrt((abs(x_hat_k(1)-thermo(1,1)))^2+(abs(x_hat_k(2)-thermo(1,2)))^2);
                              sqrt((abs(x_hat_k(1)-thermo(2,1)))^2+(abs(x_hat_k(2)-thermo(2,2)))^2);
                              sqrt((abs(x_hat_k(1)-thermo(3,1)))^2+(abs(x_hat_k(2)-thermo(3,2)))^2);
                              sqrt((abs(x_hat_k(1)-thermo(4,1)))^2+(abs(x_hat_k(2)-thermo(4,2)))^2);
                              sqrt((abs(x_hat_k(1)-thermo(5,1)))^2+(abs(x_hat_k(2)-thermo(5,2)))^2)];
                          
         h=[7.4E-4*x_hat_k(1)^2 - 0.0932*x_hat_k(1) + 2.96;
             2.18E-4*x_hat_k(1)^2 - 0.0527*x_hat_k(1) + 3.23;
             4.55E-3*x_hat_k(2)^2 - 0.236*x_hat_k(2) + 3.34;
             4.55E-3*x_hat_k(2)^2 - 0.236*x_hat_k(2) + 3.34;
             0.0006*thermo_dist(1)^2 - 0.0094*thermo_dist(1) + 1.4499;
             0.0006*thermo_dist(2)^2 - 0.0094*thermo_dist(2) + 1.4499;
             0.0006*thermo_dist(3)^2 - 0.0094*thermo_dist(3) + 1.4499;
             0.0006*thermo_dist(4)^2 - 0.0094*thermo_dist(4) + 1.4499;
             0.0006*thermo_dist(5)^2 - 0.0094*thermo_dist(5) + 1.4499];
         
         %update equation
         H_k=[1.48E-3*x_hat_k(1)-0.0932, 0;
                 4.36E-4*x_hat_k(1)-0.0527, 0;
                 0, 9.1E-3*x_hat_k(2)-0.236;
                 0, 9.1E-3*x_hat_k(2)-0.236;
                 -0.0012*abs(x_hat_k(1)-thermo(1,1))+0.0094*abs(x_hat_k(1)-thermo(1,1))/thermo_dist(1), -0.0012*abs(x_hat_k(2)-thermo(1,2))+0.0094*abs(x_hat_k(2)-thermo(1,2))/thermo_dist(1);
                 -0.0012*abs(x_hat_k(1)-thermo(2,1))+0.0094*abs(x_hat_k(1)-thermo(2,1))/thermo_dist(2), -0.0012*abs(x_hat_k(2)-thermo(2,2))+0.0094*abs(x_hat_k(2)-thermo(2,2))/thermo_dist(2);
                 -0.0012*abs(x_hat_k(1)-thermo(3,1))+0.0094*abs(x_hat_k(1)-thermo(3,1))/thermo_dist(3), -0.0012*abs(x_hat_k(2)-thermo(3,2))+0.0094*abs(x_hat_k(2)-thermo(3,2))/thermo_dist(3);
                 -0.0012*abs(x_hat_k(1)-thermo(4,1))+0.0094*abs(x_hat_k(1)-thermo(4,1))/thermo_dist(4), -0.0012*abs(x_hat_k(2)-thermo(4,2))+0.0094*abs(x_hat_k(2)-thermo(4,2))/thermo_dist(4);
                 -0.0012*abs(x_hat_k(1)-thermo(5,1))+0.0094*abs(x_hat_k(1)-thermo(5,1))/thermo_dist(5), -0.0012*abs(x_hat_k(2)-thermo(5,2))+0.0094*abs(x_hat_k(2)-thermo(5,2))/thermo_dist(5)];
    end
    
    K=P_k*H_k'*inv(H_k*P_k*H_k'+R); %kalman gain
    
    x_hat=x_hat_k+K*(y(:,k)-h);
    P_0=(eye(2)-K*H_k)*P_k;
    
    %Store estimates
    x_S(:,k) = x_hat;

    if test_case==1
        y_hat(:,k) = [7.4E-4*x_hat(1)^2 - 0.0932*x_hat(1) + 2.96;
                        2.18E-4*x_hat(1)^2 - 0.0527*x_hat(1) + 3.23;
                        4.55E-3*x_hat(2)^2 - 0.236*x_hat(2) + 3.34;
                        4.55E-3*x_hat(2)^2 - 0.236*x_hat(2) + 3.34];
    elseif test_case==2
        thermo_dist=[sqrt((abs(x_hat(1)-thermo(1,1)))^2+(abs(x_hat(2)-thermo(1,2)))^2);
                              sqrt((abs(x_hat(1)-thermo(2,1)))^2+(abs(x_hat(2)-thermo(2,2)))^2);
                              sqrt((abs(x_hat(1)-thermo(3,1)))^2+(abs(x_hat(2)-thermo(3,2)))^2);
                              sqrt((abs(x_hat(1)-thermo(4,1)))^2+(abs(x_hat(2)-thermo(4,2)))^2);
                              sqrt((abs(x_hat(1)-thermo(5,1)))^2+(abs(x_hat(2)-thermo(5,2)))^2)]; 
        y_hat(:,k) = [0.0006*thermo_dist(1)^2 - 0.0094*thermo_dist(1) + 1.4499;
                            0.0006*thermo_dist(2)^2 - 0.0094*thermo_dist(2) + 1.4499;
                            0.0006*thermo_dist(3)^2 - 0.0094*thermo_dist(3) + 1.4499;
                            0.0006*thermo_dist(4)^2 - 0.0094*thermo_dist(4) + 1.4499;
                            0.0006*thermo_dist(5)^2 - 0.0094*thermo_dist(5) + 1.4499];
    elseif test_case==3
        thermo_dist=[sqrt((abs(x_hat(1)-thermo(1,1)))^2+(abs(x_hat(2)-thermo(1,2)))^2);
                              sqrt((abs(x_hat(1)-thermo(2,1)))^2+(abs(x_hat(2)-thermo(2,2)))^2);
                              sqrt((abs(x_hat(1)-thermo(3,1)))^2+(abs(x_hat(2)-thermo(3,2)))^2);
                              sqrt((abs(x_hat(1)-thermo(4,1)))^2+(abs(x_hat(2)-thermo(4,2)))^2);
                              sqrt((abs(x_hat(1)-thermo(5,1)))^2+(abs(x_hat(2)-thermo(5,2)))^2)]; 
        y_hat(:,k) = [7.4E-4*x_hat(1)^2 - 0.0932*x_hat(1) + 2.96;
                            2.18E-4*x_hat(1)^2 - 0.0527*x_hat(1) + 3.23;
                            4.55E-3*x_hat(2)^2 - 0.236*x_hat(2) + 3.34;
                            4.55E-3*x_hat(2)^2 - 0.236*x_hat(2) + 3.34;
                            0.0006*thermo_dist(1)^2 - 0.0094*thermo_dist(1) + 1.4499;
                            0.0006*thermo_dist(2)^2 - 0.0094*thermo_dist(2) + 1.4499;
                            0.0006*thermo_dist(3)^2 - 0.0094*thermo_dist(3) + 1.4499;
                            0.0006*thermo_dist(4)^2 - 0.0094*thermo_dist(4) + 1.4499;
                            0.0006*thermo_dist(5)^2 - 0.0094*thermo_dist(5) + 1.4499];
    end
end

for k=1 : length(time_data)
    if position==1
        if test_case==1
            short_sensor_position=mean(abs(log(nw_voltage_lab(1,k)/1.9837)/-0.0007))-6;
            mid_sensor_position=mean(abs(log(nw_voltage_lab(2,k)/1.494)/0.0004))-12;
        elseif test_case==2
            thermo1=mean(abs(log(nw_voltage_lab(1,k)/coeff1(1))/coeff1(2)));
            thermo2=mean(abs(log(nw_voltage_lab(2,k)/coeff2(1))/coeff2(2)));
            thermo3=mean(abs(log(nw_voltage_lab(3,k)/coeff3(1))/coeff3(2)));
            thermo4=mean(abs(log(nw_voltage_lab(4,k)/coeff4(1))/coeff4(2)));
            thermo5=mean(abs(log(nw_voltage_lab(5,k)/coeff5(1))/coeff5(2)));
        elseif test_case==3
            thermo1=mean(abs(log(nw_voltage_lab(1,k)/coeff1(1))/coeff1(2)));
            thermo2=mean(abs(log(nw_voltage_lab(2,k)/coeff2(1))/coeff2(2)));
            thermo3=mean(abs(log(nw_voltage_lab(3,k)/coeff3(1))/coeff3(2)));
            thermo4=mean(abs(log(nw_voltage_lab(4,k)/coeff4(1))/coeff4(2)));
            thermo5=mean(abs(log(nw_voltage_lab(5,k)/coeff5(1))/coeff5(2)));
            short_sensor_position=mean(abs(log(nw_voltage_lab(6,k)/1.9837)/-0.0007))-6;
            mid_sensor_position=mean(abs(log(nw_voltage_lab(7,k)/1.494)/0.0004))-12;
        end
    elseif position==2
        if test_case==1
            short_sensor_position=mean(abs(log(se_voltage_lab(1,k)/1.0786)/0.00008))-6;
            long_sensor_position=mean(abs(log(se_voltage_lab(2,k)/2.6675)/-0.0005))-22;
        elseif test_case==2
            thermo1=mean(abs(log(se_voltage_lab(1,k)/coeff1(1))/coeff1(2)));
            thermo2=mean(abs(log(se_voltage_lab(2,k)/coeff2(1))/coeff2(2)));
            thermo3=mean(abs(log(se_voltage_lab(3,k)/coeff3(1))/coeff3(2)));
            thermo4=mean(abs(log(se_voltage_lab(4,k)/coeff4(1))/coeff4(2)));
            thermo5=mean(abs(log(se_voltage_lab(5,k)/coeff5(1))/coeff5(2)));
        elseif test_case==3
            thermo1=mean(abs(log(se_voltage_lab(1,k)/coeff1(1))/coeff1(2)));
            thermo2=mean(abs(log(se_voltage_lab(2,k)/coeff2(1))/coeff2(2)));
            thermo3=mean(abs(log(se_voltage_lab(3,k)/coeff3(1))/coeff3(2)));
            thermo4=mean(abs(log(se_voltage_lab(4,k)/coeff4(1))/coeff4(2)));
            thermo5=mean(abs(log(se_voltage_lab(5,k)/coeff5(1))/coeff5(2)));
            short_sensor_position=mean(abs(log(se_voltage_lab(6,k)/1.0786)/0.00008))-6;
            long_sensor_position=mean(abs(log(se_voltage_lab(7,k)/2.6675)/-0.0005))-22;
        end
    end
end

rectangle('Position',[x_actual(1)-1.5 x_actual(2)-1.5 3 3], 'EdgeColor', 'b', 'LineWidth', 2)
rectangle('Position',[x_hat(1)-1.5 x_hat(2)-1.5 3 3], 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
rectangle('Position',[mid_sensor_position-1.5 short_sensor_position-1.5 3 3], 'EdgeColor', 'black', 'LineWidth',2);
axis([0 10 0 10])
if test_case==1
    title('Localization using IR Sensors', 'FontSize', 18);
elseif test_case==2
    title('Localization using Thermocouples', 'FontSize',18);
else
    title('Localization using IR Sensors & Thermocouples', 'FontSize', 18);
end
xlabel('x (cm)', 'FontSize',16);
ylabel('y (cm)', 'FontSize', 16);
actual = line(NaN,NaN,'LineWidth',2,'Color','b');
final_estimate = line(NaN,NaN,'LineWidth',2,'LineStyle','--','Color','r');
ir = line(NaN,NaN,'LineWidth',2, 'Color', 'black');
grid on
pbaspect([1 1 1]);
set(gcf,'color','w');
legend('Actual','EKF', 'IR Sensors and Thermocouples');