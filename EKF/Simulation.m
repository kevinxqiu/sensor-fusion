close all;
clear all;
clc;

dt=0.2; %time step
T_f=5; %final time
T=0:dt:T_f;
test_case = 3; %1=IR, 2=thermocouple, 3=both

x_actual=[1.5;8.5];
x_pred=[3;2]; 
thermo=[1, 9;
               1, 1;
               5, 5;
               9, 9;
               9, 1]; %thermocouple position    

A=[1 0;0 1];
Q=[0.01 0;0 0.01]; %uncertainty related to system model

%uncertainty related to measurement model
if test_case==1
    R=[3.49e-6, 0, 0, 0; 0, 5.67e-5, 0, 0; 0, 0, 5.75e-6, 0; 0, 0, 0, 5.75e-6]; 
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
        H_k=[-0.0012*abs(x_hat_k(1)-thermo(1,1))+0.0094*abs(x_hat_k(1)-thermo(1,1))/thermo_dist(1),
                  -0.0012*abs(x_hat_k(2)-thermo(1,2))+0.0094*abs(x_hat_k(2)-thermo(1,2))/thermo_dist(1);
                  -0.0012*abs(x_hat_k(1)-thermo(2,1))+0.0094*abs(x_hat_k(1)-thermo(2,1))/thermo_dist(2),
                  -0.0012*abs(x_hat_k(2)-thermo(2,2))+0.0094*abs(x_hat_k(2)-thermo(2,2))/thermo_dist(2);
                  -0.0012*abs(x_hat_k(1)-thermo(3,1))+0.0094*abs(x_hat_k(1)-thermo(3,1))/thermo_dist(3),
                  -0.0012*abs(x_hat_k(2)-thermo(3,2))+0.0094*abs(x_hat_k(2)-thermo(3,2))/thermo_dist(3);
                  -0.0012*abs(x_hat_k(1)-thermo(4,1))+0.0094*abs(x_hat_k(1)-thermo(4,1))/thermo_dist(4),
                  -0.0012*abs(x_hat_k(2)-thermo(4,2))+0.0094*abs(x_hat_k(2)-thermo(4,2))/thermo_dist(4);
                  -0.0012*abs(x_hat_k(1)-thermo(5,1))+0.0094*abs(x_hat_k(1)-thermo(5,1))/thermo_dist(5),
                  -0.0012*abs(x_hat_k(2)-thermo(5,2))+0.0094*abs(x_hat_k(2)-thermo(5,2))/thermo_dist(5)];
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

%Plot full trajectory results
%figure;
%subplot(2,2,1)
%hold on
%plot(T,x(1,2:end)) %State
%plot(T,x_S(1,:)) %Estimation
%plot(T,xhat_S(1,:)) %Prediction
%title('X Position State and Estimates', 'FontSize', 18)
%xlabel('Time (s)', 'FontSize', 16);
%ylabel('Position (cm)', 'FontSize', 16);
%legend('State', 'Estimate');

%subplot(2,2,2)
%hold on;
%plot(T,x(2,2:end)) %State
%plot(T,x_S(2,:)) %Estimate
%plot(T,xhat_S(2,:)) %Prediction
%title('Y Position State and Estimates', 'FontSize', 18)
%xlabel('Time (s)', 'FontSize', 16);
%ylabel('Position (cm)', 'FontSize', 16);
%legend('State', 'Estimate');


rectangle('Position',[x_actual(1)-1.5 x_actual(2)-1.5 3 3], 'EdgeColor', 'b', 'LineWidth', 2)
rectangle('Position',[x_pred(1)-1.5 x_pred(2)-1.5 3 3], 'LineWidth', 2);
rectangle('Position',[x_S(1,1)-1.5 x_S(2,1)-1.5 3 3], 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', ':');
rectangle('Position',[x_hat(1)-1.5 x_hat(2)-1.5 3 3], 'EdgeColor', 'r', 'LineWidth', 2, 'LineStyle', '--');
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
predicted = line(NaN,NaN,'LineWidth',2, 'Color', 'black');
initial_estimate = line(NaN,NaN,'LineWidth',2,'LineStyle',':', 'Color', 'r');
final_estimate = line(NaN,NaN,'LineWidth',2,'LineStyle','--','Color','r');
grid on
pbaspect([1 1 1]);
set(gcf,'color','w');
legend('Actual','Predicted', 'Initial Estimate','Final Estimate');