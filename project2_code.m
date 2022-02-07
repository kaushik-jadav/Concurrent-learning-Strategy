clear all;
clc;
init = 0.001;
last = 10;
m_init = 1;
m_last = 3;
l_init = 0.25;
l_last = 0.75;
N = 1;
tf = 100;
ts = 1/100;
t = [0:ts:tf];
g=9.8;

% Gains
alpha_NL = 25;
beta_NL = 1;
Gamma_NL = eye(5)*1e-9;
phi_des_bar = [pi/8, pi/4];
f_phi_des = 0.2;
b_phi_des = [pi/2, pi/4];
a_phi_des = pi/2;

 %INITIAL CONDITIONS
phi_dotdot_IC = [0; 0];
phi_dot_IC = [0; 0];

cntnotry = 1;
cont = 1;
for i = 1:10
    % Theta 
%     l = [l_init+ (l_last - l_init)*rand(1); l_init+ (l_last - l_init)*rand(1)];
%     m = [m_init+ (m_last - m_init)*rand(1); m_init+ (m_last - m_init)*rand(1)];
    l1 = l_init+ (l_last - l_init)*rand(1); 
    l2 = l_init+ (l_last - l_init)*rand(1);
    m1 = m_init+ (m_last - m_init)*rand(1);
    m2 = m_init+ (m_last - m_init)*rand(1);
%     l1 = l_init + (l_last+l_init)*rand(N,1);
%     l2 = l_init + (l_last+l_init)*rand(N,1);
%     m1 =l_init + (l_last+l_init)*rand(N,1);
%     m2 = l_init + (l_last+l_init)*rand(N,1);
    theta1 = m1*l1*l1 + m2*l1*l1 + m2*l2*l2; 
    theta2 = m2*l1*l2;
    theta3 = m2*l2*l2;
    theta4 = (m1 + m2)*l1*g;
    theta5 = m2*l2*g;
    theta = [theta1, theta2, theta3, theta4, theta5];

    l1 = l_init+ (l_last - l_init)*rand(1);
    l2 = l_init+ (l_last - l_init)*rand(1);
    m1 = m_init+ (m_last - m_init)*rand(1);
    m2 = m_init+ (m_last - m_init)*rand(1);
    theta1 = m1*l1*l1 + m2*l1*l1 + m2*l2*l2; 
    theta2 = m2*l1*l2;
    theta3 = m2*l2*l2;
    theta4 = (m1 + m2)*l1*g;
    theta5 = m2*l2*g;
    new_theta = [theta1; theta2; theta3; theta4; theta5];

    THD_IC = new_theta;
    
    % Theta_hat
%     ll = [l_init+ (l_last - l_init)*rand(1); l_init+ (l_last - l_init)*rand(1)];
%     mm = [m_init+ (m_last - m_init)*rand(1); m_init+ (m_last - m_init)*rand(1)];
%     thetaH1 = mm(1)*ll(1)*ll(1) + mm(2)*ll(1)*ll(1) + mm(2)*ll(2)*ll(2);
%     thetaH2 = mm(2)*ll(1)*ll(2);
%     thetaH3 = mm(2)*ll(2)*ll(2);
%     thetaH4 = (mm(1) + mm(2))*ll(1)*g;
%     thetaH5 = mm(2)*ll(2)*g;
%     THD_IC = [thetaH1, thetaH2, thetaH3, thetaH4, thetaH5];
    try
        out = sim('Proj2_INLC2')
%         out = sim('INLC_P2')
    %     figure(1)
    %     plot(out.e_nl)
    % %     ylim([-50,50]);
    %     figure(2)
    %     plot(out.r_nl)
    % %     ylim([-80,0]);
        enm_x(i,:) = out.e_nl.Data(1,:,:);
        enm_y(i,:) = out.e_nl.Data(2,:,:);
        rnm_x(i,:) = out.r_nl.Data(1,:,:);
        rnm_y(i,:) = out.r_nl.Data(2,:,:);
    
        Tnm_x(i,:) = out.Tau.Data(1,:,:);
        Tnm_y(i,:) = out.Tau.Data(2,:,:);
    
    
        ffnm_x(i,:) = out.Feedfwd.Data(1,:,:);
        ffnm_y(i,:) = out.Feedfwd.Data(2,:,:);
        fbnm_x(i,:) = out.Feedbck.Data(1,:,:);
        fbnm_y(i,:) = out.Feedbck.Data(2,:,:);
    
        
        t_t_x1(i,:) = out.theta_tilde.Data(:,1)';
        t_t_x2(i,:) = out.theta_tilde.Data(:,2)';
        t_t_x3(i,:) = out.theta_tilde.Data(:,3)';
        t_t_x4(i,:) = out.theta_tilde.Data(:,4)';
        t_t_x5(i,:) = out.theta_tilde.Data(:,5)';
        ttt_x(i,:)=[ t_t_x1(i,:),  t_t_x2(i,:), t_t_x3(i,:),  t_t_x4(i,:), t_t_x5(i,:)];
        catch
            norm(i)=1000;
            i;
            cont = cont +1;
    end
    
end
for jj = 1:size(enm_x,2)
%      Nex(jj) =  sqrt(sum(enm_x(jj,:).^2));
%      Ney(jj) = sqrt(sum(enm_y(jj,:).^2));
%      Nrx(jj) = sqrt(sum(rnm_x(jj,:).^2));
%      Nry(jj) = sqrt(sum(rnm_y(jj,:).^2));

     Nex(jj) =  sqrt(sum(enm_x(:,jj).^2));
     Ney(jj) = sqrt(sum(enm_y(:,jj).^2));
     Nrx(jj) = sqrt(sum(rnm_x(:,jj).^2));
     Nry(jj) = sqrt(sum(rnm_y(:,jj).^2));

     Tau_xx(jj) = sqrt(sum(Tnm_x(:,jj).^2));
     Tau_yy(jj) = sqrt(sum(Tnm_y(:,jj).^2));
    
     ff_xx(jj) = sum(ffnm_x(:,jj).^2)^0.5;
     ff_yy(jj) = sum(ffnm_y(:,jj).^2)^0.5;
     fb_xx(jj) = sum(fbnm_x(:,jj).^2)^0.5;
     fb_yy(jj) = sum(fbnm_y(:,jj).^2)^0.5;

%      ttx1(jj) = sum(t_t_x1(:,jj).^2)^0.5;
%      ttx2(jj) = sum(t_t_x2(:,jj).^2)^0.5;
%      ttx3(jj) = sum(t_t_x3(:,jj).^2)^0.5;
%      ttx4(jj) = sum(t_t_x4(:,jj).^2)^0.5;
%      ttx5(jj) = sum(t_t_x5(:,jj).^2)^0.5;
% 
%      tax(jj,:) = [ttx1(jj) ttx2(jj) ttx3(jj) ttx4(jj) ttx5(jj)];
end
for jj = 1:5

     ttx1(jj) = sum(t_t_x1(:,jj).^2)^0.5;
     ttx2(jj) = sum(t_t_x2(:,jj).^2)^0.5;
     ttx3(jj) = sum(t_t_x3(:,jj).^2)^0.5;
     ttx4(jj) = sum(t_t_x4(:,jj).^2)^0.5;
     ttx5(jj) = sum(t_t_x5(:,jj).^2)^0.5;

     tax(jj,:) = [ttx1(jj) ttx2(jj) ttx3(jj) ttx4(jj) ttx5(jj)];

     
end    
     
%      Tau_xx = sqrt(sum(Tnm_x.^2));
%      Tau_yy = sqrt(sum(Tnm_y.^2));
%     
%      ff_xx = sqrt(sum(ffnm_x.^2));
%      ff_yy = sqrt(sum(ffnm_y.^2));
%      fb_xx = sqrt(sum(fbnm_x.^2));
%      fb_yy = sqrt(sum(fbnm_y.^2));
    
%      total_error = [Nex, Ney];
     time = linspace(0,100,size(Nex,2));
         figure(1)
         hold on
         plot(time, Nex)
         xlabel('Time', fontsize=12)
         ylabel('Norm of error (x)', fontsize=12)
         
         plot(time,Ney)
         xlabel('Time', fontsize=12)
         ylabel('Norm of error (Y)', fontsize=12)
         legend('Error Norm - X', 'Error Norm - Y')
         ylim([-50,50]);

         figure(2)
         hold on
         plot(time,Nrx)
         xlabel('Time', fontsize=12)
         ylabel('Norm of filtered Tracking error (x)', fontsize=12)
         plot(time,Nry)
         xlabel('Time', fontsize=12)
         ylabel('Norm of filtered tracking error (y)', fontsize=12)
         xlim([0,5]);
         legend('filtered tracking error - X', 'filtered tracking error - Y')
         grid on
    
         figure(3)
         hold on
         plot(time,Tau_xx)
         plot(time,Tau_yy)
         xlabel('Time', fontsize=12)
         ylabel('Tau', fontsize=12)
    
    
         %xlim([0,3]);
         %ylim([0,1]);
         grid on
    
         figure(4)
         hold on
         plot(time,ff_xx)
         plot(time,ff_yy)
         xlabel('Time', fontsize=12)
         ylabel('feedforward', fontsize=12)
         legend('feedforward x', 'feedforward y')
         %xlim([0,3]);
         %ylim([0,7e+6]);
         grid on
    
         figure(5)
         hold on
         plot(time,fb_xx)
         plot(time,fb_yy)
         xlabel('Time', fontsize=12)
         ylabel('feedback', fontsize=12)
         legend('feedback for x', 'feedback for y');
         grid on

         figure(6)
         total_ff = [ff_xx, ff_yy];
         plot(total_ff)

         figure(7)
         total_fb = [fb_xx, fb_yy];
         plot(total_fb)

         figure(8)
         plot(enm_x)
         plot(enm_y)
         title('Error')

         figure(9)
         plot(rnm_x)
         plot(rnm_y)
         title('Filtered Error')

         figure(10)
         plot(out.theta_tilde)
         title('Parametric error')
%          figure(8)
%          plot(tax)
% 