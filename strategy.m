function [Tau, Feedbck, Feedfwd, Ymstar, C, G, e, r, YtauSum, YYsum, theta_hat_dot] = fcn(alpha, beta, Gamma, g, theta, phi_des, phi_dot_des, phi_dotdot_des, phi, phi_dot, phi_dotdot, theta_hat, sumterm)

c1 = cos(phi(1));
c2 = cos(phi(2));
s1 = sin(phi(1));
s2 = sin(phi(2));
c12 = cos(phi(1) + phi(2));
s12 = sin(phi(1) + phi(2));

Ymstar = [theta(1) + 2*theta(2)*c2,theta(3)+theta(2)*c2; theta(3)+theta(2)*c2,theta(3)];
e = phi_des - phi;
e_dot = phi_dot_des - phi_dot;
% e_dotdot = phi_dotdot_des - phi_dotdot;
r = e_dot + alpha*e;
% r_dot = e_dotdot + alpha*e_dot;
eta = phi_dotdot_des + alpha*e_dot;
Y_M = [eta(1), (2*c2*eta(1)+ c2*eta(2)), eta(2), 0, 0; 0, (c2*eta(1)), (eta(1)+eta(2)), 0, 0];
M = Y_M*theta;
Y_C = [0, -(2*s2*phi_dot(1)*phi_dot(2)+s2*phi_dot(2)*phi_dot(2)), 0, 0, 0; 0, (s2*phi_dot(1)*phi_dot(1)), 0, 0, 0];
C = Y_C*theta;
Y_G = [0, 0, 0, (c1), (c12); 0, 0, 0, 0, (c12)];
G = Y_G*theta;
Y_MD = [0, -(2*s2*phi_dot(2)*r(1) + s2*phi_dot(2)*r(2)), 0, 0, 0; 0, (-s2*phi_dot(2)*r(1)), 0, 0, 0];
Y = Y_M + Y_C + Y_G + 0.5*Y_MD;
% theta_hat_dot = Gamma*(Y')*r;
Tau = Y*theta_hat + e + beta*r;
Feedbck = e + beta*r;
Feedfwd = Y*theta_hat;

% CL
% Ytau = Y'*Tau;
% YtYThetaHat = Y'*Y;

Ybuff = zeros(2,5,10);
YYsum = zeros(5,5);
YtauSum = zeros(5,1);
YYsumMinEig = 0.0;
lambdaCL = 0.1;
lambdaCLMet = false;
YYminDiff = 0.1;
if lambdaCLMet == false
    YSV = svd(Y);
    if (min(YSV) > YYminDiff && norm(Tau) > YYminDiff)
        minDiff = 100.0;
        for Yi= 1:size(Ybuff,1)
                    YYdiffi = norm(Yi-Y)
                    if YYdiffi < minDiff
                        minDiff = YYdiffi;
                    end
        end            
        if minDiff > YYminDiff
                    Ybuff = circshift(Ybuff, 1); Ybuff(:,:,1) = Y;
                    YY = Y'*Y;
                    Ytau = Y'*Tau;
                    YYsum = YYsum + YY;
                    YtauSum = YtauSum + Ytau;
                    YYsumMinEig = min(eig(YYsum));
                    if YYsumMinEig > lambdaCL
%                         TCL = t
                        lambdaCLMet = true;
                    end
        end
    end
end    

kcl = 0.1;
theta_hat_dot = Gamma*(Y')*r + Gamma*sumterm;









