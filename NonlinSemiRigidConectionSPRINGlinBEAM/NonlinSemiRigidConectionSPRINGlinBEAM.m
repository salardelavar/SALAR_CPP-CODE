%***********************************************************%
%                >> IN THE NAME OF GOD <<                   %
% Analysis of 1st order Linear and Nonlinear Semi Rigid     %
% Conection beam subjected to Pushover lateral load         %
% Checking the analysis by Force Control                    %
%-----------------------------------------------------------%
%     This program is written by salar delavar ghashghaei   %  
%          E-mail:salar.d.ghashghaei@gmail.com              %
%-----------------------------------------------------------%
%                                     Force                 %
%                                       |                   %
%                                       v                   %
%                  ||--------------------                   %
%                   |<-      L,EI      ->|                  %
%***********************************************************%
clear all;close all;clc
%% Define Parameters in mm,kN
P3=0; % [kN.mm]
P4=-1; % [kN]
P5=.075; % [kN]
P6=0; % [kN.mm]
L = 1000; % [mm]
EI= 200*100^4/12; % [kN.mm^2]
EA = 200*10000; % [kN]
m = 100; % number of calculation
itermax = 50;% maximum number of iterations
tolerance = 1e-12; % specified tolerance for convergence
lanX=1;
lanY=0;
% Nonlinear Rotational Spring
ty=.015; % Yeild Rotation (rad)
My=15e+3; % Yeild Moment (kN.mm)
tu=.025;  % Ultimate Rotation (rad)
Mu=18e+3;  % Ultimate Moment (kN.mm)
k1 = My/ty;
k2 =(Mu-My)/(tu-ty);
u = zeros(4,1);% Define intial displacement u
%% Element stifness coefficient
A=4*EI/L;B=6*EI/L^2;C=2*EI/L;D=12*EI/L^3;G=EA/L;
%%% monitor cpu time
starttime = cputime;
%% Gradually increase the applied load
for i=1:m   
it = 0; % initialize iteration count
residual = 100; % initialize residual
% Define the applied load
 F = [P3;P4;P5*i;P6];
 Ma = F(3)*L;
 if abs(Ma) >= Mu
        disp('      ## Moment reached to Ultimate Moment ##')
        break
        end
while (residual > tolerance)
        % assemble global K matrix
        k = ((k2-k1)/tu)*u(1)/2 +k1;
        K=[A+k B*lanY -B*lanX C; % evaluate 1 function
           B*lanY (G*lanX^2+D*lanY^2) (G-D)*lanX*lanY B*lanY;% evaluate 2 function
          -B*lanX (G-D)*lanX*lanY (G*lanY^2+D*lanX^2) -B*lanX;% evaluate 3 function
           C B*lanY -B*lanX A]; % evaluate 4 function
        f=K*u-F;       
        %Jacobian matrix 
        J=[A+((k2-k1)/tu)*u(1)+k1 B*lanY -B*lanX C;
           B*lanY (G*lanX^2+D*lanY^2) (G-D)*lanX*lanY B*lanY;
           -B*lanX (G-D)*lanX*lanY (G*lanY^2+D*lanX^2) -B*lanX;
           C B*lanY -B*lanX A];
        %calculate du1 & du2 & du3 & du4
        du = J^-1 *(-f);
        %Calculate the residual (internal-external force deviation)
        residual = max(abs(du)); % evaluate residual
        it = it + 1; % increment iteration count
        if it == itermax
          fprintf('(-)For increment %1.0f trail iteration reached to Ultimate %1.0f\n',i,it)
             disp('    ## The solution for this step is not converged ##')
            break
        end
        u = u+du; % update u
end
Mi(i) =[A B*lanY -B*lanX C]*u;
   if it < itermax
   fprintf('(+)It is converged in %1.0f iterations for increment %1.0f\n',it,i)
   end
        if abs(u(1)) >= tu
        disp('      ## spring reached to Ultimate Rotation ##')
        break
        end
% Force and Dispalcement for each increment
F3(i) = F(3);
    U1(i) = u(1);
    U2(i) = u(2);
    U3(i) = u(3);
    U4(i) = u(4);
end 
disp('=====================================================================');
disp('rotation(D3)   X-displacement(D4)   Y-displacement(D5)   rotation(D6)');
disp('---------------------------------------------------------------------')
disp([U1' U2' U3' U4'])
disp('=====================================================================');
D1=[0;U3'];F1=[0;F3'];
%% Linear Analysis
for ii=1:i 
% calculate Force
 F = [P4;P5*ii;P6];
% assemble global K matrix
 K=[(G*lanX^2+D*lanY^2) (G-D)*lanX*lanY B*lanY;% evaluate 2 function
    (G-D)*lanX*lanY (G*lanY^2+D*lanX^2) -B*lanX;% evaluate 3 function
     B*lanY -B*lanX A];% evaluate 4 function
    % evaluate 4 function
        U = K^-1 *F;
% Force and Dispalcement for each increment
F2(ii) = F(2);
U1(ii) = U(1);
U2(ii) = U(2);
U3(ii) = U(3);
end
D2=[0;U2'];F2=[0;F2'];
%% imaging
figure (1)
IMAGE=imread('NonlinSemiRigidConectionSPRINGlinBEAM.jpg');
image(IMAGE);axis image;axis off;
figure(2)
p1=plot(D1,F1,'r',D2,F2);grid on;set(p1,'LineWidth',3);legend('Non-Linear', 'Linear','Location','NorthEastOutside');
xlabel('Deflection (mm)');ylabel('Force (kN)');
title('Deflection diagram of nonlinear rotational spring with linear element ','color','b');
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s)= %7.4f \n',totaltime)