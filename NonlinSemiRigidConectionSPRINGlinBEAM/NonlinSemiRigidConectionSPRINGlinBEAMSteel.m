%***********************************************************%
%                >> IN THE NAME OF GOD <<                   %
% Analysis of 1st order Linear and Nonlinear Semi Rigid     %
% Conection beam subjected to Pushover lateral load         %
% Checking the analysis by Force Control                    %
% Newton-Raphson Method : Initial stiffness procedure       %
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
clear all;close all;clc;
% Define Parameters in mm,kN
P4=0; % [kN]
P5=.05; % [kN]
P6=0; % [kN.mm]
L = 4000; % [mm]
EI= 200*100^4/12; % [kN.mm^2]
EA = 200*10000; % [kN]
m = 150; % number of calculation
itermax = 200;% maximum number of iterations
tolerance = 1e-6; % specified tolerance for convergence
n = 2; % Momen-rotation shape parameter
u = zeros(4,1);% initial guess values
lanX=1;lanY=0;
%%% monitor cpu time
starttime = cputime;
% Element stifness coefficient
A=4*EI/L;B=6*EI/L^2;C=2*EI/L;D=12*EI/L^3;G=EA/L;
% Nonlinear Rotational Spring
t0=.05; % Yeild rotaion
Mu=10e+3; % Yeild moment
tu=.1;  % ultimate rotation
Rki=Mu/t0;

% Gradually increase the applied load
for i=1:m     
it = 0; % initialize iteration count
residual = 100; % initialize residual
% calculate Force
 F = [0;P4;P5*i;P6];
 a1=(Rki)/((1+(abs((Rki*u(1))/Mu))^n)^(1/n));
%% First-order Nonlinear Analysis
 Kinit=[A+a1 B*lanY -B*lanX C;
           B*lanY (G*lanX^2+D*lanY^2) (G-D)*lanX*lanY B*lanY;
          -B*lanX (G-D)*lanX*lanY (G*lanY^2+D*lanX^2) -B*lanX;
           C B*lanY -B*lanX A];
while (residual > tolerance)
        % assemble global K matrix
        K=[A+((Rki)/((1+(abs((Rki*u(1))/Mu))^n)^(1/n))) B*lanY -B*lanX C; 
           B*lanY (G*lanX^2+D*lanY^2) (G-D)*lanX*lanY B*lanY;
          -B*lanX (G-D)*lanX*lanY (G*lanY^2+D*lanX^2) -B*lanX;
           C B*lanY -B*lanX A];
        f=F-K*u;
        %calculate du1 & du2 & du3 & du4
        du = Kinit^-1 *(f);
        %Calculate the residual (internal-external force deviation)
        residual = max(abs(du)); % evaluate residual
        it = it + 1; % increment iteration count
        if it == itermax
          fprintf('(-)For increment %1.0f trail iteration reached to Ultimate %1.0f\n',i,it)
             disp('    ## The solution for this step is not converged ##')
        end
        u = u+du; % update u
end
   if it < itermax
   fprintf('(+)It is converged in %1.0f iterations for increment %1.0f\n',it,i)
   end
        if abs(u(1)) >= tu
        disp('      ## spring reached to Ultimate Rotation ##')
        break
        end
 Mi(i) = [A+((Rki)/((1+(abs((Rki*u(1))/Mu))^n)^(1/n))) B*lanY -B*lanX C]*u;       
% Internal element force           
         % Displacement Transformation Matrix
        T = [lanX lanY 0 0 0 0;
            -lanY lanX 0 0 0 0;
                   0 0 1 0 0 0;
             0 0 0 lanX lanY 0;
            0 0 0 -lanY lanX 0;
                   0 0 0 0 0 1];
        % Stiffness Matrix for each element
        Kele = [G 0 0 -G 0 0;
             0 D B 0 -D B;
             0 B A 0 -B C;
             -G 0 0 G 0 0;
           0 -D -B 0 D -B;
             0 B C 0 -B A]; 
        Fele = Kele*T*[0;0;u(1);u(2);u(3);u(4)];    
% Force and Dispalcement for each increment
    F3(i) = F(3);
    U1(i) = u(1);
    U2(i) = u(2);
    U3(i) = u(3);
    U4(i) = u(4);
    INT_N_f2(i) = roundn(Fele(2),-3);% Interal force of base Shear
    INT_N_f3(i) = roundn(Fele(3),-3);
    INT_N_f4(i) = roundn(Fele(4),-3);
    INT_N_f5(i) = roundn(Fele(5),-3);
    INT_N_f6(i) = roundn(Fele(6),-3);
end 
disp('=====================================================================');
disp('rotation(D3)   X-displacement(D4)   Y-displacement(D5)   rotation(D6)');
disp('---------------------------------------------------------------------')
disp([U1' U2' U3' U4'])
disp('=====================================================================');
D1=[0;U3'];F1=[0;F3'];BS1N=[0;INT_N_f2'];M1=[0;INT_N_f3'];te1=[0;U1'];
%% First-order Linear Analysis
for ii=1:i
% calculate Force
 F = [0;P4;P5*ii;P6];
% assemble global K matrix
        K=[A+Rki*1000000 B*lanY -B*lanX C;
           B*lanY (G*lanX^2+D*lanY^2) (G-D)*lanX*lanY B*lanY
          -B*lanX (G-D)*lanX*lanY (G*lanY^2+D*lanX^2) -B*lanX;
           C B*lanY -B*lanX A];
        u = K^-1 *F;
  Mi(ii) = [A+Rki*1000 B*lanY -B*lanX C]*u;      
% Internal element force          
         % Displacement Transformation Matrix
        T = [lanX lanY 0 0 0 0;
            -lanY lanX 0 0 0 0;
                   0 0 1 0 0 0;
             0 0 0 lanX lanY 0;
            0 0 0 -lanY lanX 0;
                   0 0 0 0 0 1];
        % Stiffness Matrix for each element
        Kele = [G 0 0 -G 0 0;
             0 D B 0 -D B;
             0 B A 0 -B C;
             -G 0 0 G 0 0;
           0 -D -B 0 D -B;
             0 B C 0 -B A]; 
        Fele = Kele*T*[0;0;u(1);u(2);u(3);u(4)];
        
% Force and Dispalcement for each increment
    F3(ii) = F(3);
    U1(ii) = u(1);
    U2(ii) = u(2);
    U3(ii) = u(3);
    U4(ii) = u(4);
    INT_L_f2(ii) = roundn(Fele(2),-3);
    INT_L_f3(ii) = roundn(Fele(3),-3);
    INT_L_f4(ii) = roundn(Fele(4),-3);
    INT_L_f5(ii) = roundn(Fele(5),-3);
    INT_L_f6(ii) = roundn(Fele(6),-3); 
end
D2=[0;U3'];F2=[0;F3'];BS1L=[0;INT_L_f2'];M2=[0;INT_L_f3'];te2=[0;U1'];
disp('+ ========== First-order Nonlinear ========== +');
disp('     (f3)       (f4)     (f5)      (f6)     ');
disp('----------------------------------------------');
disp([INT_N_f3' INT_N_f4' INT_N_f5' INT_N_f6']);
disp('+ ========== First-order Linear ============= +');
disp('     (f3)       (f4)     (f5)      (f6)     ');
disp('----------------------------------------------');
disp([INT_L_f3' INT_L_f4' INT_L_f5' INT_L_f6']); 
disp('==============================================');
%% imaging
figure (1)
IMAGE=imread('NonlinSemiRigidConectionSPRINGlinBEAMSteel.jpg');
image(IMAGE);axis image;axis off;
figure(2)
p1=plot(te1,M1,te2,M2,'r--');grid on;set(p1,'LineWidth',3);legend('Semi-Rigid','Rigid','Location','NorthEastOutside');
xlabel('Rotation (rad)');ylabel('Moment (kN.mm)');
title('Moment-Rotation diagram of Rigid and Semi-Rigid spring with linear element ','color','b');
figure(3)
p2=plot(D1,F1,D2,F2,'r--');grid on;set(p2,'LineWidth',3);legend('Semi-Rigid','Rigid','Location','NorthEastOutside');
xlabel('Displacement (mm)');ylabel('Force (kN)');
title('Force-Displacement diagram of Rigid and Semi-Rigid spring with linear element','color','b');
%%%  print time of computation
totaltime = cputime - starttime;
fprintf('\nTotal time (s)= %7.4f \n',totaltime)