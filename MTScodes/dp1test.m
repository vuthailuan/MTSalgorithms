% driver for Prothero Robinson test problem:
%      u' = gamma*u + epn*v-gamma*cos(t) - epn*cos(w*t) - sin(t),
%      v' = epn*u + v -epn*cos(t) - epn*cos(w*t) - sin(t)-w*sin(w*t),
% where u(0) = 1.0, v(0) = 1.0, with parameters w=100,
% gamma = -1/w and epn = 1/w.  We evaluate over the time interval [0,2]. 
% A multiple time stepping method is used 
%
% Here the goal is to determine a good time step for a high order ERK
% method that results in high accuracy 
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% Spring 2018

clear
close all

% Set functional parameters
fn    = 'f_p1';                                    % brusselator function
Jn    = 'J_p1';
gn    = 'g_p1';                                    % non-stiff part 
An    = 'A_p1';                                    % stiff part
yn    = 'y_p1';


% Set problem parameters
global Pdata;
Pdata.w     = 100; 
Pdata.gamma = -1/Pdata.w; 
Pdata.epsilon    = 1/Pdata.w;

% Set time parameters for 
Tf   = 2;                                        % end time 
n    = 3;
tout = linspace(0,Tf,n);                          % immediate times for solution
h    = 1e-2;

% Initialize space and true solution
Ytrue = feval(yn,tout);
Y2 = zeros(2,n);
Y2(:,1) = Ytrue(:,1);

% run with a high order explicit RK method
mname = 'Fehlberg-8-7-ERK';
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
tic
for i = 2:length(tout)
    Y0 = Y2(:,i-1);
    [t,Yout,ns] = solve_ERKfast(fn,[tout(i-1),tout(i)], Y0, B,h);
    Y2(:,i) = Yout(:,2);
end
toc
% Error calculation 
err_max = max(max(abs(Y2-Ytrue)));
err_rms = sqrt(sum(sum((Y2-Ytrue).^2))/numel(Y2));
fprintf('Accuracy/Work Results:\n')
fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);

