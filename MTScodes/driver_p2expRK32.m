% driver for stiff brusselator test problem (expRK32):
%      u' = a - (w+1)*u + u^2*v,
%      v' = w*u - u^2*v,
%      w' = (b-w)/ep - u*w,
% where u(0) = 1.2, v(0) = 3.1 and w(0) = 3, with parameters a=1,
% b=3.5 and ep=1e-2.  We evaluate over the time interval [0,10]. 
% A multiple time stepping method is used 
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% Spring 2018

clear
close all

% Set functional parameters
fn    = 'f_p2';                                    % brusselator function
Jn    = 'J_p2';
gn    = 'g_p2';                                    % non-stiff part 
An    = 'A_p2';                                    % stiff part
mname = 'ERK-3-3';                                 % explicit RK method

% Set problem parameters
global Pdata;
Pdata.a = 1; 
Pdata.b = 3.5; 
Pdata.ep = 1e-2;

% Set time parameters for 
Tf   = 2;                                        % end time 
n    = 3;
tout = linspace(0,Tf,n);                          % immediate times for solution
h    = 1e-1*0.5.^(0:5);                           % large time step
m    = 20;                                        % divisor for smaller time step
c2   = 1/3;

% Parameters for built-in MATLAB ode solver
hmin = 1e-7;                                                                
hmax = 1.0;
rtol = 1e-3;
atol = 1e-15*ones(3,1);

% Initialize problem variables/ allocate space
u0       = 1.2;
v0       = 3.1;
w0       = 3;
Y        = zeros(3,n);
Y(:,1)   = [u0; v0; w0];
ns       = zeros(1,n);
err_max  = zeros(1,length(h));
err_rms  = zeros(1,length(h));
max_rate = zeros(1,length(h)-1);
rms_rate = zeros(1,length(h)-1);
time     = zeros(1,length(h));

% Compute "true" solution using MATLAB solver
fprintf('Running MATLAB ode15s: \n');
tic
opts = odeset('RelTol',1e-13, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y(:,1), opts);
toc

% Compute solution using MTS solver
fprintf('Running multiple time stepping integrator: \n');
fprintf('Coefficients used: m = %i, c2 = %g, (fixed) c3 = %g.\n',...
    m,c2,2/3);
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with inner ERK integrator: %s (order = %i)\n',mname,B(s+1,1))

for j = 1:length(h)
    tstart = tic;
    for i = 2:n 
        Y0     = Y(:,i-1);
        Yout   = expRK32(An, gn, mname,Y0,m,[tout(i-1),tout(i)],c2, h(j));
        Y(:,i) = Yout;
    end
    telapsed = toc(tstart);
    time(j) = telapsed;
    
    % Error calculation
    err_max(j) = max(max(abs(Y'- Ytrue)));
       
    fprintf('Accuracy/Work Results, h = %g:\n',h(j))
    fprintf('   maxerr = %.5e,   ',err_max(j));
    
    % Rate of convergence calculation
    if j > 1
        max_rate(j-1) = log(err_max(j-1)/err_max(j))/log(h(j-1)/h(j));
        fprintf('   maxrate = %.5e \n',max_rate(j-1));      
    end
    fprintf('Time elapsed = %g.\n',telapsed);
end

% Convergence plot
figure
loglog(h,err_max,'b','LineWidth',2);
title([mname,' Error'],'Fontsize',14),xlabel('h','FontSize',12),ylabel('Error','FontSize',12)
legend('absolute','Location','Best')
print('-dpng',['expRK32Error(',mname,')'])

% Time plot
figure
loglog(time,err_max,'b','LineWidth',2);
title([mname,' Time Elapsed'],'Fontsize',14),xlabel('time','Fontsize',12),ylabel('Error','Fontsize',12)
print('-dpng',['expRK32efficiency(',mname,')'])

%     Zall = [[Z_oil_xmin,Z_oilCurrent',Z_oilCurrent(N)];[Z_gas_xmin,Z_gasCurrent',Z_gasCurrent(N)];...
%     [Z_wat_xmin,Z_watCurrent',Z_watCurrent(N)]];
%     filename = ['ZResultscase',num2str(testcase),num2str(current_time),num2str(N),'.txt'];
%     fileID = fopen(filename,'w');
%     formatSpec1 = '%16s %16s %16s\n';
%     formatSpec2 = '%16.16f %16.16f %16.16f\n';
%     fprintf(fileID,formatSpec1,'Z_oil','Z_gas','Z_wat');
%     fprintf(fileID,formatSpec2,Zall);
%     fclose(fileID);

% end of script
