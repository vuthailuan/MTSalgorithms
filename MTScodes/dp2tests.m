% driver for stiff brusselator test problem (expRK4356): (test driver) 
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
gn    = 'g_p2';                                    % non-stiff part 
An    = 'A_p2';                                    % stiff part
Jn    = 'J_p2';
mname = 'ERK-4-4';                                 % explicit RK method

% Set problem parameters
global Pdata;
Pdata.a = 1; 
Pdata.b = 3.5; 
Pdata.ep = 1e-2;
Pdata.ffast = 0;
Pdata.fslow = 0;

% Set time parameters for 
Tf   = 2;                                         % end time 
n    = 3;
tout = linspace(0,Tf,n);                          % immediate times for solution
h    = 1e-1*0.5.^(0:5);                           % large time step
m    = 20;                                        % divisor for smaller time step
%c    = [1/2,1/2,1/3,1/2,1];
c    = [1/2,1/2,1/3,5/6,1/3];

% Parameters for high order ERK
h_erk = 1e-4;

% Parameters for high order IRK
hmin = 1e-7;
hmax = 1.0; 
rtol = 1e-14;
atol = 1e-14*ones(3,1);
err_irk    = zeros(1,length(h));
irk_rate   = zeros(1,length(h)-1);

% Initialize problem variables/ allocate space
u0         = 1.2;
v0         = 3.1;
w0         = 3;
Y          = zeros(3,n);
Ytrue      = zeros(3,n);
Y(:,1)     = [u0; v0; w0];
Ytrue(:,1) = [u0;v0;w0];
nffast     = zeros(1,length(h));
nfslow     = zeros(1,length(h));
err_max    = zeros(1,length(h));
max_rate   = zeros(1,length(h)-1);
time       = zeros(1,length(h));

% % Compute "true" solution with a high order explicit RK method
% highERK = 'Fehlberg-8-7-ERK';
% D = butcher(highERK);  s = numel(D(1,:))-1;
% fprintf('\nRunning with high order ERK integrator: %s (order = %i)\n',highERK,D(s+1,1))
% 
% tic
% for i = 2:length(tout)
%     Y0 = Ytrue(:,i-1);
%     [t,Yout,nflocal] = solve_ERKfast(fn,[tout(i-1),tout(i)], Y0, D,h_erk);
%     Ytrue(:,i) = Yout(:,2);
% end
% toc

% Compute "true" solutionwith a fully-implicit RK method
highirk = 'Gauss-6-12-IRK';
D = butcher(highirk);  s = numel(D(1,:))-1;
fprintf('\nRunning with IRK integrator: %s (order = %i)\n',highirk,D(s+1,1))
tic
[t,Yirk,ns,nl] = solve_IRK(fn, Jn, tout, [u0;v0;w0], D, rtol, atol, hmin, hmax);
toc

% Compute solution using MTS solver
fprintf('Running multiple time stepping integrator: \n');
fprintf('Coefficients used: m = %i, c2 = %g, c3 = %g, c4 = %g, c5 = %g, c6 = %g. \n',...
    m,c(1),c(2),c(3),c(4),c(5));
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with inner ERK integrator: %s (order = %i)\n',mname,B(s+1,1))

for j = 1:length(h)
    tstart = tic;
    for i = 2:n 
        Y0     = Y(:,i-1);
        Yout   = expRK43s6(An, gn, mname,Y0,m,[tout(i-1),tout(i)],c, h(j));
        Y(:,i) = Yout;
    end
    telapsed = toc(tstart);
    time(j) = telapsed;
    
    % Function calls 
    nffast(j) = Pdata.ffast;
    nfslow(j) = Pdata.fslow;
    fprintf('nffast = %g ,  nfslow = %g .\n',nffast(j),nfslow(j));
    
    % Error calculation
%     err_max(j) = max(max(abs(Y - Ytrue)));
    err_irk(j) = max(max(abs(Y - Yirk)));

    fprintf('Accuracy/Work Results, h = %g:\n',h(j))
%     fprintf('   maxerr = %.5e,    ',err_max(j));
    fprintf('   irkerr = %.5e,    ',err_irk(j));

    % Rate of convergence calculation
    if j > 1
%         max_rate(j-1) = log(err_max(j-1)/err_max(j))/log(h(j-1)/h(j));
        irk_rate(j-1) = log(err_irk(j-1)/err_irk(j))/log(h(j-1)/h(j));
%         fprintf('   maxrate = %.5e \n',max_rate(j-1));  
        fprintf('   irkrate = %.5e \n',irk_rate(j-1));
    end
    fprintf('Time elapsed = %g.\n',telapsed);
    
        
    % Reset global variables
    Pdata.ffast = 0;
    Pdata.fslow = 0;
end

% Convergence plot
figure
loglog(h,err_irk,'b','LineWidth',1.5);
title([mname,' Error'],'Fontsize',14),xlabel('h','FontSize',12),ylabel('Error','FontSize',12)
legend('absolute','Location','Best')
print('-dpng',['p2_expRK43s6Error(',mname,')'])

% Time plot
figure
loglog(time,err_irk,'b','LineWidth',1.5);
title([mname,' Time Elapsed'],'Fontsize',14),xlabel('time','Fontsize',12),ylabel('Error','Fontsize',12)
print('-dpng',['p2_expRK43s6efficiency (',mname,')'])

% Function evaluation plot
figure
loglog(nffast+nfslow,err_irk,'b','LineWidth',1.5);
title([mname,' Function Calls'],'Fontsize',14),xlabel('nf','Fontsize',12),ylabel('Error','Fontsize',12)
print('-dpng',['p2_expRK43s6fcalls(',mname,')'])

% % Store results
% filename = ['p2expRK43s6m',num2str(m),'.txt'];
% fileID   = fopen(filename,'w');
% %formatSpec1  =       %for h values 
% %formatSpec2  =       %for errorvalues
% %formatSpec3  =       %for time values
% %formatSpec4  =       %for function evaluation totals


% end of script
