function driver2tests(mname,m)
% driver for stiff brusselator test problem:
%      u' = a - (w+1)*u + u^2*v,
%      v' = w*u - u^2*v,
%      w' = (b-w)/ep - u*w,
% where u(0) = 1.2, v(0) = 3.1 and w(0) = 3, with parameters a=1,
% b=3.5 and ep=1e-2.  We evaluate over the time interval [0,10]. 
% A multiple time stepping method is used.
% 
% INPUTS: 
% mname - name of Explicit RK method
% m     - divisor for smaller time step
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% Spring 2018

% Set functional parameters
fn    = 'f_p2';                                    % brusselator function
gn    = 'g_p2';                                    % non-stiff part 
An    = 'A_p2';                                    % stiff part
                                        % explicit RK method

% Set problem parameters
global Pdata;
Pdata.a = 1; 
Pdata.b = 3.5; 
Pdata.ep = 1e-2;

% Set time parameters for 
Tf   = 10;                                        % end time 
n    = 101;
tout = linspace(0,Tf,n);                          % immediate times for solution
h    = 1e-2*[1,0.5,0.25,0.125];                           % large time step

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

% Compute "true" solution using MATLAB solver
fprintf('Running MATLAB ode15s: \n');
tic
opts = odeset('RelTol',1e-13, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y(:,1), opts);
toc

% Compute solution using MTS solver
fprintf('Running multiple time stepping integrator: \n');
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with inner ERK integrator: %s (order = %i)\n',mname,B(s+1,1))

for mval = m
    fprintf('m value = %i\n',mval);
    for j = 1:length(h)
        tic
        for i = 2:n 
            Y0     = Y(:,i-1);
            Yout   = solve_MTS_ETD(An, gn, mname,Y0,mval,[tout(i-1),tout(i)], h(j));
            Y(:,i) = Yout;
        end
        toc
        % Error calculation
        err_max(j) = max(max(abs(Y'- Ytrue)));
        err_rms(j) = sqrt(sum(sum((Y'- Ytrue).^2))/numel(Y));
       
        fprintf('Accuracy/Work Results, h = %g:\n',h(j))
        fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max(j), err_rms(j));
    
        % Rate of convergence calculation
        if j > 1
            max_rate(j-1) = log(err_max(j-1)/err_max(j))/log(h(j-1)/h(j));
            rms_rate(j-1) = log(err_rms(j-1)/err_rms(j))/log(h(j-1)/h(j));
            fprintf('   maxrate = %.5e,  rmsrate = %.5e\n',max_rate(j-1), rms_rate(j-1));      
        end
    end
end

% Convergence plot
% loglog(h,err_max,'b',h,err_rms,'r','LineWidth',1.2);
% title('Error','Fontsize',14),xlabel('h','FontSize',12),ylabel('Error','FontSize',12)
% legend('absolute','rms','Location','Best')
% print('-dpng',mname)

% end of script
end