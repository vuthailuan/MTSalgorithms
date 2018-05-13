% driver for stiff brusselator test problem:
%      u' = a - (w+1)*u + u^2*v,
%      v' = w*u - u^2*v,
%      w' = (b-w)/ep - u*w,
% where u(0) = 1.2, v(0) = 3.1 and w(0) = 3, with prameters a=1,
% b=3.5 and ep=1e-3.  We evaluate over the time interval [0,10].  
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% Spring 2018
% Adjusted from the code by D. Reynolds ("driver_p2.m", Aug 2012)

clear
close all

% set problem parameters
fn = 'f_p2';
Jn = 'J_p2';
mname = 'Cooper6-ERK';
Tf = 10;
n  = 101;
tout = linspace(0,Tf,n);
hmin = 1e-7;
hmax = 1.0;
h    = 1e-2*[1,.5,.25,0.125,.0625,0.03125];
rtol = 1e-3;
atol = 1e-15*ones(3,1);
global Pdata;
Pdata.a = 1; 
Pdata.b = 3.5; 
Pdata.ep = 1e-2;
u0 = 1.2;
v0 = 3.1;
w0 = 3;
Y0 = [u0; v0; w0];
Y  = zeros(3,n);
Y(:,1) = Y0;
ns = zeros(1,n);
err_max = zeros(1,length(h));
err_rms = zeros(1,length(h));

% plot "true" solution 
opts = odeset('RelTol',1e-13, 'AbsTol',atol,'InitialStep',hmin/10, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
figure()
plot(t,Ytrue)
xlabel('t','FontSize',12), ylabel('y','FontSize',12)
title('Brusselator ODE test','FontSize',14)
set(gca,'FontSize',12)
print('-dpng','fastbrusselator')


% run with an explicit RK method
B = butcher(mname);  s = numel(B(1,:))-1;
fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
for j = 1:length(h)
    tic
    for i = 2:n
     Y0 = Y(:,i-1);
     [t,Yout,ns(i)] = solve_ERKfast(fn,[tout(i-1),tout(i)], Y0,B,h(j));
     Y(:,i) = Yout(:,2);
    end
    err_max(j) = max(max(abs(Y'-Ytrue)));
    err_rms(j) = sqrt(sum(sum((Y'-Ytrue).^2))/numel(Y));
    fprintf('Accuracy/Work Results, h = %g:\n',h(j))
    fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max(j), err_rms(j));
    toc
end 

% Convergence plot
loglog(h,err_max,'b',h,err_rms,'r','LineWidth',1.2);
title('Error','Fontsize',14),xlabel('h','FontSize',12),ylabel('Error','FontSize',12)
legend('absolute','rms','Location','Best')
%print('-dpng',mname)

% Rate of convergence
max_rate = zeros(1,length(h)-1);
rms_rate = zeros(1,length(h)-1);

for j = 2:length(h)
    max_rate(j-1) = log(err_max(j-1)/err_max(j))/log(h(j-1)/h(j));
    rms_rate(j-1) = log(err_rms(j-1)/err_rms(j))/log(h(j-1)/h(j));
end
max_rate
rms_rate
% end of script
