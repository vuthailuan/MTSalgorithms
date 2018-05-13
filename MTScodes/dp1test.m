% driver for Prothero Robinson test problem:
%      u' = gamma*u + epn*v-gamma*cos(t) - epn*cos(w*t) - sin(t),
%      v' = epn*u + v -epn*cos(t) - epn*cos(w*t) - sin(t)-w*sin(w*t),
% where u(0) = 1.0, v(0) = 1.0, with parameters w=100,
% gamma = -1/w and epn = 1/w.  We evaluate over the time interval [0,2]. 
% A multiple time stepping method is used 
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
mname = 'ERK-3-3';                                 % explicit RK method

% Set problem parameters
global Pdata;
Pdata.w     = 100; 
Pdata.gamma = -1/Pdata.w; 
Pdata.epsilon    = 1/Pdata.w;

% Set time parameters for 
Tf   = 2;                                        % end time 
n    = 3;
tout = linspace(0,Tf,n);                          % immediate times for solution
h    = 1e-2*0.5.^(0:5);                           % large time step
herk = 1e-6;
m    = 10;                                        % divisor for smaller time step
c2   = 1/3;




% % Compute solution using MTS solver
% fprintf('Running multiple time stepping integrator: \n');
% fprintf('Coefficients used: m = %i, c2 = %g, (fixed) c3 = %g.\n',...
%     m,c2,2/3);
% B = butcher(mname);  s = numel(B(1,:))-1;
% fprintf('\nRunning with inner ERK integrator: %s (order = %i)\n',mname,B(s+1,1))
% 
% for j = 1:length(h)
%     tstart = tic;
%     for i = 2:n 
%         Y0     = Y(:,i-1);
%         Ytrue(:,i) = feval(yn,tout(i));
%         Yout   = expRK32(An, gn, mname,Y0,m,[tout(i-1),tout(i)],c2, h(j));
%         Y(:,i) = Yout;
%         
%     end
%     telapsed = toc(tstart);
%     time(j) = telapsed;  
%     % Error calculation
%     err_max(j) = max(max(abs(Y- Ytrue)));
%        
%     fprintf('Accuracy/Work Results, h = %g:\n',h(j))
%     fprintf('   maxerr = %.5e,   ',err_max(j));
%     
%     % Rate of convergence calculation
%     if j > 1
%         max_rate(j-1) = log(err_max(j-1)/err_max(j))/log(h(j-1)/h(j));
%         fprintf('   maxrate = %.5e \n',max_rate(j-1));      
%     end
%     fprintf('Time elapsed = %g.\n',telapsed);
%     
%     
% end
% run with an explicit RK method
    Ytrue = feval(yn,tout);
    Y2 = zeros(2,n);
    Y2(:,1) = Ytrue(:,1);
    mname2 = 'Fehlberg-8-7-ERK';
    B2 = butcher(mname2);  s2 = numel(B2(1,:))-1;
    fprintf('\nRunning with ERK integrator: %s (order = %i)\n',mname2,B2(s2+1,1))
    tic
    for i = 2:length(tout)
        Y0 = Y2(:,i-1);
        [t,Yout,ns] = solve_ERKfast(fn,[tout(i-1),tout(i)], Y0, B2,herk);
        Y2(:,i) = Yout(:,2);
    end
    toc
    %keyboard
    err_max = max(max(abs(Y2-Ytrue)));
    err_rms = sqrt(sum(sum((Y2-Ytrue).^2))/numel(Y2));
    fprintf('Accuracy/Work Results:\n')
    fprintf('   maxerr = %.5e,  rmserr = %.5e\n',err_max, err_rms);

% % Convergence plot
% figure
% loglog(h,err_max,'b','LineWidth',2);
% title([mname,' Error'],'Fontsize',14),xlabel('h','FontSize',12),ylabel('Error','FontSize',12)
% legend('absolute','Location','Best')
% print('-dpng',['p1_expRK32 Error (',mname,')'])
% 
% % Time plot
% figure
% loglog(time,err_max,'b','LineWidth',2);
% title([mname,' Time Elapsed'],'Fontsize',14),xlabel('h','Fontsize',12),ylabel('Time','Fontsize',12)
% print('-dpng',['p1_expRK32 efficiency (',mname,')'])