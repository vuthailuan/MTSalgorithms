function [u_n] = expRK43s6(An, gn, mname, u0, m, tvals,c, h)
% usage: [u_n] = expRK43s6(An, g, mname, u0, m, T,c, h)
%
% Multiple time stepping solver for the vector-valued ODE problem 
%     u' = Au(t_ + g(t,u(t)), t in tvals, y in R^m,
%     u(t0) = [u1(t0), u2(t0), ..., um(t0)]'.
%
% Inputs:
%     An     = string holding function name for A
%     gn     = string holding function name for g(t,u(t))
%     mname  = string holding function name for explicit RK method
%              Butcher table formatting 
%                 B = [c A;
%                      q b]
%                 Here, c is a vector of stage time fractions (s-by-1),
%                    A is a matrix of Butcher coefficients (s-by-s),
%                    q is an integer denoting the method order of accuracy,
%                    b is a vector of solution weights (1-by-s),
%     tvals  = [t0,tN] initial and final time
%     u0     = initial value array (column vector of length m)
%     c      = c values used
%     h      = large time step
%
% Outputs: 
%     u_n    = column vector with solution at final time
%
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% Spring 2018

global Pdata;
Pdata.ffast;
Pdata.fslow;
% Set problem parameters
c_2      = c(1);
c_3      = c(2);
c_4      = c(3);
c_5      = c(4);
c_6      = c(5);
n        = 0;
u_n      = u0;
t_n      = tvals(1);
ONEMSM   = 1-sqrt(eps);  

% Info from butcher table for solve_ERKfast()
B = butcher(mname); 


while t_n < tvals(2)*ONEMSM
    
    % Set initial condition
    Y0 = u_n;
    
    % Compute U_{n,2}
    p_n2 = feval(gn,t_n,u_n);
    fcn = @(t,y) feval(An)*y + p_n2;
    h_fast = c_2*h/ceil(c_2*m);
    [~,Y,nflocal] = solve_ERKfast(fcn,[0,c_2*h],Y0,B,h_fast);
    Pdata.ffast   = Pdata.ffast + nflocal + 1;     % one from p_n2 evaluation?
    U_n2 = Y(:,2);
    
    % Define intermediate function
    D_ni   = @(t,c,U) feval(gn,t + c*h,U) - p_n2;
    
    % Intermediate function D_n2
    D_n2 = D_ni(t_n,c_2,U_n2);
    
    % p_{n,3} and p_{n,4} are the same, use same polynomial
    p_n3 = @(t) p_n2 + (t/(h*c_2))*D_n2;
    fcn = @(t,y) feval(An)*y + p_n3(t);
    
    % Determine time step, and times to evaluate at 
    h_fast = c_3*h/ceil(c_3*m);
    tsteps = 0:h_fast:c_3*h;
    
    % Determine the largest index so that c_4*h < tsteps(index)
    index  = floor(c_4*h/h_fast);
    
    % Set an intermediate Y value
    [~,Y,nflocal] = solve_ERKfast(fcn,[0,tsteps(index)],Y0,B,h_fast);
    Yim = Y(:,2);
    Pdata.ffast   = Pdata.ffast + nflocal + 1;     % one from D_n2 evaluation?
    
    % Determine step leftover 
    h_leftover = c_4*h - tsteps(index);
    
    % Solve for U_{n,4}
    [~,Y,nflocal] = solve_ERKfast(fcn,[tsteps(index),c_4*h],Yim,B,h_leftover);
    U_n4 = Y(:,2);
    Pdata.ffast   = Pdata.ffast + nflocal; 
    
    %Solve for U_{n,3}
    [~,Y,nflocal] = solve_ERKfast(fcn,[tsteps(index),c_3*h],Yim,B,h_fast);
    U_n3 = Y(:,2);
    Pdata.ffast   = Pdata.ffast + nflocal;
    
    % Intermediate functions D_n3 and D_n4
    D_n3 = D_ni(t_n,c_3,U_n3);
    D_n4 = D_ni(t_n,c_4,U_n4);
    
    % p_{n,5} and p_{n,5} are the same, use same polynomial
    p_n5 = @(t) p_n2 + (t/h)*(-c_4/(c_3*(c_3-c_4)))*D_n3 +...
        (t/h)*(c_3/(c_4*(c_3-c_4)))*D_n4 +...
        (t^2/(2*h^2))*(2/(c_3*(c_3-c_4)))*D_n3 - ...
        (t^2/(2*h^2))*(2/(c_4*(c_3-c_4)))*D_n4;
    fcn = @(t,y) feval(An)*y + p_n5(t);
    
    % Time step, times to evaluate at
    h_fast = c_5*h/ceil(c_5*m);
    tsteps = 0:h_fast:c_5*h;

    % Determine the largest index so that c_6*h < tsteps(index)
    index  = floor(c_6*h/h_fast);
    
    % Set an intermediate Y value
    [~,Y,nflocal] = solve_ERKfast(fcn,[0,tsteps(index)],Y0,B,h_fast);
    Yim = Y(:,2);
    Pdata.ffast   = Pdata.ffast + nflocal + 2;     % two from D_n3/D_n4 evaluation?

    % Determine step leftover 
    h_leftover = c_6*h - tsteps(index);
    
    % Solve for U_{n,6}
    [~,Y,nflocal] = solve_ERKfast(fcn,[tsteps(index),c_6*h],Yim,B,h_leftover);
    U_n6 = Y(:,2);
    Pdata.ffast   = Pdata.ffast + nflocal;
    
    %Solve for U_{n,5}
    [~,Y,nflocal] = solve_ERKfast(fcn,[tsteps(index),c_5*h],Yim,B,h_fast);
    U_n5 = Y(:,2);
    Pdata.ffast   = Pdata.ffast + nflocal;
    
    % Intermediate functions D_n5 and D_n6
    D_n5 = D_ni(t_n,c_5,U_n5);
    D_n6 = D_ni(t_n,c_6,U_n6);
    
    % Compute u_{n+1}
    q_n6 = @(t) p_n2 + (t/h)*(-c_6/(c_5*(c_5-c_6)))*D_n5 +...
        (t/h)*(c_5/(c_6*(c_5-c_6)))*D_n6 +...
        (t^2/(2*h^2))*(2/(c_5*(c_5-c_6)))*D_n5 - ...
        (t^2/(2*h^2))*(2/(c_6*(c_5-c_6)))*D_n6; 

%     % Alternative q
%     q_n6 = @(t) p_n2 + (t/h)*(-c_4/(c_3*(c_3-c_4)))*D_n3 +...
%         (t/h)*(c_3/(c_4*(c_3-c_4)))*D_n4 +...
%         (t^2/(2*h^2))*(2/(c_3*(c_3-c_4)))*D_n3 - ...
%         (t^2/(2*h^2))*(2/(c_4*(c_3-c_4)))*D_n4;
    
    fcn = @(t,y) feval(An)*y + q_n6(t);
    h_fast = h/m;
    [~,Y,nflocal] = solve_ERKfast(fcn,[0,h],Y0,B,h_fast);
    u_np1 = Y(:,2);
    Pdata.fslow = Pdata.fslow + nflocal + 2; 
    
    % Update time step
    t_n = t_n + h;
    n = n+1;
    
    % Update u value
    u_n = u_np1;   

end
end