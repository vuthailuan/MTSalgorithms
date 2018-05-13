function [u_n] = expRK32s3(An, gn, mname, u0, m, tvals,c_2, h)
% usage: [u_n] = expRK32s3(An, g, mname, u0, m, T,c, h)
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
%     c_2    = free variable but equal to 2/3
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

% Set problem parameters
c_3      = 2/3;
n        = 0;
u_n      = u0;
t_n      = tvals(1);
ONEMSM   = 1-sqrt(eps);  

% Info from butcher table for solve_ERKfast()
B = butcher(mname); 

while t_n < tvals(2)*ONEMSM
    
    % Set initial condition
    Y0 = u_n;
    
    % Determine h_fast
    h_fast = c_2*h/ceil(c_2*m);
    
    % Find  U_n2
    p_n2 = feval(gn,t_n,u_n);
    fcn = @(t,y) feval(An)*y + p_n2;
    [~,Y,~] = solve_ERKfast(fcn,[0,c_2*h],Y0,B,h_fast);
    U_n2 = Y(:,2); 
    
    % Define intermediate function
    D_ni   = @(t,c,U) feval(gn,t+c*h,U) - p_n2;
    
    % Determine h_fast
    h_fast = c_3*h/ceil(c_3*m);
    
    % Find U_{n,3}
    D_n2 = D_ni(t_n,c_2,U_n2);
    p_n3 = @(t) p_n2 + 4*t/(9*c_2*c_3^2*h)*D_n2;
    fcn = @(t,y) feval(An)*y + p_n3(t);
    [~,Y,~] = solve_ERKfast(fcn,[0,c_3*h],Y0,B,h_fast);
    U_n3 = Y(:,2);
    
    % Find q_{n,3}(tau) as in 3.19b
    D_n3 = D_ni(t_n,c_3,U_n3);
    q_n3 = @(t) p_n2 + (t/h)*(-2/(3*c_2*(c_2-c_3)))*D_n2 + ...
        (t/h)*(c_2/(c_3*(c_2-c_3)))*D_n3 + ...
        (t^2/(2*h^2))*(2/(c_2*(c_2-c_3)))*D_n2 - ...
        (t^2/(2*h^2))*(2/(c_3*(c_2-c_3)))*D_n3;
   
    % Determine h_fast
    h_fast = h/m;
    
    % Set right hand side for u_{n+1} ODE
    fcn = @(t,y) feval(An)*y + q_n3(t);
    
    % Solve for u_{n+1} on [0,h] 
    [~,Y,~] = solve_ERKfast(fcn,[0,h],Y0,B,h_fast);
    u_np1 = Y(:,2);
    
    % Update time step
    t_n = t_n + h;
    n = n+1;
    
    % Update u value
    u_n = u_np1;   

end
end