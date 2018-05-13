function [u_n] = expRK32(An, gn, mname, u0, m, tvals,c_2, h)
% usage: [u_n] = expRK32(An, g, mname, u0, m, T,c_2, h)
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
%     c_2    = free variable not 2/3
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
n        = 0;
u_n      = u0;
t_n      = tvals(1);
ONEMSM   = 1-sqrt(eps);  

% Info from butcher table for solve_ERKfast()
B = butcher(mname); 



while t_n < tvals(2)*ONEMSM
    
    % Set initial condition
    Y0 = u_n;
    
    % Find p_{n,2}(tau) as in 3.19a
    p_n2 = feval(gn,t_n,u_n);
    
    % Right hand side for U_{n,2} ODE
    fcn = @(t,y) feval(An)*y + p_n2;
        
    % Determine h_fast
    h_fast = c_2*h/ceil(c_2*m);
    
    % Solve for U_n2 using solve_ERKfast (3.20)
    [~,Y,~] = solve_ERKfast(fcn,[0,c_2*h],Y0,B,h_fast);
    U_n2 = Y(:,2);
    
    % Find p_{n,3}(tau) using U_{n,1} and U_{n,2}
    p_n3 = @(t) 1/(c_2*h)*t*feval(gn,t_n + c_2*h,U_n2) + (1-t/(c_2*h))*p_n2;
      
    % Right hand side for U_{n,3} ODE
     fcn = @(t,y) feval(An)*y + p_n3(t);
    
    % Determine h_fast
    h_fast = (2/3*h)/ceil(2/3*m);
    
    % Solve for U_{n,3} using solve_ERKfast (3.20)
    [~,Y,~] = solve_ERKfast(fcn,[0,2*h/3],Y0,B,h_fast);
    U_n3 = Y(:,2);
    
    % Find q_{n,3}(tau) as in 3.19b
    q_n3 = @(t) 3*t/(2*h)*feval(gn,t_n + 2*h/3,U_n3) + (1 - 3*t/(2*h))*p_n2;
    
    % Set right hand side for u_{n+1} ODE
    fcn = @(t,y) feval(An)*y + q_n3(t);   
    
    % Determine h_fast
    h_fast = h/m;
    
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