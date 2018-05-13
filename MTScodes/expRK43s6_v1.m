function [u_n] = expRK43s6_v1(An, gn, mname, u0, m, tvals,c, h)
% usage: [u_n] = expRK43s6_v1(An, g, mname, u0, m, T,c, h)
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
    [~,Y,~] = solve_ERKfast(fcn,[0,c_2*h],Y0,B,h_fast);
    U_n2 = Y(:,2);
    
    % Define intermediate function
    D_ni   = @(t,c,U) feval(gn,t + c*h,U) - p_n2;
    
    % Intermediate function D_n2
    D_n2 = D_ni(t_n,c_2,U_n2);
    
    % Compute U_{n,3} 
    p_n3 = @(t) p_n2 + (t/(h*c_2))*D_n2;
    fcn = @(t,y) feval(An)*y + p_n3(t);
    h_fast = c_3*h/ceil(c_3*m);
    [~,Y,~] = solve_ERKfast(fcn,[0,c_3*h],Y0,B,h_fast);
    U_n3 = Y(:,2);
    
    % Compute U_{n,4}
    p_n4 = @(t) p_n2 + (t/(h*c_2))*D_n2;
    fcn = @(t,y) feval(An)*y + p_n4(t);
    h_fast = c_4*h/ceil(c_4*m);
    [~,Y,~] = solve_ERKfast(fcn,[0,c_4*h],Y0,B,h_fast);
    U_n4 = Y(:,2);
    
    % Intermediate functions D_n3 and D_n4
    D_n3 = D_ni(t_n,c_3,U_n3);
    D_n4 = D_ni(t_n,c_4,U_n4);
    
    % Compute U_{n,5}
    p_n5 = @(t) p_n2 + (t/h)*(-c_4/(c_3*(c_3-c_4)))*D_n3 +...
        (t/h)*(c_3/(c_4*(c_3-c_4)))*D_n4 +...
        (t^2/(2*h^2))*(2/(c_3*(c_3-c_4)))*D_n3 - ...
        (t^2/(2*h^2))*(2/(c_4*(c_3-c_4)))*D_n4;
    fcn = @(t,y) feval(An)*y + p_n5(t);
    h_fast = c_5*h/ceil(c_5*m);
    [~,Y,~] = solve_ERKfast(fcn,[0,c_5*h],Y0,B,h_fast);
    U_n5 = Y(:,2);
    
    % Compute U_{n,6}
    p_n6 = @(t) p_n2 + (t/h)*(-c_4/(c_3*(c_3-c_4)))*D_n3 +...
        (t/h)*(c_3/(c_4*(c_3-c_4)))*D_n4 +...
        (t^2/(2*h^2))*(2/(c_3*(c_3-c_4)))*D_n3 - ...
        (t^2/(2*h^2))*(2/(c_4*(c_3-c_4)))*D_n4;
    fcn = @(t,y) feval(An)*y + p_n6(t);
    h_fast = c_6*h/ceil(c_6*m);
    [~,Y,~] = solve_ERKfast(fcn,[0,c_6*h],Y0,B,h_fast);
    U_n6 = Y(:,2);
    
    % Intermediate functions D_n5 and D_n6
    D_n5 = D_ni(t_n,c_5,U_n5);
    D_n6 = D_ni(t_n,c_6,U_n6);
    
    % Compute u_{n+1}
    q_n6 = @(t) p_n2 + (t/h)*(-c_6/(c_5*(c_5-c_6)))*D_n5 +...
        (t/h)*(c_5/(c_6*(c_5-c_6)))*D_n6 +...
        (t^2/(2*h^2))*(2/(c_5*(c_5-c_6)))*D_n5 - ...
        (t^2/(2*h^2))*(2/(c_6*(c_5-c_6)))*D_n6; 

%     % Alternate q
%     q_n6 = @(t) p_n2 + (t/h)*(-c_4/(c_3*(c_3-c_4)))*D_n3 +...
%         (t/h)*(c_3/(c_4*(c_3-c_4)))*D_n4 +...
%         (t^2/(2*h^2))*(2/(c_3*(c_3-c_4)))*D_n3 - ...
%         (t^2/(2*h^2))*(2/(c_4*(c_3-c_4)))*D_n4;
    
    fcn = @(t,y) feval(An)*y + q_n6(t);
    h_fast = h/m;
    [~,Y,~] = solve_ERKfast(fcn,[0,h],Y0,B,h_fast);
    u_np1 = Y(:,2);
    
    % Update time step
    t_n = t_n + h;
    n = n+1;
    
    % Update u value
    u_n = u_np1;   

end
end