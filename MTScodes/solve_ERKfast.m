function [tvals,Y,nf] = solve_ERKfast(fcn,tvals,Y0,B,h)
% usage: [tvals,Y,nf] = solve_ERKfast(fcn,tvals,Y0,B,h)
%
% Explicit Runge-Kutta solver for the vector-valued ODE problem 
%     y' = F(t,Y), t in tvals, y in R^m,
%     Y(t0) = [y1(t0), y2(t0), ..., ym(t0)]'.
%
% Inputs:
%     fcn    = string holding function name for F(t,Y)
%     tvals  = vector of start time and final time [t0,tf]
%     Y0     = initial value array (column vector of length m)
%     B      = Butcher matrix for ERK coefficients, of the form
%                 B = [c A;
%                      q b]
%              Here, c is a vector of stage time fractions (s-by-1),
%                    A is a matrix of Butcher coefficients (s-by-s),
%                    q is an integer denoting the method order of accuracy,
%                    b is a vector of solution weights (1-by-s),
%     h      = time step
%
% Outputs: 
%     tvals  = the same as the input array tvals
%     y      = [y(t0),y(tN)], where each
%               y(t*) is a column vector of length m.
%     nf = number of function calls
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% Spring 2018
% Adjusted from the code by D. Reynolds ("solve_ERK.m", Aug 2012)


   
% extract ERK method information from B
[~, Bcols] = size(B);
s = Bcols - 1;        % number of stages
c = B(1:s,1);         % stage time fraction array
A = B(1:s,2:s+1);     % RK coefficients

% initialize output arrays
m = length(Y0);
Y = zeros(m,1);
Y(:,1) = Y0;

% set the solver parameters
ONEMSM   = 1-sqrt(eps);  

% initialize temporary variables
t = tvals(1);
Ynew = Y0;
      
% initialize work & function call counter
nsteps = 0;
tstep = 2; 

% create Fdata structure for evaluating solution
Fdata.fname = fcn;    % ODE RHS function name
Fdata.B     = B;      % Butcher table 
Fdata.nf    = 0;

% loop over internal time steps to get to desired output time
while (t < tvals(tstep)*ONEMSM) 
    
    Fdata.h = h;
    Fdata.yold = Y0;

    % set Fdata values for this step
    Fdata.t    = t;    % time step

    % initialize data storage for multiple stages
    z = zeros(m,s);
    
    % loop over stages
    for stage=1:s
        % construct stage solution
	    %    zi = y_n + h*sum_{j=1}^{i-1} (A(i,j)*f(zj))
	    z(:,stage) = Y0;
        for j=1:stage-1
            z(:,stage) = z(:,stage) + h*A(stage,j)*feval(fcn,t+h*c(j),z(:,j));
            Fdata.nf = Fdata.nf + 1;
        end
    end

    % increment number of internal time steps taken
    nsteps = nsteps + 1;

    % compute new solution (and embedding if available)
    Ynew = Y_ERK(z,Fdata);
    
    % update solution and time for last successful step
    Y0 = Ynew;
    t  = t + h;
           
end  % end while loop attempting to solve steps to next output time

% store updated solution in output array
Y(:,tstep) = Ynew;
nf = Fdata.nf;
% end solve_ERK function
end

function y = Y_ERK(z, Fdata)
% usage: y = Y_ERK(z, Fdata)
%
% Inputs:
%    z     = stage solutions [z1, ..., zs]
%    Fdata = structure containing extra problem information
%
% Outputs: 
%    y     = step solution built from the z values


% extract method information from Fdata
B = Fdata.B;
[~, Bcols] = size(B);
s = Bcols - 1;
c = B(1:s,1);
b = (B(s+1,2:s+1))';

% get some problem information
[zrows,zcols] = size(z);
nvar = zrows;
if (zcols ~= s)
   error('Y_ERK error: z has incorrect number of stages');
end

% call RHS at our stages
f = zeros(nvar,s);
for is=1:s
   t = Fdata.t + Fdata.h*c(is);
   f(:,is) = feval(Fdata.fname, t, z(:,is));
   Fdata.nf = Fdata.nf + 1;
end

% form the solutions
y  = Fdata.yold + Fdata.h*f*b;

% end of function
end
