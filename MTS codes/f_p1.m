function dy = f_p1(t, y)
% usage: dy = f_p1(t, y)
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% August 2012
% All Rights Reserved

% model parameters
global Pdata;
w  = Pdata.w; 
gamma  = Pdata.gamma; 
epsilon = Pdata.epsilon;

% form the ODE RHS
dy = [gamma*y(1)+epsilon*y(2)-gamma*cos(t)-epsilon*cos(w*t)-sin(t);...
    epsilon*y(1)-y(2)-epsilon*cos(t)-cos(w*t)-w*sin(w*t)];
end
% end function
