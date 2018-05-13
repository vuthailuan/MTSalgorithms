function dy = g_p1(t, y)
% usage: dy = g_p1(t, y)
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% March 2018

% model parameters
global Pdata;
w  = Pdata.w; 
gamma  = Pdata.gamma; 
epsilon = Pdata.epsilon;

% form the ODE RHS
dy(1,1) = -gamma*cos(t)-epsilon*cos(w*t)-sin(t);
dy(2,1) = -epsilon*cos(t)-cos(w*t)-w*sin(w*t);

% end function
