function dy = g_p2(t, y)
% usage: dy = g_p2(t, y)
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% March 2018

% model parameters
global Pdata;
a  = Pdata.a; 
b  = Pdata.b; 
ep = Pdata.ep;

% form the ODE RHS
dy = [a - (y(3)+1)*y(1) + y(1)*y(1)*y(2);
      y(3)*y(1) - y(1)*y(1)*y(2);
      b/ep - y(3)*y(1)];

% end function
