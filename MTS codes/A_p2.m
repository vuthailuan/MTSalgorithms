function A = A_p2()
% usage: A = A_p2()
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% March 2018


% model parameters
global Pdata;
ep = Pdata.ep;

% form the ODE RHS
A = [0,0,0;0,0,0;0,0,-1/ep];

% end function
