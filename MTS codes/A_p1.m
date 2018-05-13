function A = A_p1()
% usage: A = A_p1()
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% March 2018


% model parameters
global Pdata;
epsilon = Pdata.epsilon;
gamma = Pdata.gamma;

% form the ODE RHS
A = [gamma,epsilon;epsilon,1];

% end function
