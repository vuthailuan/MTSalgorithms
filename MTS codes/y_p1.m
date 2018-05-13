function y_true = y_p1(t)
% usage: dy = f_p1(t)
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% April 2018

% model parameters
global Pdata;
w  = Pdata.w; 
gamma  = Pdata.gamma; 
epsilon = Pdata.epsilon;

y_true = [cos(t);cos(w*t)];
end