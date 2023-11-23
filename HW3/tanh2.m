%(Junlin Chen, 505947745, tanh2.m)
function Z = tanh2(X)
% Tanh applies the tanh activation function on the input to get a nonlinear output.
% Inputs:
% X: A M x N matrix representing the output of the neurons, which serves as the input of the ...
%activation function. M is the number of neurons, and N is the number of examples
% Outputs:
% Z: a M x N matrix representing the output after the tanh activation function

Z = (exp(X)-exp(-X))./(exp(X)+exp(-X)); %element wise division
end