%(Junlin Chen, 505947745, softmax.m)
function Z = softmax(X)
% softmax applies the softmax function on the input to get the corresponding probability distribution along all ...
%possible classes.
% Inputs:
% X: A K x N matrix. K is the number of all classes, and N is the number of examples.
% Outputs:
% Z: A K x N matrix representing the probability of the dis
Z = exp(X)./sum(exp(X)); %softmax formula, element wise division
end