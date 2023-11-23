%(Junlin Chen, 505947745, compute_cost.m)
function cost = compute_cost(AL, Y)
% The given function calculates the cross-entropy loss, ...
%also known as the log loss, given the predicted values AL and the true values Y
% Inputs:
% AL: final predicted values, a K x N matrix. K is the number of all classes, and N is the number of examples.
% Y: ground truth labels, a K x N matrix. K is the number of all classes, and N is the number of examples.
% Output:
% cost: the entropy loss

    cost = -sum(Y .* log(AL)); %cross entropy loss formula
end