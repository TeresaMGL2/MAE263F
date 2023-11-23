%(Junlin Chen, 505947745, update_parameters.m)
function parameters = update_parameters(parameters, gradients, learning_rate)
% Updates the parameters of a feedforward neural network using gradient descent
% parameters: learned parameters, a struct containing W1, b1, W2, b2, etc.
% gradients: gradients of the cost with respect to each parameter, a struct containing dW1, db1, dW2, db2, etc.
% learning_rate: learning rate for gradient descent
% returns: updated parameters, a struct containing W1, b1, W2, b2, etc.
for i=1:length(parameters)
    parameters{i}.W = parameters{i}.W-learning_rate*gradients{i}.dW; %using gradient descent algorithm
    parameters{i}.b = parameters{i}.b-learning_rate*gradients{i}.db;
end
end