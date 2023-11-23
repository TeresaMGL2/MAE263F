%(Junlin Chen, 505947745, forward_propagation.m)
function activations = forward_propagation(X, parameters)
% Computes the output of a feedforward neural network given input data and learned parameters
% X: input data, shape (input size, number of examples)
% parameters: learned parameters, a struct containing W1, b1, W2, b2, etc.
% returns: array of activations at each layer, including input and output layers
activations = cell(1, length(parameters)+1); %stores all activation value
activations{1} = X;
input = X;
for i = 1:length(parameters)
    z = parameters{i}.W*input+parameters{i}.b; %applying weight and bias
    if i==length(parameters)
        a = softmax(z); %if at output layer, use softmax
    else
        a = tanh2(z);%else, use activation
    end
    activations{i+1} = a;
    input = a;
end
end