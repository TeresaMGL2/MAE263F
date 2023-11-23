%(Junlin Chen, 505947745, initialize_parameters.m)
function parameters = initialize_parameters(layer_dims)
% Initializes the weights and biases of the feedforward neural network
% Inputs:
% layer dims: array of layer dimensions, including input and output layers
% Output:
% parameters: a struct containing W1, b1, W2, b2, etc.
% (Hint: parameters{1}.W should be a ( layer dims(2) x 784) matrix, which is the initialized weights connecting ...
%the input layer and the first hidden layer. parameters{1}.b should be a ( layer dims(2) x 1 ) ...
%array, which is the initialized biases connecting the input layer and the first hidden layer. )

L = length(layer_dims);
parameters = cell(1, L-1);
for i = 1:(L-1) %omit ouput layer
        parameters{i}.W  = randn(layer_dims(i+1), layer_dims(i)); %{} directly access the value inside cell, can get the value by parameters{i}.W
        parameters{i}.b = zeros(layer_dims(i+1), 1); %essentially stored both W and b inside one cell, access by p{i}.W/b
end

end