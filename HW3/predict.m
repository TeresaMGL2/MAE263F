%(Junlin Chen, 505947745, predict.m)
function Y_pred = predict(X, parameters)
% Returns the predicted classes of the inference images using parameters containing the weights and biases of the neural network.
% Inputs:
%         X: Inference images with sizes 784 x N. N is the number of images. 
%         parameters: a struct containing the weights and biases (W1, b1,
%         W2, b2, etc.)
% Output: 
%         Y_pred: a 10 x N array containing the predicted labels of each input image.
%
% (Hint: The predicted class is the class which has the highest probability.  )
a = forward_propagation(X,parameters);
a4 = a{end}; %predicted label distribution
Y_pred = zeros(10,size(X,2));
for i=1:size(X,2)
[~,Yl] = max(a4(:,i)); %find the most likely index
Y_pred(Yl,i) = 1;
end
end