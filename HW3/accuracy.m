%(Junlin Chen, 505947745, accuracy.m)
function acc = accuracy(Y_pred, Y)
% Returns the accuracy, which is the percentage of correct predictions.
% Inputs:
%         Y_pred: a 10 x N array containing the predicted labels of each input image.
%         Y: a 10 x N array containing the actual labels of each input image.
% Output: 
%         acc: The percentage of accurately predicted samples.
count = 0;
for i=1:size(Y,2) %loop all label
[~,Yl] = max(Y(:,i)); %actual index
[~,Ypl] = max(Y_pred(:,i)); %predicted index
if Yl==Ypl %if the same index
    count=count+1;
end
end
acc=count/size(Y,2);%calculate % accuracy
end