function coefficients = ar(history, order)
% Disposable modified-covariance fit for local WEC-Sim diagnostics.
% The published model's System Identification Toolbox is unavailable here.
history = double(history(:));
n = numel(history);
if order < 1 || n <= order
    error('ar:invalidHistory', 'Need more samples than AR coefficients');
end
count = n - order;
forward = zeros(count, order);
backward = zeros(count, order);
for j = 1:count
    forward(j, :) = history(order + j - 1:-1:j).';
    backward(j, :) = history(j + 1:order + j).';
end
targets = [history(order + 1:end); history(1:count)];
coefficients = [forward; backward] \ targets;
end
