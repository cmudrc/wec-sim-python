function future = forecast(coefficients, history, horizon)
% Disposable forecast companion for the numeric coefficient vector above.
history = double(history(:));
coefficients = double(coefficients(:));
order = numel(coefficients);
future = zeros(horizon, 1);
for j = 1:horizon
    next = coefficients.' * history(end:-1:end - order + 1);
    future(j) = next;
    history(end + 1, 1) = next;
end
end
