function x_bar = predict_epidemic_model(control_bar, MODEL_PARAM, PREV_STATE)

x_1 = PREV_STATE.x_1;

T     = size(control_bar, 1);
x_bar = zeros(T+1, length(x_1));

x_bar(1,:) = x_1;
for k = 1 : T
    x_bar(k+1,:) = epidemic_model(x_bar(k,:), control_bar(k, :), MODEL_PARAM)';
end

end