function [anomaly] = anomalyR2a(x,x_pred,thrshl)
    anomaly = zeros(length(x),1);
    for i = 1:length(x) % Iterate through the elements of the vector
        if(abs(x_pred(i)-x(i)) > thrshl) % if difference is abovce threshold
            anomaly(i) = 1; % Set anomaly vector to true
        end
    end
end

