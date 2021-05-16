function [anomaly] = anomalyR2a(x,x_pred,thrshl)
    
    anomaly = zeros(length(x),1);
    for i = 1:length(x)
        %if((x_pred(i)-x(i))/mean(abs(x)) > epsl)
        %   anomaly(i) = 1; 
        %end
        if(abs(x_pred(i)-x(i)) > thrshl)
            anomaly(i) = 1;
        end
    end

end

