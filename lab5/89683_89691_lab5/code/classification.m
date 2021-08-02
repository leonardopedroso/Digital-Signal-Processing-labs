function isWater = classification(I,sigma_water,sigma_ice)
    % likelihoods
    water_likelihood = raylpdf(I,sigma_water);
    ice_likelihood = raylpdf(I,sigma_ice);

    % classification
    isWater = zeros(size(I));
    % for each pixel in the image
    for i=1:size(I,1)
        for j=1:size(I,2)
            % if its intensity is more likely to correspond to water
            if water_likelihood(i,j) >= ice_likelihood(i,j)
                % signals it should correspond to water
                isWater(i,j) = 1;
            end
        end
    end
end

