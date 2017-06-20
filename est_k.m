function [k] = est_k(Y, disttype)
    % Structure of Y: feature x observations
    % disttype: euclidean, pearson 
    if(strcmp(disttype, 'euclid'))
        [R, dmax] = similarity_euclid(Y');
    elseif(strcmp(disttype, 'pearson'))
        R = similarity_pearson(Y');
    end
    idx = apcluster(R,median(R(:)));
    k = length(unique(idx));
end

