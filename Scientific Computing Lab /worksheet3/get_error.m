
function e = get_error(A, B)
% GET_ERROR calculates the root-mean-squared distance between two vectors 
% or matrices

    e = sqrt(sum(sum((A-B).^2))/numel(A));

end

