function is_pos_def = isPositiveDefinite(A)
% Description: this function returns true if A is positive definite, 
% false otherwise.

try
    chol(A);
    is_pos_def = true;  % if chol does not throw an error, A is positive definite.
catch ME
    is_pos_def = false;  % if chol throws an error, A is not positive definite.
end

end
