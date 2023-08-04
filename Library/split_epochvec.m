function [year, month, day, hour, min, sec] = split_epochvec(epoch_vec)
% Description: this function recieves an epoch vector and returns each
% component as a single variable.

year = epoch_vec(1);
month = epoch_vec(2);
day = epoch_vec(3);
hour = epoch_vec(4);
min = epoch_vec(5);
sec = epoch_vec(6);

end