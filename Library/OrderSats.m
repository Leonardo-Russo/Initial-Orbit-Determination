function [NewSats, indices] = OrderSats(Sats)
% Description: this function takes as input the Satellites Data Structures
% and yields a new version in which the observable times are linear in
% time.

N = length(Sats);

indices = zeros(N, 1);

for i = 1 : N

    indices(i) = find(Sats(i).Observable(:), 1, "first");   % index of first observation

end

[indices, sort_indices] = sort(indices,'ascend');

NewSats = Sats(sort_indices);

end