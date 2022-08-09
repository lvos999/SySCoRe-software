function [state2act] = NonDeterministicLabelling(outputs,Polytopes,rel)
%NonDeterministicLabelling Map the states to the labelling number
% TODO: These compuations should be done over the output space. 
epsil = rel.epsilon;
Polytope_array_large = [];
Polytope_array_small = [];
Polytope_containment_large = [];
Polytope_containtment_small = [];

n_p = length(Polytopes); % number of polytopes
for p = 1:n_p
[PolytopeLarge,PolytopeSmall] = IncreaseDecreasePolytope(Polytopes(p), epsil);

Polytope_array_large = [Polytope_array_large;PolytopeLarge.computeVRep];
Polytope_array_small = [Polytope_array_small;PolytopeSmall.computeVRep];

end

% for all states compute whether they are contained in the polytopes
Polytope_containment_large = Polytope_array_large.contains(outputs);
Polytope_containment_small = Polytope_array_small.contains(outputs);
state2act = [];
bin_vals = dec2bin(0:2^n_p-1);
for index = 1:2^n_p
    el = bin_vals(index, :);
    x_true = [];
    for p_index = 1:n_p
        if el(p_index)=='0'
            % could polytope p not hold?
            % check this for all states b looking at the smaller version of
            % polytope p
            x_true = [x_true; 1-Polytope_containment_small(p_index,:)];
        else
            x_true = [x_true; Polytope_containment_large(p_index,:)];

        end
    end
 
    if size(x_true,1)>1
    state2act = [state2act; all(x_true)];
    else
    state2act = [state2act; x_true];
    end
 
end
end
