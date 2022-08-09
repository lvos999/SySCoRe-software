function [sccs] = StronglyConnectedComponents(dfa)

DFA_STATES = length(dfa.S);

s = [];
t = [];

for i = 1:DFA_STATES
    for j = unique(dfa.trans(i, :))
        if j ~=  0 && j ~= i
            s = [s i];
            t = [t j];
        end
    end
end

G = digraph(s, t);
sccs = conncomp(G);

end