
function [linMod, R] = ReduceX(linMod, inputset,  output_region, type, horizon, varargin)
% REDUCEX reduce the polytopic statespace set based on an output set and a
% property of either time-bounded invariance or reachability (not
% implemented)

% check whether conditions under which this code holds are satisfied
if ~strcmp(type, 'invariance')
    error(['The type = ', type, ' is not yet known or implemented'])

end
R_init = Polyhedron(output_region.A*linMod.C, output_region.b);



system = LTISystem('A',  linMod.A, 'B', linMod.B);
system.x.min = min(linMod.X.V, [], 1);
system.x.max = max(linMod.X.V, [], 1);
system.u.min = inputset(:,1);
system.u.max = inputset(:,2);

if linMod.U.Dim ~=1
warning('works for 1 inputs only')
end
U = Polyhedron(inputset');

%% Debugging figures

figure(1)
subplot(2,1,1)
hold off
R_init.plot() 
title(['backward invariance set i = ', num2str(1)])

%% init reach set
R= cell(horizon,1);

R{1}= R_init;
V = [ ];
for i=1:horizon


R_ = system.reachableSet('X', R{i}, 'U', U, 'N', 1, ...
                         'direction', 'backward');
R{i+1} = R_ & R_init;
V = [V;R{i+1}.V];
%% Debugging figures
% figure(i+1)
% subplot(2,1,1)
% hold off
% R{i+1}.plot() 
% title(['backward invariance set i = ', num2str(i+1)])
% xl = xlim;yl = ylim;
% subplot(2,1,2)
% xlim(xl)
% ylim(yl)
% hold on
% scatter(R{i+1}.V(:,1), R{i+1}.V(:,2), 60, [204 204 255]/256, 'filled')

end

%% Debugging figures
% figure(i+1)
% subplot(2,1,1)
% hold on
% R{i+1}.plot() 
X = Polyhedron(V);
X.minVRep
linMod.X = X;

end
