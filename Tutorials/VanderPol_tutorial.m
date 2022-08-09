clear all
clc
close all
%% Van der Pol tutorial

disp(' Van der Pol Oscillator tutorial ')

% Load model into sysNonLin
Vanderpol
theta_dev = 0.3;
Theta = [1-theta_dev, 1+theta_dev];


% Bounds on state space 
x1l = -4;   % Lowerbound x1
x1u = 4;   % Upperbound x1
x2l = -4;   % Lowerbound x2
x2u = 4;   % Upperbound x2
sysNonLin.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');

% Bounds on  and input space
ul = [-1];   % Lowerbound input u
uu = [1];     % Upperbound input u
sysNonLin.U = Polyhedron([ul(1),uu(1)]');



% Stay in bounded region
P1 = Polyhedron(combvec([-3,3],[-3,3])');

% Reach in bounded region
P2 = Polyhedron(combvec([2,3],[-1,1])');

sysNonLin.regions = [P1;P2]; % regions that get specific atomic propositions
sysNonLin.AP = {'p1', 'p2'}; % with the corresponding atomic propositions


Plot_sysLTI(sysNonLin)


%% 1. Synthesize scLTL formula (or input DFA yourself)
%%% use LTL2BA and check if determinstic and accepting state with loop with 1.
% input: (sc)LTL formula and atomic propositions (see readme in folder
% LTL2BA)
% output: struct DFA containing (among other) the transitions

formula = '(p1 U p2)';  % p1 = safe region, p2 = target region
% formula should use atomic propositions in sysLTI.AP. 

% Make sure your current folder is the main SySCoRe folder
[DFA] = TranslateSpec(formula,sysNonLin.AP);



%% 2. Construct abstract model by gridding it
disp('start gridding');tic
% input: sysLTI, sigma, space bounds for input (controller) and state space
% output: matrix P with P(i,j,k) the probability of going from state i to state j with
% input uhat(:,k)
lu = 10;
uhat = linspace(ul(1),uu(1),lu);
l = [100, 100];  % number of grid cells in x1- and x2-direction
tol=10^-8;

% Grid system with efficient nonlinear gridding
sysAbs = GridSpace_nonlin_tensor(sysNonLin,uhat,l,tol);


% Save some extra system parameters into struct
sysAbs.orig = sysNonLin;

label = zeros(1,prod(l));
[label(1:prod(l))] = deal(1);
index_labels = sysNonLin.regions.contains(sysAbs.states);
[label(index_labels(1,:))] = deal(3);
[label(index_labels(2,:))] = deal(4);
sysAbs.labels = label;
 
toc
disp('---> finished gridding')

%% 3. Compute delta based on epsilon
disp('start computing eps delta');tic;



disp('1. compute max deviation (d) of f as a function of theta, x, u')
dx1 = max(abs(1-sysNonLin.X.V(:,1).^2));
dx2 = max(abs(sysNonLin.X.V(:,2)));
dtheta = max(abs(sysNonLin.param.theta-Theta));
d = dtheta * ...
    max([max(abs(1-sysNonLin.X.V(:,1).^2)),1]) * ...
    dx2;



warning('  -------------- NOT IMPLEMENTED : compute delta and epsilon --------------------- ')

% TODO: change gridSize from scalar to vector
epsilon = 0.0;     
delta = 0.001;



rel = SimRel(epsilon,delta,eye(2));

disp(['delta = ', num2str(delta), ', epsilon = ', num2str(epsilon) ])

toc
rel.NonDetLabels  = NonDeterministicLabelling(sysAbs.outputs, sysNonLin.regions, rel);

disp('---> finished computing eps delta')


%% Synthesize controller
disp('start computing robust controller')

% input: DFA, P, epsilon, delta, N
% output: robust satisfaction probability satProp and robust control policy
% pol

% NOT INCLUDED:
% - epsilon (so we use epsilon = 0)
% - refinement, so we get abstract policy.
% - we have converged value function, not satisfaction probability (so no
% q0_bar)

N = 30;     % time horizon

[satProp,pol] = SynthesizeRobustController(sysAbs,DFA, rel, N, true);

disp('---> finished computing robust controller')

[X1hat, X2hat] = ndgrid(sysAbs.hx{1},sysAbs.hx{2});

figure;
surf(X1hat, X2hat,reshape(satProp,l(1),l(2)),'EdgeColor','interp')
xlabel('x_1')
ylabel('x_2')
title('Value function')