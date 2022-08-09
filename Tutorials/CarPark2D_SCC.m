clc
clear
close all


%% Specify system parameters and regions
% LTI systems of the form
% x(t+1) = Ax(t) + Bu(t) + Bw w(t)
% y(t) = Cx(t) + Du(t)
% Define all system parameters (incl Bw) into a struct
 A = 0.9*eye(2);
 B = 0.7*eye(2);
 C = eye(2);
 D = zeros(2);
 Bw = eye(2);
 dim = length(A);

% Specify mean and variance of disturbance w(t) --> should be 2x1 and 2x2!
 mu = 0; % mean of disturbance
 sigma = 1; % variance of disturbance

 
sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);
 
% Bounds on state space 
    x1l = -10;   % Lowerbound x1
    x1u = 10;   % Upperbound x1
    x2l = -10;   % Lowerbound x2
    x2u = 10;   % Upperbound x2
sysLTI.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');

% Bounds on  and input space
    ul = [-1;-1];   % Lowerbound input u
    uu = [1;1];     % Upperbound input u
sysLTI.U = Polyhedron(combvec([ul(1),uu(1)],[ul(2),uu(2)])');

% Specify regions for the specification
    p1x = [4 4 10 10 4];    % x1-coordinates
    p1y = [0 -4 -4 0 0];    % x2-coordinates
    p1 = [p1x; p1y];         % parking region
P1 = Polyhedron(p1');

    p2x = [4 4 10 10 4];    % x1-coordinates
    p2y = [0 4 4 0 0];      % x2-coordinates
    p2 = [p2x; p2y];        % avoid region
P2 = Polyhedron(p2');

sysLTI.regions = [P1;P2]; % regions that get specific atomic propositions
sysLTI.AP = {'p1', 'p2'}; % with the corresponding atomic propositions

Plot_sysLTI(sysLTI)


%% Synthesize scLTL formula (or input DFA yourself)
%%% use LTL2BA and check if determinstic and accepting state with loop with 1.
% input: (sc)LTL formula and atomic propositions (see readme in folder
% LTL2BA)
% output: struct DFA containing (among other) the transitions

formula = '(!p2 U p1) & (F p2)';  % p1 = parking, p2 = avoid region
AP = {'p1', 'p2'};

% Make sure your current folder is the main SySCoRe folder
[DFA] = TranslateSpec(formula,AP);
DFA.P = [P1;P2];
%% Construct abstract model'
disp('start gridding')
% input: dim, sys, sigma, space bounds for input (controller) and state space
% output: matrix P with P(i,j,k) the probability of going from state i to state j with
% input uhat(:,k)

lu = 3;
uhat = combvec(linspace(ul(1),uu(1),lu),linspace(ul(2),uu(2),lu));
tic
l = 50;  % number of grid cells in x1- and x2-direction
tol=10^-6;
% TODO: add sink state to P?
[sysAbs] = Gridding(sysLTI,uhat,l,tol);
toc
%%

% determine labels for abstract states
label = zeros(1,l^2);
[label(1:l^2)] = deal(1);
inP1 = inpolygon(sysAbs.states(1,:),sysAbs.states(2,:),p1x,p1y);
inP2 = inpolygon(sysAbs.states(1,:),sysAbs.states(2,:),p2x,p2y);
[label(inP1)] = deal(3);
[label(inP2)] = deal(2);
sysAbs.labels = label;


disp('---> finished gridding')

%% Compute delta based on epsilon
disp('start computing eps delta')
tic
% TODO: change gridSize from scalar to vector
epsilon = 1.005;    % should be larger than vector beta!

% can be skipped if it takes too long
% including u = uhat + K(x-xhat). 
% Assuming that uuf=-ulf and that uuf(1)=uuf(2)!
% [delta, D_m, K] = ComputeDelta(epsilon,sysLTI,sysLTI.mu,sysLTI.sigma,sysAbs.beta);

delta = 0.0163; % this value corresponds to l = 100 or gridSize=0.2 AND K
K= zeros(2,2);
D_m=eye(2);
toc
disp('---> finished computing eps delta')
%%
rel = SimRel(epsilon,delta,D_m);
rel.NonDetLabels  = NonDeterministicLabelling(sysAbs.outputs, sysLTI.regions, rel);

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

N = 100000;     % time horizon

[satProp,pol] = SynthesizeRobustController_SCC(sysAbs,DFA, rel, N, true);
% [satProp,pol] = SynthesizeRobustController(sysAbs,DFA, rel, N, true);


VisualizePolicy(sysAbs, pol, DFA,l, sysLTI)

%% Start simulation
x0 = [-2;-4];
xsim = x0;
% find initial abstract state
diff = abs(x0.*ones(size(sysAbs.states))-sysAbs.states);
dis = sqrt(diff(1,:).^2+diff(2,:).^2);
[~,j] = min(dis);

xhat0 = sysAbs.states(:,j);
xhatsim = [xhat0];
uhat = pol(:,j);
indexing = 1:length(sysAbs.states);
reached_p1 = false;
u = uhat+K*(x0-xhat0);
for i = 1:N
    w1 = mu + sqrt(sigma).*randn(1);
    w2 = mu + sqrt(sigma).*randn(1);
    w = [w1;w2];
    % compute next state
    xnext = sysLTI.A*xsim(:,i)+sysLTI.B*u+sysLTI.Bw*w;
    xsim = [xsim, xnext];

    % stop if you have reached the parking region
    if inpolygon(xnext(1),xnext(2),p1x,p1y)
        disp('Reached parking area - p1');
        reached_p1 = true;
    end

    if reached_p1 && inpolygon(xnext(1),xnext(2),p2x,p2y)
        disp('Reached parking area - p2');
        break;
    end

    % find next abstract state, by looking at the coupled update (option
    % 1), or (option 2) looking at the maximizer in R wrt the value
    % function. Use option 2:
    
    inR = rel.inR(xsim(:,end),sysAbs.states);
    %  XhatSpace(:,inR)
    indices_valid = indexing(inR);
    [value_max, index_aux] = max(satProp(inR)); 
    j = indices_valid(index_aux); % find maximizing index of abstract state
    xhatnext = sysAbs.states(:,j);

    % compute next input
    uhat = [pol(:,j)];
    u = uhat+K*(xnext-xhatnext);
    
    figure(10)
    quiver(xsim(1,end),xsim(2,end),u(1),u(2))
    hold on
        

end

%% Show results
figure;
hx1 = sysAbs.hx;
hx2 = sysAbs.hx;

[X2hat, X1hat] = ndgrid(hx2,hx1);
surf(X1hat,X2hat,reshape(satProp,l,l),'EdgeColor','interp')
xlabel('x_1')
ylabel('x_2')
title('Value function')

figure(10);
plot(p1x,p1y,'g+-','LineWidth',2)
axis equal
hold on
plot(p2x,p2y,'r+-','LineWidth',2)
xlim([x1l x1u])
ylim([x2l x2u])
xlabel('x_1')
ylabel('x_2')
plot(xsim(1,:),xsim(2,:))
plot(xsim(1,:),xsim(2,:),'rx')
