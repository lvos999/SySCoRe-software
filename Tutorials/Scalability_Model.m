clear all
close all


%% Specify system parameters and regions
% LTI systems of the form
% x(t+1) = Ax(t) + Bu(t) + Bw w(t)
% y(t) = Cx(t) + Du(t)

% define parameters
Ns = .1;
n = 1:2;
Ns_sequence = [1 Ns.^(n)./factorial(n)];
c_sequence = Ns_sequence.*0;
c_sequence(1)=1;

% define dynamics
dim =length(Ns_sequence)-1;
A = toeplitz(c_sequence(1:end-1),Ns_sequence(1:end-1));
B = Ns_sequence(end-1:-1:1)';
C = eye(dim);
D = zeros(dim,1);
Bw = .1*eye(dim);

% Specify mean and variance of disturbance w(t) --> should be 2x1 and 2x2!
mu = 0;
sigma = 1;


%% 
% save all system parameters (incl Bw) into a struct
sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);

%% Spaces and sets

% Take state space equal to safe set
x1l = -10;   % Lowerbound x1
x1u = 10;   % Upperbound x1
x2l = -10;   % Lowerbound x2
x2u = 10;   % Upperbound x2
sysLTI.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');

% single input 
ul = [-1];   % Lowerbound input u
uu = [1];     % Upperbound input u
sysLTI.U = Polyhedron([ul(1),uu(1)]');



% Specify regions for the specification
% (using same format as used in function inpolygon)
p1 = combvec([-10, 10], [-10, 10]);
P1 = Polyhedron(p1');

p2 = combvec([-8, 8], [-8, 8]);        % avoid region
% p2x = p2(1,:);    % x1-coordinates
% p2y = p2(2,:);      % x2-coordinates
P2 = Polyhedron(p2');

% Definte the regions and the atomic propositions
sysLTI.regions = [P1;P2];
sysLTI.AP = {'p1','p2'};

Plot_sysLTI(sysLTI)



%% Synthesize scLTL formula  
%%% use LTL2BA and check if determinstic and accepting state with loop with 1.
% input: (sc)LTL formula and atomic propositions (see readme in folder
% LTL2BA)
% output: struct DFA containing (among other) the transitions
N=5; % Horizon of specification

formula = '( p1 & X p1 & X X p2) ';  % p1 = parking, p2 = avoid region

[DFA] = TranslateSpec(formula, sysLTI.AP);

 
%% Construct abstract model
disp('start gridding')
% input: dim, sys, sigma, space bounds for input (controller) and state space
% output: matrix P with P(i,j,k) the probability of going from state i to state j with
% input uhat(:,k)
% Specify division of input space for actuation and feedback
% ula = 0.6*ul;   % part of input for actuation (lowerbound)
% uua = 0.6*uu;
ula = 1*ul;   % part of input for actuation (lowerbound)
uua = 1*uu;
ulf = ul-ula;   % part of input for feedback (lowerbound)
uuf = uu-uua;


lu = 20;
uhat = combvec(linspace(ula(1),uua(1),lu));

l = 200;  % number of grid cells in x1- and x2-direction
tol=10^-4;


sysAbs = GridSpace_nonlin_tensor(sysLTI,uhat,l,tol);


% determine labels for abstract states
label = zeros(1,l^dim);
[label(1:l^dim)] = deal(1);
 index_labels = sysLTI.regions.contains(sysAbs.states);
[label(index_labels(1,:))] = deal(3);
[label(index_labels(2,:))] = deal(4);
sysAbs.labels = label;

disp(['---> finished gridding in ', num2str(toc), ' seconds.'])
%% Compute delta based on epsilon
disp('start computing eps delta');tic;


% TODO: epsilons delat doesnt work

epsilon = max(sysAbs.beta*10);    % should be larger than vector beta!

% Compute beta as a set
beta = sysAbs.beta';
beta = Polyhedron(combvec([-beta(1),beta(1)],[-beta(2),beta(2)])');
sysAbs.beta = beta;

% can be skipped if it takes too long
% including u = uhat + K(x-xhat) is not working yet, add bounds on K!
% [delta, D_m, K] = ComputeDelta(epsilon,sysLTI,mu,sigma,gridSize,dim,uuf);
 [delta, D_m, K] = ComputeDelta(epsilon,sysLTI,sysLTI.mu,sysLTI.sigma,sysAbs.beta);

%  delta = 0.0163; % this value corresponds to l = 100 or gridSize=0.2 AND K
%= zeros(2,2);
% D_m=eye(2);

disp(['delta = ', num2str(delta), ', epsilon = ', num2str(epsilon) ])
rel = SimRel(epsilon,delta,D_m);



disp(['---> finished computing eps delta in ', num2str(toc), ' seconds'] )

%% Synthesize controller
% input: DFA, P, epsilon, delta, N
% output: robust satisfaction probability satProp and robust control policy
% pol

% NOT INCLUDED:
% - epsilon (so we use epsilon = 0)
% - refinement, so we get abstract policy.
% - we have converged value function, not satisfaction probability (so no
% q0_bar)
disp('Start computing Robust policy')
tic

rel.NonDetLabels  = NonDeterministicLabelling(sysAbs.states, sysLTI.regions, rel);


N = 30;     % time horizon
[V,pol]  = SynthesizeRobustController(sysAbs,DFA,rel,N,uhat,false);

disp(['---> finished computing robust policy in ', num2str(toc), ' seconds.'] )

% check resulting policy
VisualizePolicy(sysAbs, pol, DFA, sysLTI)


%% Start simulation
x0 = [7.5;-7.5];
xsim = x0;
q = DFA.S0;
i=1;
u = uhat;
    if P1.contains(xsim(:,i))
        if P2.contains(xsim(:,i))
            label = 4;
        else
            label = 3; 
        end
    else
        label =1;
    end
q_next = DFA.trans(q,label);

% find initial abstract state
diff = D_m^.5*(abs(x0.*ones(size(sysAbs.states))-sysAbs.states));
dis = sqrt(diff(1,:).^2+diff(2,:).^2);
[~,j] = min(dis);

xhat0 = sysAbs.states(:,j);
xhatsim = [xhat0];
uhat = pol(:,j, q_next);
indexing = 1:length(sysAbs.states);
disp([' Satisfaction probability at xhat0 = ', num2str(V(DFA.S0,j))])
% interface function (CHANGE!)   
% todo: add dfa dynamics for specifications beyond reach-(avoid)

for i = 1:N
    % stop if you have reached the parking region


    disp(['x = ', mat2str(xsim(:,i)), ',    q^+ = ',num2str(q_next)])    % compute next continuous state
   
    if q_next ==DFA.F
        disp("reached target")
        break;
    elseif q_next == DFA.sink  
          disp("failed")
        break;

    end
    w1 = mu + sqrt(sigma).*randn(1);
    w2 = mu + sqrt(sigma).*randn(1);
    w = [w1;w2];
    
    u = uhat;



    xnext = sysLTI.A*xsim(:,i)+sysLTI.B*u+sysLTI.Bw*w;

    xsim = [xsim, xnext];

    if P1.contains(xsim(:,i))
        if P2.contains(xsim(:,i))
            label = 4;
        else
            label = 3; 
        end
    else
        label =1;
    end
    q_next = DFA.trans(q,label);
    

    % find next abstract state, by looking at the coupled update (option
    % 1), or (option 2) looking at the maximizer in R wrt the value
    % function. Use option 2:
    
    diff = abs(xsim(:,end).*ones(size(sysAbs.states))-sysAbs.states);
%     inR = diag(diff'*D_m*diff)<= epsilon^2;
    inR = (([1 1]*((D_m^.5*diff).^2)).^.5)<=epsilon;
    %  XhatSpace(:,inR)

    indices_valid = indexing(inR);
    V_q = V(q_next,:);
    [value_max, index_aux] = max(V_q(inR)); 
    j = indices_valid(index_aux); % find maximizing index of abstract state
    
    xhatnext = sysAbs.states(:,j);
    
    uhat = [pol(:,j,q_next)];
 

      q = q_next;

end

%% Show results
figure;
 [X1hat, X2hat] = ndgrid(sysAbs.hx{1},sysAbs.hx{2});

surf(X1hat,X2hat,reshape(V(DFA.S0,:),l,l),'EdgeColor','interp')
xlabel('x_1')
ylabel('x_2')
title('Value function')

figure(10);
plot(P1,'Color', 'blue')
hold on
plot(P2,'Color', 'green')
xlim([x1l x1u])
ylim([x2l x2u])
xlabel('x_1')
ylabel('x_2')
plot(xsim(1,:),xsim(2,:),'rx','LineWidth', 2)



 
