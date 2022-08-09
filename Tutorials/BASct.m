%%% Building_automation case study
% sysLTI -> original model
% sysLTIr -> reduced order model
%
%% Case study: Model with large number of continuous variables
% author: Birgit
% author of 7D model: Nathalie Cauchi
% -------------------------------------------------------
% x_c[k+1] = A_cx_c[k] + B_cu_c[k] +F_cd_c[k] + Q_c
% y_c[k]   = [1 0 0 0 0 0 0]
% x_c = [T_z1 T_z2 T_w5 T_w6 T_w2 T_w3 T_w7]^T
% u_c = T_sa
% d_c =[T_out T_hall CO2_1 CO2_2 T_rw,r1 T_rw,r2]^T
% -------------------------------------------------------
% -------------------------------------------------------

clc; clear;close all;

tStart = tic;

% Load parameters needed to build model
loadBASmodel

%% Specify parameters of 7-dimensional model
A = Z1m.A;
B = Z1m.B;
C = Z1m.C;
D = zeros(1,1);
Bw = Z1m.F(:,1:6);
dim = Z1m.dim;

mu = [9;15;500;500;35;35]; % mean of disturbance
sigma = [1, zeros(1,5); 0, 1, zeros(1,4); 0, 0 , 100, zeros(1,3);
                0,0,0,100,zeros(1,2); zeros(1,4), 5, 0;
                zeros(1,5), 5];% variance of disturbance

% TODO: program this as a model transformation (method for class)
% Transform model, such that w comes from Gaussian distribution with mean 0
% and variance eye().
Bw = Bw*sigma^(1/2);
Qc = Qc+Bw*mu;
mu = zeros(6,1);
sigma = eye(6);

% input bounds
ul = 15;
uu = 30;    % originally uu = 30!

% Work with system with state x-xss and output y-C*xss
% sysLTI = AffineModel(A,B,C,D,Bw,Qc,mu,sigma);
sysLTI = LinModel(A,B,C,D,Bw,mu,sigma);

% Remove Qc by looking at steady state solution
uss = (uu-ul)/2+ul;
uss = 15;
xss = -inv(A-eye(dim))*(B*uss+Qc);

% New bounds on input space (u-uss)
ul = ul-uss;     % Lowerbound input u
uu = uu-uss;     % Upperbound input u
sysLTI.U = Polyhedron(combvec([ul(1),uu(1)])');

% Specify division of input space for actuation and feedback
ula = -0.6*((uu-ul)/2)+((uu-ul)/2)+ul;   % part of input for actuation R*uhat, with R=1
uua = 0.6*((uu-ul)/2)+((uu-ul)/2)+ul;
ulf = -0.175*((uu-ul)/2);   % part of input for feedback K(x-P*xhat)
uuf = 0.175*((uu-ul)/2);
%%% NOT INCLUDED YET!
ulq = -0.225*((uu-ul)/2);   % part of input for Q*xhat
uuq = 0.225*((uu-ul)/2);

% TODO: prune this, I dont think it is needed
% New bounds on state space (x-xss)
x1l = 19.5-xss(1);   % Lowerbound x1
x1u = 20.5-xss(1);   % Upperbound x1
x2l = 19-xss(2);   % Lowerbound x2
x2u = 22-xss(2);   % Upperbound x2
x3l = 18-xss(3);   % Lowerbound x3
x3u = 22-xss(3);   % Upperbound x3
x4l = 18-xss(4);   % Lowerbound x4
x4u = 22-xss(4);   % Upperbound x4
x5l = 18-xss(5);   % Lowerbound x5
x5u = 22-xss(5);   % Upperbound x5
x6l = 18-xss(6);   % Lowerbound x6
x6u = 22-xss(6);   % Upperbound x6
x7l = 18-xss(7);   % Lowerbound x7
x7u = 22-xss(7);   % Upperbound x7
Xbound = [x1l, x1u; x2l, x2u; x3l, x3u; x4l, x4u; x5l, x5u; x6l, x6u; x7l x7u];
sysLTI.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u],[x3l,x3u],[x4l,x4u],[x5l,x5u],[x6l,x6u],[x7l,x7u])');

% Specify regions for the specification
P1 = sysLTI.X;
P1 = Polyhedron(combvec([x1l,x1u],[-1000,1000],[-1000,1000],[-1000,1000],[-1000,1000],[-1000,1000],[-1000,1000])');

%% Reduced order model
% TODO: See whether this can become less hard coded
f = 0.098;
dimr = 2;
[sysLTIr,F] = ModelReduction(sysLTI,dimr,f);

% Compute P and Q
[P, Q, R] = ComputeProjection(sysLTI,sysLTIr);

% Bounds on state space % x_c = [T_z1 ???]^T
x1l = x1l;   % Lowerbound x1
x1u = x1u;   % Upperbound x1


% Compute state space bounds (HARD CODED?)
Ax = [];
bx = [];
for i = 1:2^dim
    xv = sysLTI.X.V(i,:)';

    for j = 1:dim
        if P(j,2)>=0
            Axi(2*(j-1)+1,:) = [0 -P(j,2)];
            Axi(2*j,:) = [P(j,2) 0];
            bxi(2*(j-1)+1,1) = [-xv(j)+P(j,1)*x1u];
            bxi(2*j,1) = [xv(j)-P(j,1)*x1l];
        else
            Axi(2*(j-1)+1,:) = [0 P(j,2)];
            Axi(2*j,:) = [-P(j,2) 0];
            bxi(2*(j-1)+1,1) = [xv(j)-P(j,1)*x1u];
            bxi(2*j,1) = [-xv(j)+P(j,1)*x1l];
        end
    end

    Ax = [Ax; Axi];
    bx = [bx; bxi];
end

xhatBound = linprog([-1;1], Ax,bx);
x2l = xhatBound(1);   % Lowerbound x2
x2u = xhatBound(2);   % Upperbound x2

sysLTIr.X = Polyhedron(combvec([x1l,x1u],[x2l,x2u])');

sysLTIr.U = sysLTI.U;

%% Specify output region
% Specify regions for the specification
TSP = 20-xss(1);               % temperature set point
p1x = [19.5-xss(1) 20.5-xss(1)  20.5-xss(1) 19.5-xss(1) 19.5-xss(1)];    % x1-coordinates
p1y = [-1000 -1000 1000 1000 -1000];    % x2-coordinates
p1 = [p1x; p1y];        % goal region
P1r = Polyhedron(p1x(1:2)');

sysLTIr.regions = [P1r]; % regions that get specific atomic propositions
sysLTIr.AP = {'p1'}; % with the corresponding atomic propositions

%Plot_sysLTI(sysLTI)

%% Synthesize scLTL formula (or input DFA yourself)
N = 6;
formula = '(p1 & X p1 & X X p1 & X X X p1 & X X X p1 & X X X X p1 & X X X X X p1)';
[DFA] = TranslateSpec(formula, sysLTIr.AP);


%% Reduce State space
[~, output_set] = IncreaseDecreasePolytope(P1r, 0.1);
[sysLTIr,~] = ReduceX(sysLTIr, [ula(1),uua(1)], output_set, 'invariance', 5);


%% Construct abstract model by gridding it
disp('start gridding');tic

lu = 3;
uhat = combvec(linspace(ula(1),uua(1),lu));

l = [3000*3000];  % number of grid cells in x1- and x2-direction
tol=10^-6;
sysAbs = GridSpace_nonlin_tensor_v2(sysLTIr,uhat,l,tol);
 l = sysAbs.l;% load updated l
% Save some extra system parameters into struct
sysAbs.orig = sysLTIr;

label = zeros(1,prod(l));
[label(1:prod(l))] = deal(1);
inP1 = inpolygon(sysAbs.states(1,:),sysAbs.states(2,:),p1x,p1y);
[label(inP1)] = deal(2);
sysAbs.labels = label;

toc
disp(['---> finished gridding in ', num2str(toc), ' seconds.'])
%% Compute delta based on epsilon

% Compute polytope of (B*R-P*Br)uhat+(Bw-P*Brw)w
beta = sysAbs.beta;
Uhat = Polyhedron(sysAbs.inputs');

Wlb = sysLTIr.mu-3*sum(sysLTIr.sigma,2);
Wub = sysLTIr.mu+3*sum(sysLTIr.sigma,2);
Wset = Polyhedron('lb', Wlb, 'ub', Wub);

% Compute additional error, by truncating the disturbance!
onemindel = mvncdf(Wlb,Wub,mu,sigma);
del_trunc = 1-onemindel;

%Z = P*beta+(sysLTI.B*R-P*sysLTIr.B)*Uhat+(sysLTI.Bw-P*sysLTIr.Bw)*Wset;
Z = (sysLTI.B*R-P*sysLTIr.B)*Uhat+(sysLTI.Bw-P*sysLTIr.Bw)*Wset;
Zred = Z;
Zred = Z.minVRep();
%Zred = Zred.outerApprox;
% [ADD] bounding box of Zred, with a specified number of vertices!

% [OPTIONAL] This function computes the bounds for epsilon (takes a lot of time)
% ARCH: (epsilon,delta) = (0.28592;0.01)
% epsilonBounds = [eps_max,eps_min]
%[epsilonBounds] = ComputeEpsilonBounds2(sysLTI,sysLTIr,mu,sigma,uuf,Zred,P)
% computed for l= [3000, 3000]; and P
%= [1 0; -0.41 -0.38; 12.25, 3.38; 11.18, 3.08; 12.25, 3.38; 11.18, 3.08; 10.17, 2.85]
%epsilonBounds = [0.1510, 0.0134]; % f=0.02
%epsilonBounds = [0.0956, 0.0257]; % f=0.09, uss=22.5
%epsilonBounds = [0.2060, 0.0320]; % f=0.09, uss=20
epsilonBounds = [0.2413, 0.0520]; 

%epsilon = epsilonBounds(2)+(epsilonBounds(1)-epsilonBounds(2))/2;
%epsilon_1 = 0.1510;
epsilon_1 = epsilonBounds(1);
if epsilon_1>=epsilonBounds(2) && epsilon_1<=epsilonBounds(1)
    disp('Feasible epsilon chosen')
elseif epsilon_1>epsilonBounds(1)
    disp('Epsilon is larger then necessary')
elseif epsilon_1<=epsilonBounds(2)
    disp('Infeasible epsilon chosen')
end

disp('start computing simulation relation');tic

% Compute MOR simulation relation
[delta_1, D_1, K_1, F_1] = ComputeDelta_intPQRK2(epsilon_1,sysLTI,sysLTIr,mu,sigma,uuf,Zred,P);
%delta_1 = 0;
del = delta_1;
%D_1 = [5.1188   -0.2477   -0.6695    0.1937   -0.6731    0.1914   -0.4921;
%  -0.2477    0.0233    0.0410   -0.0154    0.0414   -0.0151    0.0260;
%  -0.6695    0.0410    0.1165   -0.0399    0.1144   -0.0376    0.0798;
%   0.1937   -0.0154   -0.0399    0.0205   -0.0384    0.0174   -0.0228;
%  -0.6731    0.0414    0.1144   -0.0384    0.1181   -0.0396    0.0803;
%   0.1914   -0.0151   -0.0376    0.0174   -0.0396    0.0201   -0.0225;
%  -0.4921    0.0260    0.0798   -0.0228    0.0803   -0.0225    0.0647];

%K_1 = [ -1.2534   -0.0733   -0.3829   -0.0327 -0.3829   -0.0327   -0.3922];

%F_1 = [ 0.000037308676334  -0.000072348569933  -0.000186013608414  -0.000037356752355  -0.000185299814691  -0.000036564339144  -0.000200435305736;
%        0.000025867352865  -0.000063449015250  -0.000160687437544  -0.000032035622611 -0.000161355281380  -0.000032685421257  -0.000174521367980;
%        -0.000046052034194   0.000002407948595   0.000000212426544   0.000000522482792 0.000000182575528   0.000000478061814   0.000000716516353;
%        0.000000002443962  -0.000000012662394   0.000000004384287  -0.000000016552246    0.000000004529298  -0.000000016430980  -0.000000004869327;
%        -0.001031290876769   0.000053923685707   0.000004757087473   0.000011700498111 0.000004088602770   0.000010705733173   0.000016045692442;
%        0.000001412094297  -0.000007316191383   0.000002533192453  -0.000009563705307    0.000002616978604  -0.000009493639153  -0.000002813443718];

delta_1 = delta_1+del_trunc;

% Compute gridding simulation relation
epsilon_2 = max(0.35-epsilon_1,0);
[delta_2, D_2, K_2] = ComputeDelta(epsilon_2,sysLTIr,sysLTIr.mu,sysLTIr.sigma,beta);
% delta_2 = 0;
% D_2 = [6561.9, 1770.5; 1770.5 479.9];

delta = delta_1+delta_2;
epsilon = epsilon_1+epsilon_2;

rel_1 = SimRel(epsilon_1,delta_1,D_1);
rel_2 = SimRel(epsilon_2,delta_2,D_2);
rel = rel_1.Combine(rel_2,sysLTIr.X);
disp(['delta = ', num2str(rel.delta), ', epsilon = ', num2str(rel.epsilon) ])

toc
disp('---> finished computing simulation relation')

rel.NonDetLabels  = NonDeterministicLabelling(sysAbs.outputs, sysLTIr.regions, rel);

%% Synthesize controller
disp('---> Start Value iteration')
[V,pol]  = SynthesizeRobustController(sysAbs,DFA,rel,N,false);

q = DFA.S0;
% Fix Value at initial q0 based on labeling
for i = 1:size(sysAbs.labels,2)
    if sysAbs.labels(i) == 2
        V(q,i) = V(q,i);
    else
        V(q,i) = 0;
    end
end

% Determine computation time
tEnd = toc(tStart);

X1 = reshape(sysAbs.states(1,:),l)+xss(1);
X2 = reshape(sysAbs.states(2,:),l)+xss(2);

VC = reshape(V(q,:),l);

% Plot satisfaction probability
figure;
surf(X1(1:15:end,1:15:end),X2(1:15:end,1:15:end),VC(1:15:end,1:15:end),'EdgeColor','interp')
xlabel('x_{r1}', 'FontSize', 16)
ylabel('x_{r2}', 'FontSize', 16)
%title('Satisfaction probability')
%colorbar
view(2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

% for i = 1:size(C,1)
%     for j= 1:size(C,2)
%         if C(i,j) == 0
%             C(i,j) = 3;
%         else %if C(i,j) ==  max(V(q,:))
%             C(i,j) = 4;
%         end
%     end
% end
%            
% figure;
% h = pcolor(X1(1:100:end,1:100:end),X2(1:100:end,1:100:end),C(1:100:end,1:100:end))
% xlabel('x_{r1}', 'FontSize', 16)
% ylabel('x_{r2}', 'FontSize', 16)
% set(h, 'EdgeColor', 'none')
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName','Times','fontsize',16)

%VisualizePolicy(sysAbs,pol,DFA)

%% Start simulation
YSIM = [];
USIM = [];
% Based on reduced-order dynamics instead of original 7D dynamics
[i,j] = max(V(DFA.S0,:));
xr0 = sysAbs.states(:,j);
x0 = P*xr0;

for l = 1:6 % run the simulation 6 times
    xsim = [x0];
    xrsim = [xr0];

    q = DFA.S0;
        if P1r.contains(sysLTIr.C*xr0)
            label = 2;
        else
            label = 1;
        end
    q_next = DFA.trans(q,label);

    % find initial abstract state
    % [opt 1] look at x0
    %diff = abs(x0.*ones(1,size(sysAbs.states,2))-P*sysAbs.states);
    %dis = sqrt(sum(diff.^2,1));
    %[~,j] = min(dis);

    % [opt 2] look at xr0
    diff = abs(xr0.*ones(1,size(sysAbs.states,2))-sysAbs.states);
    dis = sqrt(sum(diff.^2,1));
    [~,j] = min(dis);

    xhat0 = sysAbs.states(:,j);
    xhatsim = [xhat0];

    % satisfaction probability of this initial state
    SatProp = V(q,j);
    disp(['Value function at x0 equals ', mat2str(SatProp,4),])

    uhat = pol(:,j, q_next);
    indexing = 1:length(sysAbs.states);

    ursim = [];
    usim = [];
    q = q_next;
    for i = 1:N
        disp(['x = ', mat2str(xsim(:,i)+xss,4), ',    q^+ = ',num2str(q_next)])    % compute next continuous state

        if q ==DFA.F
            disp("satisfied specification")
            break;
        elseif q == DFA.sink
            disp("failed specification")
            break;
        end
        w = mvnrnd(sysLTIr.mu,sysLTIr.sigma)';

        u = R*uhat + Q*xrsim(:,i) + K_1*(xsim(:,i)-P*xrsim(:,i));
        ur = uhat;
        ursim = [ursim, ur];
        usim = [usim, u];

        xnext = sysLTI.A*xsim(:,i)+sysLTI.B*u+sysLTI.Bw*w;
        xsim = [xsim, xnext];

        wr = w+F_1*(xsim(:,i)-P*xrsim(:,i)); %wr = w+gamma
        xrnext = sysLTIr.A*xrsim(:,i)+sysLTIr.B*ur+sysLTIr.Bw*wr;
        xrsim = [xrsim, xrnext];

        if P1.contains(xsim(:,end)) % for some reason this does not always work!?! Because -1000<<<-4 ?
            label = 2;
        elseif xsim(1,end)>=19.5-xss(1) && xsim(1,end)<=20.5-xss(1)
            label = 2;
        else
            label = 1;
        end
        q_next = DFA.trans(q,label);

        % find next abstract state, by looking at the coupled update (option
        % 1), or (option 2) looking at the maximizer in R wrt the value
        % function. Use option 2:

        diff2 = abs(xrsim(:,end).*ones(1,size(sysAbs.states,2))-sysAbs.states);
        inR2 = (([1 1]*((D_2^.5*diff2).^2)).^.5)<=epsilon_2;
        %  XhatSpace(:,inR)
        inR = inR2;
        indices_valid = indexing(inR);
        V_q = V(q_next,:);
        [value_max, index_aux] = max(V_q(inR));
        j = indices_valid(index_aux); % find maximizing index of abstract state
        xhatnext = sysAbs.states(:,j);
        xhatsim = [xhatsim, xhatnext];

        uhat = [pol(:,j,q_next)];

        q = q_next;

    end
    
    % compute output
    y = C*xsim+C*xss;
    usim_ss = usim+uss;
    
    YSIM = [YSIM; y];
    USIM = [USIM; usim_ss];
    
end

%% plot results

figure;
subplot(2,1,1)
plot(y(1,:))
hold on
plot(y(1,:),'bo')
plot(1:N,19.5*ones(1,N),'r--')
plot(1:N,20.5*ones(1,N),'r--')
title('State evolution')

subplot(2,1,2)
plot(usim+uss);
hold on
plot(usim+uss,'bo')
title('Input')

k = 0:5;
figure;
plot(k,YSIM')
hold on
plot(0:N-1,19.5*ones(1,N),'r--')
plot(0:N-1,20.5*ones(1,N),'r--')
title('State evolution')

tEnd =toc;
%% Display execution time 
% end time is before simulation part
disp(['Total runtime = ', mat2str(tEnd)])

%% Verify input bound
Qx = Q*sysLTIr.X;
Qx = Qx.minVRep();

u_lb = ula+ulf+min(Qx.V);
u_ub = uua+uuf+max(Qx.V);

% use only (abstract) states with V>0
[~,j] = find(V(DFA.S0,:)>0);
XX = sysAbs.states(:,j);

u_lb2 = ula+ulf+min(Q*XX);
u_ub2 = uua+uuf+max(Q*XX);

disp(['Input is between ', mat2str(u_lb2+uss,4), ' and ', mat2str(u_ub2+uss,4)])