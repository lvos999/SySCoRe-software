function [sysAbs] = GridSpace_2d(sys,uhat,l,tol)
% [P,X1hat,X2hat,XhatSpace,gridSize] = GridSpace(dim,sys,Bw,sigma,cf,uhat,l) 
%
% Inputs: 
% sys   = LTI systems with fields A and B 
% sigma = variance on noise signal
% cf    = Upper and lower bounds on x variables. See example: cf = [x1l,x1u,x2l,x2u]';
% uhat  = Finite set of inputs,  example uhat = combvec(linspace(ul(1),uu(1),3),linspace(ul(2),uu(2),3));
% l     = Number of finite states in each dimension 
%
% 
% Outputs:
% P = matrix describing the transition probabilities. 
% P(i,j,k) is the probability of going from state i to state j with
% input uhat(:,k)
% output 2: abstract state space
%
% Recent updates
% - May 19: updated the XhatSpace matrix. 
%
% TODO: 
% - Add Bw
% - Allow non-square space (but rectangular)
% - Allow delta as a vector (different grid size in x-direction than in
%   y-direction)
% - Add a sink state, such that the columns of P add up to 1
%
%
%% 
A = sys.A;
B = sys.B;
dim = length(A);
try 
    dim = sys.dim;
catch
    error('sys.dim does not exist')
end
try
    Bw = sys.Bw;
catch
    error('sys.Bw does not exist')
end
try
    sigma = sys.sigma;
catch
    error('sys.sigma does not exist')
end


if sys.Bw ~= eye(dim) 
    error('Bw different from identity has not been implemented yet') 
end

if dim ~= 2
    error('Bw different from identity has not been implemented yet') 
end
    
 
if sum(sys.mu) ~= 0
    error('mu different from 0 has not been implemented yet') 
end


try
X = sys.X;
catch 
     error('sys.X does not exist or is not a Polyhedron')
end

X.computeVRep;
Xl = min(X.V);
Xu = max(X.V);


%===================== Marginals for uniform grid==========================
% width of single slot
gridSize = (Xu-Xl)/l;   

% selection of representative points in each dimension
hx = Xl+gridSize/2:gridSize:Xu;

% intervals 
mx = Xl:gridSize:Xl+l*gridSize ;        

P = zeros(l^dim,l^dim,length(uhat)); 
% finite location i1,2 per dimension maps to global index: i1, i2 ->  (i1-1)*l+i2


for k = 1:length(uhat)
    for i1=1:l % current x1
        for i2=1:l % current x2
            t1 = normcdf(mx, A(1,:)*[hx(i1);hx(i2)]+B(1,:)*uhat(:,k), sigma);  % Add Bw
            t2 = normcdf(mx, A(2,:)*[hx(i1);hx(i2)]+B(2,:)*uhat(:,k), sigma);

            cp1 = diff(t1); % equivalent to (t1(2:length(t1))-t1(1:length(t1)-1));
            cp2 = diff(t2); % equivalent to (t2(2:length(t2))-t2(1:length(t2)-1));

            cp1(cp1<tol) = 0; % equivalent to (t1(2:length(t1))-t1(1:length(t1)-1));
            cp2(cp2<tol) = 0; % equivalent to (t2(2:length(t2))-t2(1:length(t2)-1));

            p12 = reshape(cp2'*cp1,1,[]);
            P((i1-1)*l+i2,:,k) = p12;
            
        end
    end
end
P = TransitionProbability(P);
beta = Polyhedron((diag(gridSize)*(ff2n(dim)-0.5)')');
% make meshgrid of abstract states
[X1hat, X2hat] = meshgrid(hx,hx);
XhatSpace = combvec(hx,hx);
states = XhatSpace(end:-1:1,:);  % TODO dont change order of combvec, but change order in the above forloop
sysAbs = MDP_model(P,hx,states,beta, sys);

sysAbs.orig = sys;
sysAbs.P = P;
sysAbs.states = states;
sysAbs.outputs = sys.C*states;
sysAbs.inputs = uhat; 

%==========================================================================