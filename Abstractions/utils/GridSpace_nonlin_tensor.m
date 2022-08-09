function [sysAbs] = GridSpace_nonlin_tensor(sys, uhat,l, tol)
% [P,X1hat,X2hat,XhatSpace,gridSize] = GridSpace(sys, uhat,l) 
%
% Inputs: 
% sys   = LTI systems with fields A and B 
% uhat  = Finite set of inputs,  example uhat = combvec(linspace(ul(1),uu(1),3),linspace(ul(2),uu(2),3));
% l     = Number of finite states in each dimension  [l1 l2 l3 ...]
%
% 
% Outputs:
% P = matrix describing the transition probabilities. 
% P(i,j,k) is the probability of going from state i to state j with
% input uhat(:,k)
% output 2: abstract state spaceh
%
% Recent updates
% - May 19: updated the XatSpace matrix. 
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


try 
    dim = sys.dim;
catch
    error('sys.dim does not exist')
end

try
    mu = sys.mu;
catch
    mu = zeros(dim,1);
end

try
    sys.X.computeVRep;
    Xl = min(sys.X.V);
    Xu = max(sys.X.V);
catch 
     error('sys.X does not exist or is not a Polyhedron')
end

% cf    = Upper and lower bounds on x variables. See example: cf = [x1l,x1u,x2l,x2u]';

 
if Bw ~= eye(dim) 
    error('Bw different from identity has not been implemented yet') 
end

if dim ~= 2
    error('Dimension different from 2 has not been implemented yet') 
end
    
if length(l)==1
    l = l*ones(1,dim);
end
if length(sigma)==1
    sigma = sigma*ones(1,dim);
end
if length(mu)==1
    mu = mu*ones(1,dim);
end

%===================== Compute uniform grid==========================
% width of single slot
gridSize = (Xu-Xl)./l;   

% selection of representative points in each dimension 
% using a cell function
% hx = Xl+gridSize/2:gridSize:Xu;
hx  = arrayfun( @(xl,xu,gridsize)  xl+gridsize/2:gridsize:xu, Xl,Xu,gridSize,'UniformOutput',false);
XhatSpace = combvec(hx{:});


%===================== Compute deterministic probability matrix==========================
% P_det is going to be a sparse matrix. 
index_total = 1:length(XhatSpace);
j_indices = [];
i_indices = [];
for k = 1:length(uhat) % for each control action
    x_n = sys.f_det(XhatSpace,uhat(:,k)*ones(1,size(XhatSpace,2)));
    x_n_ind = ones(dim,1)+floor(diag(gridSize.^-1)*(x_n-Xl'));
    indices = min([[1;1]<=x_n_ind;l'>=x_n_ind],[],1);
    xi_indices= arrayfun(@(i)x_n_ind(i,indices),1:dim,'UniformOutput',false);
    i_indices = [i_indices, sub2ind(l,xi_indices{:})];
    j_indices = [j_indices, index_total(indices)+(k-1)*length(XhatSpace)];
end
    P_det = sparse(i_indices,j_indices, ones(size(i_indices)),length(XhatSpace),length(XhatSpace)*length(uhat));

%===================== Compute stochastic probability matrix==========================
% P_stoch
% compute extended grid (2x the size) around the origin. 
mx2  = arrayfun( @(gridsize,l_i)  -gridsize*0.5:gridsize:gridsize*(l_i-0.5),gridSize,l,'UniformOutput',false);


P = cell(dim,1); 
for k = 1:length(uhat)
            for d_index = 1:dim
               cp =  diff(normcdf(mx2{d_index},  mu(d_index), sigma(d_index)));
               cp(cp<tol) = 0; % equivalent to (t1(2:length(t1))-t1(1:length(t1)-1));
               P{d_index} = toeplitz(cp);
            end
end

P = TensorTransitionProbability(l,P_det,P{:});
beta = 2*gridSize; 
states = XhatSpace;
sysAbs = MDP_model(P,hx,states,beta, sys);
sysAbs.inputs= uhat;  




%==========================================================================