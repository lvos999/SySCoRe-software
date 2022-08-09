function [sysAbs] = Gridding(varargin)
%GRIDDING This function grids the state space and input space of dynamic
%systems
% sysAbs = Gridding(sys, uhat, l, tol)
% 
% Inputs: 
% sys   = LTI systems with fields A and B 
% uhat  = Finite set of inputs,  example uhat = combvec(linspace(ul(1),uu(1),3),linspace(ul(2),uu(2),3));
% l     = Number of finite states in each dimension  [l1 l2 l3 ...]
% tol  =  tolerance for truncating to 0
%  
% Options:
% 'TensorComputation' = false/true 

sys = varargin{1};
uhat = varargin{2};
l = varargin{3};
tol = varargin{4};

TensorComputation = false;
for i = 5:length(varargin)
    % try to find 'TensorComputation'
    if strcmp(varargin{i},'TensorComputation')
        TensorComputation = boolean(varargin{i+1});
    end
    
end


% 'TensorComputation' is False
if ~TensorComputation 
    if sys.dim == 2
        [sysAbs] = GridSpace_2d(sys, uhat, l, tol);
    else
        [sysAbs] = GridSpace_nd(sys, uhat, l, tol);
    end
end

% 'TensorComputation' is True
if TensorComputation
    if sys.dim == 2
        [sysAbs] = GridSpace_nonlin_tensor_v2(sys, uhat,l, tol);
    else 
        Error('Tensor computation is only implemente for 2D')
    end
end

end

