classdef LinModel
    %LINMODEL Class of LTI systems with noise on the transitions. 
    % This model class can be used to grid the continuous dynamics to a
    % finite MDP
    
    properties
        type = "LTI";
        A % A matrix of LTI dynamics
        B % B Matrix of LTI dynamics
        C  % C matrix of LTI matrix towards the output dimensions
        D % D matrix of LTI matrix (Should be equal to zero)
        Bw % B matrix for noise
        mu % noise mean
        sigma % noise variance
        dim  % the dimension of the state space
        
        X % state space
        U % input space 
        regions % regions labelled with atomic propositions
        AP % atomic propositions
        P % the probability transition matrix. (Class: TransitionProbability)
        
        
    end
    
    methods
        function obj = LinModel(A,B,C,D,Bw,mu,sigma)
            %LINMODEL Construct an instance of this class
            %   Load all values A,B,C,D, ...
            obj.A = A;
            obj.B = B;
            obj.D = D;
            obj.C = C;
            obj.Bw = Bw;
            obj.mu = mu;
            obj.sigma = sigma;
            obj.dim = size(Bw, 1);
        end

        function x_n = f_det(obj,x,u)
            %F compute the next states undisturbed by noise based on the deterministic/nominal dynamics
            x_n = obj.A*x+obj.B*u;
        end
        function prop_number = label(obj,x,u)
            %F compute the next states undisturbed by noise based on the deterministic/nominal dynamics
            error('Not implemented')
        end
        
        
    end
end
