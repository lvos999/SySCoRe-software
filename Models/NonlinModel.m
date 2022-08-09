classdef NonlinModel
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type = "NonLin";
        f % the deterministic function for the state update
        C % C matrix to map from state to to output space
        Bw % B matrix for noise
        mu % noise mean
        sigma % noise variance
        dim  % the dimension of the state space

        param % any additional matrices to be shipped
        
        X
        U
        regions
        AP
        P
        
        
    end
    
    methods
        function obj = NonlinModel(f,C,Bw,mu,sigma,param)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.f = f;
            obj.C =  C;
            obj.Bw = Bw;
            obj.mu = mu;
            obj.sigma = sigma;
            obj.dim = size(Bw, 1);
            obj.param = param;
        end
        
                
        function x_n = f_det(obj,x,u)
            %F compute the next states undisturbed by noise based on the deterministic/nominal dynamics
            x_n =   obj.f(x,u,obj);
        end
        
    end
end
