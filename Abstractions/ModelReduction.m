function [sysLTIr,F] = ModelReduction(sysLTI,dimr,f)
% This function creates a reduced order model sysLTIr of dimension dimr based on the original
% model sysLTI by using balanced truncations on a closed loop system. 

% Inputs: 
% - sysLTI = original affine system
% - dimr = desired dimension of reduced-order model
% - f = constant (min xCCx + ufu) used to construct feedback matrix F

% Outputs:
% - sysLTIr = reduced order affine system

% get a decent guess for the feedback matrix
[~,~, F]=dare(sysLTI.A,sysLTI.B,sysLTI.C'*sysLTI.C,f);

sysclosed=ss(sysLTI.A-sysLTI.B*F,[sysLTI.B,sysLTI.Bw],sysLTI.C,sysLTI.D,-1); %(ignore disturbance)
%sysclosed=ss(sysLTI.A,[sysLTI.B,sysLTI.Bw],sysLTI.C,sysLTI.D,-1); %(ignore disturbance)
sysred=balred(sysclosed,dimr);
sysred=ss(tf(sysred));

%%% HARD CODED!
% Make sure that the reduced model has C = [1 0]
if sysred.c(1)~= 1
    T=eye(length(sysred.c));
    T(length(sysred.c),length(sysred.c))= sysred.c(end);
    sysred=ss2ss(sysred,T);

    % flip outputs
    T  = [0 1; 1 0];
    sysred=ss2ss(sysred,T);
end

display('Reduced order model obtained' )

%% Parameters of reduced order model
Ar = sysred.A;
Br = sysred.B(:,1);
Cr = sysred.C;
Dr = zeros(1,1);
Bwr = sysred.B(:,2:end);

mur = sysLTI.mu; % mean of disturbance
sigmar = sysLTI.sigma;% variance of disturbance

sysLTIr = LinModel(Ar,Br,Cr,Dr,Bwr,mur,sigmar);

% optional: verify if behaviour of original and reduced-order model is
% similar
%sys = ss(sysLTI.A-sysLTI.B*F,sysLTI.B,sysLTI.C,sysLTI.D);
%sysr = ss(sysLTIr.A,sysLTIr.B,sysLTIr.C,sysLTIr.D);
%bodeplot(sys,sysr,'r--')


end

