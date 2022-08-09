function [P,Q,R] = ComputeProjection(sysLTI,sysLTIr)
% This fun computes the projection matrices P and Q

% Inputs:
% - sysLTI = original model
% - sysLTIr = reduced-order model

% Outputs:
% - matrix P 
% - matrix Q
% - R = 1

P = sym('P', [sysLTI.dim,sysLTIr.dim]);
Q = sym('Q', [1,sysLTIr.dim]);

A = sysLTI.A;
B = sysLTI.B;
C = sysLTI.C;

Ar = sysLTIr.A;
Cr = sysLTIr.C;

S = solve([A*P+B*Q-P*Ar;Cr-C*P] == [zeros(sysLTI.dim,1);0]);

P = double(vpa([   S.P1_1, S.P1_2; S.P2_1, S.P2_2; S.P3_1, S.P3_2; S.P4_1, S.P4_2;
       S.P5_1, S.P5_2; S.P6_1, S.P6_2; S.P7_1, S.P7_2],4));

Q = double(vpa([S.Q1,S.Q2],4));
R = 1;

%[ADD] bound on Q*xhat --> code below gives infeasible!
% P = sdpvar(sysLTI.dim,sysLTIr.dim,'full');
% Q = sdpvar(1,sysLTIr.dim,'full');
%
% obj = [A*P+B*Q-P*Ar];
%
% LMI = [(A*P+B*Q-P*Ar == 0); (Cr-C*P == 0)];
% Xhatr = sysLTIr.X;
% for i = 1:size(Xhatr.V,1)
%     LMI2 = [(-Q*Xhatr.V(i,:)'+uuq >= 0); (Q*Xhatr.V(i,:)'-ulq >= 0)];
%     LMI = [LMI; LMI2];
% end
%
% sol = optimize(LMI,obj);
%
% P = double(value(P));
% Q = double(value(Q));

end

