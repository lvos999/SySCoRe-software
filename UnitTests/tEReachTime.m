%% Test Class Definition
classdef tEReachTime < matlab.unittest.TestCase

    properties
        P
        MDP
        DFA
    end

    methods  (TestClassSetup)
        function createArguments(testCase)
            % Define model
            GW = createGridWorld(2,2);
            GW.ObstacleStates = '[1,2]';
            GW.TerminalStates = '[2,2]';
            Ps = permute(GW.T, [2,1,3]);
            testCase.P = reshape(Ps, [4, 4,4]);
            hx1 = [1:2];
            hx2 = [1:2];
            testCase.MDP = MDP_model(testCase.P,{hx1,hx2});
            testCase.MDP.states = GW.States;
            testCase.MDP.inputs = 1:length(GW.Actions);
            
            % Define specification
            formula = '( !p1 U p2)';  % p1 = safe region, p2 = target region
            testCase.DFA = TranslateSpec(formula,{'p1','p2'});

            % Labelling
            label = zeros(1,length(testCase.MDP.states));
            for idx = 1:length(testCase.MDP.states)
                 label(idx) = 1 + 1* contains(GW.TerminalStates, GW.States(idx))...
                           + 2* contains(GW.ObstacleStates, GW.States(idx));
            end

            testCase.MDP.labels = label;


        end
    end
    %% Test Method Block
    methods (Test)
        %% Test basic functinality of EReachTime
        function testBasicLoad(testCase)
            createArguments(testCase)
        end

        function testControl(testCase)
            createArguments(testCase)

         [satProp,pol] = SynthesizeController(testCase.MDP,testCase.DFA,...
             3,.1);
         pol
         
         
        end
    end
end