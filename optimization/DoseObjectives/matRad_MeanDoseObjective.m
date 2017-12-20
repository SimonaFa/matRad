classdef matRad_MeanDoseObjective < matRad_DoseObjective
    %MATRAD_MEANDOSE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        name = 'Mean Dose'
    end
    
    properties
        parameters = {'reference mean'; 0}
        penalty = 1
    end
    
    methods 
        %% Calculates the Objective Function value
        function fDose = computeDoseObjectiveFunction(obj,dose)
            fDose = obj.penalty * abs(mean(dose(:)) - obj.parameters{2,1});
        end
        
        %% Calculates the Objective Function gradient
        function fDoseGrad   = computeDoseObjectiveGradient(obj,dose)
            fDoseGrad = (obj.penalty/numel(dose))*sign(dose - obj.parameters{2,1});
        end
    end
    
end

