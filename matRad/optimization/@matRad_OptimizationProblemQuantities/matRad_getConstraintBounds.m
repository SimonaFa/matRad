function [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
% matRad IPOPT get constraint bounds wrapper function
% 
% call
%   [cl,cu] = matRad_getConstraintBounds(optiProb,cst)
%
% input
%   cst:            matRad cst struct
%
% output
%   cl: lower bounds on constraints
%   cu: lower bounds on constraints
%
% References
%   -
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 the matRad development team. 
% 
% This file is part of the matRad project. It is subject to the license 
% terms in the LICENSE file found in the top-level directory of this 
% distribution and at https://github.com/e0404/matRad/LICENSE.md. No part 
% of the matRad project, including this file, may be copied, modified, 
% propagated, or distributed except according to the terms contained in the 
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BPtype = class(optiProb.BP);
%isEffectBP = strcmp(BPtype,'matRad_EffectProjection');

% Initialize bounds
cl = [];
cu = [];

% compute objective function for every VOI.
for  i = 1:size(cst,1)

    % Only take OAR or target VOI.
    if ~any(cellfun(@isempty,cst{i,4})) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )

        % loop over the number of constraints for the current VOI
        for j = 1:numel(cst{i,6})
            
            constraint = cst{i,6}{j};
            
            % only perform computations for constraints
%{
            if ~isempty(strfind(cst{i,6}(j).type,'constraint'))

                if isequal(options.quantityOpt,'effect')
                    param = cst{i,5}.alphaX .* cst{i,6}(j).dose + cst{i,5}.betaX .* cst{i,6}(j).dose.^2; 
                else 
                    param = cst{i,6}(j).dose;
                end

                if strcmp(cst{i,6}(j).robustness,'none') || strcmp(cst{i,6}(j).robustness,'probabilistic') || strcmp(cst{i,6}(j).robustness,'VWWC') ||...
                   strcmp(cst{i,6}(j).robustness,'COWC') || strcmp(cst{i,6}(j).robustness,'VWWC_CONF')|| strcmp(cst{i,6}(j).robustness,'OWC')


                    [clTmp,cuTmp] = matRad_getConstBounds(cst{i,6}(j),param);
%}
            %if ~isempty(strfind(cst{i,6}{j}.type,'constraint'))
            
            if isa(constraint,'DoseConstraints.matRad_DoseConstraint') || isa(constraint, 'OmegaConstraints.matRad_VarianceConstraint')
                
                %if isEffectBP
                   
                %    doses  = optiFunc.getDoseParameters();
                %    effect = cst{i,5}.alphaX*doses + cst{i,5}.betaX*doses.^2;
                    
                %    optiFunc = optiFunc.setDoseParameters(effect);
                %end
                quantityConstrained = constraint.quantity;
                quantityNames = cellfun(@(x) x.quantityName,optiProb.BP.quantities, 'UniformOutput',false);
                quantityConstrainedInstance = optiProb.BP.quantities{strcmp(quantityConstrained,quantityNames)};

                constraint = quantityConstrainedInstance.setBiologicalDosePrescriptions(constraint,cst{i,5}.alphaX,cst{i,5}.betaX);
  
                 cl = [cl;constraint.lowerBounds(numel(cst{i,4}{1}))];
                 cu = [cu;constraint.upperBounds(numel(cst{i,4}{1}))];
                    
                %end
            end

        end % over all objectives of structure

    end % if structure not empty and target or oar

end % over all structures
   

