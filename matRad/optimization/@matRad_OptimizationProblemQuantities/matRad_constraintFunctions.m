function c = matRad_constraintFunctions(optiProb,w,dij,cst)
% matRad IPOPT callback: constraint function for inverse planning 
% supporting max dose constraint, min dose constraint, min mean dose constraint, 
% max mean dose constraint, min EUD constraint, max EUD constraint, 
% max DVH constraint, min DVH constraint 
% 
% call
%   c = matRad_constraintFunctions(optiProb,w,dij,cst)
%
% input
%   optiProb:   option struct defining the type of optimization
%   w:          bixel weight vector
%   dij:        dose influence matrix
%   cst:        matRad cst struct
%
% output
%   c:          value of constraints
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


% get current dose / effect / RBExDose vector
optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

% get the used scenarios
useScen  = optiProb.BP.scenarios;
scenProb = optiProb.BP.scenarioProb;

% retrieve matching 4D scenarios
fullScen = cell(ndims(d),1);
[fullScen{:}] = ind2sub(size(d),useScen);
contourScen = fullScen{1};

% Initializes constraints
c = [];

% compute objective function for every VOI.
for  i = 1:size(cst,1)
   
   % Only take OAR or target VOI.
   if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
      
      % loop over the number of constraints for the current VOI
      for j = 1:numel(cst{i,6})
         
         constraint = cst{i,6}{j};
         
         % retrieve the robustness type
         robustness = constraint.robustness;
         
         % only perform computations for constraints
         if isa(constraint,'DoseConstraints.matRad_DoseConstraint')
            quantityConstrained = constraint.quantity;

            quantityNames = cellfun(@(x) x.quantityName,optiProb.BP.quantities, 'UniformOutput',false);
            quantityConstrainedInstance = optiProb.BP.quantities{strcmp(quantityConstrained,quantityNames)};
            % rescale dose parameters to biological optimization quantity if required
            constraint = quantityConstrainedInstance.setBiologicalDosePrescriptions(constraint,cst{i,5}.alphaX,cst{i,5}.betaX);

            switch robustness
                case 'none' % if conventional opt: just sum objectives of nominal dose
                  
                  if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                      d_i = d.(quantityConstrained){1}(cst{i,4}{1});
                  elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                      d_i = d.(quantityConstrained){i};
                  end
                  
                  c = [c; constraint.computeDoseConstraintFunction(d_i)];
                  
               case 'PROB' % if prob opt: sum up expectation value of objectives
                  
                   if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                      d_i = d.(quantityConstrained){1}(cst{i,4}{1});
                  elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                      d_i = d.(quantityConstrained){i};
                  end

                  c = [c; constraint.computeDoseConstraintFunction(d_i)];
                  
               case 'VWWC'  % voxel-wise worst case - takes minimum dose in TARGET and maximum in OAR
                  contourIx = unique(contourScen);
                  if ~isscalar(contourIx)
                     % voxels need to be tracked through the 4D CT,
                     % not yet implemented
                     matRad_cfg.dispError('4D VWWC optimization is currently not supported');
                  end
                  
                  % prepare min/max dose vector
                  if ~exist('d_tmp','var')
                     d_tmp = [d{useScen}];
                  end
                  
                  d_Scen = d_tmp(cst{i,4}{contourIx},:);
                  
                  d_max = max(d_Scen,[],2);
                  d_min = min(d_Scen,[],2);
                  
                  if isequal(cst{i,3},'OAR')
                     d_i = d_max;
                  elseif isequal(cst{i,3},'TARGET')
                     d_i = d_min;
                  end
                  
                  c = [c; constraint.computeDoseConstraintFunction(d_i)];
                  
               case 'VWWC_INV'  %inverse voxel-wise conformitiy - takes maximum dose in TARGET and minimum in OAR
                  contourIx = unique(contourScen);
                  if ~isscalar(contourIx)
                     % voxels need to be tracked through the 4D CT,
                     % not yet implemented
                     matRad_cfg.dispError('4D inverted VWWC optimization is currently not supported');
                  end
                  
                  % prepare min/max dose vector
                  if ~exist('d_tmp','var')
                     d_tmp = [d{:}];
                  end
                  
                  d_Scen = d_tmp(cst{i,4}{contourIx},:);
                  d_max = max(d_Scen,[],2);
                  d_min = min(d_Scen,[],2);
                  
                  if isequal(cst{i,3},'OAR')
                     d_i = d_min;
                  elseif isequal(cst{i,3},'TARGET')
                     d_i = d_max;
                  end
                  
                  c = [c; constraint.computeDoseConstraintFunction(d_i)];   
               otherwise
                  matRad_cfg.dispError('Robustness setting %s not yet supported!',constraint.robustness);
            end
            
            
         elseif isa(constraint, 'OmegaConstraints.matRad_VarianceConstraint')

            quantityConstrained = constraint.quantity;
            quantityNames = cellfun(@(x) x.quantityName,optiProb.BP.quantities, 'UniformOutput',false);
            quantityConstrainedInstance = optiProb.BP.quantities{strcmp(quantityConstrained,quantityNames)};
            % rescale dose parameters to biological optimization quantity if required
            constraint = quantityConstrainedInstance.setBiologicalDosePrescriptions(constraint,cst{i,5}.alphaX,cst{i,5}.betaX);
            
            switch robustness

                case 'PROB'
                    allVoxels = arrayfun(@(scenStruct) scenStruct{1}, cst{i,4}, 'UniformOutput',false);
                    nVoxels = numel(unique([allVoxels{:}]));
                    d_i = d.(quantityConstrained){i};
                    c = [c; constraint.computeVarianceConstraintFunction(d_i, nVoxels)];
            end
           
         end
         
      end % if we are a constraint
      
   end % over all defined constraints & objectives
   
end % if structure not empty and oar or target


end % over all structures
