function jacob = matRad_constraintJacobian(optiProb,w,dij,cst)
% matRad IPOPT callback: jacobian function for inverse planning 
% supporting max dose constraint, min dose constraint, min mean dose constraint, 
% max mean dose constraint, min EUD constraint, max EUD constraint, max DVH 
% constraint, min DVH constraint 
% 
% call
%   jacob = matRad_jacobFunc(optiProb,w,dij,cst)
%
% input
%   optiProb: option struct defining the type of optimization
%   w:    bixel weight vector
%   dij:  dose influence matrix
%   cst:  matRad cst struct
%
% output
%   jacob: jacobian of constraint function
%
% References
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%d = matRad_backProjection(w,dij,optiProb);
%d = optiProb.matRad_backProjection(w,dij);
optiProb.BP.compute(dij,w);
d = optiProb.BP.GetResult();

% initialize jacobian (only single scenario supported in optimization)
% jacob = sparse([]);
% 
% % initialize projection matrices and id containers
% DoseProjection{1}          = sparse([]);
% mAlphaDoseProjection{1}    = sparse([]);
% mSqrtBetaDoseProjection{1} = sparse([]);
% voxelID                     = [];
% constraintID                = [];
% 
% % get the used scenarios
% useScen  = optiProb.BP.scenarios;
% scenProb = optiProb.BP.scenarioProb;

% retrieve matching 4D scenarios
% fullScen      = cell(ndims(d),1);
% [fullScen{:}] = ind2sub(size(d),useScen);
% contourScen   = fullScen{1};

% compute objective function for every VOI.
for i = 1:size(cst,1)
   
   % Only take OAR or target VOI.
   if ~isempty(cst{i,4}{1}) && ( isequal(cst{i,3},'OAR') || isequal(cst{i,3},'TARGET') )
      
      % loop over the number of constraints for the current VOI
      for j = 1:numel(cst{i,6})
         
         constraint = cst{i,6}{j}; %Get the Optimization Object
         
         % retrieve the robustness type
         robustness = constraint.robustness;
        
         % only perform computations for constraints
         if isa(constraint,'DoseConstraints.matRad_DoseConstraint')
    
             quantityConstrained = constraint.quantity;

             quantityNames = cellfun(@(x) x.quantityName,optiProb.BP.quantities, 'UniformOutput',false);
             quantityConstrainedInstance = optiProb.BP.quantities{strcmp(quantityConstrained,quantityNames)};
             % rescale dose parameters to biological optimization quantity if required
             constraint = quantityConstrainedInstance.setBiologicalDosePrescriptions(constraint,cst{i,5}.alphaX,cst{i,5}.betaX);
            
            
            % rescale dose parameters to biological optimization quantity if required
            % constraint = optiProb.BP.setBiologicalDosePrescriptions(constraint,cst{i,5}.alphaX,cst{i,5}.betaX);
            
            if ~exist('fJacob', 'var') || ~isfield(fJacob,quantityConstrained)
                if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                    fJacob.(quantityConstrained) = cell(1);%cell(size(d.(quantityConstrained)));
                    fJacob.(quantityConstrained) = {zeros(dij.doseGrid.numOfVoxels,1)};
                elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                    fJacob.(quantityConstrained)= cell(size(cst,1),1);
                    %fJacob.(quantityConstrained)(:) = {};           
                end
            
            end
            
            switch robustness
               
                case 'none' % if conventional opt: just sum objectiveectives of nominal dose
                  if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                      d_i = d.(quantityConstrained){1}(cst{i,4}{1});
                  elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                      d_i = d.(quantityConstrained){i};
                  end

                  jacobSub = constraint.computeDoseConstraintJacobian(d_i);
               case 'PROB' % if prob opt: sum up expectation value of objectives
                  if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                      d_i = d.(quantityConstrained){1}(cst{i,4}{1});
                  elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                      d_i = d.(quantityConstrained){i};
                  end
                  
                  jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                  
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
                  jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                  
               case 'VWWC_INV'  %inverse voxel-wise conformitiy - takes maximum dose in TARGET and minimum in OAR
                  contourIx = unique(contourScen);
                  if ~isscalar(contourIx)
                     % voxels need to be tracked through the 4D CT,
                     % not yet implemented
                     matRad_cfg.dispError('4D inverted VWWC optimization is currently not supported');
                  end
                  
                  % prepare min/max dose vector
                  if ~exist('d_tmp','var')
                     d_tmp = [d{useScen}];
                  end
                  
                  d_Scen = d_tmp(cst{i,4}{contourIx},:);
                  
                  d_max = max(d_Scen,[],2);
                  d_min = min(d_Scen,[],2);
                  
                  if isequal(cst{i,3},'OAR')
                     d_i = d_min;
                  elseif isequal(cst{i,3},'TARGET')
                     d_i = d_max;
                  end
                  jacobSub = constraint.computeDoseConstraintJacobian(d_i);
                  
               otherwise
                  matRad_cfg.dispError('Robustness setting %s not yet supported!',constraint.robustness);
            end
            
            nConst = size(jacobSub,2);


            if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                startIx = size(fJacob.(quantityConstrained){1},2) + 1;

                fJacob.(quantityConstrained) = {[fJacob.(quantityConstrained){1},sparse(dij.doseGrid.numOfVoxels,nConst)]}; 
                fJacob.(quantityConstrained){1}(cst{i,4}{1},startIx:end) = jacobSub;
            elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                startIx = size(fJacob.(quantityConstrained){i},2) + 1;

                fJacob.(quantityConstrained)(i) = {[fJacob.(quantityConstrained){i},sparse(1,nConst)]}; 
                fJacob.(quantityConstrained){i}(1,startIx:end) = jacobSub;
            end
            % %Iterate through columns of the sub-jacobian
            % if isa(optiProb.BP,'matRad_DoseProjection') && ~isempty(jacobSub) || isa(optiProb.BP,'matRad_ConstantRBEProjection')
            % 
            %    startIx = size(DoseProjection{1},2) + 1;
            %    %First append the Projection matrix with sparse zeros
            %    DoseProjection{1}          = [DoseProjection{1},sparse(dij.doseGrid.numOfVoxels,nConst)];
            % 
            %    %Now directly write the jacobian in there
            %    DoseProjection{1}(cst{i,4}{1},startIx:end) = jacobSub;
            % 
            % elseif isa(optiProb.BP,'matRad_EffectProjection') && ~isempty(jacobSub)
            % 
            %    if isa(optiProb.BP,'matRad_VariableRBEProjection')
            %       scaledEffect = (dij.gamma(cst{i,4}{1}) + d_i);
            %       jacobSub     = jacobSub./(2*dij.bx(cst{i,4}{1}) .* scaledEffect);
            %    end
            % 
            %    startIx = size(mAlphaDoseProjection{1},2) + 1;
            % 
            %    %First append the alphaDose matrix with sparse
            %    %zeros then insert
            %    mAlphaDoseProjection{1}    = [mAlphaDoseProjection{1},sparse(dij.doseGrid.numOfVoxels,nConst)];
            %    mAlphaDoseProjection{1}(cst{i,4}{1},startIx:end) = jacobSub;
            % 
            %    %The betadose has a different structure due to the
            %    %quadratic transformation, but in principle the
            %    %same as above
            %    mSqrtBetaDoseProjection{1} =  [mSqrtBetaDoseProjection{1}, sparse(repmat(cst{i,4}{1},nConst,1),repmat(1:numel(cst{i,4}{1}),1,nConst),2*reshape(jacobSub',[],1),dij.doseGrid.numOfVoxels,nConst*numel(cst{i,4}{1}))];
            % 
            %    if isempty(constraintID)
            %       newID = 1;
            %    else
            %       newID = constraintID(end)+1;
            %    end
            % 
            %    voxelID = [voxelID;repmat(cst{i,4}{1},nConst,1)];                         %Keep track of voxels for organizing the sqrt(beta)Dose projection later
            %    constraintID = [constraintID, ...
            %       reshape(ones(numel(cst{i,4}{1}),1)*[newID:newID+nConst-1],[1 nConst*numel(cst{i,4}{1})])];  %Keep track of constraints for organizing the sqrt(beta)Dose projection later
            % end
            
        
         elseif isa(constraint, 'OmegaConstraints.matRad_VarianceConstraint')
            quantityConstrained = constraint.quantity;
            quantityNames = cellfun(@(x) x.quantityName,optiProb.BP.quantities, 'UniformOutput',false);
            quantityConstrainedInstance = optiProb.BP.quantities{strcmp(quantityConstrained,quantityNames)};
            % rescale dose parameters to biological optimization quantity if required
            constraint = quantityConstrainedInstance.setBiologicalDosePrescriptions(constraint,cst{i,5}.alphaX,cst{i,5}.betaX);
            
            if ~exist('fJacob', 'var') || ~isfield(fJacob,quantityConstrained)
                if isa(quantityConstrainedInstance, 'matRad_DistributionQuantity')
                    fJacob.(quantityConstrained) = cell(1);
                elseif isa(quantityConstrainedInstance, 'matRad_ScalarQuantity')
                    %fJacob.(quantityConstrained) = cell(1);
                    fJacob.(quantityConstrained) = cell(size(cst,1),1);           
                end
            
            end


            allVoxels = arrayfun(@(scenStruct) scenStruct{1}, cst{i,4}, 'UniformOutput',false);
            nVoxels = numel(unique([allVoxels{:}]));

            switch robustness

                case 'PROB'
                  d_i = d.(quantityConstrained){i};
                  jacobSub = constraint.computeVarianceConstraintJacobian(d_i, nVoxels);    
            end

            nConst = size(jacobSub,2);

            startIx = size(fJacob.(quantityConstrained){1},2) + 1;
            fJacob.(quantityConstrained){i} = [fJacob.(quantityConstrained){i},jacobSub];
         end
         
      end
      
   end
   
end


optiProb.BP.computeConstraintJacobian(dij,fJacob,w);
j = optiProb.BP.wJacob;

jacob = sparse([]);
for qtIdx=optiProb.BP.constrainedQuantities
    nScensOrStructs   = find(cellfun(@(x) ~isempty(x), j.(qtIdx{1})))';
    for elementIdx=nScensOrStructs
        jacob = [jacob; j.(qtIdx{1}){elementIdx}];
    end
end
% scenario = 1;
% % enter if statement also for protons using a constant RBE
% if isa(optiProb.BP,'matRad_DoseProjection')
% 
%    if ~isempty(DoseProjection{scenario})
%       jacob = DoseProjection{scenario}' * dij.physicalDose{scenario};
%    end
% 
% elseif isa(optiProb.BP,'matRad_ConstantRBEProjection')
% 
%    if ~isempty(DoseProjection{scenario})
%       jacob = DoseProjection{scenario}' * dij.RBE * dij.physicalDose{scenario};
%    end
% 
% elseif isa(optiProb.BP,'matRad_EffectProjection')
% 
%    if ~isempty(mSqrtBetaDoseProjection{scenario}) && ~isempty(mAlphaDoseProjection{scenario})
%       mSqrtBetaDoseProjection{scenario} = mSqrtBetaDoseProjection{scenario}' * dij.mSqrtBetaDose{scenario} * w;
%       mSqrtBetaDoseProjection{scenario} = sparse(voxelID,constraintID,mSqrtBetaDoseProjection{scenario},...
%          size(mAlphaDoseProjection{scenario},1),size(mAlphaDoseProjection{scenario},2));
% 
%       jacob   = mAlphaDoseProjection{scenario}' * dij.mAlphaDose{scenario} +...
%          mSqrtBetaDoseProjection{scenario}' * dij.mSqrtBetaDose{scenario};
% 
%    end

gradientChecker = 0;
if gradientChecker == 1
    f =  matRad_constraintFunctions(optiProb,w,dij,cst);
    epsilon = 1e-5;

    ix = randi([dij.totalNumOfBixels],numel(f),5);

    for jIx=1:numel(f)
        for i=ix(jIx,:)
    
            wInit = w;
            wInit(i) = wInit(i) + epsilon;
            fDel = matRad_constraintFunctions(optiProb,wInit,dij,cst);
            numGrad = (fDel(jIx) - f(jIx))/epsilon;
            diff = (numGrad/jacob(jIx,i) - 1)*100;
            fprintf(['grad val #' num2str(i) '- rel diff numerical and analytical gradient = ' num2str(diff) '\n']);
            %fprintf([' any nan or zero for photons' num2str(sum(isnan(glog{1}))) ',' num2str(sum(~logical(glog{1}))) ' for protons: ' num2str(sum(isnan(glog{2}))) ',' num2str(sum(~logical(glog{2}))) '\n']);
        end
    end
end

end
