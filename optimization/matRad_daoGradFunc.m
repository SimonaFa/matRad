function g = matRad_daoGradFunc(apertureInfoVec,dij,cst,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad IPOPT callback: gradient function for direct aperture optimization
%
% call
%   g = matRad_daoGradFunc(apertureInfoVec,apertureInfo,dij,cst,type)
%
% input
%   apertureInfoVec:  aperture info in form of vector
%   dij:              matRad dij struct as generated by bixel-based dose calculation
%   cst:              matRad cst struct
%   options:          option struct defining the type of optimization
%
% output
%   g: gradient
%
% References
%   [1] http://dx.doi.org/10.1118/1.4914863
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2015 the matRad development team.
%
% This file is part of the matRad project. It is subject to the license
% terms in the LICENSE file found in the top-level directory of this
% distribution and at https://github.com/e0404/matRad/LICENSES.txt. No part
% of the matRad project, including this file, may be copied, modified,
% propagated, or distributed except according to the terms contained in the
% LICENSE file.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read in the global apertureInfo and apertureVector variables
global matRad_global_apertureInfo;
% update apertureInfo from the global variable
apertureInfo = matRad_global_apertureInfo;

% update apertureInfo, bixel weight vector an mapping of leafes to bixels
if ~isequal(apertureInfoVec,apertureInfo.apertureVector)
    apertureInfo = matRad_daoVec2ApertureInfo(apertureInfo,apertureInfoVec);
    matRad_global_apertureInfo = apertureInfo;
end

% bixel based gradient calculation
bixelG = matRad_gradFuncWrapper(apertureInfo.bixelWeights,dij,cst,options);

% allocate gradient vector for aperture weights and leaf positions
%changed from NaNs to zeros
g = zeros(size(apertureInfoVec,1),1);

if apertureInfo.runVMAT
    %we're doing VMAT
    offset = 1;
    DAOBeams = find([apertureInfo.propVMAT.beam.DAOBeam]);
    
    % 1. calculate aperatureGrad
    % 2. find corresponding bixel to the leaf Positions and aperture
    % weights to calculate the gradient
    
    % loop over all beams
    
    for i = 1:numel(apertureInfo.beam)
        
        % get used bixels in beam
        ix = ~isnan(apertureInfo.beam(i).bixelIndMap);
        
        if apertureInfo.propVMAT.beam(i).DAOBeam
            %optimized beam, do regular gradient
            %must always add to existing gradient, since gradient comes
            %from optimized and interpolated beams
            g(offset) = g(offset)+apertureInfo.beam(i).shape(1).shapeMap(ix)' ...
                * bixelG(apertureInfo.beam(i).bixelIndMap(ix)) ./ apertureInfo.beam(i).shape(1).jacobiScale;
            
            
            %gradient wrt leaf positions
            indInOptVec = apertureInfo.beam(i).shape(1).vectorOffset-1+[(1:apertureInfo.beam(i).numOfActiveLeafPairs) apertureInfo.totalNumOfLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)];
            indInBixVec = apertureInfo.beam(i).bixOffset-1+[(1:apertureInfo.beam(i).numOfActiveLeafPairs) apertureInfo.doseTotalNumOfLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)];
            
            g(indInOptVec) = g(indInOptVec)+apertureInfo.beam(i).shape(1).weight*bixelG(apertureInfo.bixelIndices(indInBixVec)) / apertureInfo.bixelWidth;
            
            %increment offset
            offset = offset+1;
        else
            %not optimized beam, aperture weight is interpolated between
            %previous and next optimized weights
            
            %give fraction of gradient to previous optimized beam
            
            %first weight
            lastDAOInd = find(DAOBeams == apertureInfo.propVMAT.beam(i).lastDAOIndex,1);
            g(lastDAOInd) = g(lastDAOInd)+apertureInfo.propVMAT.beam(i).fracFromLastDAO.*(apertureInfo.beam(i).time./apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).time) ...
                *apertureInfo.beam(i).shape(1).shapeMap(ix)' * bixelG(apertureInfo.beam(i).bixelIndMap(ix))./apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape(1).jacobiScale;
            %g(lastOptInd) = g(lastOptInd)+(apertureInfo.beam(i).fracFromLastDAO*apertureInfo.beam(i).doseAngleBordersDiff*apertureInfo.beam(apertureInfo.beam(i).lastDAOIndex).gantryRot ...
            %./(apertureInfo.beam(apertureInfo.beam(i).lastDAOIndex).doseAngleBordersDiff*apertureInfo.beam(i).gantryRot))*apertureInfo.beam(i).shape(1).shapeMap(ix)' ...
            %* bixelG(apertureInfo.beam(i).bixelIndMap(ix)) ./ apertureInfo.beam(apertureInfo.beam(i).lastDAOIndex).shape(1).jacobiScale;
            
            %now leaf pos
            indInOptVec = apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape(1).vectorOffset-1+[(1:apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).numOfActiveLeafPairs) apertureInfo.totalNumOfLeafPairs+(1:apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).numOfActiveLeafPairs)];
            indInBixVec = apertureInfo.beam(i).bixOffset-1+[(1:apertureInfo.beam(i).numOfActiveLeafPairs) apertureInfo.doseTotalNumOfLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)];
            
            g(indInOptVec) = g(indInOptVec)+apertureInfo.propVMAT.beam(i).fracFromLastDAO*apertureInfo.beam(i).shape(1).weight*bixelG(apertureInfo.bixelIndices(indInBixVec)) / apertureInfo.bixelWidth;
            
            %now time
            lastDAOIndTime = lastDAOInd+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
            g(lastDAOIndTime) = g(lastDAOIndTime)+(apertureInfo.propVMAT.beam(i).doseAngleBordersDiff.*apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).timeFacCurr) ...
                .*(-apertureInfo.propVMAT.beam(i).fracFromLastDAO.*apertureInfo.propVMAT.beam(i).timeFracFromNextDAO.*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape(1).weight./apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff).*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).time./apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).time.^2) ...
                +(1-apertureInfo.propVMAT.beam(i).fracFromLastDAO).*apertureInfo.propVMAT.beam(i).timeFracFromLastDAO.*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape(1).weight./apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff).*(1./apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).time)) ...
                * apertureInfo.beam(i).shape(1).shapeMap(ix)' * bixelG(apertureInfo.beam(i).bixelIndMap(ix));
            
            %give the other fraction to next optimized beam
            
            %first weight
            nextDAOInd = find(DAOBeams == apertureInfo.propVMAT.beam(i).nextDAOIndex,1);
            g(nextDAOInd) = g(nextDAOInd)+(1-apertureInfo.propVMAT.beam(i).fracFromLastDAO).*(apertureInfo.beam(i).time./apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).time) ...
                *apertureInfo.beam(i).shape(1).shapeMap(ix)' * bixelG(apertureInfo.beam(i).bixelIndMap(ix))./apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape(1).jacobiScale;
            %g(nextOptInd) = g(nextOptInd)+((1-apertureInfo.beam(i).fracFromLastDAO)*apertureInfo.beam(i).doseAngleBordersDiff*apertureInfo.beam(apertureInfo.beam(i).nextDAOIndex).gantryRot ...
            %./(apertureInfo.beam(apertureInfo.beam(i).nextDAOIndex).doseAngleBordersDiff*apertureInfo.beam(i).gantryRot))*apertureInfo.beam(i).shape(1).shapeMap(ix)' ...
            %* bixelG(apertureInfo.beam(i).bixelIndMap(ix)) ./ apertureInfo.beam(apertureInfo.beam(i).nextDAOIndex).shape(1).jacobiScale;
            
            %now leaf pos
            indInOptVec = apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape(1).vectorOffset-1+[(1:apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).numOfActiveLeafPairs) apertureInfo.totalNumOfLeafPairs+(1:apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).numOfActiveLeafPairs)];
            indInBixVec = apertureInfo.beam(i).bixOffset-1+[(1:apertureInfo.beam(i).numOfActiveLeafPairs) apertureInfo.doseTotalNumOfLeafPairs+(1:apertureInfo.beam(i).numOfActiveLeafPairs)];
            
            g(indInOptVec) = g(indInOptVec)+(1-apertureInfo.propVMAT.beam(i).fracFromLastDAO)*apertureInfo.beam(i).shape(1).weight*bixelG(apertureInfo.bixelIndices(indInBixVec)) / apertureInfo.bixelWidth;
            
            %now time
            nextDAOIndTime = nextDAOInd+apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2;
            g(nextDAOIndTime) = g(nextDAOIndTime)+(apertureInfo.propVMAT.beam(i).doseAngleBordersDiff.*apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).timeFacCurr) ...
                .*(apertureInfo.propVMAT.beam(i).fracFromLastDAO.*apertureInfo.propVMAT.beam(i).timeFracFromNextDAO.*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).shape(1).weight./apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).doseAngleBordersDiff).*(1./apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).time) ...
                -(1-apertureInfo.propVMAT.beam(i).fracFromLastDAO).*apertureInfo.propVMAT.beam(i).timeFracFromLastDAO.*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).shape(1).weight./apertureInfo.propVMAT.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).doseAngleBordersDiff).*(apertureInfo.beam(apertureInfo.propVMAT.beam(i).lastDAOIndex).time./apertureInfo.beam(apertureInfo.propVMAT.beam(i).nextDAOIndex).time.^2)) ...
                * apertureInfo.beam(i).shape(1).shapeMap(ix)' * bixelG(apertureInfo.beam(i).bixelIndMap(ix));
            
        end
    end
else
    %we're not doing VMAT
    
    % 1. calculate aperatureGrad
    % loop over all beams
    offset = 0;
    for i = 1:numel(apertureInfo.beam)
        
        % get used bixels in beam
        ix = ~isnan(apertureInfo.beam(i).bixelIndMap);
        
        % loop over all shapes and add up the gradients x openingFrac for this shape
        for j = 1:apertureInfo.beam(i).numOfShapes
            g(j+offset) = apertureInfo.beam(i).shape(j).shapeMap(ix)' ...
                * bixelG(apertureInfo.beam(i).bixelIndMap(ix))./apertureInfo.beam(i).shape(j).jacobiScale;
        end
        
        % increment offset
        offset = offset + apertureInfo.beam(i).numOfShapes;
        
    end
    
    % 2. find corresponding bixel to the leaf Positions and aperture
    % weights to calculate the gradient
    g(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2) = ...
        apertureInfoVec(apertureInfo.mappingMx(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2)) ...
        .* bixelG(apertureInfo.bixelIndices(1:apertureInfo.totalNumOfLeafPairs*2)) ./ ...
        (apertureInfo.bixelWidth.*apertureInfo.jacobiScale(apertureInfo.mappingMx(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs*2,2)));
    
end




% correct the sign for the left leaf positions
g(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs) = ...
    -g(apertureInfo.totalNumOfShapes+1:apertureInfo.totalNumOfShapes+apertureInfo.totalNumOfLeafPairs);

