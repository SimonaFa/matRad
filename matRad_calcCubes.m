function resultGUI = matRad_calcCubes(w,dij,cst,scenNum)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matRad computation of all cubes for the resultGUI struct which is used
% as result container and for visualization in matRad's GUI
%
% call
%   resultGUI = matRad_calcCubes(w,dij,cst)
%
% input
%   w:       bixel weight vector
%   dij:     dose influence matrix
%   cst:     matRad cst struct
%   scenNum: optional: number of scenario to calculated (default 1)
%
% output
%   resultGUI: matRad result struct
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

if nargin < 4
    scenNum = 1;
end

resultGUI.w = w;

% calc dose and reshape from 1D vector to 2D array
numOfBeams = max(unique(dij.beamNum));

for i = 1:(numOfBeams + 1)
    if i == (numOfBeams + 1)
      beamInfo(i).suffix = '';
      beamInfo(i).isCurrBeam = ones(size(dij.beamNum));
    else
      beamInfo(i).suffix = ['Beam', num2str(i)];
      beamInfo(i).isCurrBeam = (dij.beamNum == i);      
    end
end

for i = 1:length(beamInfo)
    resultGUI.(['physicalDose', beamInfo(i).suffix]) = reshape(full(dij.physicalDose{scenNum} * (resultGUI.w .* beamInfo(i).isCurrBeam)),dij.dimensions);
end

% consider RBE for protons
if isfield(dij,'RBE')
   fprintf(['matRad: applying a constant RBE of ' num2str(dij.RBE) ' \n']);
   for i = 1:length(beamInfo)
        resultGUI.(['RBExDose', beamInfo(i).suffix]) = resultGUI.(['physicalDose', beamInfo(i).suffix]) * dij.RBE;
   end
end

% consider VOI priorities
[cst,resultGUI.overlapCube]  = matRad_setOverlapPriorities(cst,dij.dimensions);

if isfield(dij,'mLETDose')
    for i = 1:length(beamInfo)
        LETDoseCube                      = dij.mLETDose{scenNum} * (resultGUI.w .* beamInfo(i).isCurrBeam);
        resultGUI.(['LET', beamInfo(i).suffix]) = zeros(dij.dimensions);
        ix                               = resultGUI.(['physicalDose', beamInfo(i).suffix]) > 0;
        resultGUI.(['LET', beamInfo(i).suffix])(ix) = LETDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);
    end
end

% consider biological optimization for carbon ions
if isfield(dij,'mAlphaDose') && isfield(dij,'mSqrtBetaDose')
    for i = 1:length(beamInfo)

    a_x = zeros(size(resultGUI.(['physicalDose', beamInfo(i).suffix])));
    b_x = zeros(size(resultGUI.(['physicalDose', beamInfo(i).suffix])));

    for j = 1:size(cst,1)
        % Only take OAR or target VOI.
        if isequal(cst{j,3},'OAR') || isequal(cst{j,3},'TARGET') 
            a_x(cst{j,4}{scenNum}) = cst{j,5}.alphaX;
            b_x(cst{j,4}{scenNum}) = cst{j,5}.betaX;
        end
    end
    
    ix = b_x~=0;
    
    resultGUI.(['effect', beamInfo(i).suffix]) = full(dij.mAlphaDose{scenNum}*(resultGUI.w .* beamInfo(i).isCurrBeam)+(dij.mSqrtBetaDose{scenNum}*(resultGUI.w .* beamInfo(i).isCurrBeam)).^2);
    resultGUI.(['effect', beamInfo(i).suffix]) = reshape(resultGUI.(['effect', beamInfo(i).suffix]),dij.dimensions);
    
    resultGUI.(['RBExDose', beamInfo(i).suffix])     = zeros(size(resultGUI.(['effect', beamInfo(i).suffix])));
    resultGUI.(['RBExDose', beamInfo(i).suffix])(ix) = ((sqrt(a_x(ix).^2 + 4 .* b_x(ix) .* resultGUI.(['effect', beamInfo(i).suffix])(ix)) - a_x(ix))./(2.*b_x(ix)));
                                 
    resultGUI.(['RBE', beamInfo(i).suffix])          = resultGUI.(['RBExDose', beamInfo(i).suffix])./resultGUI.(['physicalDose', beamInfo(i).suffix]);
   
    resultGUI.(['alpha', beamInfo(i).suffix])     = zeros(size(resultGUI.(['effect', beamInfo(i).suffix])));
    resultGUI.(['beta', beamInfo(i).suffix])      = zeros(size(resultGUI.(['effect', beamInfo(i).suffix])));
    AlphaDoseCube       = full(dij.mAlphaDose{scenNum} * (resultGUI.w .* beamInfo(i).isCurrBeam));
    resultGUI.(['alpha', beamInfo(i).suffix])(ix) = AlphaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix);
    SqrtBetaDoseCube    = full(dij.mSqrtBetaDose{scenNum} * (resultGUI.w .* beamInfo(i).isCurrBeam));
    resultGUI.(['beta', beamInfo(i).suffix])(ix)  = (SqrtBetaDoseCube(ix)./resultGUI.(['physicalDose', beamInfo(i).suffix])(ix)).^2;
  end
end
