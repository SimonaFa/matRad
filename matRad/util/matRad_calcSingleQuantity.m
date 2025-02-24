function cubeQuantity = matRad_calcSingleQuantity(dij, w, quantityStruct, scenNum)
    
    % Dummy way to access quantity class dijField
    obj = feval(quantityStruct.className);
    obj.useScenarios(1)= scenNum;
    %idxDijField = find(strcmp(fieldnames(dij), obj.dijField));

    %resultGUI = obj.matRad_calcQuantityCubes(dij, resultGUI, scenNum, beamInfo);
    try
        cubeQuantity = obj.computeQuantity(dij, scenNum, w);
        %cubeQuantity = obj.getResult(dij,w);
    catch
    end
end