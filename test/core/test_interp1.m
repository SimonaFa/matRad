function test_suite = test_interp1

test_functions=localfunctions();

initTestSuite;

% Test basic values
function test_matrad_interp1_values
    %R = realmax; % For R = realmax, the test may often fail because of Inf in matRad_interp1: should first divide and then multiply.
    R = 10^100;

    % Pick x1 and y1 values
    x1el1 = (2*rand - 1)*R;
    if x1el1 < 0
        x1el2 = x1el1 + rand*(R);
    else
        x1el2 = x1el1 + rand*(R - x1el1);
    end
    
    if x1el1<=x1el2
        x1 = [x1el1; x1el2];
    else
        x1 = [x1el2; x1el1];
    end

    y1el1 = (2*rand - 1)*R;
    if y1el1 < 0
        y1el2 = y1el1 + rand*(R);
    else
        y1el2 = y1el1 + rand*(R - y1el1);
    end
    
    y1 = [y1el1; y1el2];

    x2 = rand*(x1(2)-x1(1)) + x1(1);

    y2 = matRad_interp1(x1, y1, x2);
    expectedy2 = y1(1) + (x2 - x1(1))*((y1(2) - y1(1))/(x1(2) - x1(1)));

    if ~(y2 == expectedy2)
        ciao = 1;
    end

    assertTrue(~isnan(y2));
    assertElementsAlmostEqual(y2, expectedy2);

    y1 = flip(y1);
    y2 = matRad_interp1(x1, y1, x2);
    expectedy2 = y1(1) + (x2 - x1(1))*(y1(2) - y1(1))/(x1(2) - x1(1));
    assertTrue(~isnan(y2));
    assertElementsAlmostEqual(y2, expectedy2);

    x2 = x1(2) + rand*(R - x1(2));
    if x1(2) == R
        assertTrue(~isnan(matRad_interp1(x1, y1, x2)));
    else
        assertTrue(isnan(matRad_interp1(x1, y1, x2)));
    end

    x2 = x1(1) + rand*(-R - x1(1));
    if x1(1) == -R
        assertTrue(~isnan(matRad_interp1(x1, y1, x2)));
    else
        assertTrue(isnan(matRad_interp1(x1, y1, x2)));
    end

    
% Test if conditions: extrapolation
function test_matRad_interp1_extrapolation
    R = 10^100;
    realExtrap = R*(2*rand-1);
   
    % Pick x1 and y1 values
    x1el1 = (2*rand - 1)*R;
    if x1el1 < 0
        x1el2 = x1el1 + rand*(R);
    else
        x1el2 = x1el1 + rand*(R - x1el1);
    end
    if x1el1<=x1el2
        x1 = [x1el1; x1el2];
    else
        x1 = [x1el2; x1el1];
    end

    y1el1 = (2*rand - 1)*R;
    if y1el1 < 0
        y1el2 = y1el1 + rand*(R);
    else
        y1el2 = y1el1 + rand*(R - y1el1);
    end
    
    y1 = [y1el1; y1el2]; 

    % Pick 2 values outside boundaries
    xLow = x1(1) + rand*(-R - x1(1));
    xUp = x1(2) + rand*(R - x1(2));
    x2 = rand*(x1(2)-x1(1)) + x1(1);
    x2 = [xLow; x2; xUp];
    
    outIdx = find(x2<x1(1) | x2>x1(2));
    % real number
    y2 = matRad_interp1(x1, y1, x2, realExtrap);
    assertElementsAlmostEqual(y2(outIdx), realExtrap.*ones(size(y2(outIdx))));
    % NaN & 'none'
    y2 = matRad_interp1(x1, y1, x2, NaN);
    y3 = matRad_interp1(x1, y1, x2, 'none');
    assertTrue(~any(~isnan(y2(outIdx))));
    assertElementsAlmostEqual(y2, y3);
    % linear & extrap
    y2 = matRad_interp1(x1, y1, x2, 'extrap');
    y3 = matRad_interp1(x1, y1, x2, 'linear');
    assertElementsAlmostEqual(y2, y3);
    
    
function test_matRad_interp1_multiple1D
    R = 10^100;
    x1 = (2.*rand(1000, 1)- 1).*R ;
    x1 = unique(x1);
    x1 = sort(x1);
    y1 = (2.*rand(size(x1))- 1).*R;



function test_matRad_interp1_errors
    x1 = [1; 1];
    y1 = [2; 3];
    x2 = [2.5, 6];
    try
        y2 = matRad_interp1(x1, y1, x2);
    catch
        assertTrue(true);
    end

    % If x1 has repetitions, and x2 is a vector, error is correct
    % (sample points must be unique)
    % if x2 is a scalar, the result is NaN
    y2 = matRad_interp1(x1, y1, 1);
    assertTrue(isnan(y2));
    y2 = matRad_interp1(x1, y1, -20);
    assertTrue(isnan(y2));



