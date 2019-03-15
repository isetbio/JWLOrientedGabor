function [cMosaic, cparams] = getConeMosaic(eccen, expParams, sparams, spatialDensity)

% ----- CONE MOSAIC -----------------------------------------
    % Update ratio of cone types in mosaic if requested:
    if nargin < 4
        spatialDensity = [0 0.6 0.3 0.1]; % Relative density of cone types, K, L, M, S. 
    else
        if expParams.verbose; fprintf('(%s): LMS ratio was reset to %1.1f:%1.1f:%1.1f\n', mfilename, spatialDensity(2),spatialDensity(3),spatialDensity(4)); end
    end

    % Make CONE MOSAIC for a given eccentricity and polar angle
    whichEye = 'left';
    
    % Specify retinal location where stimulus is presented
    cparams.eccentricity      = eccen;             % Visual angle of stimulus center, in deg
    cparams.polarAngle        = expParams.polarAngle; % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
    
    % Cone mosaic field of view in degrees
    cparams.cmFOV             = sparams.sceneFOV; % degrees
    
    % Compute x,y position in m of center of retinal patch from ecc and angle
    [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
    x = x * expParams.deg2m;  y = y * expParams.deg2m;
    
    % Create coneMosaic for particular location and eye
    cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye, 'spatialDensity', spatialDensity);
    
    % Set the field of view (degrees)
    cMosaic.setSizeToFOV(cparams.cmFOV);
    
    % Add photon noise
    cMosaic.noiseFlag = expParams.cparams.noise; % 'random' 'frozen' 'none'
    
    % CURRENT: Set outer segment to be computed with linear filters
    cMosaic.os = osLinear;
    
return

