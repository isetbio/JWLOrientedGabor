%% s_calculateDensityFromConeMosaic

whichEye          = 'left';
polarAngleDeg     = [0,90,180,270];
cparams.cmFOV     =  2; % degrees
eccentricities     = 0:40;

% Convertion deg to m
deg2m  = 1/3 * 0.001; % 3 deg per mm, .001 mm per meter

% Predefine density vector
allDensity = size(eccentricities);

for pa = polarAngleDeg
    for eccen = eccentricities
        
        % Specify retinal location where stimulus is presented
        cparams.eccentricity = eccen;             % Visual angle of stimulus center, in deg
        cparams.polarAngle   = deg2rad(pa);   % Polar angle (radians): 0 is right, pi/2 is superior, pi is left, 3*pi/2 inferior
        
        % Compute x,y position in m of center of retinal patch from ecc and angle
        [x, y] = pol2cart(cparams.polarAngle, cparams.eccentricity);
        x = x * deg2m;  y = y * deg2m;
        
        cMosaic = coneMosaic('center', [x, y], 'whichEye', whichEye);
        
        % Set the field of view (degrees)
        cMosaic.setSizeToFOV(cparams.cmFOV);
        
        allDensity(eccen+1,pa==polarAngleDeg,:) = eccen2density(cMosaic, 'm');
        
    end
end

%% Visualize

figure('Color','w'); 
plot(repmat(eccentricities,4,1)',allDensity, 'o-', 'LineWidth',2,'MarkerSize',10);
xlabel('Eccentricity (deg)'); ylabel('Cone Density per m2')
set(gca,'TickDir', 'out','TickLength',[0.015 0.015]); box off;
legend({'Polar Angle 0','Polar Angle 90','Polar Angle 180','Polar Angle 270'})
