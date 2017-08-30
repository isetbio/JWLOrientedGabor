% Create Thibos virtual Zernike coeffs for 100 virtual subjects


eccen = 15;
pupilsize = 4.5;

Zcoeffs = VirtualEyes(100,pupilsize);


for n = 1:100
    
        % For each eccentricity (81 total)
        for ec = 1:eccen
            
            % Create human default wvf
            wvf = wvfCreate;
            
            % Set pupil size to the one reported in Jaeken and Artal (2012)
            wvf = wvfSet(wvf, 'measuredpupilsize', pupilsize);
            wvf = wvfSet(wvf, 'calcpupilsize', pupilsize);
            
            % Create ZCoefs for all 15 zernike coeffs at once
            thisZcoeffs = Zcoeffs(:,n);
            
            wvf =  wvfSet(wvf,'zcoeffs',Zcoeffs);
            
            % compute the PSF
            wvf = wvfComputePSF(wvf);
            psf(:,:,ec,eye,p) = wvfGet(wvf, 'psf');
            
            % Get psf support
            psfSupport = wvfGet(wvf, 'spatial Support', 'min');
            
            % Plot observers PSF
            if verbose
                set(0, 'currentfigure', fh1);
                subplot(4,21,ec);
                xl = [min(psfSupport) max(psfSupport)];
                
                imagesc(xl,xl, psf(:,:,ec,eye,p)); axis square; colormap gray
                xlabel('X (minutes)')
                ylabel('Y (minutes)')
                %             title(sprintf('PSF piston coeff: 1, eccen: %1.2f', eccen(ii)))
                title(sprintf('eccen: %2.0f', eccenLabel(ec)))
            end
            
            % Get otf
            otfSupport = wvfGet(wvf, 'otfSupport', 'deg');
            
            % Get the the OTF
            otf(:,:,ec,eye,p) = wvfGet(wvf,'otf');
            
            if verbose
                % Plot observers OTF
                set(0, 'currentfigure', fh2);
                subplot(4,21,ec);
                xl = [min(otfSupport) max(otfSupport)];
                imagesc(xl,xl, abs(otf(:,:,ec,p))); axis square; colormap gray
                xlabel('X (c/d)')
                ylabel('Y (c/d)')
                %     title(sprintf('OTF piston coeff: 1, eccen: %1.2f', eccen(ii)))
                title(sprintf('eccen: %2.0f', eccenLabel(ec)))
            end
            
        end
        
        %         hgexport(fh1, fullfile(proj_path, sprintf('psf_patient%d_eye%d.eps',p,eye)))
        %         hgexport(fh2, fullfile(proj_path, sprintf('otf_patient%d_eye%d.eps',p,eye)))
        
        
    end
end

% PLOT MEAN OTF EYE1
figure(3); clf; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% mean OTf across observers
otf_mn = mean(otf, 5);

for ec = 1:eccen
    subplot(4,21,ec);
    [xSfGridCyclesDeg, ySfGridCyclesDeg] = meshgrid(otfSupport);
    xl = [min(otfSupport) max(otfSupport)];
    imagesc(xl, xl, abs(otf_mn(:,:,ec,1))); axis square; colormap gray;
    xlabel('X (c/d)')
    ylabel('Y (c/d)')
    title(sprintf('eccen: %2.0f', eccenLabel(ec)))
    
end