

for ii = 2
    
    clear wvf oi
    
    fH1 = figure(ii); set(fH1, 'Position', [1215, 288, 1092, 1057]); clf; hold all;
    fH2 = figure(ii+10); set(fH2, 'Position', [1215, 288, 1092, 1057]);clf; hold all;
    
    subplotIdx = 1;
    

        
        for y = [-5,0,5]
            
            for x = [-5,5]
            
            location = [x,y];
            

            [wvf, oi] = wvfLoadWavefrontOpticsData('wvfZcoefsSource', 'Polans2015', 'jIndex', 3:14, 'whichEye', 'right', 'eccentricity', location, 'whichGroup', ii, 'verbose', false);
            
            
            psfSupport = wvfGet(wvf, 'spatial Support', 'deg');
            
            set(0, 'CurrentFigure', fH1)
            subplot(3,3,subplotIdx)
            xl = [min(psfSupport) max(psfSupport)];
            imagesc(xl,xl, wvf.psf{1}); axis square; colormap gray
            xlabel('X (minutes)')
            ylabel('Y (minutes)')
            title(sprintf('Location X: %d, Y: %d', x,y))
            
            % Get otf support
            otfSupport = wvfGet(wvf, 'otfSupport', 'deg');
            
            % Get the  OTF
            otf = wvfGet(wvf,'otf');
            
            set(0, 'CurrentFigure', fH2)
            % % Plot observers OTF
            subplot(3,3,subplotIdx)
            xl = [min(otfSupport) max(otfSupport)];
            imagesc(xl,xl, abs(fftshift(otf))); axis square; colormap gray
            xlabel('X (c/d)')
            ylabel('Y (c/d)')
            title(sprintf('Location X: %d, Y: %d', x,y))
            
            subplotIdx = subplotIdx+1;
            
        end
    end
    
   hgexport(fH1, sprintf('~/Desktop/subject%d_Polans_PSF.eps',ii))
   hgexport(fH2, sprintf('~/Desktop/subject%d_Polans_OTF.eps',ii))
 
    
end
