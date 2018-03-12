function oi = oiDefocus(defocus)
% OIDEFOCUS - Create an human wvf oi with some amount of defocus
%
% Input:
%   Zernike defocus coefficient (in microns)
% 
% Output:
%  oi - Human wvf structure with the defocus set
%
% See also:  wvfDefocusDioptersToMicrons and wvfDefocusDioptersToMicrons,
%            oiCreate
%
% JW/EK ISETBIO Team, 2017

pupilMM = 3;
zCoefs = wvfLoadThibosVirtualEyes(pupilMM);
wave = 400:10:700; wave = wave(:);

wvfP = wvfCreate('wave',wave,'zcoeffs',zCoefs,'name',sprintf('human-%d',pupilMM));
wvfP = wvfSet(wvfP,'calc pupil size',pupilMM);
wvfP = wvfSet(wvfP,'zcoeffs',defocus,{'defocus'});
wvfP = wvfComputePSF(wvfP);


% % Look at PSF and OTF 
% % Plot the PSF
% [u,p,f] = wvfPlot(wvfP,'2d psf space','um',550);
% set(gca,'xlim',[-20 20],'ylim',[-20 20],'zlim',[0 0.03],'FontSize',16,'LineWidth',2);
% 
%  % Get the support and OTF
% otfSupport = wvfGet(wvfP, 'otfSupport', 'mm',550);
% otf = wvfGet(wvfP,'otf');
% % set(gcf,'LineWidth',5)
% 
%  % Plot observers OTF
% wvfPlot(wvfP,'2d OTF','mm',550);
% set(gca,'xlim',[-400 400] ,'ylim',[-400 400], 'zlim', [0 1.2],'FontSize',16,'LineWidth',2);


oi = wvf2oi(wvfP);
oi = oiSet(oi,'optics lens',Lens);
oi = oiSet(oi,'optics model', 'shiftInvariant');
oi = oiSet(oi, 'oi name', sprintf('human-wvf-%d',round(10*defocus)));

end
