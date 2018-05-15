% s_createAbsorptionVideo.m

% Get video name and create a video object
vname = '~/Desktop/absorptions110_C0.10';
vObj = VideoWriter(vname);

% Define frame rate
vObj.FrameRate = 30;

% Open video object
open(vObj);

% Plot absorptions (one frame = one time step) and write into video object
hFig = figure;
maxRate = max(absorptions(:));
colormap(gray(maxRate));
for ii=1:size(absorptions,4)
    
    imagesc(squeeze(absorptions(1,:,:,ii,1))); set(gca,'CLim',[0 maxRate]);
    axis image; axis off; box off; drawnow;

    F = getframe; writeVideo(vObj,F);
end

% Save and close the video file
writeVideo(vObj,F);
close(vObj);
