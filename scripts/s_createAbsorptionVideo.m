vname = '~/Desktop/absorptions220';
vObj = VideoWriter(vname);
vObj.FrameRate = 30;
open(vObj);


hFig = figure;
colormap(gray(max(absorptions(:))));
for ii=1:size(absorptions,4)
    
    imagesc(squeeze(absorptions(1,:,:,ii,1))); axis image; axis off; box off; drawnow;
    
    
    
    F = getframe; writeVideo(vObj,F);
end

% Write the video object if save is true

writeVideo(vObj,F);
close(vObj);
