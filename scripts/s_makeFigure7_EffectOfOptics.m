%% s_makeFigure7_EffectOfOptics

% Script to plot figure 6 (effect of defocus on computational observer
% performance


% Requires computational model to be ran with 'defocus' condition, 
% absorptions to be classified and averaged across iterations.


%% Visualize different levels of defocus (Figure 7A)
s_plotMTF_defocusLevels


%% Plot effect of eye movements on psychometric function (Figure 7B)
plotPsychometricFunctions('defocus', false)