%% s_makeFigure6_EffectOfEyemovements

% Script to plot figure 6 (effect of eyemovements on computational observer
% model)

% Requires computational model to be ran with 'eyemov' condition, 
% absorptions to be classified and averaged across iterations.


%% Visualize eye movement paths of 5 trials (Figure 6A)
s_VisualizeEyeMovements


%% Plot effect of eye movements on psychometric function (Figure 6b)
plotPsychometricFunctions('eyemov', false)