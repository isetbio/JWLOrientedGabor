function [expName, saveFolder, seed] = prepHPC(taskID)

% function to run different experiments with different seeds on NYU's HPC



switch taskID
    
    case {1, 2, 3, 4, 5}
        expName   = 'defocus';
        saveFolder = sprintf('run%d', taskID);
        seed      = taskID;
        
    case {6, 7, 8, 9, 10}
        expName   = 'conedensity';
        saveFolder = sprintf('run%d', taskID-5);
        seed      = taskID-5;
        
    case {11, 12, 13, 14, 15}
        expName   = 'eyemov';
        saveFolder = sprintf('run%d', taskID-10);
        seed      = taskID-10;
        
    case {16, 17, 18, 19, 20}
        expName   = 'conetypes';
        saveFolder = sprintf('run%d', taskID-15);
        seed      = taskID-15;
        
    case {21, 22, 23, 24, 25}
        expName   = 'conetypesmixed';
        saveFolder = sprintf('run%d', taskID-20);
        seed      = taskID-20;
        
    case {26, 27, 28, 29, 30}
        expName   = 'defaultnophaseshift';
        saveFolder = sprintf('run%d', taskID-25);
        seed      = taskID-25;
        
    case {31, 32, 33, 34, 35}
        expName   = 'eyemovnophaseshift';
        saveFolder = sprintf('run%d', taskID-30);
        seed      = taskID-30;
end



return