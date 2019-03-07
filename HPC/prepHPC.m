function [expName, subFolder, seed] = prepHPC(taskID)

% function to run different experiments with different seeds on NYU's HPC



switch taskID
    
    case {1, 2, 3, 4, 5}
        expName   = 'defocus';
        subFolder = sprintf('run%d', taskID);
        seed      = taskID;
        
    case {6, 7, 8, 9, 10}
        expName   = 'conedensity';
        subFolder = sprintf('run%d', taskID-5);
        seed      = taskID-5;
        
    case {11, 12, 13, 14, 15}
        expName   = 'eyemov';
        subFolder = sprintf('run%d', taskID-10);
        seed      = taskID-10;
        
    case {16, 17, 19, 19, 20}
        expName   = 'conetypes';
        subFolder = sprintf('run%d', taskID-15);
        seed      = taskID-15;
        
    case {21, 22, 23, 24, 25}
        expName   = 'conetypesMixed';
        subFolder = sprintf('run%d', taskID-20);
        seed      = taskID-20;
end



return