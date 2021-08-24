function expParams = hpcArrayID2eccen(hpcArrayID, expParams)
% function to use hpcArrayID to select single eccentricity to simulate
% this is helpful when wanting to run multiple jobs in parallel on the High
% Performance Cluster (HPC)
    expParams.eccentricities = expParams.eccentricities(hpcArrayID);

end

