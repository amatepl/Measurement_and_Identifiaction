function [umat,ymat]=ReadDataLab2(N,Nrep,Drep,FileName)
% Starting from a measurement with multiple repetitions, it combines all
% the measured input signals into one matrix, and all the measured output
% signals into another matrix. 

% Example: [umat,ymat] = ReadDataLab2(8000,40,32,'Group1_Output1.mat')

% Input values:
    % N: Length of one repetition (Ex: length(Su))
    % Nrep: How many repetitions were measured
    % Drep: How many repetitions you want to retain (it selects the repetitions 'backwards', to avoid retaining those with transient)
    % Filename: Name of the matfile where the measurement is stored
    
% Output values:
    % umat: Matrix with the selected repetitions of the input. Size:[NxDrep]
    % ymat: Matrix with the selected repetitions of the output. Size:[NxDrep]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outData = load([FileName]);

umat_allrep = zeros(N,Nrep-1);
ymat_allrep = zeros(N,Nrep-1);

for iRep=0:Nrep-1
    if iRep == 0 % The first repetition is named with no index after Su/Sy
        NameVar_u = 'Su';
        NameVar_y = 'Sy';
    else % The second repetition is named Su1/Sy1, etc..
        NameVar_u = sprintf('Su%d',iRep);
        NameVar_y = sprintf('Sy%d',iRep);
    end
    umat_allrep(:,iRep+1) = getfield(outData, NameVar_u);    
    ymat_allrep(:,iRep+1) = getfield(outData, NameVar_y);
end



umat = umat_allrep(:,Nrep-Drep:end);
ymat = ymat_allrep(:,Nrep-Drep:end);

end
