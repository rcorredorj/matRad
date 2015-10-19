function indices = matRad_calcFluenceMapsQuality( w, stf, visBool )
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate some indices describing the modulation complexity of the
% fluences maps. 
%    - Modulation index (MI): as described by Webb [1] and extended to 2D 
%      by Giorgia et al. [2] 
% 
% call
%   indices = matRad_calcFluenceMapsQuality(resultSequencing, visBool)
%
% input
%   w:                  weights
%   stf:                matRad steering information struct                    
%   visBool:            toggle on/off visualization (optional)
%
% output
%   indices             modulation indices
%
% References
%   [1] http://www.ncbi.nlm.nih.gov/pubmed/11474945
%   [2] http://www.ro-journal.com/content/2/1/42
% 
% Developer: Ricardo A. Corredor @rcorredorj
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indices = {};

% if visBool not set toogle off visualization
if nargin < 2
    visBool = 0;
end

numOfBeams = numel(stf);

offset = 0;

for i = 1:numOfBeams
    
    numOfRaysPerBeam = stf(i).numOfRays; 
    
    % get relevant weights for current beam
    wOfCurrBeams = w(1+offset:numOfRaysPerBeam+offset);
    
    X = ones(numOfRaysPerBeam,1)*NaN;
    Z = ones(numOfRaysPerBeam,1)*NaN;
        
    for j=1:stf(i).numOfRays
      X(j) = stf(i).ray(j).rayPos_bev(:,1);
      Z(j) = stf(i).ray(j).rayPos_bev(:,3);
    end
        
    % sort bixels into matrix
    minX = min(X);
    maxX = max(X);
    minZ = min(Z);
    maxZ = max(Z);
    
    dimOfFluenceMxX = (maxX-minX)/stf(i).bixelWidth + 1;
    dimOfFluenceMxZ = (maxZ-minZ)/stf(i).bixelWidth + 1;
    
    %Create the fluence matrix.
    fluenceMx = zeros(dimOfFluenceMxZ,dimOfFluenceMxX);
    
    % Calculate X and Z position of every fluence's matrix spot
    % z axis = axis of leaf movement!
    xPos = (X-minX)/stf(i).bixelWidth+1;
    zPos = (Z-minZ)/stf(i).bixelWidth+1;
    
    % Make subscripts for fluence matrix
    indInFluenceMx = zPos + (xPos-1)*dimOfFluenceMxZ;
    
    %Save weights in fluence matrix.
    fluenceMx(indInFluenceMx) = wOfCurrBeams;
    m = size(fluenceMx,1);
    n = size(fluenceMx,2);
    
    sd = std(fluenceMx(:));
    
    % prepare sequencer
    % calFac = max(fluenceMx(:));
    
    diffX = conv2(fluenceMx,[1 -1],'same');
    diffY = conv2(fluenceMx,[-1 1]','same'); % considering that Y increments towards the top of the screen
    diffXY = conv2(fluenceMx,[0 -1; 1 0],'same'); % considering that Y increments towards the top of the screen
    
    diffX = abs(diffX(1:m-1,1:n-1));
    diffY = abs(diffY(1:m-1,1:n-1));
    diffXY = abs(diffXY(1:m-1,1:n-1));



    res=0.01;
    fXVals = 0:res:2;
    std_F = sd * fXVals;

    zX = zeros(numel(fXVals),1);
    zY = zeros(numel(fXVals),1);
    zXY = zeros(numel(fXVals),1);
    for k = 1:numel(std_F)
        zX(k) = (1/((m-1)*(n-1)))*sum(diffX(:) > std_F(k));
        zY(k) = (1/((m-1)*(n-1)))*sum(diffY(:) > std_F(k));
        zXY(k) = (1/((m-1)*(n-1)))*sum(diffXY(:) > std_F(k));
    end

    Zf = (zX + zY + zXY)/3;

    factors = [0.1, 0.3, 0.5, 0.6, 0.8, 1.0];
    MI = zeros(numel(factors),1);
    for k = 1:numel(factors)
        halfPos = round(factors(k)/res) + 1 ;
        MI(k) = sum(Zf(1:halfPos));
    end
    indices.MI.beam{i} = MI;
    
    if visBool
       figure;
       colormap jet;
       subplot(1,2,1);
       imagesc(fluenceMx);
       title( ['Beam No. ' num2str(i) ' - MI (0.5)=' num2str(MI(3))]);
       subplot(1,2,2);
       plot(factors,MI,'-o');
       title('MI at different factors');
    end
     
end
if visBool
    
    figure;
    vals = [indices.MI.beam{:}]; 
    vals = reshape(vals,numel(indices.MI.beam),numel(vals)/numel(indices.MI.beam));
    boxplot(vals,'plotstyle','compact')
    
end


end

