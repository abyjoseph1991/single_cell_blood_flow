% Velocity_Step4_Radon
% Author: Aby Joseph, University of Rochester
% License: GPL-3.0-or-later
% Last modified: 01-24-2019

clear;close all;
f_No=0;
fontsize1=18;

% Variables typically not changed
moveFileAfterProcessing=1;
radonMode=1; % Mode #1 - Use this. Conventional "peak variance of radon" approach. Mode #2 - "peak intensity of radon" approach.
shuffleON=1; % toggle button for subtracting "shuffled" radon-variance from original radon to improve slope detection
calculateSNR1=1; % Recommended value=1. Switch to 0 if you want to speed up computation
calculateSNR2=1; % Recommended value=1. Switch to 0 if you want to speed up computation
includeOtherQuadrantON=1; % includes the entire 180 degress in the radon angle search space. The density of angles in quadrant without streaks is set low (typically 1 degree step size).
circularMaskingON=1;
radonClippingON=0;
featureDetectionON=0; % if OFF, conventional array of equal sized kernels will be used
findQuadrant=0;

% Variables typically changed
videoDurationToAnalyze_Step4=10; % in seconds
fractionalOverlapTime=4; % If  equal to 4, 75% overlap in gliding boxes. If equal to 2, 50% overlap etc. Must be integer.
fractionalOverlapSpace=4;
dxMicrons=15; % in microns % REVISIT
nShuffles=3;
% quadrant=1;
Vstep=0.1; % in mm/s
Vmin=0; % in mm/s
Vmax=100; % in mm/s

Step5folder='E:\PhD\TRBF\Experimental DATA\Step 5 MAT files\';
Step4folder='E:\PhD\TRBF\Experimental DATA\Step 4 MAT files\';
Step3folder='E:\PhD\TRBF\Experimental DATA\Step 3 MAT files\';
Step2folder_Processed='E:\PhD\TRBF\Experimental DATA\Step 2 MAT files\PROCESSED\';
RawDir_Processed='E:\PhD\TRBF\Experimental DATA\RAW MAT files\PROCESSED\';
cd(Step3folder);

% Append current date and time to filenames of results
currentDateAndTime=clock;
year=num2str(currentDateAndTime(1));
month=num2str(currentDateAndTime(2),'%02d');
day=num2str(currentDateAndTime(3),'%02d');
hour=num2str(currentDateAndTime(4),'%02d');
minute=num2str(currentDateAndTime(5),'%02d');
appendToResults=['_',year,month,day,hour,minute];

namecontains='*_Step3*.mat*';
[dirFilenames(1).name,filePath,~]=uigetfile(namecontains,...
    'Select file to process - click cancel to process all files in parent folder');
if filePath==0
    dirFilenames=dir([namecontains]);
    nFiles=numel(dirFilenames);
    processDir=Step3folder;
else
    processDir=[filePath];
    nFiles=1;
end
screenSize=get(0, 'MonitorPositions');

cd(processDir);
for fileNo=1:nFiles
    tic;
    f_No=0;
    fileName=dirFilenames(fileNo).name;
    display(fileName);
    
    % Loading RAW.mat
    namecontains2=[fileName(1:25),'*RAW*.mat'];
    dirFilenames2=dir([RawDir_Processed,namecontains2]);
    nFiles2=numel(dirFilenames2);
    if nFiles2==1
        load([RawDir_Processed,dirFilenames2(1).name],'desinusoidON','postFPGA',...
            'linescanDir','nFrames','height','width','FOV',...
            'mic_pix','freq');
    else
        [RawFileName,RawFilePath,~]=uigetfile([RawDir_Processed,...
            namecontains2],'Select RAW MAT file');
        load([RawFilePath,RawFileName],'desinusoidON','postFPGA',...
            'linescanDir','nFrames','height','width','FOV',...
            'mic_pix','freq');
    end

    % Loading Step2.mat
    namecontains3=[fileName(1:25),'*Step2*.mat'];
    dirFilenames3=dir([Step2folder_Processed,namecontains3]);
    nFiles3=numel(dirFilenames3);
    if nFiles3==1
        load([Step2folder_Processed,dirFilenames3(1).name],...
            'insertNaNsWhenPoorRegistration',...
        'crossCorrelationCoeffcientThreshold','referenceLine',...
        'timeWindowAveraging','sizeStripRegistration',...
        'dataRegisteredWithStaticLines','line_avg0','line_std0',...
        'line_avg1','line_std1','timeAxisForMotionTrace','motionTrace',...
        'data','analysisBoundaries','videoDurationToAnalyze_Step2','quadrant');
    else
        [Step2FileName,Step2FilePath,~]=uigetfile([Step2folder_Processed,...
            namecontains3],'Select Step2 MAT file');
        load([Step2FilePath,Step2FileName],'insertNaNsWhenPoorRegistration',...
        'crossCorrelationCoeffcientThreshold','referenceLine',...
        'timeWindowAveraging','sizeStripRegistration',...
        'dataRegisteredWithStaticLines','line_avg0','line_std0',...
        'line_avg1','line_std1','timeAxisForMotionTrace','motionTrace',...
        'data','analysisBoundaries','videoDurationToAnalyze_Step2','quadrant');
    end
    
    % Loading Step3.mat
    load([processDir,fileName],'vesselAngle');
    
    if isnan(videoDurationToAnalyze_Step4)
        videoDurationToAnalyze_Step4=size(data,2)/freq;
    end    

    data=data(:,1:floor(freq*videoDurationToAnalyze_Step4));
    data=single(data);
    data=data./255;
    dataWidth_Step4=size(data,2);
    dataHeight_Step4=size(data,1);

    NANthreshold=100; % Maximum allowed percetage occurrance of NaNs along the time dimension for a given radial coordinate to be considered to be inside the vessel lumen.
    dx=round(dxMicrons/mic_pix); %kernel size in pixels

    dx=dx+(1-mod(dx,2));
    d=(dx-1)/2; % defines half-length of lines overlayed in space-time image (to visualize calculated slopes)
    
    dtStep=round((dx-1)/fractionalOverlapTime); 
    dtStepSeconds=dtStep./freq;
    maxExpectedHeartRate=400; % in beats per minute. For anesthetized mouse
    measuredSamplesPerCardiacCycle=round((1/(maxExpectedHeartRate/60))/...
        dtStepSeconds);
    dxStep=round((dx-1)/fractionalOverlapSpace);
    dxStepMicrons=dxStep*mic_pix;
    nShufflesQD=3;

    % creating circular mask
    mask=zeros(dx);
    for k=1:dx
        for m=1:dx
            if (sqrt((k-0.5-dx/2)^2+(m-0.5-dx/2)^2)<=dx/2)
                mask(k,m)=1;
            end
        end
    end
    if circularMaskingON==0
        mask=ones(dx);
    end
    
    %Quadrant determination (QD) for Radon angle search space --> quadrant=1 for 0 to 90 degrees in Radon space, quadrant=2 for 90 to 180 degrees in Radon space (follow MATLAB convention for Radon angles)
    anglesQD=0:179;
    kernelQD=data(analysisBoundaries(1):analysisBoundaries(2),1:1+...
        diff(analysisBoundaries)); % diff(analysisBoundaries) forced to be even number in Step 3, such that kernelQD has odd number of rows and columns
    kernelSquareQD=kernelQD;
    dxQD=size(kernelQD,1);
    maskQD=zeros(size(kernelQD,1),size(kernelQD,2));
    for k=1:dxQD
        for m=1:dxQD
            if (sqrt((k-0.5-dxQD/2)^2+(m-0.5-dxQD/2)^2)<=dxQD/2)
                maskQD(k,m)=1;
            end
        end
    end
    if circularMaskingON==0
        maskQD=ones(dxQD);
    end
    
    if gaussianWindowingON==0
        kernelQD=kernelQD.*maskQD;
    else
        gaussianSigmaQD=floor((dxQD-1)/4); %sigma of gaussian filter calculated such that when filter size is chosen to be exactly equal to kernel size, the relation "FilterSize=2*ceil(2*SIGMA)+1" approximately holds
        hQD=fspecial('gaussian',dxQD,gaussianSigmaQD);
        kernelQD=kernelQD.*hQD;
        kernelQD=kernelQD-min(kernelQD(:));
        kernelQD=kernelQD./max(kernelQD(:));        
    end
    [rQD,xpQD]=radon(kernelQD,anglesQD);
    stddevQD=std(rQD,0,1);
    
    stddevShuffledQD=zeros(size(stddevQD));
    for k=1:nShufflesQD
        shuffledIndices=randperm(size(kernelQD,1));
        kernelSquareShuffledQD=kernelSquareQD(:,shuffledIndices);
        if gaussianWindowingON==0
            kernelShuffledQD=kernelSquareShuffledQD.*maskQD;
        else
            kernelShuffledQD=kernelSquareShuffledQD.*hQD;
            kernelShuffledQD=kernelShuffledQD-min(kernelShuffledQD(:));
            kernelShuffledQD=kernelShuffledQD./max(kernelShuffledQD(:));
        end
        rShuffledQD=radon(kernelShuffledQD,anglesQD);
        stddevShuffledQD=std(rShuffledQD,0,1)+stddevShuffledQD;
    end
    stddevShuffledQD=stddevShuffledQD./nShufflesQD;
    stddevCleanQD=stddevQD-stddevShuffledQD;
    stddevCleanQD_display=stddevCleanQD;
    stddevCleanQD=stddevCleanQD-min(stddevCleanQD(:));
    stddevCleanQD=stddevCleanQD./max(stddevCleanQD(:));
    [~,locsClean]=findpeaks(stddevCleanQD,'SortStr','descend','NPeaks',1);
    SNR_QD=stddevCleanQD(locsClean)./mean(stddevCleanQD(:));
    thetaQD_Variance=anglesQD(locsClean(1));

    [rQDMax,I]=max(rQD(:));
    [I_row,I_col]=ind2sub(size(rQD),I);
    thetaQD_Intensity=anglesQD(I_col);
    
    if findQuadrant==1
        if thetaQD_Variance<=90
            quadrant=1;
        else
            quadrant=2;
        end
    end
    
    if gaussianWindowingON==1
        figure(4);imagesc(hQD);
    end

    % determining angles for radon transform
    V_max_calculation_pixels=abs(cotd(thetaQD_Variance))*(1/...
        abs(cosd(vesselAngle))); %Vx in pixel space
    V_max_calculation=V_max_calculation_pixels*(1e-3)*mic_pix*freq; %Vx in mm/s
    
    Vees=fliplr(Vmin:Vstep:Vmax);
    nVees=numel(Vees);
    if Vmin==0
        angles=zeros(1,2*nVees-1);
        Vees_display=zeros(1,2*nVees-1);
    else
        angles=zeros(1,2*nVees);
        Vees_display=zeros(1,2*nVees);
    end
    angles(1:nVees)=acotd(Vees*(1e3)*abs(cosd(vesselAngle))/(mic_pix*freq)); % angles in degrees
    Vees_display(1:nVees)=-1*Vees;
    anglesComp=fliplr(180-angles(1:nVees));
    if Vmin==0
        angles(nVees:end)=anglesComp;
        Vees_display(nVees:end)=fliplr(Vees);
    else
        angles(nVees+1:end)=anglesComp;
        Vees_display(nVees+1:end)=fliplr(Vees);
    end
    anglesFullRange=angles;
    if quadrant==1
        angles=anglesFullRange(anglesFullRange<=90);
        if includeOtherQuadrantON==1
            angles=[angles 91:179];
        end
    else
        angles=anglesFullRange(anglesFullRange>=90);
        if includeOtherQuadrantON==1
            angles=[1:89 angles];
        end
    end
    
    anglesLinearlySpaced=linspace(angles(1),angles(end),numel(angles));
    centAxis=round(mean(analysisBoundaries)); %REVISIT
    ndt=round((dataWidth_Step4-dx+1)/dtStep);
    
    ndx1=floor(min([centAxis-analysisBoundaries(1) centAxis-dx/2])/dxStep);
    ndx2=floor(min([analysisBoundaries(2)-centAxis dataHeight_Step4-...
        centAxis-dx/2])/dxStep);
    ndx=2*min([ndx1 ndx2])+1;
    
    theta_max90=zeros(ndx*ndt,1);
    if featureDetectionON==0
        theta_maxClean=NaN(ndx*ndt,1);
    end
    reachedEndT=0;
    kernelStartT=1;
    t=0;
    kernelStartXfirst=round(centAxis-((ndx-1)/2)*dxStep-(dx-1)/2);
    pos90=find(abs(angles-90)<0.001);
    centerOfKernel=NaN(ndx*ndt,2);

    hWaitbar1=waitbar(0,'Radon 0 % Completed');
    kernelCounter1=0;
    kernelCounter2=0;        
    if featureDetectionON==0
        stddev=NaN(ndx*ndt,size(angles,2));
    end
    % Optional Vanzetta's Gaussian filtering technique (NOT CURRENTLY USED)
    gaussianSigma=floor((dx-1)/4); %sigma of gaussian filter calculated such that when filter size is chosen to be exactly equal to kernel size, the relation "FilterSize=2*ceil(2*SIGMA)+1" approximately holds
    h=fspecial('gaussian',dx,gaussianSigma);

    % Initialising parameters for creating space-time image with overlaid slope information
    if featureDetectionON==0
        SNRofEachKernel_Mode1=NaN(ndx*ndt,1);
        SNRofEachKernel_Mode2=NaN(ndx*ndt,1);
        lineXYCoordinates=NaN(ndt*ndx,4);
    end
    
    if featureDetectionON==0
        for t_counter=1:ndt
            kernelStartX=kernelStartXfirst;
            for x_counter=1:ndx
                waitbar(((t_counter-1)*ndx+x_counter)/(ndt*ndx),hWaitbar1,...
                    ['Radon ',...
                    num2str(single(100*((t_counter-1)*ndx+x_counter)/...
                    (ndt*ndx)),'%.0f'),' % completed']);
                kernel=data(kernelStartX:kernelStartX+dx-1,...
                    kernelStartT:kernelStartT+dx-1);

                kernelCounter1=kernelCounter1+1;                
                if prod(prod(~isnan(kernel)))~=0
                    kernelCounter2=kernelCounter2+1;
                    centerOfKernel(kernelCounter1,:)=[kernelStartX+dx/2 ...
                        kernelStartT+dx/2];
                    kernelSquare=kernel;
                    if gaussianWindowingON==0
                        kernel=kernel.*mask;
                    else
                        kernel=kernel.*h;
                    end
                    [r,xp]=radon(kernel,angles);
                    if kernelCounter2==1
                        if radonClippingON==1
                            rStart=(size(r,1)+1)/2-(dx+1)/2;
                            rEnd=(size(r,1)+1)/2+(dx+1)/2;
                        else
                            if radonClippingON==0
                                rStart=1;
                                rEnd=size(r,1);
                            end
                        end
                    end
                    r=r(rStart:rEnd,:);
                    xp=xp(rStart:rEnd);
                    if radonMode==1
                        stddev(kernelCounter1,:)=std(r,0,1);

                        % Shuffling
                        if shuffleON==1
                            stddevShuffled=zeros(size(stddev(kernelCounter1,:)));
                            for k=1:nShuffles
                                shuffledIndices=randperm(dx);
                                kernelShuffled=kernelSquare(:,shuffledIndices);
                                if gaussianWindowingON==0
                                    kernelShuffled=kernelShuffled.*mask;
                                else
                                    kernelShuffled=kernelShuffled.*h;
                                end
                                rShuffled=radon(kernelShuffled,angles);
                                rShuffled=rShuffled(rStart:rEnd,:);
                                stddevShuffled=std(rShuffled,0,1)+stddevShuffled;
                            end
                            stddevShuffled=stddevShuffled./nShuffles;
                            stddev(kernelCounter1,:)=stddev(kernelCounter1,...
                                :)-stddevShuffled;
                        end
                        
                        % normalizing
                        stddevTEMP=stddev(kernelCounter1,:);
                        stddevTEMP=stddevTEMP-min(stddevTEMP(:));
                        stddevTEMP=stddevTEMP./max(stddevTEMP(:));
                        stddev(kernelCounter1,:)=stddevTEMP;
                        % SNR calculation only if peak intensity of Radon
                        % is used to define SNR. If peak variance is used
                        % instead, the calulation of the same is done later
                        % (~58 lines of code later)
                        if calculateSNR2==1
                            [rMax,I]=max(r(:));
                            [I_row,I_col]=ind2sub(size(r),I);
                            rInterpolated=interp2(angles,xp,r,...
                                anglesLinearlySpaced,xp);
                            rMean=mean(rInterpolated(:));
                            SNRofEachKernel_Mode2(kernelCounter1)=rMax./rMean;
                        end
                    else
                        if radonMode==2
                            [rMax,I]=max(r(:));
                            [I_row,I_col]=ind2sub(size(r),I);
                            theta_maxClean(kernelCounter1)=angles(I_col);
                            rInterpolated=interp2(angles,xp,r,...
                                anglesLinearlySpaced,xp);
                            
                            rMean=mean(rInterpolated(:));
                            SNRofEachKernel_Mode2(kernelCounter1)=rMax./rMean;
                            m2=tand(theta_maxClean(kernelCounter1)+90);
                            if m2==Inf
                                lineXYCoordinates(kernelCounter1,:)=...
                                    [centerOfKernel(kernelCounter1,2) ...
                                    centerOfKernel(kernelCounter1,1)+d ...
                                    centerOfKernel(kernelCounter1,2) ...
                                    centerOfKernel(kernelCounter1,1)-d];
                            else
                                lineXYCoordinates(kernelCounter1,:)=...
                                    [centerOfKernel(kernelCounter1,2)-...
                                    d/sqrt(1+m2^2) centerOfKernel(kernelCounter1,...
                                    1)+m2*d/sqrt(1+m2^2) centerOfKernel(...
                                    kernelCounter1,2)+d/sqrt(1+m2^2) ...
                                    centerOfKernel(kernelCounter1,1)-...
                                    m2*d/sqrt(1+m2^2)];
                            end
                        end
                    end
                    kernelStartX=kernelStartX+dxStep;
                end
            end
            kernelStartT=kernelStartT+dtStep;
        end
        
        % Excluding array values which are NaN. Array values could be NaN at
        % this stage if image has NaNs in it (inserted during registration)
        arbitrary1D=centerOfKernel(:,2); % to get a 1D array
        percentKernelsProcessedStage1=sum(isfinite(arbitrary1D))./...
            numel(arbitrary1D);
        centerOfKernel(isnan(arbitrary1D),:)=[];
        stddev(isnan(arbitrary1D),:)=[];
        lineXYCoordinates(isnan(arbitrary1D),:)=[];
        theta_maxClean(isnan(arbitrary1D),:)=[];
        SNRofEachKernel_Mode1(isnan(arbitrary1D),:)=[];
        SNRofEachKernel_Mode2(isnan(arbitrary1D),:)=[];
        originalXYgrid=centerOfKernel;
        totalKernels=kernelCounter2;
        
        if radonMode==1
            hWaitbar2=waitbar(0,'FindPeaks 0 % Completed');
        end
        cmap=colormap(jet);
        close;
        kernelCounter2=0;
        x=centerOfKernel(:,2);
        y=centerOfKernel(:,1);
        if radonMode==1
            for kernelCounter2=1:size(SNRofEachKernel_Mode1,1)
                stddevTEMP=stddev(kernelCounter2,:)';
                stddevInterpolated=interp1(angles,stddevTEMP,...
                    anglesLinearlySpaced);
                [~,locsClean]=findpeaks(stddevTEMP,'SortStr','descend',...
                    'NPeaks',1);
                if numel(locsClean)~=0
                    theta_maxClean(kernelCounter2)=angles(locsClean);
                    if calculateSNR1==1
                        SNRofEachKernel_Mode1(kernelCounter2)=stddevTEMP(...
                            locsClean)/mean(stddevInterpolated(:));
                    end
                    m=tand(theta_maxClean(kernelCounter2)+90);
                    lineXYCoordinates(kernelCounter2,:)=[x(kernelCounter2)-...
                        d/sqrt(1+m^2)...
                        y(kernelCounter2)+m*d/sqrt(1+m^2) ...
                        x(kernelCounter2)+d/sqrt(1+m^2) ...
                        y(kernelCounter2)-m*d/sqrt(1+m^2)];
                end
                waitbar(kernelCounter2/size(SNRofEachKernel_Mode1,1),...
                    hWaitbar2,['FindPeaks ',...
                    num2str(single(100*(kernelCounter2/size(...
                    SNRofEachKernel_Mode1,1))),'%.0f'),' % completed']);
            end
        end
        
        % Excluding theta_maxClean values which are NaN. theta_maxClean could be NaN at
        % this stage if findpeaks did not find any peak in stddev
        percentKernelsProcessedStage2=sum(isfinite(theta_maxClean))./...
            numel(theta_maxClean);
        centerOfKernel(isnan(theta_maxClean),:)=[];
        x(isnan(theta_maxClean),:)=[];
        y(isnan(theta_maxClean),:)=[];
        stddev(isnan(theta_maxClean),:)=[];
        lineXYCoordinates(isnan(theta_maxClean),:)=[];
        SNRofEachKernel_Mode1(isnan(theta_maxClean),:)=[];
        SNRofEachKernel_Mode2(isnan(theta_maxClean),:)=[];
        theta_maxClean(isnan(theta_maxClean),:)=[];
    end
    clear data dataRegisteredWithStaticLines stddev;

    close(hWaitbar1);
    if radonMode==1
        close(hWaitbar2);
    end    
    
    % Saving copy of code (added 06-04-2018, 1.17 pm)
    codeFileFullPath=mfilename('fullpath'); %does not store .m extension
    [codeFilePath,codeFileName,~]=fileparts(codeFileFullPath);
    newCodeFileName=['code_',fileName(1:25),'_Step4',appendToResults];
    codeSource=[codeFilePath,'\',codeFileName,'.m'];
    codeDestination=[Step5folder,newCodeFileName,'.m'];
    copyfile(codeSource,codeDestination);
    
    currentFolder=cd(Step4folder);
    save([fileName(1:25),'_Step4',appendToResults]);
    cd(currentFolder);
    
    % moving processed file to separate folder
    if moveFileAfterProcessing==1
        movefile([processDir,fileName],[processDir,'PROCESSED\',fileName]);
    end
    elapsedTime=toc;disp(['Elapsed time is ',num2str(elapsedTime/60),' minutes.']);
end
beep
