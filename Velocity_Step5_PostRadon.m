% Velocity_Step5_PostRadon
% Author: Aby Joseph, University of Rochester
% License: GPL-3.0-or-later
% Last modified: 01-24-2019

clear;close all;
f_No=0;
fontsize1=18;
    
% Variables typically not changed
moveFileAfterProcessing=0;
SNRmode=2; %Mode 1 - calculate from peak variance of Radon, Mode 2 - calculate from peak intensity of Radon
SNRotsu=0; % If 1, calculates SNRthreshold based on Otsu threshold. Otherwise uses pre-determined values in code
SNRthreshold=2; %for SNRotsu==0
SNRcolorMode=2; % 1 = binary cyan/magenta, 2 = spectrum
overlayLineWidth=4;
velocityInterpolation=0;
removeAberrantVelocities=1; % use only with high density kernel sampling
removeLonelyKernels=1; % works only when removeAberrantVelocities==1

% Variables typically changed
videoDurationToAnalyze_Step5=NaN;
findFFT=1;
writeSNRimage=0;
writeAcceptRejectimage=0;
writeVelocityImage=0;
writeTimeCroppedRegisteredImage=0;
writeVelocitySquareImage=0;
computeAverageCardiacCycle=1;
makeAverageCardiacVideo=1;
makeAllCardiacVideo=0;
printWhiteLine=0;
SNRmapRangeManualOverride=1; % SNR color map, select '1' to override default clamps
minSNR=[]; % if overridden above, defines custom colormap min clamp. Insert empty to keep default
maxSNR=3.4; % if overridden above, defines custom colormap max clamp. Insert empty to keep default
manualPeakTroughEntry=0; % use manually determined cardiac peaks and troughs

thresholdMinVelocity=0.1; % in mm/s. Velocities below this are ignored.
minV=5; % in mm/s. Lower clamp for Jet Colormap
maxV=27; % in mm/s. Upper clamp for Jet Colormap
binningWindowTime=15/1000; %in seconds (Recommended: 15 ms for mouse and 90 ms for human)
cropLaminarProfile=1; % cropping laminar profile based on lane traffic analysis
laneTrafficThresholdPercent=2; % minimum required percent kernels analyzed 
                               % in a given lane for that lane to be
                               % considered in laminar profile. Percentage
                               % calculated relative to maximum number of
                               % kernels analyzed in any lane. 
                               % Typical value=5

Step5folder='E:\PhD\TRBF\Experimental DATA\Step 5 MAT files\';
Step4folder='E:\PhD\TRBF\Experimental DATA\Step 4 MAT files\';
Step3folder_Processed='E:\PhD\TRBF\Experimental DATA\Step 3 MAT files\PROCESSED\';
Step2folder_Processed='E:\PhD\TRBF\Experimental DATA\Step 2 MAT files\PROCESSED\';
RawDir_Processed='E:\PhD\TRBF\Experimental DATA\RAW MAT files\PROCESSED\';
cd(Step4folder);

% Append current date and time to filenames of results
currentDateAndTime=clock;
year=num2str(currentDateAndTime(1));
month=num2str(currentDateAndTime(2),'%02d');
day=num2str(currentDateAndTime(3),'%02d');
hour=num2str(currentDateAndTime(4),'%02d');
minute=num2str(currentDateAndTime(5),'%02d');
appendToResults=['_',year,month,day,hour,minute];

namecontains='*_Step4*.mat*';
[dirFilenames(1).name,filePath,~]=uigetfile(namecontains,...
    'Select file to process - click cancel to process all files in parent folder');
if filePath==0
    dirFilenames=dir([namecontains]);
    nFiles=numel(dirFilenames);
    processDir=Step4folder;
else
    processDir=[filePath];
    nFiles=1;
end
screenSize=get(0, 'MonitorPositions');

cd(processDir);
for fileNo=1:nFiles
    tic;
    f_No=0;
    close all;
    fileName=dirFilenames(fileNo).name;
    display(fileName);
    
    % Loading Step2.mat
    namecontains3=[fileName(1:25),'*Step2*.mat'];
    dirFilenames3=dir([Step2folder_Processed,namecontains3]);
    nFiles3=numel(dirFilenames3);
    if nFiles3==1
        load([Step2folder_Processed,dirFilenames3(1).name],'data');
    else
        [Step2FileName,Step2FilePath,~]=uigetfile([Step2folder_Processed,...
            namecontains3],'Select Step2 MAT file');
        load([Step2FilePath,Step2FileName],'data');
    end
    
    % Loading Step4.mat
    load([processDir,fileName],'freq','mic_pix',...
        'dataWidth_Step4','dataHeight_Step4','featureDetectionON',...
        'vesselAngle',...
        'x','y','timeAxisForMotionTrace','motionTrace',...
        'SNRofEachKernel_All',...
        'analysisBoundaries','theta_maxClean',...
        'SNRofEachKernel_Mode1','SNRofEachKernel_Mode2',...
        'cmap','lineXYCoordinates',...
        'overlayLineWidth','dx','theta_maxClean','originalXYgrid',...
        'videoDurationToAnalyze_Step4','quadrant','centerOfKernel',...
        'theta_max90',...
        'theta_max90','d','fractionalOverlapTime','fractionalOverlapSpace',...
        'dtStep','dxStep','percentKernelsProcessedStage1',...
        'percentKernelsProcessedStage2','ndx','ndt','Vmin','Vmax');

    if isnan(videoDurationToAnalyze_Step5)
        videoDurationToAnalyze_Step5=size(data,2)/freq;
    end

    data=data(:,1:floor(freq*videoDurationToAnalyze_Step5));
    
    if max(data(:))>125
        data=single(data);
        data=data./255;
    end
    dataWidth=size(data,2);
    dataHeight=size(data,1);    
    
    % Choosing SNR mode. Both SNR types were calculated in Step4
    if SNRmode==1
        SNRofEachKernel=SNRofEachKernel_Mode1;
    end
    if SNRmode==2
        SNRofEachKernel=SNRofEachKernel_Mode2;
    end
    
    % Cropping all variables to new analysis time window
    if videoDurationToAnalyze_Step5<videoDurationToAnalyze_Step4-0.001
        xLimit=floor(videoDurationToAnalyze_Step5*freq);
        centerOfKernel(x>xLimit,:)=[];
        lineXYCoordinates(x>xLimit,:)=[];
        SNRofEachKernel(x>xLimit,:)=[];
        theta_maxClean(x>xLimit,:)=[];
        y(x>xLimit,:)=[];
        x(x>xLimit,:)=[];
    end
    
    V_pixel=abs(cotd(theta_maxClean))*(1/abs(cosd(vesselAngle))); %Vx in pixel space
    
    %Started 08-21-2018 11 am
    %Rejecting ROIs with wrong quadrant solutions, earlier in the pipeline
    %Setting detected angle to NaN in kernels with angle in "wrong" quadrant
    if quadrant==1
        [row,column]=find(theta_maxClean>90);
    else
        [row,column]=find(theta_maxClean<90);
    end
    nFinite_pre3=sum(isfinite(V_pixel));
    for k=1:numel(row)
        V_pixel(row(k),column(k))=NaN;
    end
    nFinite_post3=sum(isfinite(V_pixel));
    percentKernelsPassedStage3=nFinite_post3./nFinite_pre3;
    SNRofEachKernel_includingWrongQuadrants=SNRofEachKernel;
    rejected_SNRoverlay_SNR=SNRofEachKernel(isnan(V_pixel));
    SNRoverlay_SNR=SNRofEachKernel(isfinite(V_pixel));
    SNRoverlay_lineXYCoordinates=lineXYCoordinates(isfinite(V_pixel),:);
    SNRoverlay_V_pixel=V_pixel(isfinite(V_pixel));
    
    % Finding SNR threshold.
    SNRofEachKernel_norm=SNRoverlay_SNR-min(SNRoverlay_SNR(:));
    SNRofEachKernel_norm=SNRofEachKernel_norm./max(...
        SNRofEachKernel_norm(:)); % "graythresh" function seems to only work if the input (in double format) is between 0 and 1
    [level,EM]=graythresh(SNRofEachKernel_norm);
    if SNRotsu==1
        SNRthreshold=min(SNRoverlay_SNR(:))+...
            level.*(max(SNRoverlay_SNR(:))-min(SNRoverlay_SNR(:)));
    end
    
    % Creating space-time image with overlaid slope information
    if SNRmapRangeManualOverride==0
        minSNR=nanmin(SNRoverlay_SNR(:));
        maxSNR=nanmax(SNRoverlay_SNR(:));
    else
        if isempty(minSNR)
            minSNR=nanmin(SNRoverlay_SNR(:));
        end
        if isempty(maxSNR)
            maxSNR=nanmax(SNRoverlay_SNR(:));
        end
    end
    
    rangeSNR=linspace(minSNR,maxSNR,size(cmap,1));
    lineColor=NaN(numel(SNRoverlay_SNR),3); % totalKernals is no. of kernels placed in non-NaN regions of image
    kernelCounter3=0;
    for kernelCounter3=1:size(lineColor,1)
        if SNRcolorMode==1
            if SNRoverlay_SNR(kernelCounter3)>=SNRthreshold
                lineColor(kernelCounter3,:)=[1 0 1]; % magenta
            else
                lineColor(kernelCounter3,:)=[0 1 1]; % cyan
            end
        else
            if SNRcolorMode==2
                cmap_SNR=colormap(hot);
                if SNRoverlay_SNR(kernelCounter3)<=maxSNR
                    lineColor(kernelCounter3,:)=cmap_SNR((find(...
                        rangeSNR>=SNRoverlay_SNR(kernelCounter3),1)),:);
                else
                    lineColor(kernelCounter3,:)=cmap_SNR(end,:);
                end
            end
        end
    end
    
    data(isnan(data))=0;
    data=data.*255;
    dataColor=uint8(cat(3,data,data,data));
    lineColor=lineColor.*255;
    if writeSNRimage==1
        dataColor=insertShape(dataColor,'Line',...
            SNRoverlay_lineXYCoordinates,'LineWidth',overlayLineWidth,...
            'Color',lineColor);
        if printWhiteLine==1
            dataColor=insertShape(dataColor,'Line',[1 ...
                analysisBoundaries(1)-(dx-1)/2 size(dataColor,2)...
                analysisBoundaries(1)-(dx-1)/2],'LineWidth',...
                overlayLineWidth*0.5,'Color',[255 255 255]);
            dataColor=insertShape(dataColor,'Line',[1 ...
                analysisBoundaries(2)+(dx-1)/2 size(dataColor,2)...
                analysisBoundaries(2)+(dx-1)/2],'LineWidth',...
                overlayLineWidth*0.5,'Color',[255 255 255]);
        end
        currentFolder=cd(Step5folder);
        imwrite(dataColor,[fileName(1:25),'_Step5_Overlay_SNR',...
            appendToResults,'.tif'],'TIFF');
        cd(currentFolder);
    end
    
    if featureDetectionON==1
        dataThresh=uint8(dataThresh).*255;
        dataColor=uint8(cat(3,dataThresh,dataThresh,dataThresh));
        currentFolder=cd(Step5folder);
        imwrite(dataColor,[fileName(1:25),'_Step5_Binary',...
            appendToResults,'.tif'],'TIFF');
        cd(currentFolder);
        gpuDataFiltered=uint8(gpuDataFiltered.*255);
        dataColor=cat(3,gpuDataFiltered,gpuDataFiltered,gpuDataFiltered);
        currentFolder=cd(Step5folder);
        imwrite(dataColor,[fileName(1:25),'_Step5_Gaussian',...
            appendToResults,'.tif'],'TIFF');
        cd(currentFolder);
        clear dataColor gpuDataFiltered;
    end
    
    % Setting velocity to NaN in kernels with below threshold SNR
    [row,column]=find(SNRofEachKernel<SNRthreshold);
    nFinite_pre4=sum(isfinite(V_pixel));
    for k=1:numel(row)
        V_pixel(row(k),column(k))=NaN;
    end
    nFinite_post4=sum(isfinite(V_pixel));
    percentKernelsPassedStage4=nFinite_post4./nFinite_pre4;
    nKernelsWithSNRBelowThreshold=k;

    % Section started on 08-01-2018 1.02 pm
    % Setting velocity to NaN in kernels with aberrant velocity in small ROIs
    if removeAberrantVelocities==1
        countAberrant=0;
        nFinite_pre5=sum(isfinite(V_pixel));
        vROIwindowTime=2.1*dx;
        vROIwindowSpace=2.1*dx;
        for k=1:numel(V_pixel)
            if not(isnan(V_pixel(k)))
                spaceStart=y(k)-vROIwindowSpace/2;
                spaceEnd=y(k)+vROIwindowSpace/2;
                timeStart=x(k)-vROIwindowTime/2;
                timeEnd=x(k)+vROIwindowTime/2;                
                Vpixel_temp=V_pixel(y>spaceStart & y<spaceEnd & ...
                    x>timeStart & x<timeEnd);
                ROI_vmean=nanmean(Vpixel_temp);
                ROI_vstd=nanstd(Vpixel_temp);
                if V_pixel(k)>(ROI_vmean+1*ROI_vstd) || V_pixel(k)...
                        <(ROI_vmean-1*ROI_vstd)
                    countAberrant=countAberrant+1;
                    V_pixel(k)=NaN;
                end
            end
        end
        nFinite_post5=sum(isfinite(V_pixel));
        percentKernelsPassedStage5=nFinite_post5./nFinite_pre5;
    end
    
    % Section started on 08-01-2018 7.30 pm
    % Setting velocity to NaN in 'lonely' kernels in space dimension
    if removeLonelyKernels==1 && removeAberrantVelocities==1
        countLonely=0;
        nFinite_pre6=sum(isfinite(V_pixel));
        lonelyROIwindowTime=vROIwindowTime;
        lonelyROIwindowSpace=vROIwindowSpace;
        for k=1:numel(V_pixel)
            if not(isnan(V_pixel(k)))
                spaceStart=y(k)-lonelyROIwindowSpace/2;
                spaceEnd=y(k)+lonelyROIwindowSpace/2;
                timeStart=x(k)-lonelyROIwindowTime/2;
                timeEnd=x(k)+lonelyROIwindowTime/2;                
                Vpixel_temp_lonely=V_pixel(y>spaceStart & y<spaceEnd &...
                    x>timeStart & x<timeEnd);
                Vpixel_temp_lonely=Vpixel_temp_lonely(~isnan(...
                    Vpixel_temp_lonely));
                if numel(Vpixel_temp_lonely)<=1
                    countLonely=countLonely+1;
                    V_pixel(k)=NaN;
                end
            end
        end
        nFinite_post6=sum(isfinite(V_pixel));
        percentKernelsPassedStage6=nFinite_post6./nFinite_pre6;
    end
    
    % Setting velocity to NaN in kernels with velocity "many times" higher
    % than mean velocity in vessel
    thresholdMaxVelocity_pixels=nanmean(V_pixel(:))+7.*nanstd(V_pixel(:));
    thresholdMaxVelocity=thresholdMaxVelocity_pixels.*(mic_pix*freq*1e-3)
    [row,column]=find(V_pixel>thresholdMaxVelocity_pixels);
    nFinite_pre7=sum(isfinite(V_pixel));
    for k=1:numel(row)
        V_pixel(row(k),column(k))=NaN;
    end
    nFinite_post7=sum(isfinite(V_pixel));
    percentKernelsPassedStage7=nFinite_post7./nFinite_pre7;
    nKernelsWithVelocityAboveThreshold=k;

    % Setting velocity to NaN in kernels with velocity that is too low
    thresholdMinVelocity_pixels=thresholdMinVelocity/(mic_pix*freq*1e-3);
    [row,column]=find(V_pixel<thresholdMinVelocity_pixels);
    nFinite_pre8=sum(isfinite(V_pixel));
    for k=1:numel(row)
        V_pixel(row(k),column(k))=NaN;
    end
    nFinite_post8=sum(isfinite(V_pixel));
    percentKernelsPassedStage8=nFinite_post8./nFinite_pre8;
    nKernelsWithVelocityBelowMinThreshold=k;    
    
    % Deleting all NaNs in variables "V_pixel", "SNRofEachKernel", "centerOfKernel",
    ..."stddevClean", "theta_max90", "theta_maxClean", "x" and "y". Section started 03-20-2017
    
    rejected_SNRofEachKernel_V=SNRofEachKernel(isnan(V_pixel));
    rejected_lineXYCoordinates_V=lineXYCoordinates(isnan(V_pixel),:);
    rejected_centerOfKernel_V=centerOfKernel(isnan(V_pixel),:);
    rejected_theta_max90_V=theta_max90(isnan(V_pixel));
    rejected_theta_maxClean_V=theta_maxClean(isnan(V_pixel));
    rejected_x_V=x(isnan(V_pixel));
    rejected_y_V=y(isnan(V_pixel));
    rejected_V_pixel_V=V_pixel(isnan(V_pixel));
    
    SNRofEachKernel_V=SNRofEachKernel(isfinite(V_pixel));
    lineXYCoordinates_V=lineXYCoordinates(isfinite(V_pixel),:);
    centerOfKernel(isnan(V_pixel),:)=[];
    theta_max90(isnan(V_pixel))=[];
    theta_maxClean(isnan(V_pixel))=[];
    x(isnan(V_pixel))=[];
    y(isnan(V_pixel))=[];
    V_pixel(isnan(V_pixel))=[];

    % 2D interpolation of velocities (started 11-24-2017)
    if velocityInterpolation==1
        F=scatteredInterpolant(x,y,V_pixel);
        V_pixel_interpolated=F(originalXYgrid(:,2),originalXYgrid(:,1));
        percentKernelsVelocityInterpolated=100*((numel(...
            V_pixel_interpolated)-numel(V_pixel))/numel(...
            V_pixel_interpolated));
        V_noInterp=V_pixel*(1e-3)*mic_pix*freq; %Vx in mm/s
        x_noInterp=x;
        y_noInterp=y;
        x=originalXYgrid(:,2);
        y=originalXYgrid(:,1);        
        theta_maxClean_Vinterp=acotd(V_pixel_interpolated*abs(...
            cosd(vesselAngle)));
        if quadrant==2
            theta_maxClean_Vinterp=180-theta_maxClean_Vinterp;
        end
        V=V_pixel_interpolated*(1e-3)*mic_pix*freq; %Vx in mm/s
        AvergeVelocity_withoutInterp=nanmean(V_noInterp(:)) % in mm/s
    else
        V=V_pixel*(1e-3)*mic_pix*freq; %Vx in mm/s        
    end
    
    AvergeVelocity=nanmean(V(:)) % in mm/s
    StdVelocity=nanstd(V(:)) % in mm/s
    AvergeVelocity_pixels=nanmean(V_pixel(:)) % in mm/s
    StdVelocity_pixels=nanstd(V_pixel(:)); % in mm/s
    isNANmatrix=isnan(V);
    NumberNANs=sum(isNANmatrix(:));
    NumberNonNANs=numel(V)-NumberNANs;
    SemVelocity=StdVelocity./sqrt(NumberNonNANs); % in mm/s
    
    rangeV=linspace(minV,maxV,size(cmap,1));
    lineColor_V=NaN(numel(V(isfinite(V))),3);
    counter6=0;
    
    if featureDetectionON==0
        nFiniteV=sum(sum(isfinite(V)));
        for kernelCounter=1:numel(V)
            if isfinite(V(kernelCounter))
                counter6=counter6+1;
                if V(kernelCounter)<=maxV
                    lineColor_V(counter6,:)=cmap((find(rangeV>=...
                        V(kernelCounter),1)),:);
                else
                    lineColor_V(counter6,:)=cmap(end,:);
                end
            end
        end
    end
    
    lineColor_V=lineColor_V.*255;
    if writeVelocityImage==1
        % Overlaying lines with color represeting velocity
        dataColor=uint8(cat(3,data,data,data));
        dataColor=insertShape(dataColor,'Line',lineXYCoordinates_V,...
            'LineWidth',overlayLineWidth,'Color',lineColor_V);
        if printWhiteLine==1
            dataColor=insertShape(dataColor,'Line',[1 ...
                analysisBoundaries(1)-(dx-1)/2 size(dataColor,2)...
                analysisBoundaries(1)-(dx-1)/2],'LineWidth',...
                overlayLineWidth*0.5,'Color',[255 255 255]);
            dataColor=insertShape(dataColor,'Line',[1 ...
                analysisBoundaries(2)+(dx-1)/2 size(dataColor,2)...
                analysisBoundaries(2)+(dx-1)/2],'LineWidth',...
                overlayLineWidth*0.5,'Color',[255 255 255]);
        end
        currentFolder=cd(Step5folder);
        imwrite(dataColor,[fileName(1:25),...
            '_Step5_Overlay_Velocity_Lines',appendToResults,'.tif'],...
            'TIFF');
        cd(currentFolder);
    end

    if writeTimeCroppedRegisteredImage==1
        dataColor=uint8(cat(3,data,data,data));
        if printWhiteLine==1
            dataColor=insertShape(dataColor,'Line',[1 ...
                analysisBoundaries(1)-(dx-1)/2 size(dataColor,2)...
                analysisBoundaries(1)-(dx-1)/2],'LineWidth',...
                overlayLineWidth*0.5,'Color',[255 255 255]);
            dataColor=insertShape(dataColor,'Line',[1 ...
                analysisBoundaries(2)+(dx-1)/2 size(dataColor,2)...
                analysisBoundaries(2)+(dx-1)/2],'LineWidth',...
                overlayLineWidth*0.5,'Color',[255 255 255]);
        end
        currentFolder=cd(Step5folder);
        imwrite(dataColor,[fileName(1:25),...
            '_Step5_Registered_TimeCropped',appendToResults,'.tif'],'TIFF');
        cd(currentFolder);
    end
    
    % Creating 2nd generation SNR overlay, where 'accepted' kernels now
    % only include ones which passed criteria of correct quadrant and
    % velocity thresholds in addition to the SNR threshold
    if writeAcceptRejectimage==1
        acceptRejectCyan=NaN(size(rejected_lineXYCoordinates_V,1),3);
        acceptRejectCyan(:,1)=0;
        acceptRejectCyan(:,2)=255;
        acceptRejectCyan(:,3)=255; %cyan color
        acceptRejectMagenta=NaN(size(lineXYCoordinates_V,1),3);
        acceptRejectMagenta(:,1)=255;
        acceptRejectMagenta(:,2)=0;
        acceptRejectMagenta(:,3)=255; %magenta color
        dataColor=uint8(cat(3,data,data,data));
        dataColor=insertShape(dataColor,'Line',...
            rejected_lineXYCoordinates_V,'LineWidth',overlayLineWidth,...
            'Color',acceptRejectCyan);
        dataColor=insertShape(dataColor,'Line',lineXYCoordinates_V,...
            'LineWidth',overlayLineWidth,'Color',acceptRejectMagenta);
        if printWhiteLine==1
            dataColor=insertShape(dataColor,'Line',[1 ...
                analysisBoundaries(1)-(dx-1)/2 size(dataColor,2)...
                analysisBoundaries(1)-(dx-1)/2],'LineWidth',...
                overlayLineWidth*0.5,'Color',[255 255 255]);
            dataColor=insertShape(dataColor,'Line',[1 ...
                analysisBoundaries(2)+(dx-1)/2 size(dataColor,2)...
                analysisBoundaries(2)+(dx-1)/2],'LineWidth',...
                overlayLineWidth*0.5,'Color',[255 255 255]);
        end
        currentFolder=cd(Step5folder);
        imwrite(dataColor,[fileName(1:25),...
            '_Step5_Overlay_AcceptReject',appendToResults,'.tif'],'TIFF');
        cd(currentFolder);
    end

    % Creating overlay velocity image without interpolation, when
    % interpolation is on
    if velocityInterpolation==1
        rangeV=linspace(minV,maxV,size(cmap,1));
        lineXYCoordinates_V_noInterp=NaN(numel(V_noInterp(...
            isfinite(V_noInterp))),4);
        lineColor_V_noInterp=NaN(numel(V_noInterp(isfinite(...
            V_noInterp))),3);
        counter6=0;
        
        if featureDetectionON==0
            nFiniteV_noInterp=sum(sum(isfinite(V_noInterp)));
            for kernelCounter=1:numel(V_noInterp)
                if isfinite(V_noInterp(kernelCounter))
                    counter6=counter6+1;
                    m=tand(theta_maxClean(kernelCounter)+90);
                    lineXYCoordinates_V_noInterp(counter6,:)=...
                        [x_noInterp(kernelCounter)-d/sqrt(1+m^2)...
                        y_noInterp(kernelCounter)+m*d/sqrt(1+m^2) ...
                        x_noInterp(kernelCounter)+d/sqrt(1+m^2)...
                        y_noInterp(kernelCounter)-m*d/sqrt(1+m^2)];
                    if V_noInterp(kernelCounter)<=maxV
                        lineColor_V_noInterp(counter6,:)=...
                            cmap((find(rangeV>=V_noInterp(...
                            kernelCounter),1)),:);
                    else
                        lineColor_V_noInterp(counter6,:)=cmap(end,:);
                    end
                end
            end
        end
        
        if writeVelocityImage==1
            % Overlaying lines with color represeting velocity
            dataColor=uint8(cat(3,data,data,data));
            dataColor=insertShape(dataColor,'Line',...
                lineXYCoordinates_V_noInterp,...
                'LineWidth',overlayLineWidth,'Color',...
                lineColor_V_noInterp.*255);
            dataColor=insertShape(dataColor,'Line',[1 ...
                analysisBoundaries(1)-(dx-1)/2 size(dataColor,2)...
                analysisBoundaries(1)-(dx-1)/2],'LineWidth',...
                overlayLineWidth*0.5,'Color',[255 255 255]);
            dataColor=insertShape(dataColor,'Line',[1 ...
                analysisBoundaries(2)+(dx-1)/2 size(dataColor,2)...
                analysisBoundaries(2)+(dx-1)/2],'LineWidth',...
                overlayLineWidth*0.5,'Color',[255 255 255]);
            currentFolder=cd(Step5folder);
            imwrite(dataColor,[fileName(1:25),...
                '_Step5_Overlay_Velocity_Lines_noInterp',...
                appendToResults,'.tif'],'TIFF');
            cd(currentFolder);
        end
    end
    
    % Creating RGB linescan image with encoded velocity-color information
    if fractionalOverlapTime==1 && fractionalOverlapSpace==1 && ...
            writeVelocitySquareImage==1
        dataColor=uint8(cat(3,data,data,data));
        for k=1:size(lineColor_V,1)
            dataColor(y(k)-d:y(k)+d,x(k)-d:x(k)+d,1)=...
                dataColor(y(k)-d:y(k)+d,x(k)-d:x(k)+d,1).*lineColor_V(k,1);
            dataColor(y(k)-d:y(k)+d,x(k)-d:x(k)+d,2)=...
                dataColor(y(k)-d:y(k)+d,x(k)-d:x(k)+d,2).*lineColor_V(k,2);
            dataColor(y(k)-d:y(k)+d,x(k)-d:x(k)+d,3)=...
                dataColor(y(k)-d:y(k)+d,x(k)-d:x(k)+d,3).*lineColor_V(k,3);
        end
        currentFolder=cd(Step5folder);
        imwrite(dataColor,[fileName(1:25),'_Step5_RGB_Velocity',...
            appendToResults,'.tif'],'TIFF');
        cd(currentFolder);
    end
    clear data dataColor dataRegisteredWithStaticLines;
    
    if featureDetectionON==0
        % defining binning window as n*dtStep such that window size is 
        % approximately equal to 'binningWindowTime'. (Recommended: 15 ms
        % for mouse and 90 ms for human)
        velocityBinningWindowInTimePixels=floor((binningWindowTime*...
            freq)/dtStep)*dtStep; 
        velocityBinningStepSizeInTimePixels=...
            velocityBinningWindowInTimePixels;
        velocityBinningStepSizeInTimeSeconds=...
            velocityBinningStepSizeInTimePixels/freq; 
        velocityBinningWindowInTimeSeconds=...
            velocityBinningWindowInTimePixels/freq; 
        
        velocityBinningWindowInSpacePixels=dxStep;
        velocityBinningStepSizeInSpacePixels=...
            velocityBinningWindowInSpacePixels;

        velocityBinningWindowInSpaceMicrons=...
            velocityBinningWindowInSpacePixels*mic_pix;
        velocityBinningStepSizeInSpaceMicrons=...
            velocityBinningStepSizeInSpacePixels*mic_pix;
        TimePixels=velocityBinningWindowInTimePixels/2:...
            velocityBinningStepSizeInTimePixels:(...
            dataWidth_Step4-velocityBinningWindowInTimePixels/2);
        Time=TimePixels./freq; %in seconds
        SpacePixels=originalXYgrid(1:ndx,1)';
        SpacePixels_offset=SpacePixels-SpacePixels(...
            round(size(SpacePixels,2)/2));
        Space=SpacePixels_offset.*mic_pix.*abs(sind(...
            vesselAngle)); %in microns
        V_t=NaN(1,numel(Time));
        V_t_std=V_t;
        V_t_sem=V_t;
        V_x=NaN(1,numel(Space));
        V_x_std=V_x;
        V_x_sem=V_x;
        
        for k=1:numel(V_t)
            tStart=TimePixels(k)-velocityBinningWindowInTimePixels/2;
            tEnd=TimePixels(k)+velocityBinningWindowInTimePixels/2;
            V_temp=V(x>=tStart & x<=tEnd);
            if ~isempty(V_temp) && numel(V_temp>=3)
                V_t(k)=mean(V_temp,1);
                V_t_std(k)=std(V_temp,0,1);
                V_t_sem(k)=V_t_std(k)./sqrt(numel(V_temp));
            end
        end
        for k=1:numel(V_x)
            xStart=SpacePixels(k)-velocityBinningWindowInSpacePixels/2;
            xEnd=SpacePixels(k)+velocityBinningWindowInSpacePixels/2;            
            V_temp=V(y>=xStart & y<=xEnd);
            if ~isempty(V_temp) && numel(V_temp>=3)
                V_x(k)=mean(V_temp,1);
                V_x_std(k)=std(V_temp,0,1);
                V_x_sem(k)=V_x_std(k)./sqrt(numel(V_temp));
            end
        end
        V_xt=NaN(numel(V_x),numel(V_t));
        V_xt_std=V_xt;
        V_xt_sem=V_xt;
        for k=1:numel(V_x)
            xStart=SpacePixels(k)-velocityBinningWindowInSpacePixels/2;
            xEnd=SpacePixels(k)+velocityBinningWindowInSpacePixels/2;            
            for m=1:numel(V_t)
                tStart=TimePixels(m)-velocityBinningWindowInTimePixels/2;
                tEnd=TimePixels(m)+velocityBinningWindowInTimePixels/2;
                V_temp=V(y>=xStart & y<=xEnd & x>=tStart & x<=tEnd);
                if ~isempty(V_temp) && numel(V_temp>=3)
                    V_xt(k,m)=mean(V_temp,1);
                    V_xt_std(k,m)=std(V_temp,0,1);
                    V_xt_sem(k,m)=V_xt_std(k,m)./sqrt(numel(V_temp));
                end
            end
        end
    end
    
    % Computing lane traffic (i.e. counting no. of kernels analyzed in
    % each 'lane' - each radial position from center of lumen.)
    laneTraffic=NaN(1,numel(Space));
    for k=1:numel(Space)
        xStart=SpacePixels(k)-velocityBinningWindowInSpacePixels/2;
        xEnd=SpacePixels(k)+velocityBinningWindowInSpacePixels/2;
        V_temp=V(y>=xStart & y<=xEnd);
        laneTraffic(1,k)=numel(V_temp);
    end
        
    % Calculating heart rate (bpm)
    if findFFT==1
        dg=velocityBinningStepSizeInTimePixels/freq;
        V_t_isNAN=isnan(V_t);
        currentFolder=cd('E:\PhD\TRBF\CODES\MOUSE');
        seq=findseq(single(V_t_isNAN));
        cd(currentFolder);
        seqZeros=seq(abs(seq(:,1)-0)<.01,:);
        [m,I]=max(seqZeros(:,4),[],1);
        V_t_bpm=V_t(1,seqZeros(I,2):seqZeros(I,3));
        
        m=size(V_t_bpm,2);
        centx=floor(m/2)+1;
        indx=1:m;
        indx=indx-centx;
        dG=1/(m*dg);
        Gx=dG*indx;Gx=Gx.';
        Gx_display=Gx*60;
        Gy=abs(fftshift(fft(V_t_bpm)));
        Gy_display=Gy;
        if prod(~isnan(Gy))
            [~,loc3]=findpeaks(Gy,'SORTSTR','descend');
            if numel(loc3)>=2
                bpm=abs(Gx(loc3(2)))*60;
                Gy_display(loc3(1))=0;
            end
        end
    end
    
    % Finding average cardiac cycle
    bpmExist=exist('bpm','var');
    if bpmExist==1 && computeAverageCardiacCycle==1 
        cardiacHalfCycleArbitraryUnits=0.5*(1/bpm)*60/...
            velocityBinningStepSizeInTimeSeconds; % finding cardiac cycle time period in same units as step size in "V_t"
        [~,cardiacPeaks]=findpeaks(V_t,'MinPeakDistance',...
            cardiacHalfCycleArbitraryUnits*1.5);
        [~,cardiacTroughs]=findpeaks(-1*V_t,'MinPeakDistance',...
            cardiacHalfCycleArbitraryUnits*1.5);
        f_No=f_No+1;
        figure(f_No);
        currentFigure=gcf;
        plot(Time,V_t,'-ok','Linewidth',3);axis auto;hold on;
        plot(Time(cardiacPeaks),V_t(cardiacPeaks),'or','Linewidth',3);
        hold on;
        plot(Time(cardiacTroughs),V_t(cardiacTroughs),'ob','Linewidth',3);
        hold off;
        ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
        currentXLim=ax.XLim;currentYLim=ax.YLim;
        currentFigure.OuterPosition(1)=250;
        currentFigure.OuterPosition(2)=450;
        currentFigure.OuterPosition(3)=1500;
        currentFigure.OuterPosition(4)=600;
        currentOuterPosition=currentFigure.OuterPosition;
        xlabel('Time (s)');
        ylabel('Average velocity (mm/s)');
        title('');
        if manualPeakTroughEntry==0
            title('Click on peaks/troughs to be removed. Press Enter when done.');
            exitFlag=0;
            while exitFlag==0
                [xTemp,~]=ginput(1);
                if numel(xTemp)==0
                    exitFlag=1;
                else
                    xTemp=xTemp./velocityBinningStepSizeInTimeSeconds; % converting to same units as step size in "V_t"
                    k=find(cardiacPeaks>xTemp-...
                        cardiacHalfCycleArbitraryUnits/8 & ...
                        cardiacPeaks<xTemp+...
                        cardiacHalfCycleArbitraryUnits/8,1);
                    if numel(k)~=0
                        cardiacPeaks(k)=[];
                    end
                    m=find(cardiacTroughs>xTemp-...
                        cardiacHalfCycleArbitraryUnits/8 & ...
                        cardiacTroughs<xTemp+...
                        cardiacHalfCycleArbitraryUnits/8,1);
                    if numel(m)~=0
                        cardiacTroughs(m)=[];
                    end
                    plot(Time,V_t,'-ok','Linewidth',3);axis auto;hold on;
                    plot(Time(cardiacPeaks),V_t(cardiacPeaks),'or',...
                        'Linewidth',3);hold on;
                    plot(Time(cardiacTroughs),V_t(cardiacTroughs),...
                        'ob','Linewidth',3);hold off;
                    ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
                    ax.XLim=currentXLim;ax.YLim=currentYLim;
                    xlabel('Time (s)');
                    ylabel('Average velocity (mm/s)');
                    title('Click on peaks/troughs to be removed. Press Enter when done.');
                end
            end
        end
        plot(Time,V_t,'-ok','Linewidth',3);axis auto;hold on;
        plot(Time(cardiacPeaks),V_t(cardiacPeaks),'or','Linewidth',3);
        hold on;
        plot(Time(cardiacTroughs),V_t(cardiacTroughs),'ob','Linewidth',3);
        hold off;
        ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
        ax.XLim=currentXLim;ax.YLim=currentYLim;
        xlabel('Time (s)');
        ylabel('Average velocity (mm/s)');
        title('Original');
        currentFolder=cd(Step5folder);
        print(f_No,[fileName(1:25),'_Step5_cardiacWavesOriginal',...
            appendToResults],'-dtiff');
        saveas(f_No,[fileName(1:25),'_Step5_cardiacWavesOriginal',...
            appendToResults,'.fig'],'fig');
        cd(currentFolder);
        
        cardiacTroughsDiff=diff(cardiacTroughs);
        avgCardiacPeriod_in_V_t_units=mean(cardiacTroughsDiff(:));
        % computing "nSamplesPerCycle" such that sampling frequency in
        % average cardiac cycle is close to "velocityBinningStepSizeInTime"
        nSamplesPerCycle=floor(avgCardiacPeriod_in_V_t_units);
        nCardiacPeaks=numel(cardiacPeaks);
        nCardiacTroughs=numel(cardiacTroughs);
        nCardiacCycles=min([nCardiacPeaks nCardiacTroughs]);
        TimeCoorsOfCardiacPhases_Pixels=NaN(nCardiacCycles,...
            nSamplesPerCycle+2);
        cardiacPeaksTEMP=cardiacPeaks;
        cardiacTroughsTEMP=cardiacTroughs;
        k=1;
        TimeCoorsOfCardiacPhases_Pixels(1,1)=TimePixels(cardiacTroughs(1));
        previousPosition=1;
        cardiacTroughsTEMP(1)=[];
        previousElement=-1; % trough=-1 and peak=1
        % Making sure that peaks and troughs are appropriately positioned in
        % relevant array. Each cardaic cycle is clamped to the two troughs 
        % at its each side.
        while k<nCardiacCycles*3
            k=k+1;
            if mod(k,3)==1
                currentPosition=1;
                if cardiacPeaksTEMP(1)>cardiacTroughsTEMP(2)
                    cardiacTroughsTEMP(1)=[];
                end
                TimeCoorsOfCardiacPhases_Pixels(ceil(k/3),...
                    currentPosition)=TimePixels(cardiacTroughsTEMP(1));
                cardiacTroughsTEMP(1)=[];
            end
            if mod(k,3)==2
                cardiacPeaksTEMP(1)=[];
            end
            if mod(k,3)==0
                currentPosition=nSamplesPerCycle+2;
                TimeCoorsOfCardiacPhases_Pixels(ceil(k/3),...
                    currentPosition)=TimePixels(cardiacTroughsTEMP(1));
            end
        end
        
        start1=1;
        end1=nSamplesPerCycle+2;
        for k=1:nCardiacCycles
            arrayTemp=linspace(TimeCoorsOfCardiacPhases_Pixels(...
                k,start1),TimeCoorsOfCardiacPhases_Pixels(k,end1),...
                nSamplesPerCycle+2);
            TimeCoorsOfCardiacPhases_Pixels(k,start1+1:end1-1)=...
                arrayTemp(2:end-1);
        end
        
        % Average cardiac cycle - Computing Velocity vs Time
        f_No=f_No+1;
        figure(f_No);
        currentFigure=gcf;
        plot(Time,V_t,'k','Linewidth',3);axis auto;hold on;
        ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
        currentXLim=ax.XLim;currentYLim=ax.YLim;
        currentFigure.OuterPosition(1)=250;
        currentFigure.OuterPosition(2)=450;
        currentFigure.OuterPosition(3)=1500;
        currentFigure.OuterPosition(4)=600;
        currentOuterPosition=currentFigure.OuterPosition;
        xlabel('Time (s)');
        ylabel('Average velocity (mm/s)');
        title('Resampled');
        TimeCoorsOfCardiacPhases_Seconds=...
            TimeCoorsOfCardiacPhases_Pixels/freq;
        V_t_AllCardiac=NaN(nCardiacCycles,nSamplesPerCycle+2);
        V_t_AllCardiac_std=V_t_AllCardiac;
        V_t_AllCardiac_sem=V_t_AllCardiac;
        for k=1:nCardiacCycles
            for m=start1:end1 % cardiac phases other than two end troughs
                tStart=TimeCoorsOfCardiacPhases_Pixels(k,m)-...
                    velocityBinningWindowInTimePixels/2;
                tEnd=TimeCoorsOfCardiacPhases_Pixels(k,m)+...
                    velocityBinningWindowInTimePixels/2;
                V_temp=V(x>=tStart & x<=tEnd);
                if ~isempty(V_temp) && numel(V_temp>=3)
                    V_t_AllCardiac(k,m)=mean(V_temp,1);
                    V_t_AllCardiac_std(k,m)=std(V_temp,0,1);
                    V_t_AllCardiac_sem(k,m)=V_t_AllCardiac_std(k,m)./...
                        sqrt(numel(V_temp));
                end
            end
            plot(TimeCoorsOfCardiacPhases_Seconds(k,start1+1:end1-1),...
                V_t_AllCardiac(k,2:end-1),'ok','Linewidth',3);hold on;
            plot(TimeCoorsOfCardiacPhases_Seconds(k,start1),...
                V_t_AllCardiac(k,start1),'ob','Linewidth',3);hold on;
            plot(TimeCoorsOfCardiacPhases_Seconds(k,end1),...
                V_t_AllCardiac(k,end1),'ob','Linewidth',3);hold on;
        end
        currentFolder=cd(Step5folder);
        print(f_No,[fileName(1:25),'_Step5_cardiacWavesResampled',...
            appendToResults],'-dtiff');
        saveas(f_No,[fileName(1:25),'_Step5_cardiacWavesResampled',...
            appendToResults,'.fig'],'fig');
        cd(currentFolder);
        V_t_AvgCardiac=nanmean(V_t_AllCardiac,1);
        V_t_AvgCardiac_std=nanstd(V_t_AllCardiac,0,1);
        V_t_AvgCardiac_sem=V_t_AvgCardiac_std./sqrt(size(V_t_AllCardiac,1));
        
        % Computing histograms in diastolic and systolic phases
        V_systolic=[];
        for k=1:nCardiacPeaks
            tStart=TimePixels(cardiacPeaks(k))-velocityBinningWindowInTimePixels/2;
            tEnd=TimePixels(cardiacPeaks(k))+velocityBinningWindowInTimePixels/2;
            V_temp=V(x>=tStart & x<=tEnd);
            if ~isempty(V_temp) && numel(V_temp>=3)
                V_systolic=vertcat(V_systolic,V_temp);
            end
        end        
        V_diastolic=[];
        for k=1:nCardiacTroughs
            tStart=TimePixels(cardiacTroughs(k))-velocityBinningWindowInTimePixels/2;
            tEnd=TimePixels(cardiacTroughs(k))+velocityBinningWindowInTimePixels/2;
            V_temp=V(x>=tStart & x<=tEnd);
            if ~isempty(V_temp) && numel(V_temp>=3)
                V_diastolic=vertcat(V_diastolic,V_temp);
            end
        end        
        
        % Compute pulsatility index (PI), resistive index (RI) and vessel wall compliance index (CI)
        nPointsAvgCardiacCycle=size(V_t_AvgCardiac,2);
        [PI_Vmax,PI_Vmax_loc]=nanmax(V_t_AvgCardiac(:));
        PI_Vmin=(V_t_AvgCardiac(1)+V_t_AvgCardiac(end))/2;      
        PI_Vmean=nanmean(V_t_AvgCardiac(:));
        PI=(PI_Vmax-PI_Vmin)/PI_Vmean;
        RI=(PI_Vmax-PI_Vmin)/PI_Vmax;
        TimeCoorsOfCardiacPhases_Seconds_NormAvg=...
            TimeCoorsOfCardiacPhases_Seconds(:,:)-...
            TimeCoorsOfCardiacPhases_Seconds(:,1);
        TimeCoorsOfAvgCardiac=mean(...
            TimeCoorsOfCardiacPhases_Seconds_NormAvg,1);
        CI_t1=TimeCoorsOfAvgCardiac(1,1);
        CI_t2=TimeCoorsOfAvgCardiac(1,PI_Vmax_loc);
        CI_t3=TimeCoorsOfAvgCardiac(1,end);
        CI_t2=CI_t2-CI_t1;
        CI_t3=CI_t3-CI_t1;
        CI_t1=CI_t1-CI_t1;
        CI=(CI_t3-CI_t2)/(CI_t2-CI_t1);
        
        % Average cardiac cycle - Displaying Velocity vs Time (along with 
        % vascular indices PI, RI and CI)
        f_No=f_No+1;
        figure(f_No);
        currentFigure=gcf;
        errorbar(TimeCoorsOfAvgCardiac(1,:)-TimeCoorsOfAvgCardiac(1,1),...
            V_t_AvgCardiac,V_t_AvgCardiac_sem,'-ok','Linewidth',3);axis auto;       
        ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
        xlabel('Time (s)');
        ylabel('Average velocity (mm/s)');
        title({'Average Cardiac Cycle',fileName(1:25)},'Interpreter','none');
        textBox=annotation('textbox');
        textBox.String={['PI = ',num2str(PI,'%.2f')],['RI = ',...
            num2str(RI,'%.2f')],['CI = ',num2str(CI,'%.2f')],['HR = ',...
            num2str((1/CI_t3)*60,'%.1f')]};
        currentFolder=cd(Step5folder);
        print(f_No,[fileName(1:25),'_Step5_Avg Cardiac Cycle',...
            appendToResults],'-dtiff');
        saveas(f_No,[fileName(1:25),'_Step5_Avg Cardiac Cycle',...
            appendToResults,'.fig'],'fig');
        cd(currentFolder);
        
        % Average cardiac cycle - Computing Velocity vs Space as function
        % of cardiac phase
        nX=numel(V_x);
        nPhases=nSamplesPerCycle+2;
        V_x_AllCardiac=NaN(nX,nSamplesPerCycle+2,nCardiacCycles);
        V_x_AllCardiac_std=V_x_AllCardiac;
        V_x_AllCardiac_sem=V_x_AllCardiac;
        V_x_AvgCardiac=NaN(nX,nSamplesPerCycle+2);
        V_x_AvgCardiac_std=V_x_AvgCardiac;        
        V_x_AvgCardiac_sem=V_x_AvgCardiac;

        for k=1:nX
            xStart=SpacePixels(k)-velocityBinningWindowInSpacePixels/2;
            xEnd=SpacePixels(k)+velocityBinningWindowInSpacePixels/2;
            for m=1:nSamplesPerCycle+2
                for p=1:nCardiacCycles
                    tStart=TimeCoorsOfCardiacPhases_Pixels(p,m)-...
                        velocityBinningWindowInTimePixels/2;
                    tEnd=TimeCoorsOfCardiacPhases_Pixels(p,m)+...
                        velocityBinningWindowInTimePixels/2;
                    V_temp=V(y>=xStart & y<xEnd & x>=tStart & x<tEnd);
                    if ~isempty(V_temp)
                        V_x_AllCardiac(k,m,p)=nanmean(V_temp,1);
                        V_x_AllCardiac_std(k,m,p)=nanstd(V_temp,0,1);
                        V_x_AllCardiac_sem(k,m,p)=...
                            V_x_AllCardiac_std(k,m,p)./numel(V_temp);
                    end
                end
            end
        end
        V_x_AvgCardiac=nanmean(V_x_AllCardiac,3);
        V_x_AvgCardiac_std=nanstd(V_x_AllCardiac_std,0,3);
        V_x_AvgCardiac_sem=V_x_AvgCardiac_std./...
            sqrt(size(V_x_AllCardiac_std,3));
    end
    
    % Cropping laminar profile based on lane traffic analysis
    laneTrafficThreshold=(laneTrafficThresholdPercent/100)*...
        nanmax(laneTraffic);
    if cropLaminarProfile==1
        [leftLumenEdge_in_dxStepUnits]=find(laneTraffic>...
            laneTrafficThreshold,1,'first');
        [rightLumenEdge_in_dxStepUnits]=find(laneTraffic>...
            laneTrafficThreshold,1,'last');
        % cropping laminar profile of each phase in avg cardiac cycle
        if computeAverageCardiacCycle==1
            V_x_AvgCardiac_uncropped=V_x_AvgCardiac;
            V_x_AvgCardiac_std_uncropped=V_x_AvgCardiac_std;
            V_x_AvgCardiac_sem_uncropped=V_x_AvgCardiac_sem;
            V_x_AvgCardiac=V_x_AvgCardiac...
                (leftLumenEdge_in_dxStepUnits:rightLumenEdge_in_dxStepUnits,:);
            V_x_AvgCardiac_std=V_x_AvgCardiac_std...
                (leftLumenEdge_in_dxStepUnits:rightLumenEdge_in_dxStepUnits,:);
            V_x_AvgCardiac_sem=V_x_AvgCardiac_sem...
                (leftLumenEdge_in_dxStepUnits:rightLumenEdge_in_dxStepUnits,:);
        end
        % cropping average laminar profile
        V_x_uncropped=V_x;
        V_x_std_uncropped=V_x_std;
        V_x_sem_uncropped=V_x_sem;
        V_x=V_x...
            (1,leftLumenEdge_in_dxStepUnits:rightLumenEdge_in_dxStepUnits);
        V_x_std=V_x_std...
            (1,leftLumenEdge_in_dxStepUnits:rightLumenEdge_in_dxStepUnits);
        V_x_sem=V_x_sem...
            (1,leftLumenEdge_in_dxStepUnits:rightLumenEdge_in_dxStepUnits);
        Space_uncropped=Space;
        SpacePixels_uncropped=SpacePixels;
        laneTraffic_uncropped=laneTraffic;
        Space=Space...
            (1,leftLumenEdge_in_dxStepUnits:rightLumenEdge_in_dxStepUnits);
        SpacePixels=SpacePixels...
            (1,leftLumenEdge_in_dxStepUnits:rightLumenEdge_in_dxStepUnits);
        laneTraffic=laneTraffic...
            (1,leftLumenEdge_in_dxStepUnits:rightLumenEdge_in_dxStepUnits);
    end
    
    % Creating video of velocity profile variation in average cardiac cycle
    if makeAverageCardiacVideo==1 && computeAverageCardiacCycle==1
        currentFolder=cd(Step5folder);
        vidObj = VideoWriter([fileName(1:25),...
            '_Step5_AverageCardiacCycle',appendToResults,'.avi'],...
            'Uncompressed AVI');
        beatsDisplayedPerSecond=1;
        if bpmExist~=0
            bpmDisplay=bpm;
        else
            bpmDisplay=300;
        end
        vidObj.FrameRate=(nSamplesPerCycle+3)*beatsDisplayedPerSecond;
        open(vidObj);
        f_No=f_No+1;
        axYLim=[0 max(max(V_x_AvgCardiac+V_x_AvgCardiac_sem))*1.1];
        axXLim=[2*min(Space(:)) 2*max(Space(:))];
        for k=1:size(V_x_AvgCardiac,2)
            figure(f_No);
            errorbar(Space,V_x_AvgCardiac(:,k),V_x_AvgCardiac_sem(:,k),...
                'b','Linewidth',3);
            ax=gca;ax.YLim=axYLim;ax.XLim=axXLim;
            ax.FontSize=fontsize1;ax.FontWeight='bold';
            timeCounter=(k-1)*(1/bpmDisplay)*60/(nSamplesPerCycle+2);
            legend({[num2str(timeCounter,'%.2f'),' s']},'Location',...
                'northeast','Box','off','TextColor',[0 1 0],'FontSize',...
                fontsize1,'FontWeight','bold');
            title({'Velocity profile';'Mean \pm 1 S.E.M'});
            xlabel('Radial position (\mum)');
            ylabel('Blood velocity (mm/s)');
            cdata=print('-RGBImage');
            writeVideo(vidObj,cdata);
        end
        close(vidObj);
        cd(currentFolder);
    end    
    currentFolder=cd(Step5folder);
    
    % Plotting average velocity vs time (with motion trace)
    f_No=f_No+1;
    figure(f_No);
    errorbar(Time,V_t,V_t_sem,'k','LineWidth',3);axis auto;hold on;
    ylabel('Average velocity (mm/s)');
    ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
    currentXLim=ax.XLim;ax.XLim=[-.01 currentXLim(2)];
    currentYLim=ax.YLim;ax.YLim=[0 currentYLim(2)];    
    currentXLim=ax.XLim;
    yyaxis right
    plot(timeAxisForMotionTrace,motionTrace.*mic_pix,'-o','MarkerSize',3,...
        'LineWidth',3);axis auto;
    ylabel('Eye motion (µm)');
    ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
    ax.XLim=[currentXLim(1) currentXLim(2)];
    bpmExist=exist('bpm','var');
    if bpmExist==0
        title(['Average velocity vs time.',' FFT peak could not be found']);
    else
        title(['Average velocity vs time.',' FFT peak at ',num2str(bpm,...
            '%.1f'),' cycles/min']);
    end        
    xlabel('Time (s)');
    currentFolder=cd(Step5folder);
    print(f_No,[fileName(1:25),'_Step5_Motion_Trace',appendToResults],...
        '-dtiff');
    saveas(f_No,[fileName(1:25),'_Step5_Motion_Trace',appendToResults,...
        '.fig'],'fig');
    cd(currentFolder);
    
    % Plotting average velocity vs time
    f_No=f_No+1;
    figure(f_No);
    ha=errorbar(Time,V_t,V_t_sem,'k','LineWidth',3);axis auto;
    ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
    currentXLim=ax.XLim;ax.XLim=[-.01 currentXLim(2)];
    currentYLim=ax.YLim;ax.YLim=[0 currentYLim(2)];    
    bpmExist=exist('bpm','var');
    if bpmExist==0
        title(['Average velocity vs time.',' FFT peak could not be found']);
    else
        title(['Average velocity vs time.',' FFT peak at ',num2str(bpm,...
            '%.1f'),' cycles/min']);
    end        
    xlabel('Time (s)');
    ylabel('Average velocity (mm/s)');
    currentFolder=cd(Step5folder);
    print(f_No,[fileName(1:25),'_Step5_Velocity_vs_Time',...
        appendToResults],'-dtiff');
    saveas(f_No,[fileName(1:25),'_Step5_Velocity_vs_Time',...
        appendToResults,'.fig'],'fig');
    cd(currentFolder);
    
    % Plotting FFT
    if bpmExist~=0
        f_No=f_No+1;
        figure(f_No);
        plot(Gx_display,Gy_display);
        ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
        xlabel('cycles/min');
        ylabel('FFT of cardiac cycle');
        axis([-1000 1000 0 max(Gy_display(:)*1.1)]);
        currentFolder=cd(Step5folder);
        print(f_No,[fileName(1:25),'_Step5_FFT',appendToResults],...
            '-dtiff');
        saveas(f_No,[fileName(1:25),'_Step5_FFT',appendToResults,...
            '.fig'],'fig');
        cd(currentFolder);
    end

    % Plotting velocity profile along with lane traffic
    f_No=f_No+1;
    figure(f_No);
    if cropLaminarProfile==1
        errorbar(Space_uncropped,V_x_uncropped,V_x_sem_uncropped,...
            'k','Linewidth',3,'LineStyle','--');axis auto;hold on;
    end
    errorbar(Space,V_x,V_x_sem,'k','Linewidth',3);axis auto;hold on;
        
    ax=gca;currentYLim=ax.YLim;ax.YLim=[0 currentYLim(2)];ax.XLim=...
        [2*min(Space(:)) 2*max(Space(:))];
    ax.FontSize=fontsize1;ax.FontWeight='bold';
    title({'Velocity profile';'Mean \pm 1 S.E.M'});
    xlabel('Radial position (\mum)');
    ylabel('Blood velocity (mm/s)');
    yyaxis right
    if cropLaminarProfile==1
        plot(Space_uncropped,laneTraffic_uncropped,...
            '--o','LineWidth',3,'MarkerSize',3,'Color',[0 .45 .74]);
    end
    plot(Space,laneTraffic,'-o','LineWidth',3,'MarkerSize',3,'Color',...
        [0 .45 .74]);
    ylabel('nKernels analyzed');
    currentFolder=cd(Step5folder);
    print(f_No,[fileName(1:25),'_Step5_Lane_Traffic',appendToResults],...
        '-dtiff');
    saveas(f_No,[fileName(1:25),'_Step5_Lane_Traffic',appendToResults,...
        '.fig'],'fig');
    cd(currentFolder);
    
    % Plotting velocity profile along with lane traffic and STD
    f_No=f_No+1;
    figure(f_No);
    if cropLaminarProfile==1
        errorbar(Space_uncropped,V_x_uncropped,V_x_std_uncropped,...
            'k','Linewidth',3,'LineStyle','--');axis auto;hold on;
    end
    errorbar(Space,V_x,V_x_std,'k','Linewidth',3);axis auto;hold on;
        
    ax=gca;currentYLim=ax.YLim;ax.YLim=[0 currentYLim(2)];ax.XLim=...
        [2*min(Space(:)) 2*max(Space(:))];
    ax.FontSize=fontsize1;ax.FontWeight='bold';
    title({'Velocity profile';'Mean \pm 1 STD'});
    xlabel('Radial position (\mum)');
    ylabel('Blood velocity (mm/s)');
    yyaxis right
    if cropLaminarProfile==1
        plot(Space_uncropped,laneTraffic_uncropped,...
            '--o','LineWidth',3,'MarkerSize',3,'Color',[0 .45 .74]);
    end
    plot(Space,laneTraffic,'-o','LineWidth',3,'MarkerSize',3,'Color',...
        [0 .45 .74]);
    ylabel('nKernels analyzed');
    currentFolder=cd(Step5folder);
    print(f_No,[fileName(1:25),'_Step5_Lane_Traffic_STD',...
        appendToResults],'-dtiff');
    saveas(f_No,[fileName(1:25),'_Step5_Lane_Traffic_STD',...
        appendToResults,'.fig'],'fig');
    cd(currentFolder);
    
    % Plotting velocity profile with sem
    f_No=f_No+1;
    figure(f_No);
    errorbar(Space,V_x,V_x_sem,'k','Linewidth',3);axis auto;hold on;
    ax=gca;currentYLim=ax.YLim;ax.YLim=[0 currentYLim(2)];
    if numel(Space)>1
        ax.XLim=[2*min(Space(:)) 2*max(Space(:))];
    end
    ax.FontSize=fontsize1;ax.FontWeight='bold';
    title({'Velocity profile';'Mean \pm 1 S.E.M'});
    xlabel('Radial position (\mum)');
    ylabel('Blood velocity (mm/s)');
    currentFolder=cd(Step5folder);
    print(f_No,[fileName(1:25),'_Step5_Velocity_vs_Space',...
        appendToResults],'-dtiff');
    saveas(f_No,[fileName(1:25),'_Step5_Velocity_vs_Space',...
        appendToResults,'.fig'],'fig');
    cd(currentFolder);

    % Plotting velocity profile with std
    f_No=f_No+1;
    figure(f_No);
    errorbar(Space,V_x,V_x_std,'k','Linewidth',3);axis auto;hold on;
    ax=gca;currentYLim=ax.YLim;ax.YLim=[0 currentYLim(2)];
    if numel(Space)>1
        ax.XLim=[2*min(Space(:)) 2*max(Space(:))];
    end
    ax.FontSize=fontsize1;ax.FontWeight='bold';
    title({'Velocity profile';'Mean \pm 1 STD'});
    xlabel('Radial position (\mum)');
    ylabel('Blood velocity (mm/s)');
    currentFolder=cd(Step5folder);
    print(f_No,[fileName(1:25),'_Step5_Velocity_vs_Space_STD',...
        appendToResults],'-dtiff');
    saveas(f_No,[fileName(1:25),'_Step5_Velocity_vs_Space_STD',...
        appendToResults,'.fig'],'fig');
    cd(currentFolder);

    f_No=f_No+1;
    figure(f_No);
    histogram(reshape(SNRoverlay_SNR,[1 numel(SNRoverlay_SNR)]),...
        1000);hold on;
    xlabel('SNR');ylabel('Count');
    if SNRotsu==1
        title(['Otsu EM = ',num2str(EM,'%.2f'),', SNRthreshold = ',...
            num2str(SNRthreshold,'%.2f'),' (Otsu)']);
    else
        title(['SNRthreshold = ',num2str(SNRthreshold,'%.2f'),...
            ' (Manual)']);
    end
    textBox=annotation('textbox');
    textBox.String={['nKernels shown = ',num2str(numel(SNRoverlay_SNR),...
        '%.0f')],...
        ['% accepted = ',num2str(percentKernelsPassedStage4*100,...
        '%.0f'),' %']};
    textBox.Position(1)=.55;textBox.Position(2)=.75;
    ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';axYLim=ax.YLim;
    plot([SNRthreshold SNRthreshold],axYLim,'-r','LineWidth',2);
    currentFolder=cd(Step5folder);
    print(f_No,[fileName(1:25),'_Step5_SNR_histogram',appendToResults],...
        '-dtiff');
    saveas(f_No,[fileName(1:25),'_Step5_SNR_histogram',appendToResults,...
        '.fig'],'fig');
    cd(currentFolder);
    f_No=f_No+1;figure(f_No);
    histogram(reshape(V,[1 numel(V)]),1000);
    title(['Average Velocity = ',num2str(AvergeVelocity,'%.1f'),' \pm ',...
        num2str(StdVelocity,'%.1f'),' mm/s']);
    xlabel('Velocity (mm/s)');ylabel('Count');
    textBox=annotation('textbox');
    textBox.String={['nKernels shown = ',num2str(numel(V),'%.0f')],...
        ['% of initial nKernels = ',...
        num2str(100*numel(V)/(ndx*ndt),'%.2f'),' %']};
    textBox.Position(1)=.5;textBox.Position(2)=.75;
    ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
    ax.XLim=[Vmin Vmax];
    currentFolder=cd(Step5folder);
    print(f_No,[fileName(1:25),'_Step5_Velocity_histogram',...
        appendToResults],'-dtiff');
    saveas(f_No,[fileName(1:25),'_Step5_Velocity_histogram',...
        appendToResults,'.fig'],'fig');
    cd(currentFolder);
    
    % Angle histogram
    f_No=f_No+1;figure(f_No);
    histogram(reshape(theta_maxClean,[1 numel(theta_maxClean)]),1000);
    title(['Average Angle = ',num2str(nanmean(theta_maxClean),'%.1f'),...
        ' \pm ',num2str(nanstd(theta_maxClean,0,1),'%.1f'),' \circ']);
    xlabel('Angle (°)');ylabel('Count');
    textBox=annotation('textbox');
    textBox.String={['nKernels shown = ',num2str(numel(theta_maxClean),...
        '%.0f')],...
        ['% of initial nKernels = ',...
        num2str(100*numel(theta_maxClean)/(ndx*ndt),'%.2f'),' %']};
    textBox.Position(1)=.2;textBox.Position(2)=.75;
    ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
    ax.XLim=[0 180];
    currentFolder=cd(Step5folder);
    print(f_No,[fileName(1:25),'_Step5_Angle_histogram',...
        appendToResults],'-dtiff');
    saveas(f_No,[fileName(1:25),'_Step5_Angle_histogram',...
        appendToResults,'.fig'],'fig');
    cd(currentFolder);

    % Creating video of all velocity profile variation with each individual cardiac cycle
    if makeAllCardiacVideo==1
        currentFolder=cd(Step5folder);
        vidObj = VideoWriter([fileName(1:25),...
            '_Step5_PulsatilityInProfile',appendToResults,'.avi'],'Uncompressed AVI');
        beatsDisplayedPerSecond=1;
        if bpmExist~=0
            bpmDisplay=bpm;
        else
            bpmDisplay=300;
        end
        vidObj.FrameRate=beatsDisplayedPerSecond*(1/(bpmDisplay/60))/...
            (velocityBinningStepSizeInTimePixels/freq);
        open(vidObj);
        f_No=f_No+1;
        axYLim=[0 nanmax(V_xt(:)).*1.2];
        axXLim=[1.5*min(Space(:)) 1.5*max(Space(:))];
        for k=1:numel(V_t)
            figure(f_No);
            errorbar(Space,V_xt(:,k),V_xt_sem(:,k),'b','Linewidth',3);
            ax=gca;ax.YLim=axYLim;ax.XLim=axXLim;
            ax.FontSize=fontsize1;ax.FontWeight='bold';
            timeCounter=(k-1)*velocityBinningStepSizeInTimePixels/freq;
            legend({[num2str(timeCounter,'%.2f'),' s']},'Location',...
                'northeast','Box','off','TextColor',[0 1 0],'FontSize',...
                fontsize1,'FontWeight','bold');
            title({'Velocity profile';'Mean \pm 1 S.E.M'});
            xlabel('Radial position (\mum)');
            ylabel('Blood velocity (mm/s)');
            cdata=print('-RGBImage');
            writeVideo(vidObj,cdata);
        end
        close(vidObj);
        cd(currentFolder);
    end
    
    % Saving copy of code (added 06-05-2018, 4 pm)
    codeFileFullPath=mfilename('fullpath'); %does not store .m extension
    [codeFilePath,codeFileName,~]=fileparts(codeFileFullPath);
    newCodeFileName=['code_',fileName(1:25),'_Step5',appendToResults];
    codeSource=[codeFilePath,'\',codeFileName,'.m'];
    codeDestination=[Step5folder,newCodeFileName,'.m'];
    copyfile(codeSource,codeDestination);
    
    currentFolder=cd(Step5folder);
    save([fileName(1:25),'_Step5',appendToResults]);
    cd(currentFolder);
    % moving processed file to separate folder
    if moveFileAfterProcessing==1
        movefile([processDir,fileName],[processDir,'PROCESSED\',fileName]);
    end
    elapsedTime=toc;disp(['Elapsed time is ',num2str(elapsedTime/60),...
        ' minutes.']);
end
beep
