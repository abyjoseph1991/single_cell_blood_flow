% Velocity_Step5_PostRadon
% Author: Aby Joseph, University of Rochester
% License: GPL-3.0-or-later
% Last modified: 01-24-2019

clear;close all;
f_No=0;
fontsize1=18;

% Variables typically not changed
moveFileAfterProcessing=0;
SNRotsu=0; % If 1, calculates SNRthreshold based on Otsu threshold. Otherwise uses pre-determined values in code
SNRthreshold=2.5; %for SNRotsu==0
SNRcolorMode=2; % 1 = binary cyan/magenta, 2 = spectrum
overlayLineWidth=4;
velocityInterpolation=0;
removeAberrantVelocities=0; % use only with high density kernel sampling
SNRmode=1; %Mode 1 - calculate from peak variance of Radon
removeLonelyKernels=0; % works only when removeAberrantVelocities==1
videoDurationToAnalyze_Step5=1;
findFFT=0;
writeSNRimage=0;
writeAcceptRejectimage=0;
writeVelocityImage=0;
writeTimeCroppedRegisteredImage=0;
writeVelocitySquareImage=0;
computeAverageCardiacCycle=0;
makeAverageCardiacVideo=0;
makeAllCardiacVideo=0;
printWhiteLine=1;
SNRmapRangeManualOverride=1; % SNR color map, select '1' to override default clamps
minSNR=[]; % if overridden above, defines custom colormap min clamp. Insert empty to keep default
maxSNR=3.4; % if overridden above, defines custom colormap max clamp. Insert empty to keep default
manualPeakTroughEntry=0; % use manually determined cardiac peaks and troughs

thresholdMinVelocity=0.01; % in mm/s. Velocities below this are ignored.
minV=0; % in mm/s. Lower clamp for Jet Colormap
maxV=15; % in mm/s. Upper clamp for Jet Colormap
binningWindowTime=15/1000; %in seconds (Recommended: 15 ms for mouse and 90 ms for human)

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
        'SNRofEachKernel_Mode1',...
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
    nFinite_pre7=sum(isfinite(V_pixel));
    thresholdMaxVelocity_pixels=nanmean(V_pixel(:))+20.*nanstd(V_pixel(:));
    thresholdMaxVelocity=thresholdMaxVelocity_pixels.*(mic_pix*freq*1e-3)
    [row,column]=find(V_pixel>thresholdMaxVelocity_pixels);
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
        
    currentFolder=cd(Step5folder);
    
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
