% Velocity_Step2_LinescanRegistration
% Author: Aby Joseph, University of Rochester
% License: GPL-3.0-or-later
% Last modified: 01-24-2019

close all;
clear;
fNo=0;
fontsize1=18;

videoDurationToAnalyze_Step2=NaN; % in seconds
registrationON=1;
backgroundSubtractionON=1;
moveFileAfterProcessing=1;
normalizeFlag=1;
insertNaNsWhenPoorRegistration=1;

sobelFiltering=0; % to apply Sobel Filtering, make value 1.
crossCorrelationCoeffcientThreshold=0.50; % For mouse, even 0.97 retains most strips in high quality mouse data
sizeStripRegistration=601; % registration strip size. 
timeWindowAveraging=101;% width (in pixels) of temporal window for summation of common information. Must be odd number. 

Step5folder='E:\PhD\TRBF\Experimental DATA\Step 5 MAT files\';
Step4folder='E:\PhD\TRBF\Experimental DATA\Step 4 MAT files\';
Step2folder='E:\PhD\TRBF\Experimental DATA\Step 2 MAT files\';
RawDir='E:\PhD\TRBF\Experimental DATA\RAW MAT files\';

% Append current date and time to filenames of results
currentDateAndTime=clock;
year=num2str(currentDateAndTime(1));
month=num2str(currentDateAndTime(2),'%02d');
day=num2str(currentDateAndTime(3),'%02d');
hour=num2str(currentDateAndTime(4),'%02d');
minute=num2str(currentDateAndTime(5),'%02d');
appendToResults=['_',year,month,day,hour,minute];

namecontains='*_RAW*.mat';
[dirFilenames(1).name,filePath,~]=uigetfile([RawDir,namecontains],'Select file to process - click cancel to process all files in starting folder');
if filePath==0
    currentFolder=cd(RawDir);
    dirFilenames=dir([namecontains]);
    cd(currentFolder);
    nFiles=numel(dirFilenames);
    processDir=RawDir;
else
    processDir=[filePath];
    nFiles=1;
end
mkdir(processDir,'PROCESSED'); % the variable processDir is usually the same as the variable RawDir, except when a single selceted file (instead of a batch) needs to be processed in this step

for fileNo=1:nFiles
    close all;
    fileName=dirFilenames(fileNo).name;
    display(fileName);
    load([processDir,fileName],'desinusoidON','postFPGA',...
            'linescanDir','nFrames','data','height','width',...
            'mic_pix','freq');
    if isnan(videoDurationToAnalyze_Step2)
        videoDurationToAnalyze_Step2=size(data,2)/freq;
    end    
    
    data=single(data);
    data=data-min(data(:));
    data=data./max(data(:));
    data=data(:,1:floor(freq*videoDurationToAnalyze_Step2));

    % registration
    dataHeight_Step2=size(data,1);
    dataWidth_Step2=size(data,2);
    dataUnregisteredWithStaticLines=data;

    screenSize=get(0, 'MonitorPositions');
    fNo=fNo+1;
    figure(fNo);
    currentFigure=gcf;

    display1=dataUnregisteredWithStaticLines;
    imagesc(display1);colormap(gray);axis off;hold on;
    currentFigure.Name=[fileName,' - Click on two spatial extents of motion'];
    currentFigure.OuterPosition(1)=70;
    currentFigure.OuterPosition(2)=screenSize(4)-1000-1;
    currentFigure.OuterPosition(3)=screenSize(3)-70-1;
    currentFigure.OuterPosition(4)=1000;
    currentOuterPosition=currentFigure.OuterPosition;
    currentAxes=gca;
    currentAxes.Position=[0 0 1 1];
    line([dataWidth_Step2/2 dataWidth_Step2/2],[1 height],'LineWidth',1,...
        'Color',[1 0 0]);
    line([dataWidth_Step2/2+sizeStripRegistration dataWidth_Step2/...
        2+sizeStripRegistration],[1 height],'LineWidth',1,'Color',[1 0 0]);
    [~,y]=ginput(2);
    y=round(y);
    if numel(y)==2
        y=sort(y);
        xTop=y(1);
        xBottom=y(2);
    else
        if numel(y)==1
            if y(1)<size(display1,1)/2
                xTop=y(1);
                xBottom=size(display1,1);
            else
                xTop=1;
                xBottom=y(1);
            end
        else
            if numel(y)==0
                xTop=1;
                xBottom=size(display1,1);
            end
        end
    end
    figure(fNo);
    display1=dataUnregisteredWithStaticLines;
    currentFigure=gcf;
    currentFigure.Name='Press Enter to continue...';
    figure(fNo);    
    imagesc(display1);colormap(gray);axis off;hold on;
    line([dataWidth_Step2/2 dataWidth_Step2/2],[1 height],'LineWidth',1,...
        'Color',[1 0 0]);
    line([dataWidth_Step2/2+sizeStripRegistration dataWidth_Step2/2+...
        sizeStripRegistration],[1 height],'LineWidth',1,'Color',[1 0 0]);
    clear display1;
    figure(fNo);
    currentFigure=gcf;
    currentFigure.Name='Click on reference strip';
    display1=dataUnregisteredWithStaticLines;
    imagesc(display1);colormap(gray);axis off;hold on;
    line([dataWidth_Step2/2 dataWidth_Step2/2],[1 height],'LineWidth',1,...
        'Color',[1 0 0]);
    line([dataWidth_Step2/2+sizeStripRegistration dataWidth_Step2/2+...
        sizeStripRegistration],[1 height],'LineWidth',1,'Color',[1 0 0]);
    clear display1;

    [x,y]=ginput(1);
    referenceLine=round(x);

    if registrationON==1
        nStrips=floor(dataWidth_Step2/sizeStripRegistration);
        halfSizeStripRegistration=(sizeStripRegistration-1)/2;
        line_avg0=mean(dataUnregisteredWithStaticLines(:,(referenceLine-...
            halfSizeStripRegistration):(referenceLine+...
            halfSizeStripRegistration)),2);
        if normalizeFlag==1
            line_avg0=line_avg0-min(line_avg0(:));
            line_avg0=line_avg0./max(line_avg0(:));
        end
        line_std0=std(dataUnregisteredWithStaticLines(:,(referenceLine-...
            halfSizeStripRegistration):(referenceLine+...
            halfSizeStripRegistration)),[],2);
        if normalizeFlag==1
            line_std0=line_std0-min(line_std0(:));
            line_std0=line_std0./max(line_std0(:));
        end
        motionTrace=zeros(1,dataWidth_Step2-sizeStripRegistration-1);
        hWaitbar2=waitbar(0,'Registration 0 % completed');
        count1=0;
        for k=(1+halfSizeStripRegistration):(dataWidth_Step2-...
                halfSizeStripRegistration)
            count1=count1+1;
            if mod(k,1000)==0
                waitbar(k/(dataWidth_Step2),hWaitbar2,['Registration ',...
                    num2str(single(100*k/(dataWidth_Step2)),'%.0f'),...
                    ' % completed']);
            end
            line_std1=std(dataUnregisteredWithStaticLines(:,(k-...
                halfSizeStripRegistration):(k+halfSizeStripRegistration)),...
                [],2);
            if normalizeFlag==1
                line_std1=line_std1-min(line_std1(:));
                line_std1=line_std1./max(line_std1(:));
            end
            [r,lags]=xcorr(line_std0(xTop:xBottom),line_std1(xTop:...
                xBottom),'coeff');
            if insertNaNsWhenPoorRegistration==1 && max(r(:))<...
                    crossCorrelationCoeffcientThreshold
                data(:,k)=NaN(dataHeight_Step2,1);
                motionTrace(count1)=NaN;
            else
                [~,locs]=max(r);
                data(:,k)=...
                    circshift(dataUnregisteredWithStaticLines(:,k),...
                    [lags(locs) 0]);
                motionTrace(count1)=lags(locs);
            end
        end
        motionTrace=horzcat(NaN(1,halfSizeStripRegistration),motionTrace,...
            NaN(1,halfSizeStripRegistration));
        [~,kExample]=nanmax(motionTrace);
        line_avg1=mean(dataUnregisteredWithStaticLines(:,(kExample-...
            halfSizeStripRegistration):kExample+halfSizeStripRegistration),2);
        if normalizeFlag==1
            line_avg1=line_avg1-min(line_avg1(:));
            line_avg1=line_avg1./max(line_avg1(:));
        end
        line_std1=std(dataUnregisteredWithStaticLines(:,(kExample-...
            halfSizeStripRegistration):kExample+halfSizeStripRegistration),[],2);
        if normalizeFlag==1
            line_std1=line_std1-min(line_std1(:));
            line_std1=line_std1./max(line_std1(:));
        end
        screenSize=get(0, 'MonitorPositions');
    else
        motionTrace=zeros(1,dataWidth_Step2-sizeStripRegistration-1);
    end
    
    dataRegisteredWithStaticLines=data;
    
    % removing the vertical lines (static features) in linescan 
    
    if backgroundSubtractionON==1
        hWaitbar1=waitbar(0,'Removing static features 0 %');
        dataSmoothed=zeros(size(data,1),size(data,2));
        for k=1:dataHeight_Step2
            waitbar(k/dataHeight_Step2,hWaitbar1,['Removing static features ',...
                num2str(single(100*k/dataHeight_Step2),'%.0f'),' % completed']);
            dataSmoothed(k,:) = smooth(data(k,:),timeWindowAveraging,...
                'moving')'; % moving average filter
        end
        
        data=data-dataSmoothed;
        clear dataSmoothed;
        data=data-min(data(:));
        data=data./max(data(:));
    else
        data=data-min(data(:));
        data=data./max(data(:));
    end
    
    % Sobel filtering (optional, not currently implemented)
    if sobelFiltering==1
        h=fspecial('sobel');
        data=imfilter(data,h');
        data=data-min(data(:));
        data=data./max(data(:));
    end
        
    displayData=vertcat(dataUnregisteredWithStaticLines,...
        dataRegisteredWithStaticLines,data);

    fNo=fNo+1;
    figure(fNo);
    currentFigure=gcf;
    imagesc(displayData);colormap(gray);axis off;hold on;
    currentFigure.OuterPosition(1)=70;
    currentFigure.OuterPosition(2)=screenSize(4)-1000-1;
    currentFigure.OuterPosition(3)=screenSize(3)-70-1;
    currentFigure.OuterPosition(4)=1000;
    currentOuterPosition=currentFigure.OuterPosition;
    currentAxes=gca;
    currentAxes.Position=[0 0 1 1];
    line([1 dataWidth_Step2],[height-0.5 height-0.5],'LineWidth',1,...
        'Color',[1 0 0]);
    line([1 dataWidth_Step2],[2*height-0.5 2*height-0.5],'LineWidth',1,...
        'Color',[1 0 0]);
    
    timeAxisForMotionTrace=(0:numel(motionTrace)-1)./freq;
    if registrationON==1
        fNo=fNo+1;
        figure(fNo);plot(timeAxisForMotionTrace,motionTrace.*mic_pix,'-o',...
            'MarkerSize',2);axis tight;
        ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
        currentXLim=ax.XLim;ax.XLim=[currentXLim(1) timeAxisForMotionTrace(end)];
        xlabel('Time (s)','FontSize',fontsize1,'FontWeight','bold');
        ylabel('Displacement (\mum)','FontSize',fontsize1,'FontWeight','bold');
        currentFigure=gcf;
        currentFigure.OuterPosition(1)=700;
        currentFigure.OuterPosition(2)=580;
        currentFigure.OuterPosition(3)=600;
        currentFigure.OuterPosition(4)=500;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        currentFolder=cd(Step5folder);
        print(fNo,[fileName(1:25),'_Step2_Motion_trace',appendToResults],'-dtiff');
        saveas(fNo,[fileName(1:25),'_Step2_Motion_trace',appendToResults],'fig');
        cd(currentFolder);
        
        fNo=fNo+1;
        spaceCoor=(1:height);spaceCoor=spaceCoor-min(spaceCoor(:));
        figure(fNo);plot(spaceCoor,line_std1,'-r');
        hold on;plot(spaceCoor,line_std0,'-k');axis tight;
        ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
        xlabel('Space (pixels)','FontSize',fontsize1,'FontWeight','bold');
        ylabel('S.D.','FontSize',fontsize1,'FontWeight','bold');
        legend({'Sample strip','Reference strip'},'FontSize',fontsize1-8);
        currentFigure=gcf;
        currentFigure.OuterPosition(1)=400;
        currentFigure.OuterPosition(2)=50;
        currentFigure.OuterPosition(3)=600;
        currentFigure.OuterPosition(4)=500;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        
        fNo=fNo+1;
        spaceCoor=(1:height);spaceCoor=spaceCoor-min(spaceCoor(:));
        figure(fNo);plot(spaceCoor,line_avg1,'-r');
        hold on;plot(spaceCoor,line_avg0,'-k');axis tight;
        ax=gca;ax.FontSize=fontsize1;ax.FontWeight='bold';
        xlabel('Space (pixels)','FontSize',fontsize1,'FontWeight','bold');
        ylabel('Intensity','FontSize',fontsize1,'FontWeight','bold');
        legend({'Sample strip','Reference strip'},'FontSize',fontsize1-8);
        currentFigure=gcf;
        currentFigure.OuterPosition(1)=1010;
        currentFigure.OuterPosition(2)=50;
        currentFigure.OuterPosition(3)=600;
        currentFigure.OuterPosition(4)=500;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        figure(2);figure(3);
        drawnow;
        
        close(hWaitbar2);drawnow;
    end
    if backgroundSubtractionON==1
        close(hWaitbar1);drawnow;
    end
    
    % Added 07-31-2018 1:37 pm
    quadrant=input('Enter Radon quadrant: ');
    
    input('Press Enter to continue...');

    % Determining analysis boundaries
    screenSize=get(0, 'MonitorPositions');
    fNo=fNo+1;
    figure(fNo);
    currentFigure=gcf;
    display1=dataRegisteredWithStaticLines;
    imagesc(display1);colormap(gray);axis off;hold on;
    currentFigure.Name='Last Step! Click on spatial extents for Radon analysis in Step 4';
    currentFigure.OuterPosition(1)=70;
    currentFigure.OuterPosition(2)=screenSize(4)-1000-1;
    currentFigure.OuterPosition(3)=screenSize(3)-70-1;
    currentFigure.OuterPosition(4)=1000;
    currentOuterPosition=currentFigure.OuterPosition;
    currentAxes=gca;
    currentAxes.Position=[0 0 1 1];
    input('');
    figure(fNo);    
    [~,y]=ginput(2);
    y=round(y);
    if numel(y)==2
        y=sort(y);
        xTop=y(1);
        xBottom=y(2);
    else
        if numel(y)==1
            if y(1)<size(display1,1)/2
                xTop=y(1);
                xBottom=size(display1,1);
            else
                xTop=1;
                xBottom=y(1);
            end
        else
            if numel(y)==0
                xTop=1;
                xBottom=size(display1,1);
            end
        end
    end
    analysisBoundaries=[xTop xBottom];
    clear display1;    
    dataWrite=data;
    data=uint8(data.*255);
    dataUnregisteredWithStaticLines=uint8(dataUnregisteredWithStaticLines.*255);
    dataRegisteredWithStaticLines=uint8(dataRegisteredWithStaticLines.*255);    

    if registrationON==1
        save([Step2folder,fileName(1:25),'_Step2',appendToResults],'processDir',...
            'fileName','insertNaNsWhenPoorRegistration',...
            'crossCorrelationCoeffcientThreshold','referenceLine',...
            'timeWindowAveraging','sizeStripRegistration',...
            'nStrips','dataUnregisteredWithStaticLines',...
            'dataRegisteredWithStaticLines','line_avg0','line_std0',...
            'line_avg1','line_std1','timeAxisForMotionTrace','motionTrace',...
            'data','analysisBoundaries','videoDurationToAnalyze_Step2',...
            'registrationON','backgroundSubtractionON','quadrant');
    else
        save([Step2folder,fileName(1:25),'_Step2',appendToResults],'processDir',...
            'fileName','insertNaNsWhenPoorRegistration',...
            'crossCorrelationCoeffcientThreshold','referenceLine',...
            'timeWindowAveraging','sizeStripRegistration',...
            'dataUnregisteredWithStaticLines',...
            'dataRegisteredWithStaticLines',...
            'timeAxisForMotionTrace','motionTrace',...
            'data','analysisBoundaries','videoDurationToAnalyze_Step2',...
            'registrationON','backgroundSubtractionON','quadrant');
    end        
    fNo=0;
    % moving processed file to separate folder
    if moveFileAfterProcessing==1
        movefile([processDir,fileName],[processDir,'PROCESSED\',fileName]);
    end

    % Writing TIFF file        
    dataWrite(isnan(dataWrite))=0; % saved form of "data" still has NaNs
    dataWrite=dataWrite.*255;
    dataColor=uint8(cat(3,dataWrite,dataWrite,dataWrite));
    currentFolder=cd(Step5folder);
    imwrite(dataColor,[fileName(1:25),'_Step2_Registered',appendToResults,...
        '.tif'],'TIFF');
    cd(currentFolder);    

    % Saving copy of code (added 07-31-2018, 1.25 pm)
    codeFileFullPath=mfilename('fullpath'); %does not store .m extension
    [codeFilePath,codeFileName,~]=fileparts(codeFileFullPath);
    newCodeFileName=['code_',fileName(1:25),'_Step2',appendToResults];
    codeSource=[codeFilePath,'\',codeFileName,'.m'];
    codeDestination=[Step5folder,newCodeFileName,'.m'];
    copyfile(codeSource,codeDestination);  
    
end
beep
