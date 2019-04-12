% Velocity_Step3_VesselOrientation
% Author: Aby Joseph, University of Rochester
% License: GPL-3.0-or-later
% Last modified: 01-24-2019

close all;
clear;
fNo=0;
fontsize1=26;

motionContrastMode=1; % Mode #1 similar to ImageJ definition (standard deviation). 
                        % Mode #2 computes mean of conventional ("ImageJ")
                        % definition and a fractional standard deviaiton

moveFileAfterProcessing=1;
mainDataFolder='\\cvsnas3.urmc-sh.rochester.edu\aria\Mouse\Aby_Data\TRBF_data\Experimental DATA\';
Step5folder='E:\PhD\TRBF\Experimental DATA\Step 5 MAT files\';
Step3folder='E:\PhD\TRBF\Experimental DATA\Step 3 MAT files\';
Step2folder='E:\PhD\TRBF\Experimental DATA\Step 2 MAT files\';
RawDir_Processed='E:\PhD\TRBF\Experimental DATA\RAW MAT files\PROCESSED\';

% Added 07-31-2018 1.56 pm
% Append current date and time to filenames of results
currentDateAndTime=clock;
year=num2str(currentDateAndTime(1));
month=num2str(currentDateAndTime(2),'%02d');
day=num2str(currentDateAndTime(3),'%02d');
hour=num2str(currentDateAndTime(4),'%02d');
minute=num2str(currentDateAndTime(5),'%02d');
appendToResults=['_',year,month,day,hour,minute];

namecontains='*_Step2*.mat';
[dirFilenames(1).name,filePath,~]=uigetfile([Step2folder,namecontains],...
    'Select file to process - click cancel to process all files in parent folder');
if filePath==0
    dirFilenames=dir([Step2folder,namecontains]);
    nFiles=numel(dirFilenames);
    processDir=Step2folder;
else
    processDir=[filePath];
    nFiles=1;
end
mkdir(processDir,'PROCESSED'); % the variable processDir is usually the same as the variable Step2folder, except when a single selceted file (instead of a batch) needs to be processed in this step

screenSize=get(0, 'MonitorPositions');

for fileNo=1:nFiles
    close all;
    fileName=dirFilenames(fileNo).name;
    display(fileName);
    vesselAngle=input('Enter vessel angle (°): ');
    
    if isempty(vesselAngle)
        angleManuallyEntered=0;

        % Loading RAW.mat
        namecontains2=[fileName(1:25),'*RAW*.mat'];
        dirFilenames2=dir([RawDir_Processed,namecontains2]);
        nFiles2=numel(dirFilenames2);
        if nFiles2==1
            load([RawDir_Processed,dirFilenames2(1).name],'mic_pix',...
                'nFrames','height','width');
        else
            [RawFileName,RawFilePath,~]=uigetfile([RawDir_Processed,...
                namecontains2],'Select RAW MAT file');
            load([RawFilePath,RawFileName],'mic_pix','nFrames',...
                'height','width');
        end

        load([processDir,fileName],'nFramesToProcess','data',...
            'referenceLine','sizeStripRegistration');
        data=single(data);
        data=data./255;
        
        % loading registered cartesian video
        halfSizeStripRegistration=(sizeStripRegistration-1)/2;
        A=exist('cartesianVideoFilePath','var'); % essentially checks if this is the first video in batch process
        if A==0
            [cartesianVideoFile,cartesianVideoFilePath,~]=uigetfile(...
                [mainDataFolder,'*.avi'],...
                ['Select registered cartesian video (same FOV as linescan) for ',...
                fileName(1:end-10)]);
        else
            [cartesianVideoFile,cartesianVideoFilePath,~]=uigetfile(...
                [cartesianVideoFilePath,'\','*.avi'],...
                ['Select registered cartesian video (same FOV as linescan) for ',...
                fileName(1:end-10)]);
        end
        display([cartesianVideoFilePath,cartesianVideoFile]);
        mov1=VideoReader([cartesianVideoFilePath,'\',cartesianVideoFile]);
        
        heightCartesian=mov1.Height;
        widthCartesian=mov1.Width;
        if widthCartesian>600
            displayMagnification=600/widthCartesian;
        else
            displayMagnification=1;
        end
        frameRate=mov1.FrameRate;
        k=1;
        while hasFrame(mov1)
            if k==1
                cartesianVideoData=readFrame(mov1);
            else
                cartesianVideoData(:,:,k)=readFrame(mov1);
            end
            k=k+1;
        end
        cartesianVideoData=single(cartesianVideoData);
        nFramesCartesian=k-1;
        minFramesforSTD=6;
        if nFramesCartesian>=minFramesforSTD
            cartesianFrame=std(cartesianVideoData(:,:,2:end),0,3); % ignoring first frame due to unknown DeMotion bug in first frame of final registered video
            % normalizing pixel values betweeen 0 and 1
            cartesianFrame=cartesianFrame-min(cartesianFrame(:));
            cartesianFrame=cartesianFrame./max(cartesianFrame(:));
            if motionContrastMode==2
                cartesianFrameGlobal=cartesianFrame;
                cartesianFrameFractional=zeros(size(cartesianFrame));
                for k=2:nFramesCartesian-1
                    cartesianFrameFractional=cartesianFrameFractional+...
                        abs(cartesianVideoData(:,:,k+1)-...
                        cartesianVideoData(:,:,k))./(cartesianVideoData(...
                        :,:,k+1)+cartesianVideoData(:,:,k));
                end
                % replacing NaNs with zeros
                cartesianFrameFractional(isnan(cartesianFrameFractional))=0;
                % normalizing pixel values betweeen 0 and 1
                cartesianFrameFractional=cartesianFrameFractional-...
                    min(cartesianFrameFractional(:));
                cartesianFrameFractional=cartesianFrameFractional./...
                    max(cartesianFrameFractional(:));
                % final "motion contrast" image is average of intensity S.D. and Fractional variance images
                cartesianFrame=cartesianFrameFractional+cartesianFrameGlobal;
                % normalizing pixel values betweeen 0 and 1
                cartesianFrame=cartesianFrame-min(cartesianFrame(:));
                cartesianFrame=cartesianFrame./max(cartesianFrame(:));
            end
        else
            if nFramesCartesian==1
                cartesianFrame=cartesianVideoData;
            else
                cartesianFrame=mean(cartesianVideoData(:,:,2:end),3);
            end
            % normalizing pixel values betweeen 0 and 1
            cartesianFrame=cartesianFrame-min(cartesianFrame(:));
            cartesianFrame=cartesianFrame./max(cartesianFrame(:));
        end
        clear cartesianVideoData mov1;
        
        % Cropping cartesianFrame to get rid of black bars (white bars in STD
        % image) due to registeration leftover (added 2017-11-21)
        fNo=fNo+1;
        figure(fNo);
        currentFigure=gcf;
        imshow(cartesianFrame,'Border','tight','InitialMagnification',170);
        colormap(gray);axis off;hold on;
        currentOuterPosition=currentFigure.OuterPosition;
        currentFigure.OuterPosition(1)=500;
        currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(1,4)+1;
        currentOuterPosition=currentFigure.OuterPosition;
        currentFigure.NumberTitle='off';
        currentFigure.Name='Click on top left and bottom right to crop';
        [cropX,cropY]=ginput(2);
        cropX=round(cropX);cropY=round(cropY);
        cartesianFrameUncropped=cartesianFrame;
        cartesianFrame=cartesianFrame(cropY(1):cropY(2),cropX(1):cropX(2));
        heightCartesian=size(cartesianFrame,1);
        widthCartesian=size(cartesianFrame,2);
        % normalizing pixel values betweeen 0 and 1
        cartesianFrame=cartesianFrame-min(cartesianFrame(:));
        cartesianFrame=cartesianFrame./max(cartesianFrame(:));
        close(fNo);fNo=fNo-1;
        
        % padding cartesianFrame symmetrically with zeros to have same width (space axis) as linescan
        widthLinescan=height;
        if widthCartesian<widthLinescan
            leftPad=floor((widthLinescan-widthCartesian)/2);
            rightPad=ceil((widthLinescan-widthCartesian)/2);
            cartesianFramePadded=padarray(cartesianFrame,[0 leftPad],'pre');
            cartesianFramePadded=padarray(cartesianFramePadded,[0 rightPad],'post');
        else
            cartesianFramePadded=cartesianFrame;
        end
        linescanFrame=data(:,referenceLine-halfSizeStripRegistration:...
            referenceLine+halfSizeStripRegistration);
        linescanFrame=rot90(linescanFrame,-1);
        stdLinescan_1D=nanstd(data,0,2);
        stdLinescan_1D=rot90(stdLinescan_1D,-1);
        stdLinescan_1D=stdLinescan_1D-min(stdLinescan_1D(:));
        stdLinescan_1D=stdLinescan_1D./max(stdLinescan_1D(:));
        clear data;
        
        % finding linescan Glavo position using cross-correlation
        if nFramesCartesian>=minFramesforSTD
            maxC=NaN(1,heightCartesian);
            for k=1:heightCartesian
                cartesianFrame_1D=cartesianFrame(k,:);
                cartesianFrame_1D=cartesianFrame_1D-min(cartesianFrame_1D(:));
                cartesianFrame_1D=cartesianFrame_1D./max(cartesianFrame_1D(:));
                C=normxcorr2(cartesianFrame_1D,stdLinescan_1D);
                maxC(k)=max(C(numel(cartesianFrame_1D)-1:end-...
                    numel(cartesianFrame_1D)+1)); % limiting search space to only lags where there is 100% overlap between target and template image
            end
            [~,linescanGalvo]=max(maxC,[],2);
        end
        
        display1=vertcat(cartesianFramePadded,linescanFrame);
        if nFramesCartesian>=minFramesforSTD
            display1=insertText(display1,[size(display1,2)/2 1],...
                ['Motion Contrast (nFrames=',num2str(nFramesCartesian-1),...
                ')'],'AnchorPoint','CenterTop','BoxOpacity',0,'FontSize',...
                fontsize1,'TextColor',[1 0 0]);
        else
            display1=insertText(display1,[size(display1,2)/2 1],...
                ['Average Intensity (nFrames=',num2str(max([1 nFramesCartesian-1])),...
                ')'],'AnchorPoint','CenterTop','BoxOpacity',0,'FontSize',...
                fontsize1,'TextColor',[1 0 0]);
        end
        display1=insertText(display1,[size(display1,2)/2 size(display1,1)-120],...
            'Click on linescan Galvo position',...
            'AnchorPoint','CenterTop','BoxOpacity',0,'FontSize',fontsize1,...
            'TextColor',[1 0 0]);
        display1=insertText(display1,[size(display1,2)/2 size(display1,1)-70],...
            'Ignore if correct',...
            'AnchorPoint','CenterTop','BoxOpacity',0,'FontSize',fontsize1,...
            'TextColor',[1 0 0]);
        if exist('linescanGalvo','var')
            display1=insertShape(display1,'Line',...
                [1 linescanGalvo size(cartesianFramePadded,2) linescanGalvo],...
                'LineWidth',2,'Color',[1 1 1],'Opacity',1);
        end
        fNo=fNo+1;
        figure(fNo);
        currentFigure=gcf;
        imshow(display1,'Border','tight','InitialMagnification',...
            displayMagnification*100);colormap(gray);axis off;hold on;
        currentOuterPosition=currentFigure.OuterPosition;
        currentFigure.OuterPosition(1)=70;
        currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(1,4)+1;
        currentOuterPosition=currentFigure.OuterPosition;
        [~,y]=ginput(1);
        if numel(y)~=0
            linescanGalvo=round(y);
        end
        
        gaussianSigma=7.5*(0.03542/.3542); % most conservative value: for 5 degree FOV
        
        cartesianFrameGaussian=imgaussfilt(cartesianFrame,gaussianSigma);
        cartesianFrameGaussian=cartesianFrameGaussian-min(cartesianFrameGaussian(:));
        cartesianFrameGaussian=cartesianFrameGaussian./max(cartesianFrameGaussian(:)).*1;
        filterSize=2*ceil(2*gaussianSigma)+1;
        
        width=size(cartesianFrameGaussian,2);
        height=size(cartesianFrameGaussian,1);
        
        fNo=fNo+1;
        figure('OuterPosition',[62 450 620 620]);
        imshow(cartesianFrame,'Border','tight','InitialMagnification',...
            displayMagnification*100);drawnow;
        currentFigure=gcf;
        currentFigure.OuterPosition(1)=62;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(fNo,4)-1;
        cartesianFrameVector=double(reshape(cartesianFrame.*255,[1 numel(cartesianFrame)]));
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        fNo=fNo+1;
        figure('OuterPosition',[currentOuterPosition(fNo-1,1) 1 currentOuterPosition(fNo-1,3) 400]);
        histogram(cartesianFrameVector,'BinWidth',1);axis tight;
        axF2=gca;axF2.YScale='linear';axF2.FontSize=18;axF2.FontWeight='bold';
        currentFigure=gcf;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        
        fNo=fNo+1;
        figure('OuterPosition',[currentOuterPosition(fNo-2,1)+...
            currentOuterPosition(fNo-2,3)+1 currentOuterPosition(fNo-2,2) ...
            currentOuterPosition(fNo-2,3) currentOuterPosition(fNo-2,4)]);
        imshow(cartesianFrameGaussian,'Border','tight','InitialMagnification',...
            displayMagnification*100);drawnow;
        currentFigure=gcf;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(fNo,4)-1;
        imgGaussian_1D=double(reshape(cartesianFrameGaussian.*255,...
            [1 numel(cartesianFrameGaussian)]));
        fNo=fNo+1;
        figure('OuterPosition',[currentOuterPosition(fNo-1,1) ...
            currentOuterPosition(fNo-2,2) currentOuterPosition(fNo-2,3) ...
            currentOuterPosition(fNo-2,4)]);
        histogram(imgGaussian_1D,'BinWidth',1);axis tight;axF4=gca;
        axF4.YScale='linear';axF4.FontSize=18;axF4.FontWeight='bold';
        currentFigure=gcf;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        if axF4.YLim(2)>axF2.YLim(2)
            axF2.YLim(2)=axF4.YLim(2);
        else
            axF4.YLim(2)=axF2.YLim(2);
        end
        
        [level,EM]=graythresh(cartesianFrameGaussian);
        line([level*255,level*255],[0 axF4.YLim(2)],'LineWidth',2,'Color',[1 0 0]);
        cartesianFrameThresh=im2bw(cartesianFrameGaussian,level);
        fNo=fNo+1;
        figure('OuterPosition',[currentOuterPosition(fNo-2,1)+...
            currentOuterPosition(fNo-2,3)+1 currentOuterPosition(fNo-2,2) ...
            currentOuterPosition(fNo-2,3) currentOuterPosition(fNo-2,4)]);
        imgGaussianText=double(cartesianFrameGaussian);
        offset=-24;
        imgGaussianText=insertText(imgGaussianText,[round(width/2),...
            offset+25],'Click on ROI center','AnchorPoint','CenterTop',...
            'TextColor','red','BoxColor','black','FontSize',fontsize1,'BoxOpacity',0);
        imgGaussianText=insertText(imgGaussianText,[round(width/2),offset+70],...
            'Along red line','AnchorPoint','CenterTop','TextColor','red',...
            'BoxColor','black','FontSize',fontsize1,'BoxOpacity',0);
        imshow(imgGaussianText,'Border','tight','InitialMagnification',...
            displayMagnification*100);hold on;drawnow;
        line([1,width],[linescanGalvo linescanGalvo],'LineWidth',6,'Color',[1 0 0]);
        currentFigure=gcf;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        currentFigure.OuterPosition=[currentOuterPosition(fNo-2,1)+...
            currentOuterPosition(fNo-2,3)+1 currentOuterPosition(fNo-2,2) ...
            currentOuterPosition(fNo-2,3) currentOuterPosition(fNo-2,4)];
        currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(fNo,4)-1;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;

        [x,y]=ginput(1);
        if not(isempty(x))
            rotateCenterX=round(x(1));
        end
        
        % % MANUAL OVERRIDE
        % rotateCenterX=244;
        
        rotateCenterY=linescanGalvo;
        
        %padding image with zeros for rotation about a specified point
        padSizePre1=max([0 2*round(height/2-rotateCenterY)]);
        padSizePre2=max([0 2*round(width/2-rotateCenterX)]);
        padSizePost1=max([0 2*round(rotateCenterY-height/2)]);
        padSizePost2=max([0 2*round(rotateCenterX-width/2)]);
        cartesianFrameThreshPadded=padarray(cartesianFrameThresh,...
            [padSizePre1 padSizePre2],'pre');
        cartesianFrameGaussianPadded=padarray(cartesianFrameGaussian,...
            [padSizePre1 padSizePre2],'pre');
        cartesianFrameThreshPadded=padarray(cartesianFrameThreshPadded,...
            [padSizePost1 padSizePost2],'post');
        cartesianFrameGaussianPadded=padarray(cartesianFrameGaussianPadded,...
            [padSizePost1 padSizePost2],'post');
        heightPadded=size(cartesianFrameThreshPadded,1);
        widthPadded=size(cartesianFrameThreshPadded,2);
        if heightPadded>widthPadded
            cartesianFrameThreshPadded=padarray(cartesianFrameThreshPadded,...
                [0 round((heightPadded-widthPadded)/2)],'both');
            cartesianFrameGaussianPadded=padarray(cartesianFrameGaussianPadded,...
                [0 round((heightPadded-widthPadded)/2)],'both');
        else
            cartesianFrameThreshPadded=padarray(cartesianFrameThreshPadded,...
                [round((widthPadded-heightPadded)/2) 0],'both');
            cartesianFrameGaussianPadded=padarray(cartesianFrameGaussianPadded,...
                [round((widthPadded-heightPadded)/2) 0],'both');
        end
        heightPadded=size(cartesianFrameThreshPadded,1);
        widthPadded=size(cartesianFrameThreshPadded,2);
        
        fNo=fNo+1;
        figure('OuterPosition',[currentOuterPosition(fNo-5,1) ...
            currentOuterPosition(fNo-5,2)-10 currentOuterPosition(fNo-5,3) ...
            currentOuterPosition(fNo-5,4)]);
        imgGaussianPaddedText=insertText(double(cartesianFrameGaussianPadded)...
            ,[round(widthPadded/2),1],'Click on ROI radius','AnchorPoint',...
            'CenterTop','TextColor','red','BoxColor','black','FontSize',...
            fontsize1,'BoxOpacity',0);
        imgGaussianPaddedText=insertText(imgGaussianPaddedText,...
            [round(widthPadded/2),40],'Along red line','AnchorPoint',...
            'CenterTop','TextColor','red','BoxColor','black','FontSize',...
            fontsize1,'BoxOpacity',0);
        imshow(imgGaussianPaddedText,'Border','tight','InitialMagnification',...
            displayMagnification*100);hold on;drawnow;
        currentFigure=gcf;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        currentFigure.OuterPosition(1)=currentOuterPosition(fNo-5,1);
        currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(fNo,4)-1;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        plot(round(widthPadded/2),round(heightPadded/2),'.','MarkerSize',...
            40,'Color',[1 0 0]);
        line([1,widthPadded],[round(heightPadded/2) round(heightPadded/2)],...
            'LineWidth',6,'Color',[1 0 0]);
        
        %drawing reference circles
        centerX=round(widthPadded/2);
        centerY=round(heightPadded/2);
        for radiusMarkers=50:50:300
            phi = 0:pi/360:2*pi;
            xunit=radiusMarkers.*cos(phi)+centerX;
            yunit=radiusMarkers.*sin(phi)+centerY;
            plot(xunit,yunit,'.','MarkerSize',10,'Color',[rand rand rand]);
        end
        plot(centerX,centerY,'.','MarkerSize',40,'Color',[1 0 0]);
        line([1,widthPadded],[centerY centerY],'LineWidth',6,'Color',[1 0 0]);
        [x,y]=ginput(1);
        radiusROI=round(sqrt((centerX-x)^2+(centerY-y)^2));
        
        % % MANUAL OVERRIDE
        % radiusROI=100;
        
        cartesianFrameThreshPaddedROI=cartesianFrameThreshPadded(centerY-...
            radiusROI:centerY+radiusROI,centerX-radiusROI:centerX+radiusROI);
        cartesianFrameGaussianPaddedROI=cartesianFrameGaussianPadded(centerY-...
            radiusROI:centerY+radiusROI,centerX-radiusROI:centerX+radiusROI);
        
        % Rotating vessel and making mask for radon
        heightROI=size(cartesianFrameThreshPaddedROI,1);
        widthROI=size(cartesianFrameThreshPaddedROI,2);
        ROIcenterX=round(widthROI/2);
        ROIcenterY=round(heightROI/2);
        mask=zeros(heightROI,widthROI);
        for k=1:heightROI
            for m=1:widthROI
                if sqrt((k-ROIcenterY)^2+(m-ROIcenterX)^2)<radiusROI
                    mask(k,m)=1;
                end
            end
        end
        mask=im2bw(mask,0.5);
        cartesianFrameThreshPaddedROI_masked=cartesianFrameThreshPaddedROI.*mask;
        cartesianFrameGaussianPaddedROI_masked=cartesianFrameGaussianPaddedROI.*mask;
        angles=0:0.1:180;
        r=radon(cartesianFrameGaussianPaddedROI_masked,angles);
        st=std(r,1);
        [~,loc]=max(st);
        vesselAngle=angles(loc);
        
        cartesianFrameThreshPaddedRotated=imrotate(cartesianFrameThreshPadded,...
            -1*vesselAngle,'crop');
        cartesianFrameGaussianPaddedRotated=imrotate(cartesianFrameGaussianPadded,...
            -1*vesselAngle,'crop');
        cartesianFrameThreshPaddedRotatedROI=cartesianFrameThreshPaddedRotated(...
            centerY-radiusROI:centerY+radiusROI,centerX-radiusROI:centerX+radiusROI);
        cartesianFrameGaussianPaddedRotatedROI=cartesianFrameGaussianPaddedRotated(...
            centerY-radiusROI:centerY+radiusROI,centerX-radiusROI:centerX+radiusROI);
        cartesianFrameGaussianPaddedRotatedROI_masked=...
            cartesianFrameGaussianPaddedRotatedROI.*mask;
        fNo=fNo+1;
        figure('OuterPosition',[currentOuterPosition(fNo-1,1)+...
            currentOuterPosition(fNo-1,3)+1 200 200 200]);
        imshow(cartesianFrameGaussianPaddedRotatedROI_masked,'Border',...
            'tight','InitialMagnification',displayMagnification*100);hold on;drawnow;
        currentFigure=gcf;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        currentFigure.OuterPosition(1)=currentOuterPosition(fNo-1,1)+...
            currentOuterPosition(fNo-1,3)+1;
        currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(fNo,4)-1;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        fNo=fNo+1;
        figure('OuterPosition',[currentOuterPosition(fNo-1,1)+...
            currentOuterPosition(fNo-1,3)+1 currentOuterPosition(fNo-1,2) ...
            currentOuterPosition(fNo-1,3) currentOuterPosition(fNo-1,4)]);
        cartesianFrameText=insertText(double(cartesianFrameGaussian),...
            [round(widthCartesian/2),1],['Vessel angle = ',num2str(...
            abs(90-vesselAngle)),' degrees'],'AnchorPoint','CenterTop',...
            'TextColor','red','BoxColor','black','FontSize',fontsize1,'BoxOpacity',0);
        cartesianFrameText=insertText(cartesianFrameText,[round(widthCartesian/2),...
            40],'w.r.t. resonant scanner','AnchorPoint','CenterTop','TextColor',...
            'red','BoxColor','black','FontSize',fontsize1,'BoxOpacity',0);
        imshow(cartesianFrameText,'Border','tight','InitialMagnification',...
            displayMagnification*100);hold on;drawnow;
        currentFigure=gcf;
        currentFigure.Name=cartesianVideoFile(1:27);
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        currentFigure.OuterPosition(1)=currentOuterPosition(fNo-1,1)+...
            currentOuterPosition(fNo-1,3)+1;
        currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(fNo,4)-1;
        currentOuterPosition(fNo,:)=currentFigure.OuterPosition;
        mSlope=tand(vesselAngle+90);
        lineXYCoordinates=[rotateCenterX-radiusROI/sqrt(1+mSlope^2) ...
            rotateCenterY+mSlope*radiusROI/sqrt(1+mSlope^2) rotateCenterX+...
            radiusROI/sqrt(1+mSlope^2) rotateCenterY-mSlope*radiusROI/...
            sqrt(1+mSlope^2)];
        line([lineXYCoordinates(1) lineXYCoordinates(3)],[lineXYCoordinates(2) ...
            lineXYCoordinates(4)],'LineWidth',2,'Color',[0 0 1]);
        phi = 0:pi/360:2*pi;
        xunit=radiusROI.*cos(phi)+rotateCenterX;
        yunit=radiusROI.*sin(phi)+rotateCenterY;
        plot(xunit,yunit,'.','MarkerSize',10,'Color',[1 0 0]);
        vesselAngle=abs(90-vesselAngle);
        set(currentFigure,'PaperPositionMode','auto');
        figure(currentFigure);
        pause(2);
        input('Press Enter to continue... ');
        print(fNo,[Step5folder,fileName(1:25),'_Step3_VesselOrientation',...
            appendToResults],'-dtiff');
        
        vesselAngle
        save([Step3folder,fileName(1:25),'_Step3',appendToResults]);
        fNo=0;
    else
        angleManuallyEntered=1;
        save([Step3folder,fileName(1:25),'_Step3',appendToResults]);
    end
    
    % moving processed file to separate folder
    if moveFileAfterProcessing==1
        movefile([processDir,fileName],[processDir,'PROCESSED\',fileName]);
    end

    % Saving copy of code (added 07-31-2018, 2.01 pm)
    codeFileFullPath=mfilename('fullpath'); %does not store .m extension
    [codeFilePath,codeFileName,~]=fileparts(codeFileFullPath);
    newCodeFileName=['code_',fileName(1:25),'_Step3',appendToResults];
    codeSource=[codeFilePath,'\',codeFileName,'.m'];
    codeDestination=[Step5folder,newCodeFileName,'.m'];
    copyfile(codeSource,codeDestination);

end
beep;
