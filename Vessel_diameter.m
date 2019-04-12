% Calculate Vessel Diameter
% Author: Aby Joseph, University of Rochester
% License: GPL-3.0-or-later
% Last modified: 01-24-2019

clear;close all;clc;
screenSize=get(0,'ScreenSize');

f_No=0;hist_bin_size=255;
startDir='\\cvsnas3.urmc-sh.rochester.edu\aria\Mouse\Aby_Data\TRBF_data\Experimental DATA\'; %data directory
[fileName,filePath,~]=uigetfile([startDir,'*.tif*'],'Select fluorescence TIFF image');
cd(filePath);

ker_size=150; %has to be even number
ker_size_b=50; %has to be even number

FOV=input('Enter FOV: ');
originalVideoWidth=input('Enter original video width in FOV direction: ');
micronsPerDegree=34;
mic_pix=micronsPerDegree*FOV./originalVideoWidth;

imgOriginal=double(imread(fileName)); %fluorescense TIFF image
imgOriginal=imgOriginal-min(imgOriginal(:));
imgOriginal=imgOriginal./max(imgOriginal(:)).*1;

gaussianSigma=7.5*(0.5/FOV);

imgGaussian=imgaussfilt(imgOriginal,gaussianSigma);
imgGaussian=imgGaussian-min(imgGaussian(:));
imgGaussian=imgGaussian./max(imgGaussian(:)).*1;
filterSize=2*ceil(2*gaussianSigma)+1;

width=size(imgGaussian,2);
height=size(imgGaussian,1);

fontsize1=18;
f_No=f_No+1;
figure('OuterPosition',[62 450 620 620]);
imshow(imgOriginal,'Border','tight','InitialMagnification',100);drawnow;
currentFigure=gcf;
currentFigure.OuterPosition(1)=62;
currentOuterPosition=currentFigure.OuterPosition;
currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(f_No,4)-1;
imgOriginal_1D=double(reshape(imgOriginal.*255,[1 numel(imgOriginal)]));
f_No=f_No+1;
figure('OuterPosition',[currentOuterPosition(1) 1 currentOuterPosition(3) 400]);
histogram(imgOriginal_1D,'BinWidth',1);axis tight;axF2=gca;axF2.YScale='linear';axF2.FontSize=fontsize1;axF2.FontWeight='bold';
currentFigure=gcf;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;

f_No=f_No+1;
figure('OuterPosition',[currentOuterPosition(f_No-2,1)+...
    currentOuterPosition(f_No-2,3)+1 currentOuterPosition(f_No-2,2) ...
    currentOuterPosition(f_No-2,3) currentOuterPosition(f_No-2,4)]);
imshow(imgGaussian,'Border','tight','InitialMagnification',100);drawnow;
currentFigure=gcf;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(f_No,4)-1;
imgGaussian_1D=double(reshape(imgGaussian.*255,[1 numel(imgGaussian)]));
f_No=f_No+1;
figure('OuterPosition',[currentOuterPosition(f_No-1,1) ...
    currentOuterPosition(f_No-2,2) currentOuterPosition(f_No-2,3) ...
    currentOuterPosition(f_No-2,4)]);
histogram(imgGaussian_1D,'BinWidth',1);axis tight;axF4=gca;axF4.YScale=...
    'linear';axF4.FontSize=fontsize1;axF4.FontWeight='bold';
currentFigure=gcf;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
if axF4.YLim(2)>axF2.YLim(2)
    axF2.YLim(2)=axF4.YLim(2);
else
    axF2.YLim(2)=axF4.YLim(2);
end

[level, EM]=graythresh(imgGaussian);
imgThresh=im2bw(imgGaussian,level);
fontsize2=22;
f_No=f_No+1;
figure('OuterPosition',[currentOuterPosition(f_No-2,1)+...
    currentOuterPosition(f_No-2,3)+1 currentOuterPosition(f_No-2,2) ...
    currentOuterPosition(f_No-2,3) currentOuterPosition(f_No-2,4)]);
imgGaussianText=insertText(double(imgGaussian),[round(width/2),1],...
    ['Otsu effectiveness metric = ',num2str(EM,3)],'AnchorPoint',...
    'CenterTop','TextColor','red','BoxColor','black','FontSize',...
    fontsize2,'BoxOpacity',0);
imshow(imgGaussianText,'Border','tight','InitialMagnification',100);
hold on;drawnow;
currentFigure=gcf;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
currentFigure.OuterPosition=[currentOuterPosition(f_No-2,1)+...
    currentOuterPosition(f_No-2,3)+1 currentOuterPosition(f_No-2,2) currentOuterPosition(f_No-2,3) currentOuterPosition(f_No-2,4)];
currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(f_No,4)-1;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;

input('Press enter to continue..');
linescanL=input('Enter Line Number: ');

% MANUAL OVERRIDE
% linescanL=175;

if not(isempty(linescanL))
    linescanL=round((linescanL/480)*height);
    offset=0;
else
    imgGaussianText=insertText(imgGaussian,[round(width/2),25],...
        'Click on position of linescan','AnchorPoint','CenterTop',...
        'TextColor','red','BoxColor','black','FontSize',fontsize2,...
        'BoxOpacity',0);
    offset=25;
    imshow(imgGaussianText,'Border','tight','InitialMagnification',100);
    hold on;drawnow;
    [x,y]=ginput(1);
    if isempty(y)
        linescanL=round(height/2);
    else
        linescanL=round(y);
    end
end
linescanL;

figure(f_No);
imgGaussianText=insertText(imgGaussian,[round(width/2),offset+25],...
    'Click on center for vessel rotation','AnchorPoint','CenterTop',...
    'TextColor','red','BoxColor','black','FontSize',fontsize2,...
    'BoxOpacity',0);
imshow(imgGaussianText,'Border','tight','InitialMagnification',100);
hold on;drawnow;
line([1,width],[linescanL linescanL],'LineWidth',6,'Color',[1 0 0]);

[x,y]=ginput(1);
if not(isempty(x))
    rotateCenterX=round(x(1));
end

% MANUAL OVERRIDE
% rotateCenterX=244;

rotateCenterY=linescanL;

%padding image with zeros for rotation about a specified point
padSizePre1=max([0 2*round(height/2-rotateCenterY)]);
padSizePre2=max([0 2*round(width/2-rotateCenterX)]);
padSizePost1=max([0 2*round(rotateCenterY-height/2)]);
padSizePost2=max([0 2*round(rotateCenterX-width/2)]);
imgThreshPadded=padarray(imgThresh,[padSizePre1 padSizePre2],'pre');
imgGaussianPadded=padarray(imgGaussian,[padSizePre1 padSizePre2],'pre');
imgThreshPadded=padarray(imgThreshPadded,[padSizePost1 padSizePost2],'post');
imgGaussianPadded=padarray(imgGaussianPadded,[padSizePost1 padSizePost2],'post');
heightPadded=size(imgThreshPadded,1);
widthPadded=size(imgThreshPadded,2);
if heightPadded>widthPadded
    imgThreshPadded=padarray(imgThreshPadded,[0 round((heightPadded-widthPadded)/2)],'both');
    imgGaussianPadded=padarray(imgGaussianPadded,[0 round((heightPadded-widthPadded)/2)],'both');
else
    imgThreshPadded=padarray(imgThreshPadded,[round((widthPadded-heightPadded)/2) 0],'both');    
    imgGaussianPadded=padarray(imgGaussianPadded,[round((widthPadded-heightPadded)/2) 0],'both');    
end
heightPadded=size(imgThreshPadded,1);
widthPadded=size(imgThreshPadded,2);

f_No=f_No+1;
figure('OuterPosition',[currentOuterPosition(f_No-5,1) ...
    currentOuterPosition(f_No-5,2)-10 currentOuterPosition(f_No-5,3) ...
    currentOuterPosition(f_No-5,4)]);
imshow(imgGaussianPadded,'Border','tight','InitialMagnification',100);
hold on;drawnow;
currentFigure=gcf;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
currentFigure.OuterPosition(1)=currentOuterPosition(f_No-5,1);
currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(f_No,4)-1;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
plot(round(widthPadded/2),round(heightPadded/2),'.','MarkerSize',40,...
    'Color',[1 0 0]);
line([1,widthPadded],[round(heightPadded/2) round(heightPadded/2)],...
    'LineWidth',6,'Color',[1 0 0]);

%drawing reference circles
imgGaussianPaddedText=insertText(double(imgGaussianPadded),[round(...
    widthPadded/2),1],'Click on maximum distance from center for diameter measurement',...
    'AnchorPoint','CenterTop','TextColor','red','BoxColor','black',...
    'FontSize',fontsize1,'BoxOpacity',0);
imshow(imgGaussianPaddedText,'Border','tight','InitialMagnification',100);
drawnow;
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

% MANUAL OVERRIDE
% radiusROI=100;

imgThreshPaddedROI=imgThreshPadded(centerY-radiusROI:centerY+radiusROI,...
    centerX-radiusROI:centerX+radiusROI);
imgGaussianPaddedROI=imgGaussianPadded(centerY-radiusROI:centerY+...
    radiusROI,centerX-radiusROI:centerX+radiusROI);

%rotating vessel
%making mask for radon
heightROI=size(imgThreshPaddedROI,1);
widthROI=size(imgThreshPaddedROI,2);
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
imgThreshPaddedROI_masked=imgThreshPaddedROI.*mask;
imgGaussianPaddedROI_masked=imgGaussianPaddedROI.*mask;
angles=0:0.1:180;
r=radon(imgThreshPaddedROI_masked,angles);
st=std(r,1);
[~,loc]=max(st);
vesselAngle=angles(loc);

imgThreshPaddedRotated=imrotate(imgThreshPadded,-1*vesselAngle,'crop');
imgGaussianPaddedRotated=imrotate(imgGaussianPadded,-1*vesselAngle,'crop');
imgThreshPaddedRotatedROI=imgThreshPaddedRotated(centerY-radiusROI:...
    centerY+radiusROI,centerX-radiusROI:centerX+radiusROI);
imgGaussianPaddedRotatedROI=imgGaussianPaddedRotated(centerY-radiusROI:...
    centerY+radiusROI,centerX-radiusROI:centerX+radiusROI);
f_No=f_No+1;
figure('OuterPosition',[currentOuterPosition(f_No-1,1)+...
    currentOuterPosition(f_No-1,3)+1 200 200 200]);
imshow(imgGaussianPaddedRotatedROI,'Border','tight','InitialMagnification',...
    100);hold on;drawnow;
currentFigure=gcf;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
currentFigure.OuterPosition(1)=currentOuterPosition(f_No-1,1)+...
    currentOuterPosition(f_No-1,3)+1;
currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(f_No,4)-1;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
f_No=f_No+1;
figure('OuterPosition',[currentOuterPosition(f_No-1,1)+...
    currentOuterPosition(f_No-1,3)+1 currentOuterPosition(f_No-1,2) currentOuterPosition(f_No-1,3) currentOuterPosition(f_No-1,4)]);
imshow(imgGaussianPaddedRotated,'Border','tight','InitialMagnification',...
    100);hold on;drawnow;
currentFigure=gcf;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
currentFigure.OuterPosition(1)=currentOuterPosition(f_No-1,1)+...
    currentOuterPosition(f_No-1,3)+1;
currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(f_No,4)-1;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
plot(centerX,centerY,'.','MarkerSize',40,'Color',[1 0 0]);
line([centerX-radiusROI centerX+radiusROI],[centerY-radiusROI ...
    centerY-radiusROI],'LineWidth',2,'Color',[1 0 0]);
line([centerX-radiusROI centerX+radiusROI],[centerY+radiusROI ...
    centerY+radiusROI],'LineWidth',2,'Color',[1 0 0]);
line([centerX-radiusROI centerX-radiusROI],[centerY-radiusROI ...
    centerY+radiusROI],'LineWidth',2,'Color',[1 0 0]);
line([centerX+radiusROI centerX+radiusROI],[centerY-radiusROI ...
    centerY+radiusROI],'LineWidth',2,'Color',[1 0 0]);

display('To furhter restrict analysis extent along length of vessel, click on desired top and bottom points')
[x,y]=ginput(2);
if isempty(y)==1
    analysisY(1)=centerY-radiusROI;
    analysisY(2)=centerY+radiusROI;
else
    analysisY(1)=y(1);
    analysisY(2)=y(2);
end

%determining strips
stripSizeMicrons=5; %strip size in microns
stripSize=round(stripSizeMicrons./mic_pix);
if mod(stripSize,2)==0
    stripSize=stripSize+1; %stripSize always an odd number
end
nStrips=2*(floor(((abs(analysisY(2)-analysisY(1))/2)-(stripSize-1)/2)/...
    stripSize))+1; %nStrips always an odd number
stripYCoors=zeros(1,nStrips+1);
widthPadded=size(imgThreshPadded,2);
for k=1:nStrips+1
    stripYCoors(k)=round((centerY-nStrips/2*stripSize)+(k-1)*stripSize);
end

stripDiameters=zeros(1,nStrips);
tlead=zeros(1,nStrips);
ttrail=zeros(1,nStrips);
intensityProfile=zeros(nStrips,2*radiusROI+1);
for k=1:nStrips
    intensityProfile(k,:)=mean(imgGaussianPaddedRotated(stripYCoors(k):...
        stripYCoors(k+1),centerX-radiusROI:centerX+radiusROI),1);
    currentFolder=cd('E:\PhD\TRBF\CODES\MOUSE');
    [stripDiameters(k),tlead(k),ttrail(k)]=fwhm(1:2*radiusROI+1,...
        intensityProfile(k,:));
    cd(currentFolder);
end

f_No=f_No+1;
figure('OuterPosition',[currentOuterPosition(f_No-1,1) ...
    currentOuterPosition(f_No-1,2) currentOuterPosition(f_No-1,3) ...
    currentOuterPosition(f_No-1,4)]);
imshow(imgGaussianPaddedRotated,'Border','tight','InitialMagnification',...
    100);hold on;drawnow;
currentFigure=gcf;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
currentFigure.OuterPosition(1)=currentOuterPosition(f_No-1,1);
currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(f_No,4)-1;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
plot(centerX,centerY,'.','MarkerSize',40,'Color',[1 0 0]);
line([centerX-radiusROI centerX+radiusROI],[centerY-radiusROI ...
    centerY-radiusROI],'LineWidth',2,'Color',[1 0 0]);
line([centerX-radiusROI centerX+radiusROI],[centerY+radiusROI ...
    centerY+radiusROI],'LineWidth',2,'Color',[1 0 0]);
line([centerX-radiusROI centerX-radiusROI],[centerY-radiusROI ...
    centerY+radiusROI],'LineWidth',2,'Color',[1 0 0]);
line([centerX+radiusROI centerX+radiusROI],[centerY-radiusROI ...
    centerY+radiusROI],'LineWidth',2,'Color',[1 0 0]);drawnow;
for k=1:nStrips+1
    line([centerX-radiusROI+2,centerX+radiusROI-2],[stripYCoors(k),...
        stripYCoors(k)],'LineWidth',2,'Color','cyan');
end

f_No=f_No+1;
figure('OuterPosition',[currentOuterPosition(f_No-1,1) ...
    currentOuterPosition(f_No-1,2)-currentOuterPosition(f_No-1,4)-1 currentOuterPosition(f_No-1,3) currentOuterPosition(f_No-1,4)]);
for k=1:nStrips
    plot(intensityProfile(k,:),'LineWidth',2,'Color',[rand rand rand]);
    hold on;
end
axis tight;ax=gca;ax.YLim=[0 1];
meanDiameter=mean(stripDiameters(:));
meanDiameterMicrons=meanDiameter*mic_pix;
stdDiameter=std(stripDiameters(:),[],1);
stdDiameterMicrons=stdDiameter*mic_pix;
tlead=tlead+centerX-radiusROI-1;
ttrail=ttrail+centerX-radiusROI-1;
meanTlead=mean(tlead(:));
meanTtrail=mean(ttrail(:));

input('Press enter to continue...');

f_No=f_No+1;
figure('OuterPosition',[currentOuterPosition(f_No-2,1)+...
    currentOuterPosition(f_No-2,3)+1 currentOuterPosition(f_No-2,2) ...
    currentOuterPosition(f_No-2,3) currentOuterPosition(f_No-2,4)]);
imshow(imgGaussianPaddedRotated,'Border','tight','InitialMagnification',...
    100);hold on;
imgGaussianPaddedRotatedText=insertText(imgGaussianPaddedRotated,...
    [round(widthPadded/2),1],['Diameter = ',num2str(meanDiameterMicrons,...
    '%.1f'),' +- ',num2str(stdDiameterMicrons,'%.1f'),' um'],...
    'AnchorPoint','CenterTop','TextColor','red','BoxColor','black',...
    'FontSize',fontsize2,'BoxOpacity',0);
imshow(imgGaussianPaddedRotatedText,'Border','tight',...
    'InitialMagnification',100);hold on;drawnow;
currentFigure=gcf;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
currentFigure.OuterPosition(1)=currentOuterPosition(f_No-2,1)+...
    currentOuterPosition(f_No-2,3)+1;
currentFigure.OuterPosition(2)=screenSize(4)-currentOuterPosition(f_No,4)-1;
currentOuterPosition(f_No,:)=currentFigure.OuterPosition;
line([meanTlead meanTlead],[analysisY(1) analysisY(2)],'LineWidth',1,...
    'Color',[1 0 0]);
line([meanTtrail meanTtrail],[analysisY(1) analysisY(2)],'LineWidth',1,...
    'Color',[1 0 0]);drawnow;

figure(7);
line([meanTlead-(centerX-radiusROI-1) meanTlead-(centerX-radiusROI-1)],...
    [analysisY(1)-(centerY-radiusROI-1) analysisY(2)-...
    (centerY-radiusROI-1)],'LineWidth',1,'Color',[1 0 0]);
line([meanTtrail-(centerX-radiusROI-1) meanTtrail-(centerX-radiusROI-1)],...
    [analysisY(1)-(centerY-radiusROI-1) analysisY(2)-(centerY-...
    radiusROI-1)],'LineWidth',1,'Color',[1 0 0]);drawnow;

drawnow;
