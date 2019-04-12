% Velocity_Step1_MakeStripAndDesinusoid
% Author: Aby Joseph, University of Rochester
% License: GPL-3.0-or-later
% Last modified: 01-24-2019

clear;
close all;

postFPGA=1; % "0" if data collected on Matrox, "1" if collected on FPGA
desinusoidMethod=2; % Method #1 = rescaling data to min and max pixel
% values after desinusoiding. Method #2 = Recommended. No rescaling. Direct
% conversion to uint8 after desinusoiding.
writeOriginalImage=1;

micronsPerDegree=34; % micron to degree FOV conversion for mouse. Schumaker and Schaeffel 2004.

startDir='\\cvsnas3.urmc-sh.rochester.edu\aria\Mouse\Aby_Data\TRBF_data\Experimental DATA\';
saveDir='E:\PhD\TRBF\Experimental DATA\RAW MAT files\';
Step5folder='E:\PhD\TRBF\Experimental DATA\Step 5 MAT files\';
linescanDir=uigetdir(startDir,'Select directory with raw linescan AVIs');
linescanDir=[linescanDir,'\'];

% Append current date and time to filenames of results
currentDateAndTime=clock;
year=num2str(currentDateAndTime(1));
month=num2str(currentDateAndTime(2),'%02d');
day=num2str(currentDateAndTime(3),'%02d');
hour=num2str(currentDateAndTime(4),'%02d');
minute=num2str(currentDateAndTime(5),'%02d');
appendToResults=['_',year,month,day,hour,minute];

namecontains='*.avi';
[dirFilenames(1).name,filePath,~]=uigetfile([linescanDir,namecontains],...
    'Select file to process - click cancel to process all files in folder');
if filePath==0
    filePath=[linescanDir,'\'];    
    currentFolder=cd(linescanDir);
    dirFilenames=dir(namecontains);
    cd(currentFolder);
    nFiles=numel(dirFilenames);
else
    nFiles=1;
end
tic

heightLinescanVideo=608;
for fileNumber=1:nFiles
    fileName=dirFilenames(fileNumber).name
    mov=VideoReader([filePath,fileName]);
    height=mov.Height;
    width=mov.Width
    if abs(postFPGA-0)<.01 %MATROX
        fprintf('\nMatrox mode\n\n');
        if abs(height-heightLinescanVideo)<0.01
            display(fileName);
            freq=15450; %Hz freq of resonant scanner pre July 2017
            %grids
            A=exist('FOVs','var'); % checks if this is the first linescan video found
            if A==0
                [gridsFile,gridsPath,~]=uigetfile('\\cvsnas3.urmc-sh.rochester.edu\aria\Mouse\Aby_Data\TRBF_data\Experimental DATA\*.mat',['Select desinusoid MAT file for ',fileName]);
                load([gridsPath,gridsFile],'vertical_fringes_desinusoid_matrix');
                mic_pix=input('Enter microns per pixel: ');
                FOVs=struct('mic_pix',mic_pix,'gridsFileName',gridsFile);
            else
                [gridsFile,gridsPath,~]=uigetfile([gridsPath,'*.mat'],['Select desinusoid MAT file for ',fileName]);
                load([gridsPath,gridsFile],'vertical_fringes_desinusoid_matrix');
                nFOVs=numel(FOVs);
                m=1;foundFlag=0;
                while m<=nFOVs && foundFlag==0
                    if strcmp(FOVs(m).gridsFileName,gridsFile)
                        foundFlag=1;
                        display(['Microns per pixel = ',num2str(FOVs(m).mic_pix)]);
                        mic_pix=input('Blank if OK, else enter new number: ');
                        if numel(mic_pix)==0
                            mic_pix=FOVs(m).mic_pix;
                        end
                    end
                    m=m+1;
                end
                if foundFlag==0
                    mic_pix=input('Enter microns per pixel: ');
                    FOVs(nFOVs+1).mic_pix=mic_pix;
                    FOVs(end).gridsFileName=gridsFile;
                end
            end
            
            frame_fs=struct('cdata',zeros(height,width,3,'uint8'),'colormap',[]);
            k=1;
            while hasFrame(mov)
                frame_fs(k).cdata=readFrame(mov);
                k=k+1;
            end
            nFrames=k-1;
            height=height-1; % to take care of first row of black pixels in every linescan frame
            data=uint8(zeros(nFrames*height,width));
            for k=1:nFrames
                tempFrame=frame_fs(k).cdata;
                data((k-1)*height+1:k*height,:)=tempFrame(2:end,:); % deleting first row of every linescan frame
            end
            
            % Desinusoiding
            if desinusoidMethod==1
                data=single(data);
                data=data*vertical_fringes_desinusoid_matrix';
                data=data-min(data(:));
                data=data./max(data(:));
                data=single(data); % data turned to "double" format in matrix multiplication 3 lines above
            else
                if desinusoidMethod==2
                    data=single(data);
                    data=uint8(data*vertical_fringes_desinusoid_matrix');
                end
            end
            data=rot90(data);
            
            width=height; % 'width' is time dimension from now on
            height=size(data,1); % 'height' is space dimension from now on
            FOV=mic_pix*height/micronsPerDegree;
            save([saveDir,fileName(1:25),'_RAW'],'postFPGA','FOVs',...
                'vertical_fringes_desinusoid_matrix','mic_pix',...
                'linescanDir','gridsFile','gridsPath','fileName',...
                'nFrames','data','height','width','desinusoidMethod',...
                'freq','-v7.3');
            if desinusoidMethod==1
                data=data.*255;
                data=uint8(data);
            end
            dataColor=cat(3,data,data,data);
            if writeOriginalImage==1
                imwrite(dataColor,[Step5folder,fileName(1:25),'_Step1_Original_Static_Lines',appendToResults,'.tif'],'TIFF');
            end
        end
    end
    
    if abs(postFPGA-1)<.01 %FPGA
        % for data collected on FPGA with live desinusoid
        fprintf('\nFPGA mode\n\n');
        FOV=input('Enter full FOV (degrees) in fast scan direction: ');
        mic_pix=micronsPerDegree*FOV./width
%         freq=input('Enter freq: ');
        freq=15063;
        frame_fs=struct('cdata',zeros(height,width,3,'uint8'),'colormap',[]);
        k=1;
        while hasFrame(mov)
            frame_fs(k).cdata=readFrame(mov);
            k=k+1;
        end
        nFrames=k-1;
        data=uint8(zeros(nFrames*height,width));
        for k=1:nFrames
            data((k-1)*height+1:k*height,:)=frame_fs(k).cdata;
        end
        data=rot90(data);
        
        width=height; % 'width' is time dimension from now on
        height=size(data,1); % 'height' is space dimension from now on
        save([saveDir,fileName(1:25),'_RAW',appendToResults],'postFPGA',...
            'linescanDir','fileName','nFrames','data','height','width','FOV',...
            'mic_pix','freq','-v7.3');
        dataColor=cat(3,data,data,data);
        if writeOriginalImage==1
            imwrite(dataColor,[Step5folder,fileName(1:25),'_Step1_Original_Static_Lines',appendToResults,'.tif'],'TIFF');
        end
    end
    
    % Saving copy of code (added 07-31-2018, 12.27 pm)
    codeFileFullPath=mfilename('fullpath'); %does not store .m extension
    [codeFilePath,codeFileName,~]=fileparts(codeFileFullPath);
    newCodeFileName=['code_',fileName(1:25),'_Step1',appendToResults];
    codeSource=[codeFilePath,'\',codeFileName,'.m'];
    codeDestination=[Step5folder,newCodeFileName,'.m'];
    copyfile(codeSource,codeDestination);    
end
elapsedTime=toc;disp(['Elapsed time is ',num2str(elapsedTime/60),' minutes.']);
beep
