% Calculate Full Width at Half Maximum of 1D Intensity Profile
% Author: Aby Joseph, University of Rochester
% License: GPL-3.0-or-later
% Last modified: 01-24-2019

function [width,tlead,ttrail] = fwhm(x,y)

y=y-min(y(:));
y=y./max(y(:));
N=numel(y);

background=[y(1:3),y(end-2:end)];
meanBackground=mean(background(:));
stdBackground=std(background(:),0,1);
levSignal=meanBackground+3*stdBackground;

lev50=0.5;

[~,centerindex]=max(y);

check1=find(y(1:centerindex)<0.5);
check2=find(y(centerindex:end)<0.5);
if numel(check1)*numel(check2)==0
    fwhmGO=0;
else
    fwhmGO=1;    
end

if lev50<levSignal
    fwhmGO=0;
end

if fwhmGO==1
    i=centerindex; %start search for leading edge from center to left edge
    while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) && (i>1))
        i=i-1;
    end
    
    if i ~= 1
        
        interp = (lev50-y(i-1)) / (y(i)-y(i-1));
        tlead = x(i-1) + interp*(x(i)-x(i-1));
        
        i=centerindex; %start search for trailing edge from center to right edge
        while ((sign(y(i)-lev50) == sign(y(i+1)-lev50)) && (i<N))
            i = i+1;
        end
        
        % i = centerindex+1;                    %start search for next crossing at center
        % while ((sign(y(i)-lev50) == sign(y(i-1)-lev50)) && (i <= N-1))
        %     i = i+1;
        % end
        
        if i ~= N
            interp = (y(i)-lev50) / (y(i)-y(i+1));
            ttrail = x(i) + interp*(x(i+1)-x(i));
            width = ttrail - tlead;
        else
            ttrail = NaN;
            width = NaN;
        end
        
    else
        tlead = NaN;
        width = NaN;
    end
else
    tlead = NaN;
    ttrail = NaN;
    width = NaN;
end

end
