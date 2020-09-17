function [angle,distance,octile,xp,yp]=CircularFitMapping(XY,figflag,ang_shift,ratio_flag)

if ratio_flag
    ratio_XY = (prctile(XY(:,1),97.5)-prctile(XY(:,1),2.5)) / (prctile(XY(:,2),97.5)-prctile(XY(:,2),2.5));
    XY_circled = [XY(:,1) XY(:,2)*ratio_XY];
else
    ratio_XY = 1;
    XY_circled = XY;
end

[x0,y0,rad] = CircleFitByPratt(XY_circled);
ang=-linspace(0,2*pi,100000)';
if nargin < 3
    ang_shift = pi/2;
end
xp=rad*cos(ang+ang_shift)+x0; 
yp=rad*sin(ang+ang_shift)+y0;

yp = yp/ratio_XY;
y0 = y0/ratio_XY;

%%% find the closest 
[idx,D]=knnsearch(1000*[xp yp],1000*XY,'K',2);

angle = (ang(idx(:,1))+ang(idx(:,2)))/2;
distance = D(:,1);

% sort by angle
[~,idx_sort] = sort(angle);
octile = ceil(double(idx_sort*8)/length(idx_sort));

if figflag
    figure
    scatter(XY(:,1)-x0,XY(:,2)-y0,10,'filled','MarkerFaceAlpha',0.2)
    hold on
    plot(xp-x0,yp-y0,'.')
    scatter(xp(1)-x0,yp(1)-y0,40,'filled')
end

