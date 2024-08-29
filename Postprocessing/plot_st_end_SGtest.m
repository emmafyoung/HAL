% plot_st_end_SGtest.m to plot particle start and end points from HAL
% Set up for SGtest simulations forced with SG NEMO

% Set plot options (0 is off, 1 is on)
ipr=0; % Print or not

nplank=90; % Length of time from particle release to analysis, in days

noutpd=24; % Number of outputs per day in HAL output file 
           % (set in model_spec.dat)

% Read in NEMO info to overlay model bathymetry on plot

fname1='/Inputfiles/SG_mesh_hgr.nc';
fname2='/Inputfiles/SG_mesh_zgr.nc';
fname3='/Inputfiles/SG_mask.nc';
x=ncread(fname1,'glamt');
y=ncread(fname1,'gphit');
mbathy=ncread(fname2,'mbathy');
e3t=ncread(fname2,'e3t_0');
e3w=ncread(fname2,'e3w_0');
tmask=ncread(fname3,'tmaskutil');

nx=size(tmask,1);
ny=size(tmask,2);
nz=size(e3w,3);
land=find(tmask == 0);

% Find corners for pcolor
% Regular grid in both x and y
dx=x(2,1)-x(1,1);
x2=x-dx/2;
dy=y(1,2)-y(1,1);
y2=y-dy/2;

% Get range for particle exclusion (edge of model grid)
maxxp=max(max(x))-dx;
minxp=min(min(x))+dx;
maxyp=max(max(y))-dy;
minyp=min(min(y))+dy;

% Extract bathymetry

beddep=zeros(size(mbathy));
for i=1:nx;
 for j=1:ny;
  if (mbathy(i,j) == 0);
   beddep(i,j)=nan;
  else
   beddep(i,j)=sum(e3t(i,j,1:mbathy(i,j)),3);
  end;
 end;
end;

% Read in particle info.

file='/Outputfiles/IBMout.nc';
tim=ncread(file,'time');
lon=ncread(file,'longitude');
lat=ncread(file,'latitude');
flag=ncread(file,'switch');
npts=size(lon,1);

% End trajectories at analysis point

for i=1:npts
% Find where particle is switched on (flag to 1)
 parton=find(flag(i,:) == 1);
 firstt=min(parton);
 bg=firstt+(nplank*noutpd)-1;
 lon2(i)=lon(i,bg);
 lat2(i)=lat(i,bg);
% Set those at boundary to NaN
 if (lon2(i) > maxxp | lon2(i) < minxp | lat2(i) > maxyp | lat2(i) < minyp);
  lon2(i)=nan;
  lat2(i)=nan;
 end; 
end

% Plot start and end pts

figure(1)
clf
m_proj('Mercator','lat',[-58 -48],'lon',[-49 -30]);
hold on
% Start position
for i=1:npts;
  h=m_plot(lon(i,1),lat(i,1),'ro');
  set(h,'MarkerFaceColor',['r'],'MarkerSize',4.0);
end;
% Overlay end position of points
for i=1:npts;
  h=m_plot(lon2(i),lat2(i),'bo');
  set(h,'MarkerFaceColor',['b'],'MarkerSize',4.0);
end;
% Overlay bathymetry
v=[250 500 1000 2500 5000 7500];
[c,h]=m_contour(x,y,beddep,v,'k');
hold off
m_grid;
title('Predicted particle start and end locations - SGtest');
m_usercoast('/Inputfiles/SScoast5.mat','patch',[.6 .6 .6]);
hold off;
drawnow;

if (ipr == 1);
 orient landscape;
 print -dtiff st_end_SGtest;
end;

