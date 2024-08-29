% plot_traject_SGtest.m to plot trajectories of particles from HAL
% Set up for SGtest simulations forced with SG NEMO

% Option to print (1 to print)
iprint=0;

nplank=180; % Length of time from particle release to end of trajectories
            % in days

% Option for difference line colours (0 for all the same colour)
linecol=1;

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

% Regular grid in both x and y
dx=x(2,1)-x(1,1);
dy=y(1,2)-y(1,1);

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

% Extract trajectories

ic=size(lon,2); % Total number of outputs
for i=1:npts
% Find where particle is switched on (flag to 1)
% Adjustment assumes output every 1/2 hour
 parton=find(flag(i,:) == 1);
 firstt=min(parton);
 bg=firstt+(nplank*48)-1;
 lon(i,bg+1:ic)=nan;
 lat(i,bg+1:ic)=nan;
end

% Define different line colours if linecol=1

if (linecol == 1);
 MAP=colormap(jet(npts));
end;

% Plot trajectories

figure(1)
clf
m_proj('Mercator','lat',[-55 -49],'lon',[-43 -30]);
for i=1:npts;
 if (linecol == 0);
  m_plot(lon(i,:),lat(i,:),'r-');
 else
  m_plot(lon(i,:),lat(i,:),'color',MAP(i,:),'linewidth',2);
 end;
 if (i == 1);
  hold on;
 end;
end;
% Overlay bathymetry
v=[250 500 1000 2500 5000 7500];
[c,h]=m_contour(x,y,beddep,v,'k');
m_grid;
title('Predicted particle trajectories - SGtest');
m_usercoast('/Inputfiles/SScoast5.mat','patch',[.6 .6 .6]);
hold off;
drawnow;

if (iprint == 1);
 orient landscape;
 print -dtiff part_traj_SGtest
end;


