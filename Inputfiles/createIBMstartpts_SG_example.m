% createIBMstartpts_SG_example.m to create particle start points for the
% particle_spec.dat file of the HAL Lagrangian model. This example is for 
% a single release of 10 particles on the northern South Georgia shelf

% ***Ensure directory paths and names of input and output files are appropriate

% Print option (1=on)
ipr=0;

% Depth range for release sites
d1=0;
d2=300;

% Depth of particles
driftdep=50;

% Number of release events
nrel=1;

% Number of particles per grid cell
np=1;

% Grid files of underpinning oceanographic model
fname2='/Inputfiles/SG_mesh_hgr.nc';
fname3='/Inputfiles/SG_mesh_zgr.nc';

x=ncread(fname2,'glamt');
y=ncread(fname2,'gphit');
wdep=ncread(fname3,'gdepw_0');
mbathy=ncread(fname3,'mbathy');
vlat=ncread(fname2,'gphiv');

nx=size(mbathy,1);
ny=size(mbathy,2);
nz=size(wdep,3);

dx=x(2,1)-x(1,1);
x2=x-dx/2;
dy=y(1,2)-y(1,1);
y2=y-dy/2;

x1d=x(:,1);
y1d=y(1,:);

for j=1:ny;
 for i=1:nx;
  D(i,j)=wdep(i,j,(mbathy(i,j)+1));
 end
end

% Specify region for release

S1x=[-36.5 -36.0];
S1y=[-53.99 -53.97];

% Find suitable cells

suitloc=find(~isnan(D) & D <= d2 & D > d1 & ...
 ((x >= S1x(1) & x <= S1x(2) & y >= S1y(1) & y <= S1y(2))));

npts=length(suitloc)

% Number of particles per release period

nprel=npts*np

ic=1;
for j=1:nrel;
 for i=1:npts;
  for k=1:np;
   lonp(ic)=x(suitloc(i));
   latp(ic)=y(suitloc(i));
% Set release depth of particles
   depp(ic)=driftdep;
   if (D(suitloc(i)) < driftdep+2);
    disp('Release point too shallow')
    disp(D(suitloc(i)))
    return
   end
   ic=ic+1;
  end;
 end;
end;

nptot=np*npts*nrel

% Plot positions

lonmin=-39.5;
lonmax=-34.0;
latmin=-55.0;
latmax=-53.5;

figure('Name','release points');
m_proj('Mercator','lat',[latmin latmax],'lon',[lonmin lonmax]);
set(gca,'fontsize',16);
h=m_plot(lonp,latp,'bo');
set(h,'MarkerFaceColor',['b'],'MarkerSize',2.0);
% Overlay bathymetry
v=[250 500 1000 2500 5000 7500];
hold on;
[c,h]=m_contour(x,y,D,v,'k');
m_gshhs_f('patch',[0 0 0]);
m_grid;
title('Particle start points');
hold off;
drawnow;

if (ipr == 1);
 orient tall;
 print -dtiff startpos_SG_example
end;

% Output particle positions in required ascii format

fileo='/Outputfiles/particle_spec.dat.SGtest';

nplank=180; % Length of time for simulating transport (days)
stdy=3; % Release Jday
enddy=stdy+nplank;
weight=10; % Weight only used for 'super-individuals' with mortality option
sthr=12; % Release particles at midday
endhr=12;

fido=fopen(fileo,'w');

ic=1;
for j=1:nrel;
% Header
 fprintf(fido,'%i\n',j);
 fprintf(fido,'%i\n',nprel);
 for i=1:nprel;
  latd=fix(latp(ic));
  if (latp(ic) > 0 | (latp(ic) < 0 & latd == 0));
   latm=(latp(ic)-latd)*60.;
  else;
   latm=(latp(ic)-latd)*(-60.);
  end
  lond=fix(lonp(ic));
  if (lonp(ic) > 0 | (lonp(ic) < 0 & lond == 0));
   lonm=(lonp(ic)-lond)*60.;
  else
   lonm=(lonp(ic)-lond)*(-60.);
  end
  fprintf(fido,'%i %f \n',latd,latm);
  fprintf(fido,'%i %f \n',lond,lonm);
  fprintf(fido,'%i %f\n',stdy,sthr);
  fprintf(fido,'%i %f\n',enddy,endhr);
  fprintf(fido,'%f\n',depp(ic));
  fprintf(fido,'%i \n',weight);
  ic=ic+1;
 end
 stdy=stdy+5; %If nrel>1, this sets number of days between releases
 enddy=enddy+5;
end

fclose(fido);

