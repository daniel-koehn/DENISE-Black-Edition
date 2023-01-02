% Plot shot gathers
%
% Daniel Koehn
%
% Kiel, the 2nd of January 2022
%
% SU matlab code partly from
% Copyright (C) 2008, Signal Analysis and Imaging Group
% http://www-geo.phys.ualberta.ca/saig/SeismicLab
% Author: M.D.Sacchi
%
% SegyMat codes from
% https://github.com/cultpenguin/segymat
% Author: Thomas Meier Hansen

clear all
close all

% Locate directories for SegyMat
Dir='segymat';
addpath(Dir);

% Define parameters
% -----------------

% Define shot number
shotno = 1;

% Define path and filename of field data
rawdata = ['DENISE_MARMOUSI_y.su.shot',int2str(shotno)];
path2data = 'su/';
rawdatafile = [path2data,rawdata];

% Define colormap
cmap = 'gray';

% Define seismic clipping value
pclip = 5e-2;

% Fontsize of text in shot gather
FSize = 20;

% define image size
imsize_x = 1900;
imsize_z = 1200;

% Import data
% -----------

% Import seismic shot gathers & headers
[D_mod,Htr,H_mod] = ReadSu(rawdatafile,'endian','l'); 

% Get number of traces ntr & number of samples ns
[ns,ntr]=size(D_mod); 

% Get time sample interval dt and define time axis t
DT=H_mod(1).dt.*1e-6;
t=1:ns;
t=t.*DT;
t=t';

% Get receiver offsets from header
for i=1:ntr
    
    offset(i) = Htr(i).offset/1000.;
    
end

% Plot shot gather 
% ----------------
Fig = figure;
figure(Fig)

imagesc(offset,t,D_mod);
colormap(cmap);

% apply image clipping
vmax = pclip * max(abs(D_mod(:)));
vmin = - vmax;
caxis([vmin vmax]);

% Add labels and title 
xlabel('offset [m]');
ylabel('time [s]');
iter_text=['Field data shot no. ',int2str(shotno)];
title(iter_text);

% Define font size of figure title, labels, linewidth
set(get(gca,'title'),'FontSize',FSize);
set(get(gca,'title'),'FontWeight','bold');
set(get(gca,'Ylabel'),'FontSize',FSize);
set(get(gca,'Ylabel'),'FontWeight','bold');
set(get(gca,'Xlabel'),'FontSize',FSize);
set(get(gca,'Xlabel'),'FontWeight','bold');
set(gca,'FontSize',FSize);
set(gca,'FontWeight','bold');
set(gca,'Box','on');
set(gca,'Linewidth',1.0);
axis ij

% define figure size
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [5 7]);
set(Fig,'position',[0 0, imsize_x imsize_z])
set(Fig,'PaperPositionMode','Auto')       

output=['Seis_shot_',int2str(shotno)];
%saveas(Fig,output,'psc2'); 
saveas(Fig,output,'png');