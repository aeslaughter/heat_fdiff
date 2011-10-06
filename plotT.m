function plotT(varargin)
%__________________________________________________________________________
% PLOTT builds a contour plot of temperatures with depth and time
% 
% SYNTAX:
%   T = plotT;
%   T = plotT(filename);
%
% DESCRIPTION:
%   T = plotT;
%   T = plotT(filename);
%__________________________________________________________________________

% LOAD/SET THE LAST USED DIRECTORY
lastdir = cd;   % Default directory
if ispref('plotT','lastdir');
    lastdir = getpref('plotT','lastdir');
end

% GATHER THE FILENAME
if nargin == 0;
    [F,P] = uigetfile('*.*','Select file...',lastdir);
    if isnumeric(F); return; end;
    filename = [P,F];
else
    filename = varargin{1};
    P = fileparts(filename);
end

% READ THE FILE
    % Open the file and read header information
    fid = fopen(filename);
    C = textscan(fid,'%s%f',4,'delimiter','=');
    n = C{2}(1);    % Number of layers
    dz = C{2}(2);   % Layer thickness
    nt = C{2}(3);   % Number of time steps
    dt = C{2}(4);   % Time step

    % Read the temperature data
    for i = 1:nt+1;
        C = textscan(fid,'%f',n);
        T(:,i) = C{1};
    end
    
    % Close file and store directory
    fclose(fid);
    setpref('plotT','lastdir',P);

% CREATE CONTOUR PLOT    
    figure;
    [~,h] = contourf(T,50);
    set(h,'LineColor','none');
    xlabel(gca,'Time (s)');
    ylabel(gca,'Depth (m)');
    cbar = colorbar;
    ylabel(cbar,'Temperature (\circC)');
    colormap(jet);
    set(gca,'YDir','reverse');
    set(findall(gcf,'-property','FontName'),'FontName','Times New Roman');
    set(findall(gcf,'-property','Interpreter'),'Interpreter','Tex');
    y = get(gca,'Ytick');
    set(gca,'YtickLabel',(y-1)*dz);
    x = get(gca,'Xtick');
    set(gca,'XtickLabel',x*dt);

