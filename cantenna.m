%
% Tutorials / bent patch antenna
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_Bent_Patch_Antenna
%
% Tested with
%  - Matlab 2011a / Octave 4.0
%  - openEMS v0.0.33
%
% (C) 2013-2015 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm

%% variable declaration
cylinder.radius = 37.5; % Circular waveguide radius
cylinder.length = 112.5; % Cicular waveguide length

monopole.length = 37.5; % monopole length
monopole.radius = 1;  % monopole raduis
monopole.dist = 37.5;    % monopole position

%setup feeding
feed.R = 50;     %feed resistance

%% setup FDTD parameter & excitation function
FDTD = InitFDTD('CoordSystem', 1);          % init a cylindrical FDTD
f0 = 2e9; % center frequency
fc = 1e9; % 20 dB corner frequency
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'PEC' 'PEC' 'PEC' 'PEC' 'PEC' 'PEC'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
% init a cylindrical mesh
CSX = InitCSX('CoordSystem',1);

substrate.epsR = 1;
substrate.kappa  = 1e-3 * 2*pi*2.45e9 * EPS0*substrate.epsR;

%% Creating air substrate to fill inside the circular waveguide
CSX = AddMaterial( CSX, 'air');
CSX = SetMaterialProperty( CSX, 'air', 'Epsilon',substrate.epsR, 'Kappa', substrate.kappa);

%% Adding circular waveguide
CSX = AddMetal(CSX, 'circ_waveguide'); % Creating the metal for circular waveguide
CSX = AddCylinder(CSX, 'circ_waveguide', 8, [0,0,0], [0,0,cylinder.length], cylinder.radius+1);

%% Fill air inside circular waveguide
CSX = AddCylinder(CSX, 'air', 10, [0,0,0], [0,0,cylinder.length], cylinder.radius);

%% Covering one side of the waveguide
CSX = AddCylinder(CSX, 'circ_waveguide', 10, [0,0,0], [0,0,-1], cylinder.radius+1);

%% Creating monopole inside circular waveguide
CSX = AddMetal(CSX, 'monopole');
start = [monopole.length-1, 0, monopole.dist];
stop = start + [-monopole.length, 0, 0];
CSX = AddBox(CSX, 'monopole', 10, start, stop);
%CSX = AddCylinder(CSX, 'monopole', 10, [monopole.length-1, 0, monopole.dist], [-1, 0, monopole.dist], monopole.radius);

%% apply the excitation & resist as a current source
start = [monopole.length, 0, monopole.dist];
stop  = start + [-1, 0, 0];
[CSX port] = AddLumpedPort(CSX, 50 ,1 ,feed.R, start, stop, [1 0 0], true);


%% finalize the mesh
% detect all edges
mesh = DetectEdges(CSX);

% generate a smooth mesh with max. cell size: lambda_min / 20
max_res = c0 / (f0+fc) / unit / 20;
mesh = DetectEdges(CSX, [], 'SetProperty','monopole');
mesh.r = [monopole.length monopole.length-1 monopole.length+1 1];
mesh = SmoothMesh(mesh, c0 / (f0+fc) / unit / 20);
mesh_temp = [0:1:cylinder.length+10 -2*pi:2*pi/30:2*pi -10:1:cylinder.length+10];
mesh.r = [mesh.r 0:0.5:cylinder.length+10 monopole.length monopole.length-1 1];
mesh.a = [mesh.a -pi:(2*pi/40):pi];
mesh.z = [mesh.z -10:0.5:cylinder.length+10 monopole.dist monopole.dist-monopole.radius monopole.dist+monopole.radius];
disp(['Num of cells: ' num2str(numel(mesh.r)*numel(mesh.a)*numel(mesh.z))]);
CSX = DefineRectGrid( CSX, unit, mesh );

%% create nf2ff, keep some distance to the boundary conditions, e.g. 8 cells pml
start = [mesh.r(3)     mesh.a(3)     mesh.z(3)];
stop  = [mesh.r(end-2) mesh.a(end-2) mesh.z(end-2)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'Directions',[1 1 1 1 1 1]);

%% prepare simulation folder & run
Sim_Path = ['tmp_' mfilename];
Sim_CSX  = [mfilename '.xml'];

[status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % create empty simulation folder

% write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

% show the structure
CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

%% postprocessing & do the plots
freq = linspace( max([1e9,f0-fc]), f0+fc, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;
P_in = 0.5*real(port.uf.tot .* conj(port.if.tot)); % antenna feed power

% plot feed point impedance
figure
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
figure
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

drawnow

%find resonance frequncy from s11
f_res_ind = find(s11==min(s11));
f_res = freq(f_res_ind);

%%
%disp('dumping resonant current distribution to vtk file, use Paraview to visualize');
%ConvertHDF5_VTK([Sim_Path '/Jt_patch.h5'],[Sim_Path '/Jf_patch'],'Frequency',f_res,'FieldName','J-Field');

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the far field at phi=0 degree
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, [-180:2:180]*pi/180, 0, 'Outfile','pattern_phi_0.h5');
% normalized directivity as polar plot
figure
polarFF(nf2ff,'xaxis','theta','param',1,'normalize',1)

% calculate the far field at phi=0 degree
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, pi/2, (-180:2:180)*pi/180, 'Outfile','pattern_theta_90.h5');
% normalized directivity as polar plot
figure
polarFF(nf2ff,'xaxis','phi','param',1,'normalize',1)

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax) ' (' num2str(10*log10(nf2ff.Dmax)) ' dBi)'] );
disp( ['efficiency: nu_rad = ' num2str(100*nf2ff.Prad./real(P_in(f_res_ind))) ' %']);

drawnow

%%
disp( 'calculating 3D far field pattern and dumping to vtk (use Paraview to visualize)...' );
thetaRange = (0:2:180);
phiRange = (0:2:360) - 180;
nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Verbose',1,'Outfile','3D_Pattern.h5');

figure
plotFF3D(nf2ff,'logscale',-20);

