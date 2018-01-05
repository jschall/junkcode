% Tutorial on time delay and signal integrity for radar 
% and UWB applications
% 
% Tested with
%  - Octave 4.0
%  - openEMS v0.0.35
% 
% Author: Georg Michel, 2016

clear;
close all;

display_geom = 0;
run_sim = 1;

physical_constants;

% --- start of configuration section ---

% In radar and ultrawideband applications it is important to know the
% delay and fidelity of RF pulses. The delay is the retardation of the
% signal from the source to the phase center of the antenna. It is
% composed out of linear delay, dispersion and minimum-phase
% delay. Dispersion due to waveguides or frequency-dependent
% permittivity and minimum-phase delay due to resonances will degrade
% the fidelity which is the normalized similarity between excitation and
% radiated signal. In this tutorial you can examine the performance of a
% simple ultrawideband (UWB) monopole. The delay and fidelity of this
% antenna are calculated and plotted. You can compare these properties
% in different channels.
% 
% The Gaussian excitation is set to the same 3dB bandwidth as the
% channels of the IEEE 802.15.4 UWB PHY. One exeption is channel4twice
% which has the double bandwidth of channel 4. It can be seen that the
% delay is larger and the fidelity is smaller in the vicinity of the
% (undesired) resonances of the antenna. Note that for a real UWB system
% the total delay and fidelity result from both the transmitting and
% receiving antenna or twice the delay and the square of the fidelity
% for monostatic radar. 
% 
% The resolution of the delay will depend on the 'Oversampling'
% parameter to InitFDTD.  See the description of DelayFidelity
% 
% In the configuration section below you can uncomment the respective
% parameter settings. As an exercise, you can examine how the permittivity 
% of the substrate influences gain, delay and fidelity.  


%suffix = "channel1";
%f_0 = 3.5e9; % center frequency of the channel
%f_c = 0.25e9 / 0.3925; % 3dB bandwidth is 0.3925 times 20dB bandwidth for Gaussian excitation

%suffix = "channel2";
%f_0 = 4.0e9; % center frequency of the channel
%f_c = 0.25e9 / 0.3925; 

%suffix = "channel3";
%f_0 = 4.5e9; % center frequency of the channel
%f_c = 0.25e9 / 0.3925; 

%  suffix = "channel4";
%  f_0 = 4.0e9; % center frequency of the channel
%  f_c = 0.5e9 / 0.3925;

%suffix = "channel5";
%f_0 = 6.5e9; % center frequency of the channel
%f_c = 0.25e9 / 0.3925; 

suffix = "channel7";
f_0 = 6.4896e9; % center frequency of the channel
f_c = 0.4992e9 / 0.3925;

%suffix = "channel4twice"; % this is just to demonstrate the degradation of the fidelity with increasing bandwidth
%f_0 = 4.0e9; % center frequency of the channel
%f_c = 1e9 / 0.3925; 

% --- end of configuration section ---

% path and filename setup
Sim_Path = 'tmp';
Sim_CSX = 'uwb.xml';

lambda0 = C0/f_0;

helix_radius = lambda0/(2*pi);
helix_turns = 1;
helix_pitch = lambda0*.23;
feed_height = 1e-3;
ground_radius = lambda0/2;
ground_height = 1e-3;
feed_R = 50;

ang = linspace(0,2*pi,200);
coil_x = helix_radius*cos(ang);
coil_y = helix_radius*sin(ang);
coil_z = ang/2/pi*helix_pitch;

helix.x=[];
helix.y=[];
helix.z=[];
zpos = feed_height;
for n=1:helix_turns
    helix.x = [helix.x coil_x];
    helix.y = [helix.y coil_y];
    helix.z = [helix.z coil_z+zpos];
    zpos = zpos + helix_pitch;
end

helix_points(1,:) = helix.x;
helix_points(2,:) = helix.y;
helix_points(3,:) = helix.z;

mesh_res = C0 / (f_0 + f_c) / 20;

CSX = InitCSX();
CSX = AddMetal(CSX, 'Ground');
CSX = AddMetal(CSX, 'Helix');
CSX = AddCylinder(CSX,'Ground',1,[0 0 -ground_height],[0 0 0],ground_radius);
CSX = AddCurve(CSX, 'Helix', 0, helix_points);
CSX = AddWire(CSX, 'Helix', 0, helix_points, 3.22e-4);

[CSX port] = AddLumpedPort(CSX, 2, 1, feed_R, [helix_radius 0 0], [helix_radius 0 feed_height], [0 0 1], true);


mesh.x = [-(2*lambda0+ground_radius) -helix_radius 0 helix_radius (2*lambda0+ground_radius)];
mesh.y = [-(2*lambda0+ground_radius) -helix_radius 0 helix_radius (2*lambda0+ground_radius)];
mesh.z = [-(2*lambda0+ground_height) 0 feed_height (2*lambda0+feed_height+helix_turns*helix_pitch)];
mesh = SmoothMesh(mesh, mesh_res/2);
CSX = DefineRectGrid(CSX, 1, mesh);


start = [mesh.x(10)     mesh.y(10)     mesh.z(10)];
stop  = [mesh.x(end-9) mesh.y(end-9) mesh.z(end-9)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop);

FDTD = InitFDTD( 'NrTs', 30000, 'EndCriteria', 1e-5, 'OverSampling', 20);
FDTD = SetGaussExcite(FDTD, f_0, f_c );
BC   = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'};
FDTD = SetBoundaryCond(FDTD, BC);

if (run_sim==1 || display_geom==1)
    confirm_recursive_rmdir(0);
    [status, message, messageid] = rmdir( Sim_Path, 's' );
    [status, message, messageid] = mkdir( Sim_Path );

    WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);

    if (display_geom==1)
        CSXGeomPlot([Sim_Path '/' Sim_CSX]);
        return;
    end

    RunOpenEMS( Sim_Path, Sim_CSX);
end



freq  = linspace(f_0-f_c, f_0+f_c, 200);
port = calcPort(port, Sim_Path, freq);

phi = [0] * pi / 180;
theta = [-180:5:180] * pi / 180;

[delay, fidelity, nf2ff] = DelayFidelity(nf2ff, port, Sim_Path, -1i, 1, theta, phi, f_0, f_c, 'Mode', 0, 'center', [0 0 0]);

figure;
plot(theta.' * 180/pi, delay*C0*1000);
title(["delay ", suffix, " (mm)"]);

figure;
plot(theta.' * 180/pi, fidelity)
title(["fidelity"]);

f_0_nearest_ind = nthargout(2, @min, abs(nf2ff.freq -f_0));
nf2ff.Dmax(f_0_nearest_ind) *= nf2ff.Prad(f_0_nearest_ind) / calcPort(port, Sim_Path, nf2ff.freq(f_0_nearest_ind)).P_inc;
figure
plotFFdB(nf2ff, 'xaxis', 'theta', 'freq_index', f_0_nearest_ind);
title(strcat("gain ",suffix," / dBi"));

% save the plots in order to compare them afer simulating the different channels
%  print(1, ["s11_", suffix, ".png"]);
%  print(2, ["farfield_", suffix, ".png"]);
%  print(3, ["delay_mm_", suffix, ".png"]);
%  print(4, ["fidelity_", suffix, ".png"]);

% calculate the far field at phi=0 degrees and at phi=90 degrees

thetaRange = unique([0:2:90 90:4:180]);
phiRange = (0:8:360) - 180;
disp( 'calculating the 3D far field...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_0, thetaRange*pi/180, phiRange*pi/180,'Mode',0,'Outfile','3D_Pattern.h5','Verbose',1);

directivity = nf2ff.P_rad{1}/nf2ff.Prad*4*pi;
directivity_CPRH = abs(nf2ff.E_cprh{1}).^2./max(nf2ff.E_norm{1}(:)).^2*nf2ff.Dmax;
directivity_CPLH = abs(nf2ff.E_cplh{1}).^2./max(nf2ff.E_norm{1}(:)).^2*nf2ff.Dmax;

DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],directivity,thetaRange,phiRange);
DumpFF2VTK([Sim_Path '/3D_Pattern_CPRH.vtk'],directivity_CPRH,thetaRange,phiRange);
DumpFF2VTK([Sim_Path '/3D_Pattern_CPLH.vtk'],directivity_CPLH,thetaRange,phiRange);
return;
