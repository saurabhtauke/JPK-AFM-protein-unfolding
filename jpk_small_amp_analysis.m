%% JPK output to compatible files

%This program converts the JPK output data files which are then compatible
%with the MATLAB code for Data Analysis of the Protein unfolding Force
%Spectroscopy experiments. MATLAB Code used for analysis...
%(https://github.com/saurabhtauke/Small-Amplitude-AFM/blob/master/analysis.m)

clear all
%% listing and adding files to array

%navigate to the folder to list.
% this assumes the standard format 'map-data-year.month.day-time.numbers.txt'

list = ls('map*.txt')

%SOMEHOW PUT THIS INTO AN ARRAY OF STRINGS PLS PLS PLS
array = cellstr(list);

%indexing for extraction of 
sens = zeros(length(array),1);
springk = zeros(length(array),1);

%% Open and delet


for i = 1: length(array)
    
    junk = fileread(list(i,:));
    retr = strfind(junk, '# segment: retract');
    junk(1:retr-1) = [];
    

    % in newdat, find '#' and replace with '%'
    hash = strfind(junk,'#');
    
    for j = 1:length(hash)
        
        junk(hash(j)) = '%';
    end
    
    sindex = strfind(junk,'% sensitivity:');
        senstring = junk(sindex+15:sindex+35);
        sens(i) = str2num(senstring);
    
    kindex = strfind(junk,'% springConstant:');
        kstring = junk(kindex+18:kindex+35);
        springk(i) = str2num(kstring);
    
    %save this into a file
    fid = fopen(list(i,:),'w+');
    fprintf(fid,'%s\r\n',junk);
    fclose(fid);
end

%%
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

for i = 1 :length(array)
    
%% Experimental constants

lockin_sens = 10;                  % in mv

cantilever_stiffness = springk(i);          % in N/m
detection_sensitivity = sens(i);

free_amplitude = 1 *10^-10;                % in m
drive_frequency = 1000*2*pi;            % in Hz

density = 2329;                     %density of material of cantilever
cantilever_length = 350e-006;      %length of cantilver SI units
cantilever_width = 35e-006;      %length of cantilver SI units
area = cantilever_length*cantilever_width;                        %surface area of cantilever

lever_damping = 10e-006;

%% Removing the approach values and rescaling the z_voltage values and converting to nanometers

% trimm all the arrays to remove approach z values.

temp= importdata(list(i,:));
%temp= a.('data');

z = temp(:,1);
totlength = length(z);


%% data allocation

z_range = temp(:,1);
x_signal = temp(:,5)*(detection_sensitivity* lockin_sens/10000);
y_signal = -temp(:,6)*(detection_sensitivity* lockin_sens/10000);

amplitude = sqrt(x_signal.^2 + y_signal.^2);
phase = atan(y_signal./x_signal) ;                        
force = temp(:,2);            

%% Stiffness and Damping calculation   


stiffness = double(cantilever_stiffness* ((free_amplitude./amplitude).*(cos(phase)) -1)); 
%stiffness = stiffness - min(stiffness);

damping = double(cantilever_stiffness * (1.0) * ((free_amplitude./(amplitude.*(drive_frequency))) .* (sin((phase))))) ;
%damping = damping- min(damping);


stiffx = double(-1)*(((0.666*cantilever_stiffness * cantilever_length).*x_signal./free_amplitude)- 0.333*density*10*area*cantilever_length*drive_frequency*drive_frequency);
stiffx = stiffx - min(stiffx);


dampy = double((0.666*cantilever_stiffness*cantilever_length.*y_signal./(free_amplitude.*drive_frequency))- (0.333*lever_damping*cantilever_length));
%dampy = double((0.666*cantilever_stiffness*cantilever_length.*y_signal./(free_amplitude.*drive_frequency)));

%dampy = dampy- min(dampy);

relaxation_time = double(dampy./(stiffx));

force = force - force(length(force));

%% Plotting

figureindex = num2str(i);
%figurename = strcat('plotting ',figureindex);

figurename = figure;

 %'please refer to code for plot_matrix labels'
 plot_matrix = zeros(length(z_range),12);
 plot_matrix(:,1) = z_range;
 plot_matrix(:,2) = amplitude;
 plot_matrix(:,3) = phase;
 plot_matrix(:,4) = x_signal;
 plot_matrix(:,5) = y_signal;
 plot_matrix(:,6) = force;
 plot_matrix(:,7) = stiffx;
 plot_matrix(:,8) = dampy;
 plot_matrix(:,9) = relaxation_time;
 plot_matrix(:,10) = stiffness;
 plot_matrix(:,11) = damping;

% subplot(4,2,1)
% plot(z_range,amplitude,'-b')
% title('Amplitude')
% xlabel('Distance(nm)')  
% ylabel('Amplitude(Å)')


% phase = phase.*180./pi;
% subplot(4,2,2)
% plot(z_range,phase,'-b')
% title('Phase(Deg.)')
% xlabel('Distance(nm)')
% ylabel('Phase(degrees)')


% subplot(4,2,3)
% plot(z_range,x_signal,'-b')
% title('X')
% xlabel('Distance(nm)')
% ylabel('X(Å)')


% subplot(4,2,4)
% plot(z_range,y_signal,'-b')
% title('Y')
% xlabel('Distance(nm)')
% ylabel('Y(Å)')

% subplot(4,2,5)
% plot(z_range,stiffx,'-b')
% title('Stiffness by Altered Boundary Condition')
% xlabel('Distance(nm)')
% ylabel('Stiffness(N/m)')
% 
% subplot(4,2,6)
% plot(z_range,dampy,'-b')
% title('Damping by Altered Boundary Condition')
% xlabel('Distance(nm)')
% ylabel('Damping(Kg/s)')

% subplot(4,2,7)
% plot(z_range,stiffness,'-b')
% title('Stiffness By Point-Mass Model')
% xlabel('Distance(nm)') 
% ylabel('Stiffness(N/m)')
% 
% 
% subplot(4,2,8)
% plot(z_range,damping,'-b')
% title('Damping by Point-Mass Model')
% xlabel('Distance(nm)')
% ylabel('Damping(Kg/s)')


% subplot(4,2,7)
% plot(z_range,force,'-b')
% title('Force')
% xlabel('Distance(nm)')
% ylabel('Force(N)')
% 
% subplot(4,2,8)
% plot(z_range,relaxation_time,'-b')
% title('relaxation time')
% xlabel('Distance(nm)')
% ylabel('relaxation time')

%% PLOT 2

subplot(2,2,1)
plot(z_range,stiffx,'-b')
title('Stiffness by Altered Boundary Condition')
xlabel('Distance(nm)')
ylabel('Stiffness(N/m)')

subplot(2,2,2)
plot(z_range,dampy,'-b')
title('Damping by Altered Boundary Condition')
xlabel('Distance(nm)')
ylabel('Damping(Kg/s)')

subplot(2,2,3)
plot(z_range,force,'-b')
title('Force')
xlabel('Distance(nm)')
ylabel('Force(N)')

subplot(2,2,4)
plot(z_range,relaxation_time,'-b')
title('relaxation time')
xlabel('Distance(nm)')
ylabel('relaxation time')

%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH


end