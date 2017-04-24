%% JPK output to compatible files

% this is made because i forgot to add the force, stiffness and dissipation
% data to the saved mat file.
%sorry 
%

clear all
%% listing and adding files to array

%navigate to the folder to list.
% this assumes the standard format 'map-data-year.month.day-time.numbers.txt'

list = ls('map*.txt')
pist = ls('*.mat')

%SOMEHOW PUT THIS INTO AN ARRAY OF STRINGS PLS PLS PLS
array = cellstr(list);
parray = cellstr(pist);

%indexing for extraction of 
sens = zeros(length(array),1);
springk = zeros(length(array),1);

%% Open and delet


for file_index = 1: length(array)
    
    junk = fileread(list(file_index,:));   

    sindex = strfind(junk,'% sensitivity:');
        senstring = junk(sindex+15:sindex+35);
        sens(file_index) = str2num(senstring);
    
    kindex = strfind(junk,'% springConstant:');
        kstring = junk(kindex+18:kindex+35);
        springk(file_index) = str2num(kstring);
    
end

%%
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH A E S T H E T I C
%HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH


drive_frequency1 = input('pls input frequency in hertz...');            % in Hz
drive_frequency = drive_frequency1*2*pi
pulling_speed = input('pulling speed...');
cantilever_length = 130 *10^(-6);      %length of cantilver SI units
cantilever_width =  32.5 *10^(-6);      %width of cantilver SI units
cantilever_thickness = 1 *10^(-6);    


% cantilever_length = input('pls input cantilever length in microns...') *10^(-6);      %length of cantilver SI units
% cantilever_width = input('pls input cantilever width in microns...') *10^(-6);      %width of cantilver SI units
% cantilever_thickness = input('pls input cantiver thickness in microns...')*10^(-6);    


% permission = input('input a number to continue...');
% file_index

for file_index = 1 :length(array)
    
%% Experimental constants

lockin_sens = 10;                  % in mv

cantilever_stiffness = springk(file_index);          % in N/m
detection_sensitivity = sens(file_index);

%free_amplitude =  1*10^-10;                % in m

list(file_index,:)

%drive_frequency = 1200*2*pi;
density = 2329;                     %density of material of cantilever

%cantilever_length = 350*10^(-6);

%cantilever_width = 35e-006;      %length of cantilver SI units

%cantilever_thickness = 2*10^(-6);

area = cantilever_thickness*cantilever_width;                        %surface area of cantilever

lever_damping = 10e-006;

%% Removing the approach values and rescaling the z_voltage values and converting to nanometers

%Import and allocate data
temp= importdata(list(file_index,:),' ',10);
temp1= temp.('data');

z = temp1(:,1);
totlength = length(z);

%trimming initial 50 points

for n=1:50
    temp1(1,:) = [];
    
end

%% data allocation

z_range = temp1(:,1);

x_signal = temp1(:,5)*(detection_sensitivity* lockin_sens/(10*1000));
x1 = mean(x_signal((length(x_signal)-100):length(x_signal)));

y_signal = -temp1(:,6)*(detection_sensitivity* lockin_sens/(10*1000));
y1 = mean(y_signal((length(y_signal)-100):length(y_signal)));


x_signal = temp1(:,5)*(detection_sensitivity* lockin_sens/(10*1000*cantilever_length));
y_signal = -temp1(:,6)*(detection_sensitivity* lockin_sens/(10*1000*cantilever_length));


amplitude = sqrt(x_signal.^2 + y_signal.^2);
free_amplitude = sqrt(x1^2 + y1^2);

phase = atan(y_signal./x_signal) ;                        
force = temp1(:,2);            



%% Stiffness and Damping calculation   


stiffness = double(cantilever_stiffness* ((free_amplitude./amplitude).*(cos(phase)) -1)); 
%stiffness = stiffness - min(stiffness);

damping = double(cantilever_stiffness * (1.0) * ((free_amplitude./(amplitude.*(drive_frequency))) .* (sin((phase))))) ;
%damping = damping- min(damping);


stiffx = double(-1)*(((0.666*cantilever_stiffness * cantilever_length).*x_signal./free_amplitude)- 0.333*density*10*area*cantilever_length*drive_frequency*drive_frequency);
%stiffx = stiffx - min(stiffx);
free_stiff = mean(stiffx((length(stiffx)-100):length(stiffx)));
stiffx = stiffx-free_stiff;

dampy = double((0.666*cantilever_stiffness*cantilever_length.*y_signal./(free_amplitude.*drive_frequency))- (0.333*lever_damping));
%dampy = double((0.666*cantilever_stiffness*cantilever_length.*y_signal./(free_amplitude.*drive_frequency)));
free_damp = mean(dampy((length(dampy)-100):length(dampy)));
dampy = dampy - free_damp;

relaxation_time = double(dampy./(stiffx));

force = force - force(length(force));

%% save plots

v = ['fsd_dat' num2str(file_index)];
rite = struct('force',force,'stiffness',stiffx,'damping',dampy,'frequency',drive_frequency1,'pulling_speed',pulling_speed,'cantilever_stiffness',cantilever_stiffness,'cantilever_length',cantilever_length,'cantilever_width',cantilever_width,'cantilever_thickness',cantilever_thickness);
eval([v ' = rite']);


findex_name =  list(file_index,:);
txt_name = findex_name(length(list(file_index,:))-6 : length(list(file_index,:))-4);
findex_name(length(list(file_index,:))-3 : length(list(file_index,:))) = '.mat';

rot = regexp(list(file_index,:),txt_name);
if not(isempty(rot))
    rot = rot(length(rot));
end

for ri = 1:length(parray)
    
    pot = regexp(pist(ri,:),txt_name);
    if not(isempty(pot))
        pot = pot(length(pot));
    end
    
    youi = rot == pot;
    
    if youi==1
        disp('analyzed')
        save (findex_name,v,'-append');
        
    elseif isempty(youi)
        disp('not analyzed')
    end
    
end

end

