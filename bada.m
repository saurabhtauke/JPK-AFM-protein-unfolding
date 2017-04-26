clear all

list = ls('avg*.mat')

index = 0;

for i =1:length(list(:,1))
    
    pop = load(list(i,:));
    fy = fieldnames(pop);
    
    agg_force(i,1) = pop.(fy{1});
    agg_stiff(i,1) = pop.(fy{3});
    agg_damp(i,1) = pop.(fy{5});
    
end


%%
frequency = [100;300;500;800;1000;1200;1400;1500;1600;1800;2000;2200;2400;2500;3000;5000;8000;10000];
speed = [50;250;400;500;800];

fid1 = figure;
semilogx(frequency,agg_force,'r-o')
title('Force Vs Frequency')

fid2 = figure;
semilogx(frequency,agg_stiff,'ro-')
title('Stiffness Vs Frequency')

fid3 = figure;
semilogx(frequency,agg_damp,'r-o')
title('Dissipation Vs Frequency')

% fid1 = figure;
% semilogx(speed,agg_force,'r-o')
% title('Force Vs pulling speed')
% 
% fid2 = figure;
% semilogx(speed,agg_stiff,'ro-')
% title('Stiffness Vs pulling speed')
% 
% fid3 = figure;
% semilogx(speed,agg_damp,'r-o')
% title('Dissipation Vs pulling speed')

