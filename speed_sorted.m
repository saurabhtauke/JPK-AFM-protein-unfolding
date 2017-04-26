 clear all 

%% stiffness

pist = ls('m*.mat')

index = 0;

for i =1:length(pist(:,1))
    
    pop = load(pist(i,:));
    fy = fieldnames(pop);
    
    peak_dat = pop.(fy{1});
    
    %index = index+1;
    jindex = 0;
    
    for j = 1:length(peak_dat.stiffness)
        jindex = jindex+1;
        
        if not(peak_dat.stiffness(jindex,2)== 0)
            index = index+1;
            stiff_peaks(index,1) = peak_dat.stiffness(jindex,2);
        end
    end
end


%% force

index = 0;

for i =1:length(pist(:,1))
    
    pop = load(pist(i,:));
    fy = fieldnames(pop);
    
    peak_dat = pop.(fy{1});
    
    %index = index+1;
    jindex = 0;
    
    for j = 1:length(peak_dat.force)
        jindex = jindex+1;
        
        if not(peak_dat.force(jindex,2)== 0)
            index = index+1;
            force_peaks(index,1) = peak_dat.force(jindex,2);
        end
    end
end


%% damping

index = 0;

for i =1:length(pist(:,1))
    
    pop = load(pist(i,:));
    fy = fieldnames(pop);
    
    peak_dat = pop.(fy{1});
  
   % index = index+1;
    jindex = 0;
    
    for j = 1:length(peak_dat.damping)
       jindex = jindex+1;
        
         if not(peak_dat.damping(jindex,2)== 0)
            index = index+1;
            damp_peaks(index,1) = peak_dat.damping(jindex,2);
         end
    end
end

%% analysis

avg_force = mean(force_peaks);
std_force = std(force_peaks);

avg_stiff = mean(stiff_peaks);
std_stiff = std(stiff_peaks);

avg_damp = mean(damp_peaks);
std_damp = std(damp_peaks);

frequency =  input('input frequency...');
pulling_speed = input('input pulling speed...');
fname = ['avg_std' num2str(pulling_speed) '.mat'];
path = 'D:\saurabh\Data analysis\TEMP DELETE\smaol\';
filen = [path num2str(frequency) fname];

save (filen,'frequency','pulling_speed','avg_force','std_force','avg_stiff','std_stiff','avg_damp','std_damp')

%% plotting

% fig1 = figure;
% hist(force_peaks,5)
% title('Force')
% 
% fig2 = figure;
% hist(stiff_peaks,5)
% title('Stiffness')
% 
% fig3  =figure;
% hist(damp_peaks,5)
% title('Dissipation')
