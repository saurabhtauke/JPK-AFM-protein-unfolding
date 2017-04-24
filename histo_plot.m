clear all 

%% stiffness

pist = ls('*.mat')

index = 0;

for i =1:length(pist(:,1))
    
    pop = load(pist(i,:));
    fy = fieldnames(pop);
    
    peak_dat = pop.(fy{1});
    
   % index = index+1;
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


%% plotting

fig1 = figure;
hist(force_peaks,12)
title('Force')

fig2 = figure;
hist(stiff_peaks,12)
title('Stiffness')

fig3  =figure;
hist(damp_peaks,12)
title('Dissipation')