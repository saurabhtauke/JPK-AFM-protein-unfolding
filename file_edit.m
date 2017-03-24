%% JPK output to compatible files

%This program converts the JPK output data files which are then compatible
%with the MATLAB code for Data Analysis of the Protein unfolding Force
%Spectroscopy experiments. MATLAB Code used for analysis...
%(https://github.com/saurabhtauke/Small-Amplitude-AFM/blob/master/analysis.m)


%% listing and adding files to array

%navigate to the folder to list.
% this assumes the standard format 'map-data-year.month.day-time.numbers.txt'

list = ls('map*.txt')

%SOMEHOW PUT THIS INTO AN ARRAY OF STRINGS PLS PLS PLS

array = cellstr(list);


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
    
    %save this into a file
    fid = fopen(list(i,:),'w+')
    fprintf(fid,'%s\r\n',junk)
    fclose(fid)
end


