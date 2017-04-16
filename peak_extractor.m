%Data Analysis Program

%AIM: to retrive data from DATA_CURSOR on plots.

%% Force


fig1 = figure;
plot (force)
title('Force')
hold on

datacursormode on

dcm_obj = datacursormode(fig1);
disp('click on the datapoint')
disp('input 1 for storing the data')
disp('input 0 for storing the data as 0')
disp('input any other number to exit')

rubbish = input('...');

jindex = 0;

while rubbish == 1 || rubbish ==0
    jindex = jindex+1;

    c_info = getCursorInfo(dcm_obj);
    yum = c_info.Position

f_dat1(jindex,1) = yum(1);
f_dat1(jindex,2) = yum(2);

if rubbish ==0
    f_dat(jindex,2) = 0;
end

plot(f_dat1(:,1),f_dat1(:,2),'ro')

rubbish = input('ZERO or ONE...batao... ')
end

jindex = 0;

rubbish = input('input 0 or 1...');

while rubbish == 1 || rubbish == 0
    jindex = jindex+1;

    c_info = getCursorInfo(dcm_obj);
    yum = c_info.Position

f_dat2(jindex,1) = yum(1);
f_dat2(jindex,2) = yum(2);

if rubbish == 0
    f_dat2(jindex,2) = 0;
end

plot(f_dat2(:,1),f_dat2(:,2),'go')

rubbish = input('ZERO or ONE...batao... ')
end
 
hold off
 
force_dat = f_dat1-f_dat2;

%% Stiffness

fig2 = figure
plot (stiffx)
title('stiffness')
hold on
datacursormode on

dcm_obj = datacursormode(fig2);
disp('click on the datapoint')

rubbish = input('input 0 or 1...');

jindex = 0;

while rubbish == 1 || rubbish == 0
    jindex = jindex+1;

    c_info = getCursorInfo(dcm_obj);
    yum = c_info.Position

stiff_dat(jindex,1) = yum(1);
stiff_dat(jindex,2) = yum(2);

if rubbish ==0
    stiff_dat (jindex,2) = 0;
end

plot(stiff_dat(:,1),stiff_dat(:,2),'ro')

rubbish = input('ZERO or ONE...batao... ')
end

hold off
%% Dissipation

fig3 = figure
plot (dampy)
title('dampy');
hold on

datacursormode on

dcm_obj = datacursormode(fig3);
disp('click on the datapoint')

rubbish = input('input 0 or 1...');

jindex = 0;

while rubbish == 1 || rubbish ==0
    jindex = jindex+1;

    c_info = getCursorInfo(dcm_obj);
    yum = c_info.Position

damp_dat(jindex,1) = yum(1);
damp_dat(jindex,2) = yum(2);

if rubbish == 0
    damp_dat(jindex,2) = 0;
end

plot(damp_dat(:,1),damp_dat(:,2),'ro')

rubbish = input('ZERO or ONE...batao... ')
end
 
hold off


%% write to file

len = [length(force_dat),length(stiff_dat),length(damp_dat)];
combolen = max(len);

%3 layers... first is force, 2nd is stiffness, 3rd is damping
combo = zeros(combolen,2,3); 

index =1;
for i = 1:size(force_dat,1)
    combo(index,1,1) = force_dat(index,1);
    combo(index,2,1) = force_dat(index,2)
    index = index+1;
end


index =1;
for i = 1:size(force_dat,1)
    combo(index,1,1) = stiff_dat(index,1);
    combo(index,2,1) = stiff_dat(index,2)
    index = index+1;
end


index =1;
for i = 1:size(force_dat,1)
    combo(index,1,1) = damp_dat(index,1);
    combo(index,2,1) = damp_dat(index,2)
    index = index+1;
end


%  fid = fopen('data1.mat','w');
%     fprintf(fid,combo);
%     fclose(fid);

%close(fig1)
%close(fig2)
%close(fig3)
