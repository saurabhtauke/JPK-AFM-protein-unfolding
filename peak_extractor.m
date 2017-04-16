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
disp('input any other NUMBER to exit')

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

stiff_dat = zeros(length(force_dat),2);

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

damp_dat = zeros(length(force_dat),2);

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

peak_dat = struct('force',force_dat,'stiffness',stiff_dat,'damping',damp_dat);

index = 1;
peak_num = zeros(length(force_dat),1);
for i = 1:length(force_dat)
    peak_num(index) = index;
    index = index+1;
end

combo = [peak_num peak_dat.force peak_dat.stiffness peak_dat.damping]

%close(fig1)
%close(fig2)
%close(fig3)
