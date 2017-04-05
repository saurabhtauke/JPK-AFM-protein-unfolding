%run main program before this atleast once

ultaforce = -(force);
ultaforce = ultaforce-min(ultaforce);
% ultaforce = detrend(ultaforce)

[pks,locs] = findpeaks(ultaforce);

%tub = zeros(length(locs),1);
for i = 1:length(locs)
tub(locs(i))= ultaforce(locs(i));
end

% figure
% plot(ultaforce,'r.-')
% hold on
% plot (tub,'o')
% hold off

sep = -diff(ultaforce)*100/8;

[pks1,locs1] = findpeaks(sep,'MinPeakHeight',4E-10);


%find out peak and peak" around the point of largest difference

modpeak = pks1;
modlocs = locs1

for m = 1:length(locs1)
    if ultaforce(locs1(m)-1) > ultaforce(locs1(m));
        modpeak(m) = ultaforce(locs1(m)-1);
        modlocs(m) = locs1(m)-1;
    end
end

%tub1 = zeros(length(locs1),1);
for i = 1:length(locs1)
tub1(modlocs(i))= ultaforce(modlocs(i));
end

%%%%%%%%%%%%%%%

modip = pks1;
diploc = locs1;

for m = 1:length(locs1)
    if ultaforce(locs1(m)+1) < ultaforce(locs1(m));
        modip(m) = ultaforce(locs1(m)+1);
        diploc(m) = locs1(m)+1;
    end
    
    if ultaforce(locs1(m)+2) < ultaforce(locs1(m)+1);
        modip(m) = ultaforce(locs1(m)+2);
        diploc(m) = locs1(m)+2;
    end
    
end

%tub2 = zeros(length(locs1),1);
for i = 1:length(locs1)
tub2(diploc(i))= ultaforce(diploc(i));
end
%%%%%%%%%%%%%%


figure
plot(ultaforce,'r.-')
hold on
plot (tub1,'bO')
plot (tub2,'gO')
%plot(sep,'g')
hold off

shatru = ultaforce(modlocs) - ultaforce(diploc)
figure
hist(shatru,100)
