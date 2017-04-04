% run main program before this atleast once

ultaforce = -(force);
ultaforce = ultaforce-min(ultaforce);
%ultaforce = detrend(ultaforce)

[pks,locs] = findpeaks(ultaforce);

for i = 1:length(locs)
tub(locs(i))= ultaforce(locs(i));
end


plot(ultaforce)
hold on
plot (tub,'o')
hold off
