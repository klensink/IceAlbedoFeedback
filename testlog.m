clear;

max = 7;
min = 2;
mags = [];
power = 10;

for i = min:max
    
    mags = [mags; i*ones(power.^(max - i),1)];
    
end

%%
n = histc(mags,[min:max]);

subplot 211
bar(n);
xticklabel = num2str([min:max]);
xticklabel(xticklabel == ' ') = [];
set(gca,'Yscale','log','XTickLabel', num2cell(xticklabel));
ylim([10^-1 10^8])
%% 
num_events = 100000;

r = randsample(mags,num_events);

nr = histc(r,[min:max]);

subplot 212
bar(nr);
set(gca,'Yscale','log','XTickLabel', num2cell(xticklabel));
ylim([10^-1 10^8])