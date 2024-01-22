city = unique(dat{:,2});
time = [];
for y = 1999:2019
time = [time; y*ones(12,1), [1:12]'];
end

out = zeros(291,252);
for i = 1:length(city)
    cond = dat.ZipCode == city(i,1);
    int = dat(cond,:);
    for j = 1:height(int)
        idx = find(time(:,1) == int.Year(j) & time(:,2) == int.Month(j));
        out(i,idx) = 1;
    end
end