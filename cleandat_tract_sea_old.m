time = [];
for y = 1999:2019
    time = [time; y*ones(12,1), [1:12]'];
end

dat = readtable("PV_heat_sea.csv");
shp = shaperead("sea_tract_poly.shp");
city =[1:length(shp)];
geo=[];
for i = 1:length(shp)
    geo = [geo,str2num(shp(i).geoid(8:end))];
end

polygonIDs = zeros(height(dat),1);
x = dat.Longitude;
y = dat.Latitude;
% Loop through each point and check which polygon it belongs to
for i = 1:height(dat)
    isInAnyPolygon = false;
    
    % Loop through each polygon
    for j = 1:length(shp)
        % Check if the point is in the current polygon
        isInPolygon = inpolygon(x(i), y(i), shp(j).X, shp(j).Y);
        
        if any(isInPolygon)
            % Set the polygon ID for the current point
            % polygonIDs(i) = city(j);
            polygonIDs(i) = j;
            isInAnyPolygon = true;
            break; % Exit the inner loop once a polygon is found
        end
    end
    
    if ~isInAnyPolygon
        % Point is outside all polygons
        polygonIDs(i) = NaN;
    end
end

out = zeros(length(shp),252);
for i = 1:length(city)
    cond = polygonIDs == city(1,i);
    int = dat(cond,:);
    for j = 1:height(int)
        idx = find(time(:,1) == year(int.CompletedDate(j)) & time(:,2) == month(int.CompletedDate(j)));
        out(i,idx) = 1;
    end
end
out = out';
