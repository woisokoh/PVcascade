time = [];
for y = 1999:2019
    time = [time; y*ones(12,1), [1:12]'];
end

dat = readtable("PV_heat_sea.csv");
shp = shaperead("tl_2019_53_tract_seat.shp");
city =[];
for i = 1:length(shp)
    city = [city,str2num(shp(i).TRACTCE)];
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
    cond = polygonIDs == i;
    int = dat(cond,:);
    for j = 1:height(int)
        idx = find(time(:,1) == year(int.CompletedDate(j)) & time(:,2) == month(int.CompletedDate(j)));
        out(i,idx) = 1;
    end
end
out = out';


census = readtable("census_sea.csv");
shp_GEOID = [];
for i = 1:length(shp)
    shp_GEOID = [shp_GEOID;str2num(shp(i).GEOID)];
end
census_GEOID = census.GEOID;
idx = zeros(length(shp_GEOID),1);
for i = 1:length(shp_GEOID)
    idx(i,1) = find(census.GEOID == shp_GEOID(i,1));
end
census_sorted = census(idx,:);

for i = 1:length(shp)
    shp(i).PopDensity = census_sorted.PopDensity(i);
    shp(i).HomeOwn = census_sorted.HomeOwn(i);
    shp(i).SingleFamily = census_sorted.SingleFamily(i);
    shp(i).Edu = census_sorted.Edu(i);
    shp(i).HomeValue = census_sorted.HomeValue(i);
    shp(i).Income = census_sorted.Income(i);
    shp(i).White = census_sorted.White(i);
    shp(i).Poverty = census_sorted.Poverty(i);
    shp(i).Gini = census_sorted.Gini(i);
end

output_shapefile = "tl_2019_53_tract_seat_census.shp";
shapewrite(shp,output_shapefile);
