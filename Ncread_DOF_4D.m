% netcdf file read master
% Shu Wu dof.shuwu@gmail.com
% file: path of data file

function [output lon lat] = Ncread_DOF_4D(file, ...
    target_varName,target_timeName,target_zName,target_lonName,target_latName,target_lon,target_lat,target_z)

target_lon = target_lon(:);
target_lat = target_lat(:);
target_lon(target_lon<0) = target_lon(target_lon<0)+360;

if target_lat(2) > target_lat(1)
    disp('Change the latitude range from South to North');
end

ncid = netcdf.open(file,'nowrite');
disp('Output from Ncread_DOF ... ')
disp('Matlab NETCDF PACKAGE VERSION')
netcdf.inqLibVers
format = netcdf.inqFormat(ncid)
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% ivar starts from zero ~
disp('NC file variables')
for ivar = 0:nvars-1
    [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,ivar);
    varname
end
disp('NC file variables')


lon = ncread(file,target_lonName);
lat = ncread(file,target_latName);
level = ncread(file,target_zName);
TIME = ncread(file,target_timeName);
lon(lon<0) = lon(lon<0) + 360;


[c LAT_B] = min(abs(target_lat(1)-lat));
[c LAT_E] = min(abs(target_lat(2)-lat));
[c LON_B] = min(abs(target_lon(1)-lon));
[c LON_E] = min(abs(target_lon(2)-lon));

[c LEVEL] = min(abs(target_z-level));


flip_lat = 0;
if LAT_E < LAT_B
    tmp = LAT_B;
    LAT_B = LAT_E;
    LAT_E = tmp;
    flip_lat = 1;
end
N_LAT = LAT_E-LAT_B+1;
lat = lat(LAT_B:LAT_E);

varID = netcdf.inqVarID(ncid,target_varName);

[varname, xtype, dimids, atts] = netcdf.inqVar(ncid,varID);

for k = 1:length(dimids)
    
    idim = dimids(k);
    [dimname, dimlen] = netcdf.inqDim(ncid,idim);
    if strcmp(dimname,target_timeName)
        dimid_time = k;
    end
    if strcmp(dimname,target_lonName)
        dimid_lon = k;
    end
    if strcmp(dimname,target_latName)
        dimid_lat = k;
    end
    
   if strcmp(dimname,target_zName)
      dimid_z = k;
   end      
    
end

START(dimid_lat) = LAT_B;
COUNT(dimid_lat) = N_LAT;
STRIDE(dimid_lat) = 1;
START(dimid_time) = 1;
START(dimid_z) = LEVEL;
COUNT(dimid_time) = length(TIME);
COUNT(dimid_z) = 1;
STRIDE(dimid_time) = 1;
STRIDE(dimid_lon) = 1;
STRIDE(dimid_z) = 1;
if LON_E < LON_B
    START(dimid_lon) = LON_B;
    COUNT(dimid_lon) = length(lon(LON_B:end));
    data1 = ncread(file,target_varName,START,COUNT,STRIDE);
    START(dimid_lon) = 1;
    COUNT(dimid_lon) = length(lon(1:LON_E));
    data2 = ncread(file,target_varName,START,COUNT,STRIDE);
    output = cat(dimid_lon,data1,data2);
    lon = [lon(LON_B:length(lon));lon(1:LON_E)];
else
    N_LON = LON_E-LON_B+1;
    START(dimid_lon) = LON_B;
    COUNT(dimid_lon) = N_LON;
    STRIDE(dimid_lon) = 1;
    output = ncread(file,target_varName,START,COUNT,STRIDE);
    lon = lon(LON_B:LON_E);
end

if flip_lat == 1
    output = flip(output,dimid_lat);
    lat = flipud(lat);
end

output = double(output);
lat = double(lat);
lon = double(lon);

disp('End of Output from Ncread_DOF')

end


