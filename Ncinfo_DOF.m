% output netcdf file info
% Shu Wu dof.shuwu@gmail.com
% file: path of data file

function [] = Ncinfo_DOF(file)

ncid = netcdf.open(file,'nowrite');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);

% ivar starts from zero ~
disp('NC file variables')
for ivar = 0:nvars-1
    ivar
    [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,ivar);
    varname
end

end


