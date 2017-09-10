% check volume estimates for Stanford talk
% Kurt Feigl 20150124

format compact
echo on

%% average rate from from Ali et al. 2016 (2016), page 117
dVdt=unit(48.2e3,'m^3/year')
convert(dVdt,'m^3/year')
convert(dVdt,'liter/s')
convert(dVdt,'gallon/minute')

%% fn0 = '../intf/In20130513_20140511/drhomaskd_utm.grd';
dV = unit(-23556.7,'m^3')
dt = unit(years(datetime(2014,05,11)-datetime(2013,05,13)),'year')
dVdt=dV/dt;
convert(dVdt,'m^3/year')
convert(dVdt,'liter/s')
convert(dVdt,'gallon/minute')

%%
dV = unit(-0.38e6,'m^3')
convert(dV,'m^3')
convert(dV,'liter')
convert(dV,'gallon')

dt = unit(2014.5 - 2004.5,'year')
dVdt2 = dV/dt
convert(dVdt2,'m^3/year')
convert(dVdt2,'liter/s')
convert(dVdt2,'gallon/minute')

% average model
dVdt3=unit(23107,'m^3/year')
convert(dVdt3,'m^3/year')
convert(dVdt3,'liter/s')
convert(dVdt3,'gallon/minute')

 
% average observation - net in field
dVdt4=unit(-0.23,'m^3/s')
convert(dVdt4,'m^3/year')
convert(dVdt4,'liter/s')
convert(dVdt4,'gallon/minute')

% volume ratio
3645/16