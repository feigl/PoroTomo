% demonstrate time tags in Matlab
% 20160407 Kurt Feigl
% 20160708 Kurt Feigl - use HH instead of hh

%% To specify a default format, type
% datetime.setDefaultFormats('default',fmt)
% where fmt is a character vector composed of the letters A-Z and a-z described for the Format property of datetime arrays, above. For example,
%datetime.setDefaultFormats('default','yyyy-MM-dd_hh:mm:SSSSSS [z = ZZZZ]')
% sets the default datetime format to include a 4-digit year, 2-digit month number, 2-digit day number, and hour, minute, and second values.
% In addition, you can specify a default format for datetimes created without time components. For example,

% datetime.setDefaultFormats('defaultdate','yyyy-MM-dd')
% sets the default date format to include a 4-digit year, 2-digit month number, and 2-digit day number.
%% To reset the both the default format and the default date-only formats to the factory defaults, type
% datetime.setDefaultFormats('reset')

echo on

%% set an example time to PoroTomo teleconference
timedate_teleconference = datetime(2016,04,08,13,0,0)
%timedate_teleconference.Format = 'yyyy-MM-dd_hh:mm:SSSSSS [z = ZZZZ]' ; % 12-hour clock
timedate_teleconference.Format = 'yyyy-MM-dd_HH:mm:SSSSSS [z = ZZZZ]'  ; % 24 hour clock

% set the time zone of the input
timedate_teleconference.TimeZone = 'America/Chicago'
%timedate_teleconference.TimeZone = 'America/Los_Angeles'
% convert to UTC time
timedate_teleconference.TimeZone = 'UTC'


% the following line is NOT necessary
%timedate_teleconference.TimeZone = 'UTCLeapSeconds'

% calculate the epoch when the teleconference will adjourn by adding one
% hour
timedate_adjourn = timedate_teleconference + duration(1,0,0)

% display in different time zones
timedate_adjourn.TimeZone = 'America/Chicago'
%timedate_adjourn.TimeZone = 'America/Los_Angeles'