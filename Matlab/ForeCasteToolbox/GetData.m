% data = Quandl.get('NSE/OIL');
%
% %mydata = Quandl.get('NSE/OIL', 'start_date','yyyy-mm-dd','end_date','yyyy-mm-dd');
% mydata = Quandl.get('NSE/OIL', 'start_date','2010-01-1','end_date','2015-01-04');
% mydata = Quandl.get('NSE/OIL', 'collapse' ,'annual');% ("weekly"|"monthly"|"quarterly"|"annual")
%
% % Transformations: `
% mydata = Quandl.get('NSE/OIL', 'transformation' ,'rdiff'); % ("diff"|"rdiff"|"normalize"|"cumulative")
%  %* Return only n number of rows: `
%  mydata = Quandl.get('NSE/OIL','rows',5);
%
% %% ## Available Data Types ##
% % There are four options for which datatype you would like your data returned as, you choose your type as follows:
% %
% % 	Quandl.get('NSE/OIL','type','ts')
% %
% % * **Timeseries (default)**: returns a timeseries if only 1 column in data, tscollection if more. `('type','ts')`
% % * **Financial timeseries** :`('type','fints')`
% % * **CSV string**: `('type','ASCII')`
% % * **DataMatrix**: `('type','data')`
% % * **Cell Strings**: `('type','cellstr')`
%
%  mydata = Quandl.get('NSE/OIL','rows',5,'type','fints');
%  %
%     data = Quandl.datatable('ZACKS/FE');
%
%     Quandl.search('crude oil');
%%
clear
clc

formatOut = 'yy.dd';
SD='2009-01-1';
ED=date;%'2017-12-31';
% Opec Oil Price
oil=Quandl.get('OPEC/ORB', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
%oil.TimeInfo.StartDate;
oil=dataset(datenum(getabstime(oil)),oil.Data,'varnames',{'Date','oil'});
% Copper
CU=Quandl.get('LME/PR_CU', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
% CU.TimeInfo.getTimeStr;
CU=dataset(CU.CashBuyer.Data,datenum(getabstime(CU)),'varnames',{'CU','Date'});
% Aluminium
Al= Quandl.get('LME/PR_AA', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);

Al=dataset(Al.CashBuyer.Data,datenum(getabstime(Al)),'varnames',{'Al','Date'});
%  Al.TimeInfo.getTimeStr;
% Coal
Co=Quandl.get('ODA/PCOALAU_USD', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
Co=dataset(Co.Data,datenum(getabstime(Co)),'varnames',{'Co','Date'});
% Co=Quandl.get('EIA/COAL', 'collapse' ,'monthly', 'start_date','2012-01-1','end_date','2016-12-31');
%Co.TimeInfo.getTimeStr ;
% getabstime
% ORE Iron
Or=Quandl.get('ODA/PIORECR_USD', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);

Or=dataset(Or.Data,datenum(getabstime(Or)),'varnames',{'Or','Date'});

% US Energy price index
%Erg=Quandl.get('FRED/USACPIENGMINMEI', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
%Erg2=dataset(Erg.Data,datenum(getabstime(Erg),'varnames',{'Energy','Date'});

% Nasdaaq Index
NsQ=Quandl.get('NASDAQOMX/NQGI', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
NsQ= dataset(datenum(getabstime(NsQ)),NsQ.IndexValue.Data,'varnames',{'Date','NsQ'});


CpiUs= Quandl.get('BLSI/CUSR0000SA0', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
CpiUs=dataset(CpiUs.Data,datenum(getabstime(CpiUs)),'varnames',{'CPI','Date'});

% Commodity index japan banks
CMdi=Quandl.get('ODA/PALLFNF_INDEX', 'authcode', 'RnDoxxE6yUUYF7xVVkEp', 'collapse' ,'monthly', 'start_date',SD,'end_date',ED);
CMdi=dataset(CMdi.Data,datenum(getabstime(CMdi)),'varnames',{'Commodity','Date'});

clear d1

d1=join(NsQ,CpiUs,'Type','outer', 'MergeKeys',true);
d1=join(d1,CMdi,'Type','outer', 'MergeKeys',true);
d1=join(d1,oil,'Type','outer', 'MergeKeys',true);
d1 = join(d1,CU,'Type','outer', 'MergeKeys',true);
d1 = join(d1,Al,'Type','outer', 'MergeKeys',true);
d1 = join(d1,Co,'Type','outer', 'MergeKeys',true);
Data = join(d1,Or,'Type','outer', 'MergeKeys',true);
%Data.DateIndex=datenum(Data.Date);
Data = sortrows(Data,'Date','ascend');
%Data.Date
export(Data,'xlsfile','Input\Mar.xlsx');
% Steel Bilet

%Quandl.get('LME/PR_FM')