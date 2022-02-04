%% Computational finance project | Group 5
% *Goal:* Evaluate the performance of a portfolio.
% *Benchmark:* a euro area bond index all maturities
% *Assets:* bond indexes for Euro area countries

%% Initialization and data cleaning
clear
clc

% working with monthly data for now
[inputData, inputText] = xlsread('dati.xlsx', 'daily_euro_flat');
% extract daily and monthly data from xls file
dm = ExtractData(inputData, inputText, 1998, 1);
clear inputData inputText;

% compute inflation
[cpi, cpiDatesRaw] = xlsread('dati.xlsx', 'CPI');
cpiDates = datetime(cpiDatesRaw(4:end,1), 'InputFormat', 'dd/MM/yyyy');
dataFilter = (cpiDates >= dm.dates(1)) & (cpiDates <= dm.dates(end));
cpiFiltered = cpi(dataFilter);
inflation = diff(cpiFiltered)./cpiFiltered(1:end-1);
clear cpi cpiDates dataFilter cpiFiltered cpiDatesRaw;

dates = dm.dates(2:end);

% remove greece
dm.prices(:,12) = [];
dm.labels(12) = [];
dm.dailyPrices(:,12) = [];
dm.dailyPricesWithNans(:,13) = [];

% monthly returns
returns = computeReturns(dm.prices);
benchReturns = computeReturns(dm.benchPrices);

labels = dm.labels;
benchLabel = 'EU';

save data.mat