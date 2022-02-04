classdef ExtractData
    % ExtractData Extracts daily and monthly prices and dates from input
    % xlsx file
    %   Inputs:
    %   * inputData and inputText, in form of values which are
    %     the output of the xlsread function over a file downloaded from
    %     EIKON
    %   * year, to keep only data starting from that year
    %   * benchCol, benchmark index column
    
    % outputs
    properties
        rawData
        rawText
        prices              % monthly prices
        dailyPrices         % daily prices
        dailyPricesWithNans
        dailyDates          % daily dates
        dates               % monthly dates
        datesWithNans
        benchDailyPrices    % benchmark daily prices
        benchPrices         % benchmark monthly prices
        labels
        benchLabel
    end
    
    methods
        function obj = ExtractData(inputData, inputText, year, benchCol)
            %ExtractData Construct an instance of this class
            %   Initializes all the properties of the class by
            %   * removing nans
            %   * keeping only common dates
            %   * computing monthly data
            
            % convert dates to proper format
            dateNum = datenum(inputText(4:end,1), 'dd/mm/yyyy');
            obj.rawText = datestr(dateNum);
            obj.rawData = inputData;

            % extract asset labels
            labLen=length(inputText(1,:))-1;
            obj.labels = cell(labLen,1);
            % start from second column
            for i=2:labLen+1
                s = char(inputText(1,i));
                obj.labels{i-1} = s(1:2);
            end
            % move bench label to proper obj param
            obj.benchLabel = obj.labels(benchCol);
            obj.labels(benchCol) = [];
            
            % store data with nans for possible visualization
            obj.datesWithNans = datetime(obj.rawText, 'InputFormat', 'dd-MM-yyyy');
            obj.dailyPricesWithNans = obj.rawData;
            
            % filter data and remove nans
            obj = obj.filterData(year);
            
            % store daily prices
            obj.benchDailyPrices = obj.dailyPrices(:, benchCol);
            obj.dailyPrices(:, benchCol) = [];
            
            % extract monthly prices
            obj = obj.monthlyPrices();
        end
        
        function obj = monthlyPrices(obj)
            % monthly from daily
            % get the nr of unique years+months
            dateStr=obj.dailyDates;
            yearsStr=num2str(year(dateStr));
            monthStr=num2str(month(dateStr));
            uniqueMonths=unique(strcat(yearsStr,monthStr), 'rows');
            % initialize empty array for unique monthly returns and dates
            % exclude last month since it terminates on 07/12/2018
            monthsNr=length(uniqueMonths)-1;
            obj.prices=zeros(monthsNr,size(obj.dailyPrices,2));
            obj.benchPrices=zeros(monthsNr,1);
            obj.dates=NaT(monthsNr,1,'Format','dd-MMM-yyyy');
            % populate
            prevMonth=month(dateStr(1,:));
            j=1;
            for i=1:length(dateStr)
                currMonth=month(dateStr(i,:));
                if currMonth~=prevMonth
                    % month has changed, store previous value
                    obj.prices(j,:)=obj.dailyPrices(i-1,:);
                    obj.benchPrices(j,:)=obj.benchDailyPrices(i-1,:);
                    obj.dates(j,1)=datetime(dateStr(i-1,:), 'InputFormat', 'dd-MM-yyyy');
                    prevMonth=currMonth;
                    j=j+1;
                end
            end
        end
        
        function obj = filterData(obj, yearLimit)
            % filterData Keeps only dates starting from input year and
            % remove nans
            
            % restrict sample to a common specific common period
            dateFilter = year(obj.rawText)>= yearLimit;
            dateStr = obj.rawText(dateFilter, :);
            data = obj.rawData(dateFilter, :);
            % remove nan rows
            nanFilter = ~any(isnan(data), 2);
            dateStr = dateStr(nanFilter, :);
            data = data(nanFilter, :);
            % convert dates to datetime type and store values
            obj.dailyDates = datetime(dateStr, 'InputFormat', 'dd-MM-yyyy');
            obj.dailyPrices = data;
        end
    end
end

