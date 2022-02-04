function []=plotEvents(positions, shift)
    % plot relevant dates
    relevantDates = {
        {'2001/01/02', '2002/01/01','Dot-com bubble'},...
        {'2007/01/01', '2010/01/01','Subprime crisis'},...
        {'2010/01/01', '2012/01/01','European debt crisis'},...
        {'2014/06/01','NW oil crisis'},...
        {'2015/01/01', 'SNB cap'},...
        {'2016/06/01', 'Brexit'},...
        {'2016/11/08', 'Trump'}};
    for i=1:length(relevantDates)
        dateStart = datetime(relevantDates{i}{1}, 'InputFormat', 'yyyy/M/dd');
        if length(relevantDates{i}) == 3
            dateEnd = datetime(relevantDates{i}{2}, 'InputFormat', 'yyyy/M/dd');
            line([dateStart dateEnd], [positions(i)-shift positions(i)-shift], 'Color','k','LineStyle','--');
            text(dateStart + caldays(30), positions(i), relevantDates{i}{3}, 'Fontsize', 8);
        else
            xline(dateStart, '--k');
            text(dateStart + caldays(30), positions(i), relevantDates{i}{2}, 'Fontsize', 8);
        end
    end
end