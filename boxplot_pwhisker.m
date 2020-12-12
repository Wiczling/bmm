function [ h ] = boxplot_pwhisker( data, varargin )
%BOXPLOT_PWHISKER
% BOXPLOT_PWHISKER(data) plots a boxplot using Statistics Toolbox where
% the whiskers are set using 5th and 95th percentiles. Outliers are turned
% off since they are not typically used for this type of plot.
%
% BOXPLOT_PWHISKER(data, args) accepts arguments, such as {'widths',0.6};
%
% BOXPLOT_PWHISKER(data, args, p1, p2) allows the lower (p1) and upper (p2)
% percentile settings to be specified.
%
%   Examples:
%      data = [1:10; 11:20]';
%      boxplot_pwhisker(data,{'widths',0.6});
%      boxplot_pwhisker(data,{},10,90);

% January 03, 2015
% Wraps code from http://compgroups.net/comp.soft-sys.matlab/5th-95th-percentile-in-boxplot/994208

% Initialize
N = size(data,2); % number of boxes

if nargin == 1      % use default 5th and 95th percentile
    p1 = 5;
    p2 = 95;
    h = boxplot(data);
elseif nargin == 2
    p1 = 5;
    p2 = 95;
    h = boxplot(data,varargin{1}{:});
elseif nargin == 4
    p1 = varargin{2};
    p2 = varargin{3};
    h = boxplot(data,varargin{1}{:});
end

p = prctile(data,[p1 p2]);
set(h(7,:),'Visible','off') % hide outliers

% Replace lower whisker(s)
h1 = flipud(findobj(gca,'Tag','Lower Whisker'));        % the vertical line
h2 = flipud(findobj(gca,'Tag','Lower Adjacent Value')); % the horizontal line
for i=1:N;
    ydata = get(h1(length(h1)-N+i),'YData');
    ydata(1) = p(1,i);
    set(h1(length(h1)-N+i),'YData',ydata);
    ydata(:) = p(1,i);
    set(h2(length(h2)-N+i),'YData',ydata);
end

% Replace upper whisker(s)
h1 = flipud(findobj(gca,'Tag','Upper Whisker'));
h2 = flipud(findobj(gca,'Tag','Upper Adjacent Value'));
for i=1:N;
    ydata = get(h1(length(h1)-N+i),'YData');
    ydata(2) = p(2,i);
    set(h1(length(h1)-N+i),'YData',ydata);
    ydata(:) = p(2,i);
    set(h2(length(h2)-N+i),'YData',ydata);
end


end