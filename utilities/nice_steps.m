function steps = nice_steps(tickvals,nticks)
% Choose nice steps
%   Roots: fix_ticks by Jacco Geul <jacco@geul.net>
% https://www.mathworks.com/matlabcentral/fileexchange/55921-fixticks-automatic-redraw-ticks-on-any-axis
% 20180313 Kurt Feigl 


% limits
lim(1) = min(tickvals);
lim(2) = max(tickvals);

% If it is less than the minimum number or exact number of ticks are
% required change the current Tick
if length(tickvals) < nticks
    % Compute unit and its order-1 then reduce to two significant digits
    u = (lim(2)-lim(1))/(nticks-1);
    o = 10^floor(log10(u/10));
    u = round(u/o)*o;
    
    % Make sure ticks start at a nice multiple of u. However if numticks is
    % even the start must be offset by a half. In the case of uneven number
    % of ticks (recommended) it is ensured that zero is part of the ticks, 
    % if present.
    s = (floor(lim(1)/u)*u + (mod(nticks,2)-1)/2*u); 
    steps = s:u:lim(2);
else
    warning('Too many steps');
    steps = tickvals;
end

end

