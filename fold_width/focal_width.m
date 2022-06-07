%% Analysis of focal width

function focal_width = focal_width(x,y)

%  Fit parabola

fit = polyfit(x,y,2);

% Get the focus and the focal distance
x_focus = -fit(2)/(2*fit(1));
y_focus = (1-(fit(2)*fit(2)-4*fit(1)*fit(3)))/(4*fit(1));

syms z
x_comparison = max(double(solve(fit(1)*z^2 + fit(2)*z + fit(3) - y_focus == 0, z)));

focal_width = abs(x_comparison - x_focus);


end