function dy = Robertson(~,y)

dy = y;

dy(1) = -0.04*y(1) + 1e4*y(2)*y(3);
dy(2) = 0.04*y(1) - 1e4*y(2)*y(3) - 3e7*y(2)^2;
dy(3) = 3e7*y(2)^2;
end