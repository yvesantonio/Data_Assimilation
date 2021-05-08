function [ sel ] = roundSelector( i, m, radius )
% simpleSelector Select from m measurements a set of measurements distant
% radius from the ith element

sel = zeros(2*radius + 1, m);

j = 1;
for h=max(i-radius,1):min(i+radius,m)
    sel(j,h) = 1;
    j = j+1;
end

for h=i-radius:0
    sel(j,m+h) = 1;
    j = j+1;
end

for h=m+1:i+radius
    sel(j,h-m) = 1;
    j = j+1;
end
end

