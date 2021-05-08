function [ sel ] = simpleSelector( i, m, radius )
% simpleSelector Select from m measurements a set of measurements distant
% radius from the ith element

if(i-radius < 1)
    sel = zeros(i+radius, m);
    j=1:i+radius;
    if(i+radius > m)
        j=1:m;
    end
    for h=j
        sel(h,h) = 1;
    end
elseif(i+radius > m)
    sel = zeros(radius + 1 + m-i, m);
    j=i-radius:m;
    cur = 1;
    for h=j
        sel(cur,h) = 1;
        cur = cur+1;
    end
else
    sel = zeros(2*radius + 1, m);
    j=i-radius:i+radius;
    cur = 1;
    for h=j
        sel(cur,h) = 1;
        cur = cur + 1;
    end
end

end

