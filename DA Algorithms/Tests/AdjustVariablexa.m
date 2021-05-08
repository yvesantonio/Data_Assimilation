x_a_1_30 = cell(L_max, 30);
x_a_31_50 = cell(L_max, 50);

for tmp0=L_min:L_max
    for tmp1=1:30
        x_a_1_30{tmp0, tmp1} = x_a{tmp0, tmp1};
    end
end

for tmp0=L_min:L_max
    for tmp1=31:50
        x_a_31_50{tmp0, tmp1} = x_a{tmp0, tmp1};
    end
end