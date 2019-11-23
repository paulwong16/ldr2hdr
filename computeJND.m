function output = computeJND(arg)
    if arg <= 127
        output = 17*(1-sqrt((arg/127)))+3;
    else
        output = (3/128)*(arg-127)+3;
    end
end