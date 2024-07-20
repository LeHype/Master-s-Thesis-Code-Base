function [out] = GetPosition(inputArg1,inputArg2)
try
    out = GetPosition_fromCrane()';
catch
    out = zeros(5,1);
    out(3) = 0.564;
end
end

