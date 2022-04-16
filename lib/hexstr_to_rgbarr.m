function arr = hexstr_to_rgbarr(str)
    arr = zeros(3, 1);
    for j = 1:3
        arr(j) = hex2dec(str(2*j-1:2*j)) / 256;
    end
end
