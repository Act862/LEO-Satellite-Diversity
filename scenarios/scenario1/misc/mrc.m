function [BER,result] = mrc(g_2,r,sample_num,data)
    %   weighted by fading gain and combine
    r_mrc = real(sum(conj(g_2).*r,3));
    %   Detection
    result = (r_mrc > 0)*2 - 1;
    BER = getBER(result,data,sample_num);
end