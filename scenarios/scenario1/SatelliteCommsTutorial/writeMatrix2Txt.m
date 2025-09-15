function writeMatrix2Txt(A,filename)
    fid = fopen(filename,'wt');
    for ii = 1:size(A,1)
        fprintf(fid,'%g\t',A(ii,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
end