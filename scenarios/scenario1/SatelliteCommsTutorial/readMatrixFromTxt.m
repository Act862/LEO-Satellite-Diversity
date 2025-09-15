function A = readMatrixFromTxt(filename)
    fid = fopen(filename,'rt');
    if fid == -1
        error('Could not open file %s for reading.', filename);
    end
    
    % Read numeric data from file (tab-delimited)
    data = textscan(fid, '%f', 'Delimiter', '\t', 'CollectOutput', true);
    fclose(fid);

    % Convert column vector into matrix form
    % Count how many columns based on first line
    fid = fopen(filename,'rt');
    firstLine = fgetl(fid);
    fclose(fid);
    
    numCols = numel(str2num(firstLine)); %#ok<ST2NM>
    A = reshape(data{1}, numCols, []).';
end
