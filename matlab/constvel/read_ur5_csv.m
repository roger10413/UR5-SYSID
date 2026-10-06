function [t, q] = read_ur5_csv(filepath)
    % 手動解析表頭 + dlmread 讀取數值，避免依賴 readtable（Octave 需額外套件）
    fid = fopen(filepath, 'r');
    header_line = fgetl(fid);
    fclose(fid);
    col_names = strsplit(strtrim(header_line), ',');

    data = dlmread(filepath, ',', 1, 0);  % 跳過表頭讀數值矩陣

    q = struct();
    for c = 1:numel(col_names)
        name = strtrim(col_names{c});
        name = strrep(name, '-', '_');  % 保險起見去掉不合法字元
        q.(name) = data(:, c);
    end
    t = q.timestamp;
end
