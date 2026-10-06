function [pos_file, neg_file] = find_csv(data_folder)
    files = dir(fullfile(data_folder, '*.csv'));
    if isempty(files)
        error('在 %s 找不到任何 CSV 檔，請確認 DATA_FOLDER 設定正確。', data_folder);
    end
    pos_file = '';
    neg_file = '';
    for k = 1:numel(files)
        fname_lower = lower(files(k).name);
        if ~isempty(strfind(fname_lower, 'pos'))
            pos_file = fullfile(files(k).folder, files(k).name);
        elseif ~isempty(strfind(fname_lower, 'neg'))
            neg_file = fullfile(files(k).folder, files(k).name);
        end
    end
    if isempty(pos_file) || isempty(neg_file)
        names = {files.name};
        error('找不到檔名含 "pos" 與 "neg" 的 CSV 各一份。目前資料夾內的 CSV: %s', ...
              strjoin(names, ', '));
    end
end
