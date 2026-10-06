function [x_cmd, y, x_act] = stack_plateaus(plateaus, direction, kt_out)
    % 把各檔平台期樣本攤平：x_cmd=指令速度（帶號）、x_act=實測 qd、y=tau
    x_cmd = []; y = []; x_act = [];
    for k = 1:numel(plateaus)
        p = plateaus(k);
        x_cmd = [x_cmd; direction * p.level * ones(numel(p.i), 1)]; %#ok<AGROW>
        x_act = [x_act; p.qd]; %#ok<AGROW>
        y = [y; p.i * kt_out]; %#ok<AGROW>
    end
end
