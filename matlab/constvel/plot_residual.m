function plot_residual(plateaus, direction, B, Tc, data_scale, model_scale, color, label_name)
    for k = 1:numel(plateaus)
        p = plateaus(k);
        pred = (B * direction * p.level + Tc * direction) * model_scale;
        resid = p.i * data_scale - pred;
        if k == 1
            scatter(p.qd, resid, 8, color, 'filled', 'DisplayName', label_name);
        else
            scatter(p.qd, resid, 8, color, 'filled', 'HandleVisibility', 'off');
        end
    end
end
