function write_level_table(fid, plateaus, n_gear)
    fprintf(fid, '  檔位[rad/s]  樣本數  平均電流[A]  電流std[A]  馬達一圈平均後std[A]  完整圈數\n');
    for k = 1:numel(plateaus)
        p = plateaus(k);
        [~, m] = motor_rev_bins(p, n_gear);
        if numel(m) >= 2
            s2 = sprintf('%.4f', std(m));
        else
            s2 = '  (圈數不足)';
        end
        fprintf(fid, '  %.3f       %5d   %+.4f     %.4f      %10s          %d\n', ...
                p.level, numel(p.i), mean(p.i), std(p.i), s2, numel(m));
    end
end
