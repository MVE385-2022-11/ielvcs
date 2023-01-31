function lambda_color(vs, v_res)
    ax = gca;

    if length(vs) == 1
        refcols = [1 0 0];
    else
        refcols = [
            27 158 119;
            217 95 2;
            117 112 179;
            231 41 138;
            102 166 30
            ]./255;
        refcols = refcols(1:length(vs),:);
    end

    for k = 2:length(ax.Children)
        v = str2num(ax.Children(k).Tag(3:end));
        ds = abs(v - vs);
        [d, i] = min(ds);

        A = refcols(i,:).*max(1 - d/(3*v_res), 0);
        B = 0.66.*min(d./(3*v_res), 1);

        ax.Children(k).Color = A + B;
    end
end
