s_winname_abbr = ["Ham" "Han" "Wel"];

for row=3:3:R
    % one_target
    %fprintf("%s & %f ", s_winname_abbr{sim_grid(row,2)}, sim_grid(row,4));
    % one_target_snr
    %fprintf("%s & %f ", s_winname_abbr{sim_grid(row,2)}, r_SNR(sim_grid(row,3)));
    % masking_v
    %fprintf("%s & %f ", s_winname_abbr{sim_grid(row,2)}, sim_grid(row,4));
    % masking_ampmask
    %fprintf("%s & %f ", s_winname_abbr{sim_grid(row,2)}, 20*log10(r_ampmask(sim_grid(row,5))));
    % separation
    fprintf("%f & %s & %f & %f", r_SNR(sim_grid(row,3)), s_winname_abbr{sim_grid(row,2)}, sim_grid(row,4), sim_grid(row,5));

    for col=1:5
        % frac print
        %suc = b_sucfrac(row,col);
        %if suc < 0.95
        %    fprintf("& \\multicolumn{2}{c}{---}");
        %else
        %    fprintf("& %f & %f", b(row,col), b_err(row,col));
        %    if suc < 0.99
        %        fprintf("**");
        %    elseif suc < 1
        %        fprintf("*");
        %    else
        %    end
        %    fprintf(" ");
        %end

        % percent print
        fprintf("& %f ", sucfrac(row,col));

        % if lambdas should be output
        if col >= 4
            fprintf("& %d", lambdas(1 + (col==5), sim_grid(row,2), sim_grid(row,3)));
        end
        
    end
    fprintf("\\\\\n");
end
