function [b, AFC, pass_err, supp_err] = synthesize_filter(pass_band, supp_band, M, plt)

    j = 0:1:M;
    w = pi * (j + 0.5) / (M + 1);

    pass_band_indices  = find(w >= pass_band(1) & w <= pass_band(2));
    trans_band_indices = find(w >  pass_band(2) & w <  supp_band(1));
    supp_band_indices  = find(w >= supp_band(1) & w <= supp_band(2));

    Kd_p = ones(1, length(pass_band_indices));
    Kd_t_0 = 0.5 * ones(1, length(trans_band_indices));
    Kd_s = zeros(1, length(supp_band_indices));
    
    Kd_t_opt = fminsearch(@(Kd_t) find_filter_err(w, M, ...
        Kd_p, Kd_t, Kd_s, pass_band, supp_band), Kd_t_0);
    
    Kd = [Kd_p Kd_t_opt Kd_s];

    b = find_b(Kd, w, M);
    
    AFC = build_AFC(b, M);
    
    pass_err = get_pass_err;
    supp_err = get_supp_err;
    
    if nargin < 4
        plt = 0;
    end
    
    if plt == 1
%         figure, hold on, grid on
%     
%         Kd0 = [Kd_p Kd_t_0 Kd_s];
%         b0 = find_b(Kd0, w, M);
%         AFC0 = build_AFC(b0, M);
%         plot(w, Kd0, 'r*--');
%         fplot(A0, [0, pi]);

        % ----- Линейная шкала ----- %
        figure, hold on, grid on

        % AЧХ
        plot(w, Kd, 'b*');
        fplot(AFC, [-pi, 2*pi], 'b', 'LineWidth', 1);

        % идеальная АЧХ
        IAFC = @(w) 1 - heaviside(w - (pass_band(2) + supp_band(1))/2);
        fplot(IAFC, [-pi, 2*pi], 'm', 'LineWidth', 1);

        title('Оптимальная АЧХ')
        xlabel('w'), ylabel('|K(w)|')
        legend('w_jопт', 'Оптимальная АЧХ', 'Идеальная АЧХ')
        
        line([pass_band(2) pass_band(2)], [-0.1 1.1])
        line([supp_band(1) supp_band(1)], [-0.1 1.1])

        % ----- Логарифмическая шкала ----- %
        figure, hold on, grid on

        w_ = -pi:pi/50:2*pi;
        AFC_ = AFC(w_);
        AFC_(AFC_ < 10^-7) = 10^-7;
        AFC_log = 20*log10(AFC_);

        % AЧХ
        plot(w_, AFC_log, 'b', 'LineWidth', 1);

        title('Оптимальная АЧХ')
        xlabel('w'), ylabel('20lg|K(w)/C|')
    end
    
end

function[b] = find_b(Kd, w, M)
    b = zeros(1, length(w));

    for k = 0:1:M
        S = 0;
        for j = 0:1:M
            S = S + Kd(j+1) * cos(k*w(j+1));
        end
        b(M-k+1) = S / (M + 1);
    end
end

function [abs_AFC, AFC] = build_AFC(b, M)

    syms w
    
    S = 0;
    for k = 1:1:M 
        S = S + b(M-k+1)*cos(w*k);
    end
     
    abs_AFC = matlabFunction(abs(b(M+1) + 2*S)); 
    AFC = matlabFunction(b(M+1) + 2*S);
    
end

function [err] = find_filter_err(w, M, ...
    Kd_p, Kd_t, Kd_s, pass_band, supp_band)

    Kd = [Kd_p Kd_t Kd_s];
    
    b = find_b(Kd, w, M);
    
    [~,AFC] = build_AFC(b, M);

    w_p = pass_band(1):(pass_band(2) - pass_band(1))/100:pass_band(2);
    w_s = supp_band(1):(supp_band(2) - supp_band(1))/100:supp_band(2);
    
    global pass_err supp_err
    
    pass_err = max(abs(AFC(w_p) - 1));
    supp_err = max(abs(AFC(w_s)));

    err = max([pass_err supp_err]);
    
    % disp(err)
end

function res = get_pass_err
    global pass_err
    res = pass_err;
end

function res = get_supp_err
    global supp_err
    res = supp_err;
end
