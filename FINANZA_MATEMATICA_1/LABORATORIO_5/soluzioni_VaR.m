function w_sol = soluzioni_VaR(mu1, mu2, sigma1, sigma2, cov12, VaR_target, p)
    % VaR_target = mu_p + sigma_p * z
    % VaR_target - mu_p = sigma_p * z

    % Calcola z (quantile inverso della distribuzione normale)
    z = norminv(1 - p);

    % sigma^2_p = (sigma1*w)^2 + ((1-w)*sigma2)^2 + 2*w*(1-w)*cov12
    %           = a_var w^2 + b_var w + c_var

    a_var = sigma1^2 + sigma2^2 - 2 * cov12;
    b_var = 2 * (cov12 - sigma2^2);
    c_var = sigma2^2;

    % moltiplico i parametri per z^2
    a = z^2 * a_var;
    b = z^2 * b_var;
    c = z^2 * c_var;
    % (VaR_target - mu_p)^2 = a w^2 + b w + c

    % mu_p = mu1*w + (1-w)*mu2 = mu2 + w*(mu1-mu2)
    % (VaR_target - mu_p)^2 = 
    % (mu1-mu2)^2 w^2 - 2*(mu1-mu2)*(VaR_target-mu2) w + (VaR_target-mu2)^2
    % differenza di mu e VaR
    delta_mu = mu1 - mu2;
    delta_VaR = VaR_target - mu2;

    % Coefficienti equazione quadratica: a' w^2 + b' w + c' = 0
    a_ = a - delta_mu^2;
    b_ = b + 2 * delta_mu * delta_VaR;
    c_ = c - delta_VaR^2;

    disc = b_^2 - 4 * a_ * c_; % Discriminante

    w_sol = [];

    if disc >= 0
        w_1 = (-b_ + sqrt(disc)) / (2 * a_);
        w_2 = (-b_ - sqrt(disc)) / (2 * a_);
        w_sol = [w_1, w_2];
    end
end
