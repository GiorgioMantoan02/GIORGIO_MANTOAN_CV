function ES = calcolo_ES(w, mu, V, p)
    % formula con rendimenti normali :
    % ES_p (t+1) = mu + sigma * f(phi^(-1) (1-p)) / p
    sigma = sqrt(w' * V * w);
    ES = w' * mu + sigma * normpdf(norminv(1-p)) / p;
end