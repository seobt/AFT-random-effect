function jpdf = density_censored(x,dij,vi,mu,sigma )

    den=normpdf(x,vi+mu,sqrt(sigma)).^dij;
    surv=(1-normcdf(x,vi+mu,sqrt(sigma))).^(1-dij);
    den=max(den,1.0e-20);
    surv=max(surv,1.0e-20);
    jpdf=prod(den.*surv);
end
