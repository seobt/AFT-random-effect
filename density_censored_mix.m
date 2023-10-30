function jpdf = density_censored_mix(x,dij,support,weight,mu,sigma )

    jpdf=0;
    for j=1:length(support)
        jpdf=jpdf+weight(j)*density_censored(x,dij,support(j),mu,sigma);
    end
end
