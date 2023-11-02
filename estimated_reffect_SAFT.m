function [r_effect,va]=estimated_reffect_SAFT(data,beta_SAFT,SAFT_sigmae,r_dist)

support=r_dist(1,:);
weight=r_dist(2,:);
[r_effect,va]=p_reffect(data,support-beta_SAFT(1),weight,beta_SAFT,SAFT_sigmae,@test_linear);

