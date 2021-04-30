function [reaction_rates] = ffs_noisy_propensity_vectorized(reactants, params, I)
%ffs_noisy_propensity_vectorized describes reaction rates for ffs with diff
%                                inputs across sims
%
%   reactants: two dimensional array of current reactant values
%   parameters: struct of parameters containing:
%               k1, k2, k3, K3, k4, A0, B0, I
%
%   reaction_rates: array of 4 reaction rates
                    %A production
                    %A degradation
                    %B production
                    %B degradation

reaction_rates = [params.k1*I.*(params.A0-reactants(:, 1)), (params.k2*reactants(:, 1).*(reactants(:, 2)/params.B0)), (params.k3*I.*((params.B0 - reactants(:, 2))./(params.K3+1-(reactants(:, 2)/params.B0)))), (params.k4*reactants(:, 2))];

end

