function [reaction_rates] = nfs_propensity_vectorized(reactants, params)
%nfs_propensity_vectorized describes reaction rates for nfs, but for
%                          gillespie_vectorized
%
%   reactants: two dimensional array of current reactant values
%   parameters: struct of parameters containing:
%               k1, k2, k3, K3, k4, K4, A0, B0, I
%
%   reaction_rates: array of 4 reaction rates
                    %A production
                    %A degradation
                    %B production
                    %B degradation
                    
reaction_rates = [params.k1*params.I*(params.A0-reactants(:, 1)), (params.k2*reactants(:, 1).*reactants(:, 2))/params.B0, (params.k3*reactants(:, 1).*(params.B0-reactants(:, 2)))./(params.A0*(params.K3+1-(reactants(:, 2)/params.B0))), (params.k4*reactants(:, 2)*params.B0)./(params.K4*params.B0+reactants(:, 2))];

end

