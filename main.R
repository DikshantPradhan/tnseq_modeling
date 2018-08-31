
mutans_grb <- GRB_mutans_model(model = 'data/mutans_model.RData')
mutans_r_coupling <- flux_coupling_raptor(mutans_grb)

mutans_rxn_coupling <- mutans_r_coupling$coupled
mutans_flux <- mutans_r_coupling$flux