

ALL_SHELL_TO_MODEL = {
    "resolutionLimitLow":  "reflns_d_resolution_low",
    "resolutionLimitHigh": "reflns_d_resolution_high",
    "rMerge":              "reflns_pdbx_Rmerge_I_obs",
    "rMeasAllIPlusIMinus": "reflns_pdbx_Rrim_I_all",        # total Rmeas ~ Rrim
    "rPimAllIPlusIMinus":  "reflns_pdbx_Rpim_I_all",
    "meanIOverSigI":       "reflns_pdbx_netI_over_sigmaI",
    "multiplicity":        "reflns_pdbx_redundancy",
    "completeness":        "reflns_percent_possible_obs",
    "nTotalObservations":  "reflns_pdbx_number_measured_all",
    "nTotalUniqueObservations": "reflns_number_obs",
    "ccHalf":              "reflns_pdbx_CC_half",
}

HIGH_SHELL_TO_MODEL = {
    "resolutionLimitLow":  "reflns_shell_d_res_low",
    "resolutionLimitHigh": "reflns_shell_d_res_high",
    "rMerge":              "reflns_shell_Rmerge_I_obs",
    "rMeasAllIPlusIMinus": "reflns_shell_pdbx_Rrim_I_all",        # total Rmeas ~ Rrim
    "rPimAllIPlusIMinus":  "reflns_shell_pdbx_Rpim_I_all",
    "meanIOverSigI":       "reflns_shell_meanI_over_sigmaI_obs",
    "multiplicity":        "reflns_shell_pdbx_redundancy",
    "completeness":        "reflns_shell_percent_possible_obs",
    "nTotalObservations":  "reflns_shell_number_measured_all",
    "nTotalUniqueObservations": "reflns_shell_number_unique_obs",
    "ccHalf":              "reflns_shell_pdbx_CC_half",
}

