from pydantic import BaseModel, Field, validator, root_validator
from typing import Union



class ReflectionStats(BaseModel):
    # ---- _reflns.* ----
    reflns_entry_id: Union[str, float] = Field("xxxx", alias="_reflns.entry_id")
    reflns_pdbx_diffrn_id: Union[str, float] = Field("1", alias="_reflns.pdbx_diffrn_id")
    reflns_pdbx_ordinal: Union[str, float] = Field("1", alias="_reflns.pdbx_ordinal")
    reflns_d_resolution_low: Union[str, float] = Field("?", alias="_reflns.d_resolution_low")
    reflns_d_resolution_high: Union[str, float] = Field("?", alias="_reflns.d_resolution_high")
    reflns_number_obs: Union[str, float] = Field("?", alias="_reflns.number_obs")
    reflns_percent_possible_obs: Union[str, float] = Field("?", alias="_reflns.percent_possible_obs")
    reflns_pdbx_Rmerge_I_obs: Union[str, float] = Field("?", alias="_reflns.pdbx_Rmerge_I_obs")
    reflns_pdbx_netI_over_sigmaI: Union[str, float] = Field("?", alias="_reflns.pdbx_netI_over_sigmaI")
    reflns_pdbx_redundancy: Union[str, float] = Field("?", alias="_reflns.pdbx_redundancy")
    reflns_pdbx_Rrim_I_all: Union[str, float] = Field("?", alias="_reflns.pdbx_Rrim_I_all")
    reflns_pdbx_Rpim_I_all: Union[str, float] = Field("?", alias="_reflns.pdbx_Rpim_I_all")
    reflns_pdbx_CC_half: Union[str, float] = Field("?", alias="_reflns.pdbx_CC_half")
    reflns_pdbx_number_measured_all: Union[str, float] = Field("?", alias="_reflns.pdbx_number_measured_all")
    reflns_pdbx_chi_squared: Union[str, float] = Field("?", alias="_reflns.pdbx_chi_squared")
    reflns_observed_criterion_sigma_I: Union[str, float] = Field("?", alias="_reflns.observed_criterion_sigma_I")
    reflns_observed_criterion_sigma_F: Union[str, float] = Field("?", alias="_reflns.observed_criterion_sigma_F")
    reflns_number_all: Union[str, float] = Field("?", alias="_reflns.number_all")
    reflns_pdbx_Rsym_value: Union[str, float] = Field("?", alias="_reflns.pdbx_Rsym_value")
    reflns_B_iso_Wilson_estimate: Union[str, float] = Field("?", alias="_reflns.B_iso_Wilson_estimate")
    reflns_pdbx_CC_star: Union[str, float] = Field("?", alias="_reflns.pdbx_CC_star")

    # ---- _reflns_shell.* ----
    reflns_shell_pdbx_diffrn_id: Union[str, float] = Field("1", alias="_reflns_shell.pdbx_diffrn_id")
    reflns_shell_pdbx_ordinal: Union[str, float] = Field("1", alias="_reflns_shell.pdbx_ordinal")
    reflns_shell_d_res_high: Union[str, float] = Field("?", alias="_reflns_shell.d_res_high")
    reflns_shell_d_res_low: Union[str, float] = Field("?", alias="_reflns_shell.d_res_low")
    reflns_shell_number_measured_all: Union[str, float] = Field("?", alias="_reflns_shell.number_measured_all")
    reflns_shell_number_unique_obs: Union[str, float] = Field("?", alias="_reflns_shell.number_unique_obs")
    reflns_shell_Rmerge_I_obs: Union[str, float] = Field("?", alias="_reflns_shell.Rmerge_I_obs")
    reflns_shell_pdbx_chi_squared: Union[str, float] = Field("?", alias="_reflns_shell.pdbx_chi_squared")
    reflns_shell_pdbx_redundancy: Union[str, float] = Field("?", alias="_reflns_shell.pdbx_redundancy")
    reflns_shell_percent_possible_obs: Union[str, float] = Field("?", alias="_reflns_shell.percent_possible_obs")
    reflns_shell_pdbx_netI_over_sigmaI_obs: Union[str, float] = Field("?", alias="_reflns_shell.pdbx_netI_over_sigmaI_obs")
    reflns_shell_pdbx_Rrim_I_all: Union[str, float] = Field("?", alias="_reflns_shell.pdbx_Rrim_I_all")
    reflns_shell_pdbx_Rpim_I_all: Union[str, float] = Field("?", alias="_reflns_shell.pdbx_Rpim_I_all")
    reflns_shell_pdbx_CC_half: Union[str, float] = Field("?", alias="_reflns_shell.pdbx_CC_half")
    reflns_shell_percent_possible_all: Union[str, float] = Field("?", alias="_reflns_shell.percent_possible_all")
    reflns_shell_pdbx_Rsym_value: Union[str, float] = Field("?", alias="_reflns_shell.pdbx_Rsym_value")
    reflns_shell_meanI_over_sigI_obs: Union[str, float] = Field("?", alias="_reflns_shell.meanI_over_sigI_obs")
    reflns_shell_number_measured_obs: Union[str, float] = Field("?", alias="_reflns_shell.number_measured_obs")
    reflns_shell_number_unique_all: Union[str, float] = Field("?", alias="_reflns_shell.number_unique_all")
    reflns_shell_pdbx_CC_star: Union[str, float] = Field("?", alias="_reflns_shell.pdbx_CC_star")

    class Config:
        allow_population_by_field_name = True