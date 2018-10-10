\ Model Ecoli_core_model
\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
  Biomass_Ecoli_core_w_GAM + I_a_ACKr + I_a_ATPS4r + I_a_CYTBD + I_a_D_LACt2
   + I_a_FBA + I_a_FBP + I_a_FUM + I_a_GLNS + I_a_GLUN + I_a_LDH_D
   + I_a_MALS + I_a_PFK + I_a_PFL + I_a_PGM + I_a_PIt2r + I_a_PTAr
   + I_a_PYK + I_a_RPE + I_a_RPI + I_a_TALA + I_a_ACALD + I_a_ACALDt
   + I_a_ACONTa + I_a_ACONTb + I_a_AKGDH + I_a_ALCD2x + I_a_CO2t
   + I_a_FORt2 + I_a_FORti + I_a_FRUpts2 + I_a_FUMt2_2 + I_a_GLCpts
   + I_a_H2Ot + I_a_MALt2_2 + I_a_NADTRHD + I_a_NH4t + I_a_O2t + I_a_PDH
   + I_a_SUCCt2_2 + I_a_THD2 + I_a_TKT1 + I_a_TKT2
Subject To
 13dpg[c]: GAPD + PGK = 0
 2pg[c]: - ENO - PGM = 0
 3pg[c]: - 1.496 Biomass_Ecoli_core_w_GAM - PGK + PGM = 0
 6pgc[c]: - GND + PGL = 0
 6pgl[c]: G6PDH2r - PGL = 0
 ac[c]: - ACKr + ACt2r = 0
 ac[e]: - ACt2r - EX_ac(e) = 0
 acald[c]: - ACALD + ACALDt + ALCD2x = 0
 acald[e]: - ACALDt - EX_acald(e) = 0
 accoa[c]: ACALD - 3.7478 Biomass_Ecoli_core_w_GAM - CS - MALS + PDH + PFL
   - PTAr = 0
 acon_C[c]: ACONTa - ACONTb = 0
 actp[c]: ACKr + PTAr = 0
 adp[c]: ACKr + 2 ADK1 + ATPM - ATPS4r + 59.81 Biomass_Ecoli_core_w_GAM
   + GLNS + GLNabc + PFK + PGK + PPCK - PYK + SUCOAS = 0
 akg[c]: - AKGDH + AKGt2r + 4.1182 Biomass_Ecoli_core_w_GAM + GLUDy - GLUSy
   + ICDHyr = 0
 akg[e]: - AKGt2r - EX_akg(e) = 0
 amp[c]: - ADK1 + PPS = 0
 atp[c]: - ACKr - ADK1 - ATPM + ATPS4r - 59.81 Biomass_Ecoli_core_w_GAM
   - GLNS - GLNabc - PFK - PGK - PPCK - PPS + PYK - SUCOAS = 0
 cit[c]: - ACONTa + CS = 0
 co2[c]: AKGDH + CO2t + GND + ICDHyr + ME1 + ME2 + PDH - PPC + PPCK = 0
 co2[e]: - CO2t - EX_co2(e) = 0
 coa[c]: - ACALD - AKGDH + 3.7478 Biomass_Ecoli_core_w_GAM + CS + MALS
   - PDH - PFL + PTAr - SUCOAS = 0
 dhap[c]: FBA - TPI = 0
 e4p[c]: - 0.361 Biomass_Ecoli_core_w_GAM + TALA - TKT2 = 0
 etoh[c]: - ALCD2x + ETOHt2r = 0
 etoh[e]: - ETOHt2r - EX_etoh(e) = 0
 f6p[c]: - 0.0709 Biomass_Ecoli_core_w_GAM + FBP + FRUpts2 - PFK + PGI
   + TALA + TKT2 = 0
 fdp[c]: - FBA - FBP + PFK = 0
 for[c]: FORt2 - FORti + PFL = 0
 for[e]: - EX_for(e) - FORt2 + FORti = 0
 fru[e]: - EX_fru(e) - FRUpts2 = 0
 fum[c]: - FRD7 - FUM + FUMt2_2 + SUCDi = 0
 fum[e]: - EX_fum(e) - FUMt2_2 = 0
 g3p[c]: - 0.129 Biomass_Ecoli_core_w_GAM + FBA - GAPD - TALA + TKT1 + TKT2
   + TPI = 0
 g6p[c]: - 0.205 Biomass_Ecoli_core_w_GAM - G6PDH2r + GLCpts - PGI = 0
 glc_D[e]: - EX_glc(e) - GLCpts = 0
 gln_L[c]: - 0.2557 Biomass_Ecoli_core_w_GAM + GLNS + GLNabc - GLUN - GLUSy
   = 0
 gln_L[e]: - EX_gln_L(e) - GLNabc = 0
 glu_L[c]: - 4.9414 Biomass_Ecoli_core_w_GAM - GLNS - GLUDy + GLUN
   + 2 GLUSy + GLUt2r = 0
 glu_L[e]: - EX_glu_L(e) - GLUt2r = 0
 glx[c]: ICL - MALS = 0
 h2o[c]: ACONTa - ACONTb - ATPM + ATPS4r - 59.81 Biomass_Ecoli_core_w_GAM
   - CS + CYTBD + ENO - FBP - FUM - GLNabc - GLUDy - GLUN + H2Ot - MALS
   - PGL - PPC - PPS = 0
 h2o[e]: - EX_h2o(e) - H2Ot = 0
 h[c]: ACALD + ACt2r + AKGt2r + ALCD2x + ATPM + 3 ATPS4r
   + 59.81 Biomass_Ecoli_core_w_GAM + CS - 2 CYTBD + D_LACt2 + ETOHt2r
   + FORt2 + 2 FUMt2_2 + G6PDH2r + GAPD + GLNS + GLNabc + GLUDy - GLUSy
   + GLUt2r + LDH_D + MALS + 2 MALt2_2 + MDH - 4 NADH16 + PFK + PGL + PIt2r
   + PPC + 2 PPS - PYK + PYRt2r + 2 SUCCt2_2 + SUCCt3 + 2 THD2 = 0
 h[e]: - ACt2r - AKGt2r - 4 ATPS4r + 2 CYTBD - D_LACt2 - ETOHt2r - EX_h(e)
   - FORt2 - 2 FUMt2_2 - GLUt2r - 2 MALt2_2 + 3 NADH16 - PIt2r - PYRt2r
   - 2 SUCCt2_2 - SUCCt3 - 2 THD2 = 0
 icit[c]: ACONTb - ICDHyr - ICL = 0
 lac_D[c]: D_LACt2 - LDH_D = 0
 lac_D[e]: - D_LACt2 - EX_lac_D(e) = 0
 mal_L[c]: FUM + MALS + MALt2_2 - MDH - ME1 - ME2 = 0
 mal_L[e]: - EX_mal_L(e) - MALt2_2 = 0
 nad[c]: - ACALD - AKGDH - ALCD2x - 3.547 Biomass_Ecoli_core_w_GAM - GAPD
   - LDH_D - MDH - ME1 + NADH16 - NADTRHD - PDH + THD2 = 0
 nadh[c]: ACALD + AKGDH + ALCD2x + 3.547 Biomass_Ecoli_core_w_GAM + GAPD
   + LDH_D + MDH + ME1 - NADH16 + NADTRHD + PDH - THD2 = 0
 nadp[c]: 13.0279 Biomass_Ecoli_core_w_GAM - G6PDH2r - GLUDy + GLUSy - GND
   - ICDHyr - ME2 + NADTRHD - THD2 = 0
 nadph[c]: - 13.0279 Biomass_Ecoli_core_w_GAM + G6PDH2r + GLUDy - GLUSy
   + GND + ICDHyr + ME2 - NADTRHD + THD2 = 0
 nh4[c]: - GLNS + GLUDy + GLUN + NH4t = 0
 nh4[e]: - EX_nh4(e) - NH4t = 0
 o2[c]: - 0.5 CYTBD + O2t = 0
 o2[e]: - EX_o2(e) - O2t = 0
 oaa[c]: - 1.7867 Biomass_Ecoli_core_w_GAM - CS + MDH + PPC - PPCK = 0
 pep[c]: - 0.5191 Biomass_Ecoli_core_w_GAM + ENO - FRUpts2 - GLCpts - PPC
   + PPCK + PPS - PYK = 0
 pi[c]: ATPM - ATPS4r + 59.81 Biomass_Ecoli_core_w_GAM + FBP - GAPD + GLNS
   + GLNabc + PIt2r + PPC + PPS - PTAr + SUCOAS = 0
 pi[e]: - EX_pi(e) - PIt2r = 0
 pyr[c]: - 2.8328 Biomass_Ecoli_core_w_GAM + FRUpts2 + GLCpts + LDH_D + ME1
   + ME2 - PDH - PFL - PPS + PYK + PYRt2r = 0
 pyr[e]: - EX_pyr(e) - PYRt2r = 0
 q8[c]: CYTBD + FRD7 - NADH16 - SUCDi = 0
 q8h2[c]: - CYTBD - FRD7 + NADH16 + SUCDi = 0
 r5p[c]: - 0.8977 Biomass_Ecoli_core_w_GAM - RPI - TKT1 = 0
 ru5p_D[c]: GND - RPE + RPI = 0
 s7p[c]: - TALA + TKT1 = 0
 succ[c]: FRD7 + ICL + SUCCt2_2 - SUCCt3 - SUCDi - SUCOAS = 0
 succ[e]: - EX_succ(e) - SUCCt2_2 + SUCCt3 = 0
 succoa[c]: AKGDH + SUCOAS = 0
 xu5p_D[c]: RPE - TKT1 - TKT2 = 0
 a_b0351: - Ex_a_b0351 - a_ACALD_fwd_conversion_1
   - a_ACALD_rev_conversion_1 = 0
 a_b1241: - Ex_a_b1241 - a_ACALD_fwd_conversion_2
   - a_ACALD_rev_conversion_2 - a_ALCD2x_fwd_conversion_3
   - a_ALCD2x_rev_conversion_3 = 0
 a_s0001: - Ex_a_s0001 - a_ACALDt_fwd_conversion_1
   - a_ACALDt_rev_conversion_1 - a_CO2t_fwd_conversion_1
   - a_CO2t_rev_conversion_1 - a_H2Ot_fwd_conversion_2
   - a_H2Ot_rev_conversion_2 - a_NH4t_fwd_conversion_1
   - a_NH4t_rev_conversion_1 - a_O2t_fwd_conversion_1
   - a_O2t_rev_conversion_1 = 0
 a_b3115: - Ex_a_b3115 - a_ACKr_fwd_conversion_1 - a_ACKr_rev_conversion_1
   = 0
 a_b2296: - Ex_a_b2296 - a_ACKr_fwd_conversion_2 - a_ACKr_rev_conversion_2
   = 0
 a_b1849: - Ex_a_b1849 - a_ACKr_fwd_conversion_3 - a_ACKr_rev_conversion_3
   = 0
 a_b0118: - Ex_a_b0118 - a_ACONTa_fwd_conversion_1
   - a_ACONTa_rev_conversion_1 - a_ACONTb_fwd_conversion_1
   - a_ACONTb_rev_conversion_1 = 0
 a_b1276: - Ex_a_b1276 - a_ACONTa_fwd_conversion_2
   - a_ACONTa_rev_conversion_2 - a_ACONTb_fwd_conversion_2
   - a_ACONTb_rev_conversion_2 = 0
 a_b0474: - ADK1 - Ex_a_b0474 = 0
 a_b0116: - Ex_a_b0116 - a_AKGDH_fwd_conversion_1
   - a_AKGDH_rev_conversion_1 - a_PDH_fwd_conversion_1
   - a_PDH_rev_conversion_1 = 0
 a_b0726: - Ex_a_b0726 - a_AKGDH_fwd_conversion_1
   - a_AKGDH_rev_conversion_1 = 0
 a_b0727: - Ex_a_b0727 - a_AKGDH_fwd_conversion_1
   - a_AKGDH_rev_conversion_1 = 0
 a_b2587: - AKGt2r - Ex_a_b2587 = 0
 a_b0356: - Ex_a_b0356 - a_ALCD2x_fwd_conversion_1
   - a_ALCD2x_rev_conversion_1 = 0
 a_b1478: - Ex_a_b1478 - a_ALCD2x_fwd_conversion_2
   - a_ALCD2x_rev_conversion_2 = 0
 a_b3736: - Ex_a_b3736 - a_ATPS4r_fwd_conversion_1
   - a_ATPS4r_rev_conversion_1 - a_ATPS4r_fwd_conversion_2
   - a_ATPS4r_rev_conversion_2 = 0
 a_b3737: - Ex_a_b3737 - a_ATPS4r_fwd_conversion_1
   - a_ATPS4r_rev_conversion_1 - a_ATPS4r_fwd_conversion_2
   - a_ATPS4r_rev_conversion_2 = 0
 a_b3738: - Ex_a_b3738 - a_ATPS4r_fwd_conversion_1
   - a_ATPS4r_rev_conversion_1 - a_ATPS4r_fwd_conversion_2
   - a_ATPS4r_rev_conversion_2 = 0
 a_b3731: - Ex_a_b3731 - a_ATPS4r_fwd_conversion_1
   - a_ATPS4r_rev_conversion_1 - a_ATPS4r_fwd_conversion_2
   - a_ATPS4r_rev_conversion_2 = 0
 a_b3732: - Ex_a_b3732 - a_ATPS4r_fwd_conversion_1
   - a_ATPS4r_rev_conversion_1 - a_ATPS4r_fwd_conversion_2
   - a_ATPS4r_rev_conversion_2 = 0
 a_b3733: - Ex_a_b3733 - a_ATPS4r_fwd_conversion_1
   - a_ATPS4r_rev_conversion_1 - a_ATPS4r_fwd_conversion_2
   - a_ATPS4r_rev_conversion_2 = 0
 a_b3734: - Ex_a_b3734 - a_ATPS4r_fwd_conversion_1
   - a_ATPS4r_rev_conversion_1 - a_ATPS4r_fwd_conversion_2
   - a_ATPS4r_rev_conversion_2 = 0
 a_b3735: - Ex_a_b3735 - a_ATPS4r_fwd_conversion_1
   - a_ATPS4r_rev_conversion_1 - a_ATPS4r_fwd_conversion_2
   - a_ATPS4r_rev_conversion_2 = 0
 a_b3739: - Ex_a_b3739 - a_ATPS4r_fwd_conversion_2
   - a_ATPS4r_rev_conversion_2 = 0
 a_b0720: - CS - Ex_a_b0720 = 0
 a_b0978: - Ex_a_b0978 - a_CYTBD_fwd_conversion_1
   - a_CYTBD_rev_conversion_1 = 0
 a_b0979: - Ex_a_b0979 - a_CYTBD_fwd_conversion_1
   - a_CYTBD_rev_conversion_1 = 0
 a_b0733: - Ex_a_b0733 - a_CYTBD_fwd_conversion_2
   - a_CYTBD_rev_conversion_2 = 0
 a_b0734: - Ex_a_b0734 - a_CYTBD_fwd_conversion_2
   - a_CYTBD_rev_conversion_2 = 0
 a_b2975: - Ex_a_b2975 - a_D_LACt2_fwd_conversion_1
   - a_D_LACt2_rev_conversion_1 = 0
 a_b3603: - Ex_a_b3603 - a_D_LACt2_fwd_conversion_2
   - a_D_LACt2_rev_conversion_2 = 0
 a_b2779: - ENO - Ex_a_b2779 = 0
 a_b2097: - Ex_a_b2097 - a_FBA_fwd_conversion_1 - a_FBA_rev_conversion_1
   = 0
 a_b1773: - Ex_a_b1773 - a_FBA_fwd_conversion_2 - a_FBA_rev_conversion_2
   = 0
 a_b2925: - Ex_a_b2925 - a_FBA_fwd_conversion_3 - a_FBA_rev_conversion_3
   = 0
 a_b3925: - Ex_a_b3925 - a_FBP_fwd_conversion_1 - a_FBP_rev_conversion_1
   = 0
 a_b4232: - Ex_a_b4232 - a_FBP_fwd_conversion_2 - a_FBP_rev_conversion_2
   = 0
 a_b0904: - Ex_a_b0904 - a_FORt2_fwd_conversion_1
   - a_FORt2_rev_conversion_1 - a_FORti_fwd_conversion_1
   - a_FORti_rev_conversion_1 = 0
 a_b2492: - Ex_a_b2492 - a_FORt2_fwd_conversion_2
   - a_FORt2_rev_conversion_2 - a_FORti_fwd_conversion_2
   - a_FORti_rev_conversion_2 = 0
 a_b4151: - FRD7 - Ex_a_b4151 = 0
 a_b4152: - FRD7 - Ex_a_b4152 = 0
 a_b4153: - FRD7 - Ex_a_b4153 = 0
 a_b4154: - FRD7 - Ex_a_b4154 = 0
 a_b1817: - Ex_a_b1817 - a_FRUpts2_fwd_conversion_1
   - a_FRUpts2_rev_conversion_1 - a_GLCpts_fwd_conversion_2
   - a_GLCpts_rev_conversion_2 = 0
 a_b1818: - Ex_a_b1818 - a_FRUpts2_fwd_conversion_1
   - a_FRUpts2_rev_conversion_1 - a_GLCpts_fwd_conversion_2
   - a_GLCpts_rev_conversion_2 = 0
 a_b1819: - Ex_a_b1819 - a_FRUpts2_fwd_conversion_1
   - a_FRUpts2_rev_conversion_1 - a_GLCpts_fwd_conversion_2
   - a_GLCpts_rev_conversion_2 = 0
 a_b2415: - Ex_a_b2415 - a_FRUpts2_fwd_conversion_1
   - a_FRUpts2_rev_conversion_1 - a_GLCpts_fwd_conversion_1
   - a_GLCpts_rev_conversion_1 - a_GLCpts_fwd_conversion_2
   - a_GLCpts_rev_conversion_2 - a_GLCpts_fwd_conversion_3
   - a_GLCpts_rev_conversion_3 = 0
 a_b2416: - Ex_a_b2416 - a_FRUpts2_fwd_conversion_1
   - a_FRUpts2_rev_conversion_1 - a_GLCpts_fwd_conversion_1
   - a_GLCpts_rev_conversion_1 - a_GLCpts_fwd_conversion_2
   - a_GLCpts_rev_conversion_2 - a_GLCpts_fwd_conversion_3
   - a_GLCpts_rev_conversion_3 = 0
 a_b1612: - Ex_a_b1612 - a_FUM_fwd_conversion_1 - a_FUM_rev_conversion_1
   = 0
 a_b4122: - Ex_a_b4122 - a_FUM_fwd_conversion_2 - a_FUM_rev_conversion_2
   = 0
 a_b1611: - Ex_a_b1611 - a_FUM_fwd_conversion_3 - a_FUM_rev_conversion_3
   = 0
 a_b3528: - Ex_a_b3528 - a_FUMt2_2_fwd_conversion_1
   - a_FUMt2_2_rev_conversion_1 - a_MALt2_2_fwd_conversion_1
   - a_MALt2_2_rev_conversion_1 - a_SUCCt2_2_fwd_conversion_1
   - a_SUCCt2_2_rev_conversion_1 = 0
 a_b1852: - G6PDH2r - Ex_a_b1852 = 0
 a_b1779: - GAPD - Ex_a_b1779 = 0
 a_b2417: - Ex_a_b2417 - a_GLCpts_fwd_conversion_1
   - a_GLCpts_rev_conversion_1 - a_GLCpts_fwd_conversion_3
   - a_GLCpts_rev_conversion_3 = 0
 a_b1101: - Ex_a_b1101 - a_GLCpts_fwd_conversion_1
   - a_GLCpts_rev_conversion_1 = 0
 a_b1621: - Ex_a_b1621 - a_GLCpts_fwd_conversion_3
   - a_GLCpts_rev_conversion_3 = 0
 a_b3870: - Ex_a_b3870 - a_GLNS_fwd_conversion_1 - a_GLNS_rev_conversion_1
   = 0
 a_b1297: - Ex_a_b1297 - a_GLNS_fwd_conversion_2 - a_GLNS_rev_conversion_2
   = 0
 a_b0811: - GLNabc - Ex_a_b0811 = 0
 a_b0810: - GLNabc - Ex_a_b0810 = 0
 a_b0809: - GLNabc - Ex_a_b0809 = 0
 a_b1761: - GLUDy - Ex_a_b1761 = 0
 a_b1812: - Ex_a_b1812 - a_GLUN_fwd_conversion_1 - a_GLUN_rev_conversion_1
   = 0
 a_b0485: - Ex_a_b0485 - a_GLUN_fwd_conversion_2 - a_GLUN_rev_conversion_2
   = 0
 a_b1524: - Ex_a_b1524 - a_GLUN_fwd_conversion_3 - a_GLUN_rev_conversion_3
   = 0
 a_b3212: - GLUSy - Ex_a_b3212 = 0
 a_b3213: - GLUSy - Ex_a_b3213 = 0
 a_b4077: - GLUt2r - Ex_a_b4077 = 0
 a_b2029: - GND - Ex_a_b2029 = 0
 a_b0875: - Ex_a_b0875 - a_H2Ot_fwd_conversion_1 - a_H2Ot_rev_conversion_1
   = 0
 a_b1136: - ICDHyr - Ex_a_b1136 = 0
 a_b4015: - ICL - Ex_a_b4015 = 0
 a_b2133: - Ex_a_b2133 - a_LDH_D_fwd_conversion_1
   - a_LDH_D_rev_conversion_1 = 0
 a_b1380: - Ex_a_b1380 - a_LDH_D_fwd_conversion_2
   - a_LDH_D_rev_conversion_2 = 0
 a_b4014: - Ex_a_b4014 - a_MALS_fwd_conversion_1 - a_MALS_rev_conversion_1
   = 0
 a_b2976: - Ex_a_b2976 - a_MALS_fwd_conversion_2 - a_MALS_rev_conversion_2
   = 0
 a_b3236: - MDH - Ex_a_b3236 = 0
 a_b1479: - ME1 - Ex_a_b1479 = 0
 a_b2463: - ME2 - Ex_a_b2463 = 0
 a_b2276: - NADH16 - Ex_a_b2276 = 0
 a_b2277: - NADH16 - Ex_a_b2277 = 0
 a_b2278: - NADH16 - Ex_a_b2278 = 0
 a_b2279: - NADH16 - Ex_a_b2279 = 0
 a_b2280: - NADH16 - Ex_a_b2280 = 0
 a_b2281: - NADH16 - Ex_a_b2281 = 0
 a_b2282: - NADH16 - Ex_a_b2282 = 0
 a_b2283: - NADH16 - Ex_a_b2283 = 0
 a_b2284: - NADH16 - Ex_a_b2284 = 0
 a_b2285: - NADH16 - Ex_a_b2285 = 0
 a_b2286: - NADH16 - Ex_a_b2286 = 0
 a_b2287: - NADH16 - Ex_a_b2287 = 0
 a_b2288: - NADH16 - Ex_a_b2288 = 0
 a_b3962: - Ex_a_b3962 - a_NADTRHD_fwd_conversion_1
   - a_NADTRHD_rev_conversion_1 = 0
 a_b1602: - Ex_a_b1602 - a_NADTRHD_fwd_conversion_2
   - a_NADTRHD_rev_conversion_2 - a_THD2_fwd_conversion_1
   - a_THD2_rev_conversion_1 = 0
 a_b1603: - Ex_a_b1603 - a_NADTRHD_fwd_conversion_2
   - a_NADTRHD_rev_conversion_2 - a_THD2_fwd_conversion_1
   - a_THD2_rev_conversion_1 = 0
 a_b0451: - Ex_a_b0451 - a_NH4t_fwd_conversion_2 - a_NH4t_rev_conversion_2
   = 0
 a_b0114: - Ex_a_b0114 - a_PDH_fwd_conversion_1 - a_PDH_rev_conversion_1
   = 0
 a_b0115: - Ex_a_b0115 - a_PDH_fwd_conversion_1 - a_PDH_rev_conversion_1
   = 0
 a_b3916: - Ex_a_b3916 - a_PFK_fwd_conversion_1 - a_PFK_rev_conversion_1
   = 0
 a_b1723: - Ex_a_b1723 - a_PFK_fwd_conversion_2 - a_PFK_rev_conversion_2
   = 0
 a_b0902: - Ex_a_b0902 - a_PFL_fwd_conversion_1 - a_PFL_rev_conversion_1
   - a_PFL_fwd_conversion_2 - a_PFL_rev_conversion_2
   - a_PFL_fwd_conversion_3 - a_PFL_rev_conversion_3 = 0
 a_b0903: - Ex_a_b0903 - a_PFL_fwd_conversion_1 - a_PFL_rev_conversion_1
   - a_PFL_fwd_conversion_2 - a_PFL_rev_conversion_2 = 0
 a_b2579: - Ex_a_b2579 - a_PFL_fwd_conversion_1 - a_PFL_rev_conversion_1
   = 0
 a_b3114: - Ex_a_b3114 - a_PFL_fwd_conversion_3 - a_PFL_rev_conversion_3
   = 0
 a_b3951: - Ex_a_b3951 - a_PFL_fwd_conversion_4 - a_PFL_rev_conversion_4
   = 0
 a_b3952: - Ex_a_b3952 - a_PFL_fwd_conversion_4 - a_PFL_rev_conversion_4
   = 0
 a_b4025: - PGI - Ex_a_b4025 = 0
 a_b2926: - PGK - Ex_a_b2926 = 0
 a_b0767: - PGL - Ex_a_b0767 = 0
 a_b3612: - Ex_a_b3612 - a_PGM_fwd_conversion_1 - a_PGM_rev_conversion_1
   = 0
 a_b4395: - Ex_a_b4395 - a_PGM_fwd_conversion_2 - a_PGM_rev_conversion_2
   = 0
 a_b0755: - Ex_a_b0755 - a_PGM_fwd_conversion_3 - a_PGM_rev_conversion_3
   = 0
 a_b2987: - Ex_a_b2987 - a_PIt2r_fwd_conversion_1
   - a_PIt2r_rev_conversion_1 = 0
 a_b3493: - Ex_a_b3493 - a_PIt2r_fwd_conversion_2
   - a_PIt2r_rev_conversion_2 = 0
 a_b3956: - PPC - Ex_a_b3956 = 0
 a_b3403: - PPCK - Ex_a_b3403 = 0
 a_b1702: - PPS - Ex_a_b1702 = 0
 a_b2297: - Ex_a_b2297 - a_PTAr_fwd_conversion_1 - a_PTAr_rev_conversion_1
   = 0
 a_b2458: - Ex_a_b2458 - a_PTAr_fwd_conversion_2 - a_PTAr_rev_conversion_2
   = 0
 a_b1854: - Ex_a_b1854 - a_PYK_fwd_conversion_1 - a_PYK_rev_conversion_1
   = 0
 a_b1676: - Ex_a_b1676 - a_PYK_fwd_conversion_2 - a_PYK_rev_conversion_2
   = 0
 a_b3386: - Ex_a_b3386 - a_RPE_fwd_conversion_1 - a_RPE_rev_conversion_1
   = 0
 a_b4301: - Ex_a_b4301 - a_RPE_fwd_conversion_2 - a_RPE_rev_conversion_2
   = 0
 a_b2914: - Ex_a_b2914 - a_RPI_fwd_conversion_1 - a_RPI_rev_conversion_1
   = 0
 a_b4090: - Ex_a_b4090 - a_RPI_fwd_conversion_2 - a_RPI_rev_conversion_2
   = 0
 a_b0721: - SUCDi - Ex_a_b0721 = 0
 a_b0722: - SUCDi - Ex_a_b0722 = 0
 a_b0723: - SUCDi - Ex_a_b0723 = 0
 a_b0724: - SUCDi - Ex_a_b0724 = 0
 a_b0728: - SUCOAS - Ex_a_b0728 = 0
 a_b0729: - SUCOAS - Ex_a_b0729 = 0
 a_b2464: - Ex_a_b2464 - a_TALA_fwd_conversion_1 - a_TALA_rev_conversion_1
   = 0
 a_b0008: - Ex_a_b0008 - a_TALA_fwd_conversion_2 - a_TALA_rev_conversion_2
   = 0
 a_b2935: - Ex_a_b2935 - a_TKT1_fwd_conversion_1 - a_TKT1_rev_conversion_1
   - a_TKT2_fwd_conversion_1 - a_TKT2_rev_conversion_1 = 0
 a_b2465: - Ex_a_b2465 - a_TKT1_fwd_conversion_2 - a_TKT1_rev_conversion_2
   - a_TKT2_fwd_conversion_2 - a_TKT2_rev_conversion_2 = 0
 a_b3919: - TPI - Ex_a_b3919 = 0
 a_ACKr: - ACKr + a_ACKr_fwd_conversion_1 - a_ACKr_rev_conversion_1
   + a_ACKr_fwd_conversion_2 - a_ACKr_rev_conversion_2
   + a_ACKr_fwd_conversion_3 - a_ACKr_rev_conversion_3 = 0
 a_ATPS4r: - ATPS4r + a_ATPS4r_fwd_conversion_1 - a_ATPS4r_rev_conversion_1
   + a_ATPS4r_fwd_conversion_2 - a_ATPS4r_rev_conversion_2 = 0
 a_CYTBD: - CYTBD + a_CYTBD_fwd_conversion_1 - a_CYTBD_rev_conversion_1
   + a_CYTBD_fwd_conversion_2 - a_CYTBD_rev_conversion_2 = 0
 a_D_LACt2: - D_LACt2 + a_D_LACt2_fwd_conversion_1
   - a_D_LACt2_rev_conversion_1 + a_D_LACt2_fwd_conversion_2
   - a_D_LACt2_rev_conversion_2 = 0
 a_FBA: - FBA + a_FBA_fwd_conversion_1 - a_FBA_rev_conversion_1
   + a_FBA_fwd_conversion_2 - a_FBA_rev_conversion_2
   + a_FBA_fwd_conversion_3 - a_FBA_rev_conversion_3 = 0
 a_FBP: - FBP + a_FBP_fwd_conversion_1 - a_FBP_rev_conversion_1
   + a_FBP_fwd_conversion_2 - a_FBP_rev_conversion_2 = 0
 a_FUM: - FUM + a_FUM_fwd_conversion_1 - a_FUM_rev_conversion_1
   + a_FUM_fwd_conversion_2 - a_FUM_rev_conversion_2
   + a_FUM_fwd_conversion_3 - a_FUM_rev_conversion_3 = 0
 a_GLNS: - GLNS + a_GLNS_fwd_conversion_1 - a_GLNS_rev_conversion_1
   + a_GLNS_fwd_conversion_2 - a_GLNS_rev_conversion_2 = 0
 a_GLUN: - GLUN + a_GLUN_fwd_conversion_1 - a_GLUN_rev_conversion_1
   + a_GLUN_fwd_conversion_2 - a_GLUN_rev_conversion_2
   + a_GLUN_fwd_conversion_3 - a_GLUN_rev_conversion_3 = 0
 a_LDH_D: - LDH_D + a_LDH_D_fwd_conversion_1 - a_LDH_D_rev_conversion_1
   + a_LDH_D_fwd_conversion_2 - a_LDH_D_rev_conversion_2 = 0
 a_MALS: - MALS + a_MALS_fwd_conversion_1 - a_MALS_rev_conversion_1
   + a_MALS_fwd_conversion_2 - a_MALS_rev_conversion_2 = 0
 a_PFK: - PFK + a_PFK_fwd_conversion_1 - a_PFK_rev_conversion_1
   + a_PFK_fwd_conversion_2 - a_PFK_rev_conversion_2 = 0
 a_PFL: - PFL + a_PFL_fwd_conversion_1 - a_PFL_rev_conversion_1
   + a_PFL_fwd_conversion_2 - a_PFL_rev_conversion_2
   + a_PFL_fwd_conversion_3 - a_PFL_rev_conversion_3
   + a_PFL_fwd_conversion_4 - a_PFL_rev_conversion_4 = 0
 a_PGM: - PGM + a_PGM_fwd_conversion_1 - a_PGM_rev_conversion_1
   + a_PGM_fwd_conversion_2 - a_PGM_rev_conversion_2
   + a_PGM_fwd_conversion_3 - a_PGM_rev_conversion_3 = 0
 a_PIt2r: - PIt2r + a_PIt2r_fwd_conversion_1 - a_PIt2r_rev_conversion_1
   + a_PIt2r_fwd_conversion_2 - a_PIt2r_rev_conversion_2 = 0
 a_PTAr: - PTAr + a_PTAr_fwd_conversion_1 - a_PTAr_rev_conversion_1
   + a_PTAr_fwd_conversion_2 - a_PTAr_rev_conversion_2 = 0
 a_PYK: - PYK + a_PYK_fwd_conversion_1 - a_PYK_rev_conversion_1
   + a_PYK_fwd_conversion_2 - a_PYK_rev_conversion_2 = 0
 a_RPE: - RPE + a_RPE_fwd_conversion_1 - a_RPE_rev_conversion_1
   + a_RPE_fwd_conversion_2 - a_RPE_rev_conversion_2 = 0
 a_RPI: - RPI + a_RPI_fwd_conversion_1 - a_RPI_rev_conversion_1
   + a_RPI_fwd_conversion_2 - a_RPI_rev_conversion_2 = 0
 a_TALA: - TALA + a_TALA_fwd_conversion_1 - a_TALA_rev_conversion_1
   + a_TALA_fwd_conversion_2 - a_TALA_rev_conversion_2 = 0
 a_ACALD: - ACALD + a_ACALD_fwd_conversion_1 - a_ACALD_rev_conversion_1
   + a_ACALD_fwd_conversion_2 - a_ACALD_rev_conversion_2 = 0
 a_ACALDt: - ACALDt + a_ACALDt_fwd_conversion_1 - a_ACALDt_rev_conversion_1
   = 0
 a_ACONTa: - ACONTa + a_ACONTa_fwd_conversion_1 - a_ACONTa_rev_conversion_1
   + a_ACONTa_fwd_conversion_2 - a_ACONTa_rev_conversion_2 = 0
 a_ACONTb: - ACONTb + a_ACONTb_fwd_conversion_1 - a_ACONTb_rev_conversion_1
   + a_ACONTb_fwd_conversion_2 - a_ACONTb_rev_conversion_2 = 0
 a_AKGDH: - AKGDH + a_AKGDH_fwd_conversion_1 - a_AKGDH_rev_conversion_1
   = 0
 a_ALCD2x: - ALCD2x + a_ALCD2x_fwd_conversion_1 - a_ALCD2x_rev_conversion_1
   + a_ALCD2x_fwd_conversion_2 - a_ALCD2x_rev_conversion_2
   + a_ALCD2x_fwd_conversion_3 - a_ALCD2x_rev_conversion_3 = 0
 a_CO2t: - CO2t + a_CO2t_fwd_conversion_1 - a_CO2t_rev_conversion_1 = 0
 a_FORt2: - FORt2 + a_FORt2_fwd_conversion_1 - a_FORt2_rev_conversion_1
   + a_FORt2_fwd_conversion_2 - a_FORt2_rev_conversion_2 = 0
 a_FORti: - FORti + a_FORti_fwd_conversion_1 - a_FORti_rev_conversion_1
   + a_FORti_fwd_conversion_2 - a_FORti_rev_conversion_2 = 0
 a_FRUpts2: - FRUpts2 + a_FRUpts2_fwd_conversion_1
   - a_FRUpts2_rev_conversion_1 = 0
 a_FUMt2_2: - FUMt2_2 + a_FUMt2_2_fwd_conversion_1
   - a_FUMt2_2_rev_conversion_1 = 0
 a_GLCpts: - GLCpts + a_GLCpts_fwd_conversion_1 - a_GLCpts_rev_conversion_1
   + a_GLCpts_fwd_conversion_2 - a_GLCpts_rev_conversion_2
   + a_GLCpts_fwd_conversion_3 - a_GLCpts_rev_conversion_3 = 0
 a_H2Ot: - H2Ot + a_H2Ot_fwd_conversion_1 - a_H2Ot_rev_conversion_1
   + a_H2Ot_fwd_conversion_2 - a_H2Ot_rev_conversion_2 = 0
 a_MALt2_2: - MALt2_2 + a_MALt2_2_fwd_conversion_1
   - a_MALt2_2_rev_conversion_1 = 0
 a_NADTRHD: - NADTRHD + a_NADTRHD_fwd_conversion_1
   - a_NADTRHD_rev_conversion_1 + a_NADTRHD_fwd_conversion_2
   - a_NADTRHD_rev_conversion_2 = 0
 a_NH4t: - NH4t + a_NH4t_fwd_conversion_1 - a_NH4t_rev_conversion_1
   + a_NH4t_fwd_conversion_2 - a_NH4t_rev_conversion_2 = 0
 a_O2t: - O2t + a_O2t_fwd_conversion_1 - a_O2t_rev_conversion_1 = 0
 a_PDH: - PDH + a_PDH_fwd_conversion_1 - a_PDH_rev_conversion_1 = 0
 a_SUCCt2_2: - SUCCt2_2 + a_SUCCt2_2_fwd_conversion_1
   - a_SUCCt2_2_rev_conversion_1 = 0
 a_THD2: - THD2 + a_THD2_fwd_conversion_1 - a_THD2_rev_conversion_1 = 0
 a_TKT1: - TKT1 + a_TKT1_fwd_conversion_1 - a_TKT1_rev_conversion_1
   + a_TKT1_fwd_conversion_2 - a_TKT1_rev_conversion_2 = 0
 a_TKT2: - TKT2 + a_TKT2_fwd_conversion_1 - a_TKT2_rev_conversion_1
   + a_TKT2_fwd_conversion_2 - a_TKT2_rev_conversion_2 = 0
 a_ACKr_fwd_conversion_1_I: - a_ACKr_fwd_conversion_1 + 1000 I_a_ACKr >= 0
 a_ACKr_rev_conversion_1_I: - a_ACKr_rev_conversion_1 - 1000 I_a_ACKr
   >= -1000
 a_ACKr_fwd_conversion_2_I: - a_ACKr_fwd_conversion_2 + 1000 I_a_ACKr >= 0
 a_ACKr_rev_conversion_2_I: - a_ACKr_rev_conversion_2 - 1000 I_a_ACKr
   >= -1000
 a_ACKr_fwd_conversion_3_I: - a_ACKr_fwd_conversion_3 + 1000 I_a_ACKr >= 0
 a_ACKr_rev_conversion_3_I: - a_ACKr_rev_conversion_3 - 1000 I_a_ACKr
   >= -1000
 a_ATPS4r_fwd_conversion_1_I: - a_ATPS4r_fwd_conversion_1 + 1000 I_a_ATPS4r
   >= 0
 a_ATPS4r_rev_conversion_1_I: - a_ATPS4r_rev_conversion_1 - 1000 I_a_ATPS4r
   >= -1000
 a_ATPS4r_fwd_conversion_2_I: - a_ATPS4r_fwd_conversion_2 + 1000 I_a_ATPS4r
   >= 0
 a_ATPS4r_rev_conversion_2_I: - a_ATPS4r_rev_conversion_2 - 1000 I_a_ATPS4r
   >= -1000
 a_CYTBD_fwd_conversion_1_I: - a_CYTBD_fwd_conversion_1 + 1000 I_a_CYTBD
   >= 0
 a_CYTBD_rev_conversion_1_I: - a_CYTBD_rev_conversion_1 - 1000 I_a_CYTBD
   >= -1000
 a_CYTBD_fwd_conversion_2_I: - a_CYTBD_fwd_conversion_2 + 1000 I_a_CYTBD
   >= 0
 a_CYTBD_rev_conversion_2_I: - a_CYTBD_rev_conversion_2 - 1000 I_a_CYTBD
   >= -1000
 a_D_LACt2_fwd_conversion_1_I: - a_D_LACt2_fwd_conversion_1
   + 1000 I_a_D_LACt2 >= 0
 a_D_LACt2_rev_conversion_1_I: - a_D_LACt2_rev_conversion_1
   - 1000 I_a_D_LACt2 >= -1000
 a_D_LACt2_fwd_conversion_2_I: - a_D_LACt2_fwd_conversion_2
   + 1000 I_a_D_LACt2 >= 0
 a_D_LACt2_rev_conversion_2_I: - a_D_LACt2_rev_conversion_2
   - 1000 I_a_D_LACt2 >= -1000
 a_FBA_fwd_conversion_1_I: - a_FBA_fwd_conversion_1 + 1000 I_a_FBA >= 0
 a_FBA_rev_conversion_1_I: - a_FBA_rev_conversion_1 - 1000 I_a_FBA
   >= -1000
 a_FBA_fwd_conversion_2_I: - a_FBA_fwd_conversion_2 + 1000 I_a_FBA >= 0
 a_FBA_rev_conversion_2_I: - a_FBA_rev_conversion_2 - 1000 I_a_FBA
   >= -1000
 a_FBA_fwd_conversion_3_I: - a_FBA_fwd_conversion_3 + 1000 I_a_FBA >= 0
 a_FBA_rev_conversion_3_I: - a_FBA_rev_conversion_3 - 1000 I_a_FBA
   >= -1000
 a_FBP_fwd_conversion_1_I: - a_FBP_fwd_conversion_1 + 1000 I_a_FBP >= 0
 a_FBP_rev_conversion_1_I: - a_FBP_rev_conversion_1 - 1000 I_a_FBP
   >= -1000
 a_FBP_fwd_conversion_2_I: - a_FBP_fwd_conversion_2 + 1000 I_a_FBP >= 0
 a_FBP_rev_conversion_2_I: - a_FBP_rev_conversion_2 - 1000 I_a_FBP
   >= -1000
 a_FUM_fwd_conversion_1_I: - a_FUM_fwd_conversion_1 + 1000 I_a_FUM >= 0
 a_FUM_rev_conversion_1_I: - a_FUM_rev_conversion_1 - 1000 I_a_FUM
   >= -1000
 a_FUM_fwd_conversion_2_I: - a_FUM_fwd_conversion_2 + 1000 I_a_FUM >= 0
 a_FUM_rev_conversion_2_I: - a_FUM_rev_conversion_2 - 1000 I_a_FUM
   >= -1000
 a_FUM_fwd_conversion_3_I: - a_FUM_fwd_conversion_3 + 1000 I_a_FUM >= 0
 a_FUM_rev_conversion_3_I: - a_FUM_rev_conversion_3 - 1000 I_a_FUM
   >= -1000
 a_GLNS_fwd_conversion_1_I: - a_GLNS_fwd_conversion_1 + 1000 I_a_GLNS >= 0
 a_GLNS_rev_conversion_1_I: - a_GLNS_rev_conversion_1 - 1000 I_a_GLNS
   >= -1000
 a_GLNS_fwd_conversion_2_I: - a_GLNS_fwd_conversion_2 + 1000 I_a_GLNS >= 0
 a_GLNS_rev_conversion_2_I: - a_GLNS_rev_conversion_2 - 1000 I_a_GLNS
   >= -1000
 a_GLUN_fwd_conversion_1_I: - a_GLUN_fwd_conversion_1 + 1000 I_a_GLUN >= 0
 a_GLUN_rev_conversion_1_I: - a_GLUN_rev_conversion_1 - 1000 I_a_GLUN
   >= -1000
 a_GLUN_fwd_conversion_2_I: - a_GLUN_fwd_conversion_2 + 1000 I_a_GLUN >= 0
 a_GLUN_rev_conversion_2_I: - a_GLUN_rev_conversion_2 - 1000 I_a_GLUN
   >= -1000
 a_GLUN_fwd_conversion_3_I: - a_GLUN_fwd_conversion_3 + 1000 I_a_GLUN >= 0
 a_GLUN_rev_conversion_3_I: - a_GLUN_rev_conversion_3 - 1000 I_a_GLUN
   >= -1000
 a_LDH_D_fwd_conversion_1_I: - a_LDH_D_fwd_conversion_1 + 1000 I_a_LDH_D
   >= 0
 a_LDH_D_rev_conversion_1_I: - a_LDH_D_rev_conversion_1 - 1000 I_a_LDH_D
   >= -1000
 a_LDH_D_fwd_conversion_2_I: - a_LDH_D_fwd_conversion_2 + 1000 I_a_LDH_D
   >= 0
 a_LDH_D_rev_conversion_2_I: - a_LDH_D_rev_conversion_2 - 1000 I_a_LDH_D
   >= -1000
 a_MALS_fwd_conversion_1_I: - a_MALS_fwd_conversion_1 + 1000 I_a_MALS >= 0
 a_MALS_rev_conversion_1_I: - a_MALS_rev_conversion_1 - 1000 I_a_MALS
   >= -1000
 a_MALS_fwd_conversion_2_I: - a_MALS_fwd_conversion_2 + 1000 I_a_MALS >= 0
 a_MALS_rev_conversion_2_I: - a_MALS_rev_conversion_2 - 1000 I_a_MALS
   >= -1000
 a_PFK_fwd_conversion_1_I: - a_PFK_fwd_conversion_1 + 1000 I_a_PFK >= 0
 a_PFK_rev_conversion_1_I: - a_PFK_rev_conversion_1 - 1000 I_a_PFK
   >= -1000
 a_PFK_fwd_conversion_2_I: - a_PFK_fwd_conversion_2 + 1000 I_a_PFK >= 0
 a_PFK_rev_conversion_2_I: - a_PFK_rev_conversion_2 - 1000 I_a_PFK
   >= -1000
 a_PFL_fwd_conversion_1_I: - a_PFL_fwd_conversion_1 + 1000 I_a_PFL >= 0
 a_PFL_rev_conversion_1_I: - a_PFL_rev_conversion_1 - 1000 I_a_PFL
   >= -1000
 a_PFL_fwd_conversion_2_I: - a_PFL_fwd_conversion_2 + 1000 I_a_PFL >= 0
 a_PFL_rev_conversion_2_I: - a_PFL_rev_conversion_2 - 1000 I_a_PFL
   >= -1000
 a_PFL_fwd_conversion_3_I: - a_PFL_fwd_conversion_3 + 1000 I_a_PFL >= 0
 a_PFL_rev_conversion_3_I: - a_PFL_rev_conversion_3 - 1000 I_a_PFL
   >= -1000
 a_PFL_fwd_conversion_4_I: - a_PFL_fwd_conversion_4 + 1000 I_a_PFL >= 0
 a_PFL_rev_conversion_4_I: - a_PFL_rev_conversion_4 - 1000 I_a_PFL
   >= -1000
 a_PGM_fwd_conversion_1_I: - a_PGM_fwd_conversion_1 + 1000 I_a_PGM >= 0
 a_PGM_rev_conversion_1_I: - a_PGM_rev_conversion_1 - 1000 I_a_PGM
   >= -1000
 a_PGM_fwd_conversion_2_I: - a_PGM_fwd_conversion_2 + 1000 I_a_PGM >= 0
 a_PGM_rev_conversion_2_I: - a_PGM_rev_conversion_2 - 1000 I_a_PGM
   >= -1000
 a_PGM_fwd_conversion_3_I: - a_PGM_fwd_conversion_3 + 1000 I_a_PGM >= 0
 a_PGM_rev_conversion_3_I: - a_PGM_rev_conversion_3 - 1000 I_a_PGM
   >= -1000
 a_PIt2r_fwd_conversion_1_I: - a_PIt2r_fwd_conversion_1 + 1000 I_a_PIt2r
   >= 0
 a_PIt2r_rev_conversion_1_I: - a_PIt2r_rev_conversion_1 - 1000 I_a_PIt2r
   >= -1000
 a_PIt2r_fwd_conversion_2_I: - a_PIt2r_fwd_conversion_2 + 1000 I_a_PIt2r
   >= 0
 a_PIt2r_rev_conversion_2_I: - a_PIt2r_rev_conversion_2 - 1000 I_a_PIt2r
   >= -1000
 a_PTAr_fwd_conversion_1_I: - a_PTAr_fwd_conversion_1 + 1000 I_a_PTAr >= 0
 a_PTAr_rev_conversion_1_I: - a_PTAr_rev_conversion_1 - 1000 I_a_PTAr
   >= -1000
 a_PTAr_fwd_conversion_2_I: - a_PTAr_fwd_conversion_2 + 1000 I_a_PTAr >= 0
 a_PTAr_rev_conversion_2_I: - a_PTAr_rev_conversion_2 - 1000 I_a_PTAr
   >= -1000
 a_PYK_fwd_conversion_1_I: - a_PYK_fwd_conversion_1 + 1000 I_a_PYK >= 0
 a_PYK_rev_conversion_1_I: - a_PYK_rev_conversion_1 - 1000 I_a_PYK
   >= -1000
 a_PYK_fwd_conversion_2_I: - a_PYK_fwd_conversion_2 + 1000 I_a_PYK >= 0
 a_PYK_rev_conversion_2_I: - a_PYK_rev_conversion_2 - 1000 I_a_PYK
   >= -1000
 a_RPE_fwd_conversion_1_I: - a_RPE_fwd_conversion_1 + 1000 I_a_RPE >= 0
 a_RPE_rev_conversion_1_I: - a_RPE_rev_conversion_1 - 1000 I_a_RPE
   >= -1000
 a_RPE_fwd_conversion_2_I: - a_RPE_fwd_conversion_2 + 1000 I_a_RPE >= 0
 a_RPE_rev_conversion_2_I: - a_RPE_rev_conversion_2 - 1000 I_a_RPE
   >= -1000
 a_RPI_fwd_conversion_1_I: - a_RPI_fwd_conversion_1 + 1000 I_a_RPI >= 0
 a_RPI_rev_conversion_1_I: - a_RPI_rev_conversion_1 - 1000 I_a_RPI
   >= -1000
 a_RPI_fwd_conversion_2_I: - a_RPI_fwd_conversion_2 + 1000 I_a_RPI >= 0
 a_RPI_rev_conversion_2_I: - a_RPI_rev_conversion_2 - 1000 I_a_RPI
   >= -1000
 a_TALA_fwd_conversion_1_I: - a_TALA_fwd_conversion_1 + 1000 I_a_TALA >= 0
 a_TALA_rev_conversion_1_I: - a_TALA_rev_conversion_1 - 1000 I_a_TALA
   >= -1000
 a_TALA_fwd_conversion_2_I: - a_TALA_fwd_conversion_2 + 1000 I_a_TALA >= 0
 a_TALA_rev_conversion_2_I: - a_TALA_rev_conversion_2 - 1000 I_a_TALA
   >= -1000
 a_ACALD_fwd_conversion_1_I: - a_ACALD_fwd_conversion_1 + 1000 I_a_ACALD
   >= 0
 a_ACALD_rev_conversion_1_I: - a_ACALD_rev_conversion_1 - 1000 I_a_ACALD
   >= -1000
 a_ACALD_fwd_conversion_2_I: - a_ACALD_fwd_conversion_2 + 1000 I_a_ACALD
   >= 0
 a_ACALD_rev_conversion_2_I: - a_ACALD_rev_conversion_2 - 1000 I_a_ACALD
   >= -1000
 a_ACALDt_fwd_conversion_1_I: - a_ACALDt_fwd_conversion_1 + 1000 I_a_ACALDt
   >= 0
 a_ACALDt_rev_conversion_1_I: - a_ACALDt_rev_conversion_1 - 1000 I_a_ACALDt
   >= -1000
 a_ACONTa_fwd_conversion_1_I: - a_ACONTa_fwd_conversion_1 + 1000 I_a_ACONTa
   >= 0
 a_ACONTa_rev_conversion_1_I: - a_ACONTa_rev_conversion_1 - 1000 I_a_ACONTa
   >= -1000
 a_ACONTa_fwd_conversion_2_I: - a_ACONTa_fwd_conversion_2 + 1000 I_a_ACONTa
   >= 0
 a_ACONTa_rev_conversion_2_I: - a_ACONTa_rev_conversion_2 - 1000 I_a_ACONTa
   >= -1000
 a_ACONTb_fwd_conversion_1_I: - a_ACONTb_fwd_conversion_1 + 1000 I_a_ACONTb
   >= 0
 a_ACONTb_rev_conversion_1_I: - a_ACONTb_rev_conversion_1 - 1000 I_a_ACONTb
   >= -1000
 a_ACONTb_fwd_conversion_2_I: - a_ACONTb_fwd_conversion_2 + 1000 I_a_ACONTb
   >= 0
 a_ACONTb_rev_conversion_2_I: - a_ACONTb_rev_conversion_2 - 1000 I_a_ACONTb
   >= -1000
 a_AKGDH_fwd_conversion_1_I: - a_AKGDH_fwd_conversion_1 + 1000 I_a_AKGDH
   >= 0
 a_AKGDH_rev_conversion_1_I: - a_AKGDH_rev_conversion_1 - 1000 I_a_AKGDH
   >= -1000
 a_ALCD2x_fwd_conversion_1_I: - a_ALCD2x_fwd_conversion_1 + 1000 I_a_ALCD2x
   >= 0
 a_ALCD2x_rev_conversion_1_I: - a_ALCD2x_rev_conversion_1 - 1000 I_a_ALCD2x
   >= -1000
 a_ALCD2x_fwd_conversion_2_I: - a_ALCD2x_fwd_conversion_2 + 1000 I_a_ALCD2x
   >= 0
 a_ALCD2x_rev_conversion_2_I: - a_ALCD2x_rev_conversion_2 - 1000 I_a_ALCD2x
   >= -1000
 a_ALCD2x_fwd_conversion_3_I: - a_ALCD2x_fwd_conversion_3 + 1000 I_a_ALCD2x
   >= 0
 a_ALCD2x_rev_conversion_3_I: - a_ALCD2x_rev_conversion_3 - 1000 I_a_ALCD2x
   >= -1000
 a_CO2t_fwd_conversion_1_I: - a_CO2t_fwd_conversion_1 + 1000 I_a_CO2t >= 0
 a_CO2t_rev_conversion_1_I: - a_CO2t_rev_conversion_1 - 1000 I_a_CO2t
   >= -1000
 a_FORt2_fwd_conversion_1_I: - a_FORt2_fwd_conversion_1 + 1000 I_a_FORt2
   >= 0
 a_FORt2_rev_conversion_1_I: - a_FORt2_rev_conversion_1 - 1000 I_a_FORt2
   >= -1000
 a_FORt2_fwd_conversion_2_I: - a_FORt2_fwd_conversion_2 + 1000 I_a_FORt2
   >= 0
 a_FORt2_rev_conversion_2_I: - a_FORt2_rev_conversion_2 - 1000 I_a_FORt2
   >= -1000
 a_FORti_fwd_conversion_1_I: - a_FORti_fwd_conversion_1 + 1000 I_a_FORti
   >= 0
 a_FORti_rev_conversion_1_I: - a_FORti_rev_conversion_1 - 1000 I_a_FORti
   >= -1000
 a_FORti_fwd_conversion_2_I: - a_FORti_fwd_conversion_2 + 1000 I_a_FORti
   >= 0
 a_FORti_rev_conversion_2_I: - a_FORti_rev_conversion_2 - 1000 I_a_FORti
   >= -1000
 a_FRUpts2_fwd_conversion_1_I: - a_FRUpts2_fwd_conversion_1
   + 1000 I_a_FRUpts2 >= 0
 a_FRUpts2_rev_conversion_1_I: - a_FRUpts2_rev_conversion_1
   - 1000 I_a_FRUpts2 >= -1000
 a_FUMt2_2_fwd_conversion_1_I: - a_FUMt2_2_fwd_conversion_1
   + 1000 I_a_FUMt2_2 >= 0
 a_FUMt2_2_rev_conversion_1_I: - a_FUMt2_2_rev_conversion_1
   - 1000 I_a_FUMt2_2 >= -1000
 a_GLCpts_fwd_conversion_1_I: - a_GLCpts_fwd_conversion_1 + 1000 I_a_GLCpts
   >= 0
 a_GLCpts_rev_conversion_1_I: - a_GLCpts_rev_conversion_1 - 1000 I_a_GLCpts
   >= -1000
 a_GLCpts_fwd_conversion_2_I: - a_GLCpts_fwd_conversion_2 + 1000 I_a_GLCpts
   >= 0
 a_GLCpts_rev_conversion_2_I: - a_GLCpts_rev_conversion_2 - 1000 I_a_GLCpts
   >= -1000
 a_GLCpts_fwd_conversion_3_I: - a_GLCpts_fwd_conversion_3 + 1000 I_a_GLCpts
   >= 0
 a_GLCpts_rev_conversion_3_I: - a_GLCpts_rev_conversion_3 - 1000 I_a_GLCpts
   >= -1000
 a_H2Ot_fwd_conversion_1_I: - a_H2Ot_fwd_conversion_1 + 1000 I_a_H2Ot >= 0
 a_H2Ot_rev_conversion_1_I: - a_H2Ot_rev_conversion_1 - 1000 I_a_H2Ot
   >= -1000
 a_H2Ot_fwd_conversion_2_I: - a_H2Ot_fwd_conversion_2 + 1000 I_a_H2Ot >= 0
 a_H2Ot_rev_conversion_2_I: - a_H2Ot_rev_conversion_2 - 1000 I_a_H2Ot
   >= -1000
 a_MALt2_2_fwd_conversion_1_I: - a_MALt2_2_fwd_conversion_1
   + 1000 I_a_MALt2_2 >= 0
 a_MALt2_2_rev_conversion_1_I: - a_MALt2_2_rev_conversion_1
   - 1000 I_a_MALt2_2 >= -1000
 a_NADTRHD_fwd_conversion_1_I: - a_NADTRHD_fwd_conversion_1
   + 1000 I_a_NADTRHD >= 0
 a_NADTRHD_rev_conversion_1_I: - a_NADTRHD_rev_conversion_1
   - 1000 I_a_NADTRHD >= -1000
 a_NADTRHD_fwd_conversion_2_I: - a_NADTRHD_fwd_conversion_2
   + 1000 I_a_NADTRHD >= 0
 a_NADTRHD_rev_conversion_2_I: - a_NADTRHD_rev_conversion_2
   - 1000 I_a_NADTRHD >= -1000
 a_NH4t_fwd_conversion_1_I: - a_NH4t_fwd_conversion_1 + 1000 I_a_NH4t >= 0
 a_NH4t_rev_conversion_1_I: - a_NH4t_rev_conversion_1 - 1000 I_a_NH4t
   >= -1000
 a_NH4t_fwd_conversion_2_I: - a_NH4t_fwd_conversion_2 + 1000 I_a_NH4t >= 0
 a_NH4t_rev_conversion_2_I: - a_NH4t_rev_conversion_2 - 1000 I_a_NH4t
   >= -1000
 a_O2t_fwd_conversion_1_I: - a_O2t_fwd_conversion_1 + 1000 I_a_O2t >= 0
 a_O2t_rev_conversion_1_I: - a_O2t_rev_conversion_1 - 1000 I_a_O2t
   >= -1000
 a_PDH_fwd_conversion_1_I: - a_PDH_fwd_conversion_1 + 1000 I_a_PDH >= 0
 a_PDH_rev_conversion_1_I: - a_PDH_rev_conversion_1 - 1000 I_a_PDH
   >= -1000
 a_SUCCt2_2_fwd_conversion_1_I: - a_SUCCt2_2_fwd_conversion_1
   + 1000 I_a_SUCCt2_2 >= 0
 a_SUCCt2_2_rev_conversion_1_I: - a_SUCCt2_2_rev_conversion_1
   - 1000 I_a_SUCCt2_2 >= -1000
 a_THD2_fwd_conversion_1_I: - a_THD2_fwd_conversion_1 + 1000 I_a_THD2 >= 0
 a_THD2_rev_conversion_1_I: - a_THD2_rev_conversion_1 - 1000 I_a_THD2
   >= -1000
 a_TKT1_fwd_conversion_1_I: - a_TKT1_fwd_conversion_1 + 1000 I_a_TKT1 >= 0
 a_TKT1_rev_conversion_1_I: - a_TKT1_rev_conversion_1 - 1000 I_a_TKT1
   >= -1000
 a_TKT1_fwd_conversion_2_I: - a_TKT1_fwd_conversion_2 + 1000 I_a_TKT1 >= 0
 a_TKT1_rev_conversion_2_I: - a_TKT1_rev_conversion_2 - 1000 I_a_TKT1
   >= -1000
 a_TKT2_fwd_conversion_1_I: - a_TKT2_fwd_conversion_1 + 1000 I_a_TKT2 >= 0
 a_TKT2_rev_conversion_1_I: - a_TKT2_rev_conversion_1 - 1000 I_a_TKT2
   >= -1000
 a_TKT2_fwd_conversion_2_I: - a_TKT2_fwd_conversion_2 + 1000 I_a_TKT2 >= 0
 a_TKT2_rev_conversion_2_I: - a_TKT2_rev_conversion_2 - 1000 I_a_TKT2
   >= -1000
Bounds
 -1000 <= ACALD <= 1000
 -1000 <= ACALDt <= 1000
 -1000 <= ACKr <= 1000
 -1000 <= ACONTa <= 1000
 -1000 <= ACONTb <= 1000
 -1000 <= ACt2r <= 1000
 -1000 <= ADK1 <= 1000
 AKGDH <= 1000
 -1000 <= AKGt2r <= 1000
 -1000 <= ALCD2x <= 1000
 8.39 <= ATPM <= 1000
 -1000 <= ATPS4r <= 1000
 Biomass_Ecoli_core_w_GAM <= 1000
 -1000 <= CO2t <= 1000
 CS <= 1000
 CYTBD <= 1000
 -1000 <= D_LACt2 <= 1000
 -1000 <= ENO <= 1000
 -1000 <= ETOHt2r <= 1000
 EX_ac(e) <= 1000
 EX_acald(e) <= 1000
 EX_akg(e) <= 1000
 -1000 <= EX_co2(e) <= 1000
 EX_etoh(e) <= 1000
 EX_for(e) <= 1000
 EX_fru(e) <= 1000
 EX_fum(e) <= 1000
 -10 <= EX_glc(e) <= 1000
 EX_gln_L(e) <= 1000
 EX_glu_L(e) <= 1000
 -1000 <= EX_h(e) <= 1000
 -1000 <= EX_h2o(e) <= 1000
 EX_lac_D(e) <= 1000
 EX_mal_L(e) <= 1000
 -1000 <= EX_nh4(e) <= 1000
 -1000 <= EX_o2(e) <= 1000
 -1000 <= EX_pi(e) <= 1000
 EX_pyr(e) <= 1000
 EX_succ(e) <= 1000
 -1000 <= FBA <= 1000
 FBP <= 1000
 FORt2 <= 1000
 FORti <= 1000
 FRD7 <= 1000
 FRUpts2 <= 1000
 -1000 <= FUM <= 1000
 FUMt2_2 <= 1000
 -1000 <= G6PDH2r <= 1000
 -1000 <= GAPD <= 1000
 GLCpts <= 1000
 GLNS <= 1000
 GLNabc <= 1000
 -1000 <= GLUDy <= 1000
 GLUN <= 1000
 GLUSy <= 1000
 -1000 <= GLUt2r <= 1000
 GND <= 1000
 -1000 <= H2Ot <= 1000
 -1000 <= ICDHyr <= 1000
 ICL <= 1000
 -1000 <= LDH_D <= 1000
 MALS <= 1000
 MALt2_2 <= 1000
 -1000 <= MDH <= 1000
 ME1 <= 1000
 ME2 <= 1000
 NADH16 <= 1000
 NADTRHD <= 1000
 -1000 <= NH4t <= 1000
 -1000 <= O2t <= 1000
 PDH <= 1000
 PFK <= 1000
 PFL <= 1000
 -1000 <= PGI <= 1000
 -1000 <= PGK <= 1000
 PGL <= 1000
 -1000 <= PGM <= 1000
 -1000 <= PIt2r <= 1000
 PPC <= 1000
 PPCK <= 1000
 PPS <= 1000
 -1000 <= PTAr <= 1000
 PYK <= 1000
 -1000 <= PYRt2r <= 1000
 -1000 <= RPE <= 1000
 -1000 <= RPI <= 1000
 SUCCt2_2 <= 1000
 SUCCt3 <= 1000
 SUCDi <= 1000
 -1000 <= SUCOAS <= 1000
 -1000 <= TALA <= 1000
 THD2 <= 1000
 -1000 <= TKT1 <= 1000
 -1000 <= TKT2 <= 1000
 -1000 <= TPI <= 1000
 -1000 <= Ex_a_b0351 <= 1000
 -1000 <= Ex_a_b1241 <= 1000
 -1000 <= Ex_a_s0001 <= 1000
 -1000 <= Ex_a_b3115 <= 1000
 -1000 <= Ex_a_b2296 <= 1000
 -1000 <= Ex_a_b1849 <= 1000
 -1000 <= Ex_a_b0118 <= 1000
 -1000 <= Ex_a_b1276 <= 1000
 -1000 <= Ex_a_b0474 <= 1000
 -1000 <= Ex_a_b0116 <= 1000
 -1000 <= Ex_a_b0726 <= 1000
 -1000 <= Ex_a_b0727 <= 1000
 -1000 <= Ex_a_b2587 <= 1000
 -1000 <= Ex_a_b0356 <= 1000
 -1000 <= Ex_a_b1478 <= 1000
 -1000 <= Ex_a_b3736 <= 1000
 -1000 <= Ex_a_b3737 <= 1000
 -1000 <= Ex_a_b3738 <= 1000
 -1000 <= Ex_a_b3731 <= 1000
 -1000 <= Ex_a_b3732 <= 1000
 -1000 <= Ex_a_b3733 <= 1000
 -1000 <= Ex_a_b3734 <= 1000
 -1000 <= Ex_a_b3735 <= 1000
 -1000 <= Ex_a_b3739 <= 1000
 -1000 <= Ex_a_b0720 <= 1000
 -1000 <= Ex_a_b0978 <= 1000
 -1000 <= Ex_a_b0979 <= 1000
 -1000 <= Ex_a_b0733 <= 1000
 -1000 <= Ex_a_b0734 <= 1000
 -1000 <= Ex_a_b2975 <= 1000
 -1000 <= Ex_a_b3603 <= 1000
 -1000 <= Ex_a_b2779 <= 1000
 -1000 <= Ex_a_b2097 <= 1000
 -1000 <= Ex_a_b1773 <= 1000
 -1000 <= Ex_a_b2925 <= 1000
 -1000 <= Ex_a_b3925 <= 1000
 -1000 <= Ex_a_b4232 <= 1000
 -1000 <= Ex_a_b0904 <= 1000
 -1000 <= Ex_a_b2492 <= 1000
 -1000 <= Ex_a_b4151 <= 1000
 -1000 <= Ex_a_b4152 <= 1000
 -1000 <= Ex_a_b4153 <= 1000
 -1000 <= Ex_a_b4154 <= 1000
 -1000 <= Ex_a_b1817 <= 1000
 -1000 <= Ex_a_b1818 <= 1000
 -1000 <= Ex_a_b1819 <= 1000
 -1000 <= Ex_a_b2415 <= 1000
 -1000 <= Ex_a_b2416 <= 1000
 -1000 <= Ex_a_b1612 <= 1000
 -1000 <= Ex_a_b4122 <= 1000
 -1000 <= Ex_a_b1611 <= 1000
 -1000 <= Ex_a_b3528 <= 1000
 -1000 <= Ex_a_b1852 <= 1000
 -1000 <= Ex_a_b1779 <= 1000
 -1000 <= Ex_a_b2417 <= 1000
 -1000 <= Ex_a_b1101 <= 1000
 -1000 <= Ex_a_b1621 <= 1000
 -1000 <= Ex_a_b3870 <= 1000
 -1000 <= Ex_a_b1297 <= 1000
 -1000 <= Ex_a_b0811 <= 1000
 -1000 <= Ex_a_b0810 <= 1000
 -1000 <= Ex_a_b0809 <= 1000
 -1000 <= Ex_a_b1761 <= 1000
 -1000 <= Ex_a_b1812 <= 1000
 -1000 <= Ex_a_b0485 <= 1000
 -1000 <= Ex_a_b1524 <= 1000
 -1000 <= Ex_a_b3212 <= 1000
 -1000 <= Ex_a_b3213 <= 1000
 -1000 <= Ex_a_b4077 <= 1000
 -1000 <= Ex_a_b2029 <= 1000
 -1000 <= Ex_a_b0875 <= 1000
 -1000 <= Ex_a_b1136 <= 1000
 -1000 <= Ex_a_b4015 <= 1000
 -1000 <= Ex_a_b2133 <= 1000
 -1000 <= Ex_a_b1380 <= 1000
 -1000 <= Ex_a_b4014 <= 1000
 -1000 <= Ex_a_b2976 <= 1000
 -1000 <= Ex_a_b3236 <= 1000
 -1000 <= Ex_a_b1479 <= 1000
 -1000 <= Ex_a_b2463 <= 1000
 -1000 <= Ex_a_b2276 <= 1000
 -1000 <= Ex_a_b2277 <= 1000
 -1000 <= Ex_a_b2278 <= 1000
 -1000 <= Ex_a_b2279 <= 1000
 -1000 <= Ex_a_b2280 <= 1000
 -1000 <= Ex_a_b2281 <= 1000
 -1000 <= Ex_a_b2282 <= 1000
 -1000 <= Ex_a_b2283 <= 1000
 -1000 <= Ex_a_b2284 <= 1000
 -1000 <= Ex_a_b2285 <= 1000
 -1000 <= Ex_a_b2286 <= 1000
 -1000 <= Ex_a_b2287 <= 1000
 -1000 <= Ex_a_b2288 <= 1000
 -1000 <= Ex_a_b3962 <= 1000
 -1000 <= Ex_a_b1602 <= 1000
 -1000 <= Ex_a_b1603 <= 1000
 -1000 <= Ex_a_b0451 <= 1000
 -1000 <= Ex_a_b0114 <= 1000
 -1000 <= Ex_a_b0115 <= 1000
 -1000 <= Ex_a_b3916 <= 1000
 -1000 <= Ex_a_b1723 <= 1000
 -1000 <= Ex_a_b0902 <= 1000
 -1000 <= Ex_a_b0903 <= 1000
 -1000 <= Ex_a_b2579 <= 1000
 -1000 <= Ex_a_b3114 <= 1000
 -1000 <= Ex_a_b3951 <= 1000
 -1000 <= Ex_a_b3952 <= 1000
 -1000 <= Ex_a_b4025 <= 1000
 -1000 <= Ex_a_b2926 <= 1000
 -1000 <= Ex_a_b0767 <= 1000
 -1000 <= Ex_a_b3612 <= 1000
 -1000 <= Ex_a_b4395 <= 1000
 -1000 <= Ex_a_b0755 <= 1000
 -1000 <= Ex_a_b2987 <= 1000
 -1000 <= Ex_a_b3493 <= 1000
 -1000 <= Ex_a_b3956 <= 1000
 -1000 <= Ex_a_b3403 <= 1000
 -1000 <= Ex_a_b1702 <= 1000
 -1000 <= Ex_a_b2297 <= 1000
 -1000 <= Ex_a_b2458 <= 1000
 -1000 <= Ex_a_b1854 <= 1000
 -1000 <= Ex_a_b1676 <= 1000
 -1000 <= Ex_a_b3386 <= 1000
 -1000 <= Ex_a_b4301 <= 1000
 -1000 <= Ex_a_b2914 <= 1000
 -1000 <= Ex_a_b4090 <= 1000
 -1000 <= Ex_a_b0721 <= 1000
 -1000 <= Ex_a_b0722 <= 1000
 -1000 <= Ex_a_b0723 <= 1000
 -1000 <= Ex_a_b0724 <= 1000
 -1000 <= Ex_a_b0728 <= 1000
 -1000 <= Ex_a_b0729 <= 1000
 -1000 <= Ex_a_b2464 <= 1000
 -1000 <= Ex_a_b0008 <= 1000
 -1000 <= Ex_a_b2935 <= 1000
 -1000 <= Ex_a_b2465 <= 1000
 -1000 <= Ex_a_b3919 <= 1000
 a_ACKr_fwd_conversion_1 <= 1000
 a_ACKr_rev_conversion_1 <= 1000
 a_ACKr_fwd_conversion_2 <= 1000
 a_ACKr_rev_conversion_2 <= 1000
 a_ACKr_fwd_conversion_3 <= 1000
 a_ACKr_rev_conversion_3 <= 1000
 a_ATPS4r_fwd_conversion_1 <= 1000
 a_ATPS4r_rev_conversion_1 <= 1000
 a_ATPS4r_fwd_conversion_2 <= 1000
 a_ATPS4r_rev_conversion_2 <= 1000
 a_CYTBD_fwd_conversion_1 <= 1000
 a_CYTBD_rev_conversion_1 <= 1000
 a_CYTBD_fwd_conversion_2 <= 1000
 a_CYTBD_rev_conversion_2 <= 1000
 a_D_LACt2_fwd_conversion_1 <= 1000
 a_D_LACt2_rev_conversion_1 <= 1000
 a_D_LACt2_fwd_conversion_2 <= 1000
 a_D_LACt2_rev_conversion_2 <= 1000
 a_FBA_fwd_conversion_1 <= 1000
 a_FBA_rev_conversion_1 <= 1000
 a_FBA_fwd_conversion_2 <= 1000
 a_FBA_rev_conversion_2 <= 1000
 a_FBA_fwd_conversion_3 <= 1000
 a_FBA_rev_conversion_3 <= 1000
 a_FBP_fwd_conversion_1 <= 1000
 a_FBP_rev_conversion_1 <= 1000
 a_FBP_fwd_conversion_2 <= 1000
 a_FBP_rev_conversion_2 <= 1000
 a_FUM_fwd_conversion_1 <= 1000
 a_FUM_rev_conversion_1 <= 1000
 a_FUM_fwd_conversion_2 <= 1000
 a_FUM_rev_conversion_2 <= 1000
 a_FUM_fwd_conversion_3 <= 1000
 a_FUM_rev_conversion_3 <= 1000
 a_GLNS_fwd_conversion_1 <= 1000
 a_GLNS_rev_conversion_1 <= 1000
 a_GLNS_fwd_conversion_2 <= 1000
 a_GLNS_rev_conversion_2 <= 1000
 a_GLUN_fwd_conversion_1 <= 1000
 a_GLUN_rev_conversion_1 <= 1000
 a_GLUN_fwd_conversion_2 <= 1000
 a_GLUN_rev_conversion_2 <= 1000
 a_GLUN_fwd_conversion_3 <= 1000
 a_GLUN_rev_conversion_3 <= 1000
 a_LDH_D_fwd_conversion_1 <= 1000
 a_LDH_D_rev_conversion_1 <= 1000
 a_LDH_D_fwd_conversion_2 <= 1000
 a_LDH_D_rev_conversion_2 <= 1000
 a_MALS_fwd_conversion_1 <= 1000
 a_MALS_rev_conversion_1 <= 1000
 a_MALS_fwd_conversion_2 <= 1000
 a_MALS_rev_conversion_2 <= 1000
 a_PFK_fwd_conversion_1 <= 1000
 a_PFK_rev_conversion_1 <= 1000
 a_PFK_fwd_conversion_2 <= 1000
 a_PFK_rev_conversion_2 <= 1000
 a_PFL_fwd_conversion_1 <= 1000
 a_PFL_rev_conversion_1 <= 1000
 a_PFL_fwd_conversion_2 <= 1000
 a_PFL_rev_conversion_2 <= 1000
 a_PFL_fwd_conversion_3 <= 1000
 a_PFL_rev_conversion_3 <= 1000
 a_PFL_fwd_conversion_4 <= 1000
 a_PFL_rev_conversion_4 <= 1000
 a_PGM_fwd_conversion_1 <= 1000
 a_PGM_rev_conversion_1 <= 1000
 a_PGM_fwd_conversion_2 <= 1000
 a_PGM_rev_conversion_2 <= 1000
 a_PGM_fwd_conversion_3 <= 1000
 a_PGM_rev_conversion_3 <= 1000
 a_PIt2r_fwd_conversion_1 <= 1000
 a_PIt2r_rev_conversion_1 <= 1000
 a_PIt2r_fwd_conversion_2 <= 1000
 a_PIt2r_rev_conversion_2 <= 1000
 a_PTAr_fwd_conversion_1 <= 1000
 a_PTAr_rev_conversion_1 <= 1000
 a_PTAr_fwd_conversion_2 <= 1000
 a_PTAr_rev_conversion_2 <= 1000
 a_PYK_fwd_conversion_1 <= 1000
 a_PYK_rev_conversion_1 <= 1000
 a_PYK_fwd_conversion_2 <= 1000
 a_PYK_rev_conversion_2 <= 1000
 a_RPE_fwd_conversion_1 <= 1000
 a_RPE_rev_conversion_1 <= 1000
 a_RPE_fwd_conversion_2 <= 1000
 a_RPE_rev_conversion_2 <= 1000
 a_RPI_fwd_conversion_1 <= 1000
 a_RPI_rev_conversion_1 <= 1000
 a_RPI_fwd_conversion_2 <= 1000
 a_RPI_rev_conversion_2 <= 1000
 a_TALA_fwd_conversion_1 <= 1000
 a_TALA_rev_conversion_1 <= 1000
 a_TALA_fwd_conversion_2 <= 1000
 a_TALA_rev_conversion_2 <= 1000
 a_ACALD_fwd_conversion_1 <= 1000
 a_ACALD_rev_conversion_1 <= 1000
 a_ACALD_fwd_conversion_2 <= 1000
 a_ACALD_rev_conversion_2 <= 1000
 a_ACALDt_fwd_conversion_1 <= 1000
 a_ACALDt_rev_conversion_1 <= 1000
 a_ACONTa_fwd_conversion_1 <= 1000
 a_ACONTa_rev_conversion_1 <= 1000
 a_ACONTa_fwd_conversion_2 <= 1000
 a_ACONTa_rev_conversion_2 <= 1000
 a_ACONTb_fwd_conversion_1 <= 1000
 a_ACONTb_rev_conversion_1 <= 1000
 a_ACONTb_fwd_conversion_2 <= 1000
 a_ACONTb_rev_conversion_2 <= 1000
 a_AKGDH_fwd_conversion_1 <= 1000
 a_AKGDH_rev_conversion_1 <= 1000
 a_ALCD2x_fwd_conversion_1 <= 1000
 a_ALCD2x_rev_conversion_1 <= 1000
 a_ALCD2x_fwd_conversion_2 <= 1000
 a_ALCD2x_rev_conversion_2 <= 1000
 a_ALCD2x_fwd_conversion_3 <= 1000
 a_ALCD2x_rev_conversion_3 <= 1000
 a_CO2t_fwd_conversion_1 <= 1000
 a_CO2t_rev_conversion_1 <= 1000
 a_FORt2_fwd_conversion_1 <= 1000
 a_FORt2_rev_conversion_1 <= 1000
 a_FORt2_fwd_conversion_2 <= 1000
 a_FORt2_rev_conversion_2 <= 1000
 a_FORti_fwd_conversion_1 <= 1000
 a_FORti_rev_conversion_1 <= 1000
 a_FORti_fwd_conversion_2 <= 1000
 a_FORti_rev_conversion_2 <= 1000
 a_FRUpts2_fwd_conversion_1 <= 1000
 a_FRUpts2_rev_conversion_1 <= 1000
 a_FUMt2_2_fwd_conversion_1 <= 1000
 a_FUMt2_2_rev_conversion_1 <= 1000
 a_GLCpts_fwd_conversion_1 <= 1000
 a_GLCpts_rev_conversion_1 <= 1000
 a_GLCpts_fwd_conversion_2 <= 1000
 a_GLCpts_rev_conversion_2 <= 1000
 a_GLCpts_fwd_conversion_3 <= 1000
 a_GLCpts_rev_conversion_3 <= 1000
 a_H2Ot_fwd_conversion_1 <= 1000
 a_H2Ot_rev_conversion_1 <= 1000
 a_H2Ot_fwd_conversion_2 <= 1000
 a_H2Ot_rev_conversion_2 <= 1000
 a_MALt2_2_fwd_conversion_1 <= 1000
 a_MALt2_2_rev_conversion_1 <= 1000
 a_NADTRHD_fwd_conversion_1 <= 1000
 a_NADTRHD_rev_conversion_1 <= 1000
 a_NADTRHD_fwd_conversion_2 <= 1000
 a_NADTRHD_rev_conversion_2 <= 1000
 a_NH4t_fwd_conversion_1 <= 1000
 a_NH4t_rev_conversion_1 <= 1000
 a_NH4t_fwd_conversion_2 <= 1000
 a_NH4t_rev_conversion_2 <= 1000
 a_O2t_fwd_conversion_1 <= 1000
 a_O2t_rev_conversion_1 <= 1000
 a_PDH_fwd_conversion_1 <= 1000
 a_PDH_rev_conversion_1 <= 1000
 a_SUCCt2_2_fwd_conversion_1 <= 1000
 a_SUCCt2_2_rev_conversion_1 <= 1000
 a_THD2_fwd_conversion_1 <= 1000
 a_THD2_rev_conversion_1 <= 1000
 a_TKT1_fwd_conversion_1 <= 1000
 a_TKT1_rev_conversion_1 <= 1000
 a_TKT1_fwd_conversion_2 <= 1000
 a_TKT1_rev_conversion_2 <= 1000
 a_TKT2_fwd_conversion_1 <= 1000
 a_TKT2_rev_conversion_1 <= 1000
 a_TKT2_fwd_conversion_2 <= 1000
 a_TKT2_rev_conversion_2 <= 1000
Binaries
 I_a_ACKr I_a_ATPS4r I_a_CYTBD I_a_D_LACt2 I_a_FBA I_a_FBP I_a_FUM I_a_GLNS
 I_a_GLUN I_a_LDH_D I_a_MALS I_a_PFK I_a_PFL I_a_PGM I_a_PIt2r I_a_PTAr
 I_a_PYK I_a_RPE I_a_RPI I_a_TALA I_a_ACALD I_a_ACALDt I_a_ACONTa
 I_a_ACONTb I_a_AKGDH I_a_ALCD2x I_a_CO2t I_a_FORt2 I_a_FORti I_a_FRUpts2
 I_a_FUMt2_2 I_a_GLCpts I_a_H2Ot I_a_MALt2_2 I_a_NADTRHD I_a_NH4t I_a_O2t
 I_a_PDH I_a_SUCCt2_2 I_a_THD2 I_a_TKT1 I_a_TKT2
End
