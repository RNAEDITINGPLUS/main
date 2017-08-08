# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""miRmap module"""

from . import seed

class mmModel(seed.mmSeed):
    def eval_score(self, model_name=None, model=None):
        """Computes the *miRmap* score(s).

           :param model_name: Model name.
           :type model_name: str
           :param model: Model with coefficients and intercept as keys.
           :type model: dict"""
        # Parameter
        if model_name is None:
            if model is None:
                self.model = Defaults.model
            else:
                self.model = model
        else:
            self.model = Defaults.models[model_name]
        # Reset
        self.scores = []
        # Compute
        for its in range(len(self.end_sites)):
            if Defaults.model_select:
                self.model = Defaults.model_select(self, its)
            score = self.model['intercept']
            for k in self.model.keys():
                if k != 'intercept':
                    score += getattr(self, k+'s')[its] * self.model[k]
            self.scores.append(score)

    @property
    def score(self):
        """*miRmap* score with default parameters."""
        if hasattr(self, 'scores') is False:
            self.eval_score()
        self._score = sum(self.scores)
        return self._score

class Defaults(object):
    models = {}
    # ---------------------------------------------------------------
    # Current models based on Grimson et al dataset
    # ---------------------------------------------------------------
    #
    # -----------------------------------------------
    # All-feature model
    models['full_seed6'] = {}
    models['full_seed6']['tgs_au'] = -0.275016235769136
    models['full_seed6']['tgs_position'] = 5.43367028065211e-06
    models['full_seed6']['tgs_pairing3p'] = -0.00233278119760994
    models['full_seed6']['dg_duplex'] = 0.00772658898496047
    models['full_seed6']['dg_binding'] = -0.00303683833660696
    models['full_seed6']['dg_duplex_seed'] = 0.0496909801533612
    models['full_seed6']['dg_binding_seed'] = -0.048931930580652
    models['full_seed6']['dg_open'] = 0.000674676164622922
    models['full_seed6']['prob_exact'] = 0.16111635592018
    models['full_seed6']['prob_binomial'] = -0.0388333740708671
    models['full_seed6']['cons_bls'] = -0.00426314077593848
    models['full_seed6']['selec_phylop'] = -0.0112455248228072
    models['full_seed6']['intercept'] = 0.148300586692704
    models['full_seed7'] = {}
    models['full_seed7']['tgs_au'] = -0.402470212080983
    models['full_seed7']['tgs_position'] = 6.89249707831041e-05
    models['full_seed7']['tgs_pairing3p'] = -0.0129891251446967
    models['full_seed7']['dg_duplex'] = 0.0141332997802509
    models['full_seed7']['dg_binding'] = -0.0132159175462755
    models['full_seed7']['dg_duplex_seed'] = -0.0814445085121904
    models['full_seed7']['dg_binding_seed'] = 0.115558118311931
    models['full_seed7']['dg_open'] = 0.00331507347139685
    models['full_seed7']['prob_exact'] = 0.792962156550929
    models['full_seed7']['prob_binomial'] = -0.22119499646323
    models['full_seed7']['cons_bls'] = -0.0355840335642203
    models['full_seed7']['selec_phylop'] = -0.0127531995991629
    models['full_seed7']['intercept'] = 0.349448109979275
    # -----------------------------------------------
    # Python-only model
    models['python_only_seed6'] = {}
    models['python_only_seed6']['slope_tgs_au'] = -0.275594504153219
    models['python_only_seed6']['slope_tgs_position'] = 9.44582844229299e-06
    models['python_only_seed6']['slope_tgs_pairing3p'] = -0.0111209267382849
    models['python_only_seed6']['slope_prob_binomial'] = 0.0701619992923641
    models['python_only_seed6']['slope_cons_bls'] = -0.00646548621345819
    models['python_only_seed6']['intercept'] = 0.121104869645859
    models['python_only_seed7'] = {}
    models['python_only_seed7']['slope_tgs_au'] = -0.443606032336791
    models['python_only_seed7']['slope_tgs_position'] = 6.34603935320321e-05
    models['python_only_seed7']['slope_tgs_pairing3p'] = -0.0207672870210752
    models['python_only_seed7']['slope_prob_binomial'] = 0.378665477250754
    models['python_only_seed7']['slope_cons_bls'] = -0.0552713344740971
    models['python_only_seed7']['intercept'] = 0.150015113841088
    # ---------------------------------------------------------------
    # Models used in the NAR manuscript
    # ---------------------------------------------------------------
    #
    # -----------------------------------------------
    # Grimson et al dataset
    models['grimson_full_seed6'] = {}
    models['grimson_full_seed6']['slope_tgs_au'] = -0.281151452127843
    models['grimson_full_seed6']['slope_tgs_position'] = 1.58364508198562e-05
    models['grimson_full_seed6']['slope_tgs_pairing3p'] = -0.00572613771427005
    models['grimson_full_seed6']['slope_dg_duplex'] = 0.00688228600570162
    models['grimson_full_seed6']['slope_dg_binding'] = -0.000679763956488013
    models['grimson_full_seed6']['slope_dg_open'] = 0.00119023516052887
    models['grimson_full_seed6']['slope_prob_exact'] = 0.0594135340475584
    models['grimson_full_seed6']['slope_prob_binomial'] = -0.00317919728661838
    models['grimson_full_seed6']['slope_cons_bls'] = 0.000982269978720738
    models['grimson_full_seed6']['slope_selec_phylop'] = -0.0112364032773455
    models['grimson_full_seed6']['intercept'] = 0.1689635839907
    models['grimson_targetscan_seed6'] = {}
    models['grimson_targetscan_seed6']['slope_tgs_au'] = -0.298007685739703
    models['grimson_targetscan_seed6']['slope_tgs_position'] = 3.93869410309151e-05
    models['grimson_targetscan_seed6']['slope_tgs_pairing3p'] = -0.0167020153729452
    models['grimson_targetscan_seed6']['intercept'] = 0.1557243611024
    models['grimson_full_seed7'] = {}
    models['grimson_full_seed7']['slope_tgs_au'] = -0.411866584343326
    models['grimson_full_seed7']['slope_tgs_position'] = 5.08848962787636e-05
    models['grimson_full_seed7']['slope_tgs_pairing3p'] = -0.0105032269896278
    models['grimson_full_seed7']['slope_dg_duplex'] = 0.00989995601197114
    models['grimson_full_seed7']['slope_dg_binding'] = -0.00291883203289749
    models['grimson_full_seed7']['slope_dg_open'] = 0.00391409605528933
    models['grimson_full_seed7']['slope_prob_exact'] = 1.0215533027696
    models['grimson_full_seed7']['slope_prob_binomial'] = -0.358205454543431
    models['grimson_full_seed7']['slope_cons_bls'] = -0.0343128641803062
    models['grimson_full_seed7']['slope_selec_phylop'] = -0.00962162512595236
    models['grimson_full_seed7']['intercept'] = 0.0808766804281268
    models['grimson_targetscan_seed7'] = {}
    models['grimson_targetscan_seed7']['slope_tgs_au'] = -0.536812485965615
    models['grimson_targetscan_seed7']['slope_tgs_position'] = 0.000123412283564819
    models['grimson_targetscan_seed7']['slope_tgs_pairing3p'] = -0.0257302675767828
    models['grimson_targetscan_seed7']['intercept'] = 0.158640215271411
    # -----------------------------------------------
    # Selbach et al dataset
    models['selbach_psilac_full_seed6'] = {}
    models['selbach_psilac_full_seed6']['slope_tgs_au'] = -0.22547790543644
    models['selbach_psilac_full_seed6']['slope_tgs_position'] = 1.19766415606843e-07
    models['selbach_psilac_full_seed6']['slope_tgs_pairing3p'] = 0.00223612059446233
    models['selbach_psilac_full_seed6']['slope_dg_duplex'] = -0.00410022873387209
    models['selbach_psilac_full_seed6']['slope_dg_binding'] = 0.0102900341645888
    models['selbach_psilac_full_seed6']['slope_dg_open'] = 0.00111905810876512
    models['selbach_psilac_full_seed6']['slope_prob_exact'] = -0.103334478930681
    models['selbach_psilac_full_seed6']['slope_prob_binomial'] = 0.158657698087026
    models['selbach_psilac_full_seed6']['slope_cons_bls'] = -0.00111756166817755
    models['selbach_psilac_full_seed6']['slope_selec_phylop'] = -0.0137389980750266
    models['selbach_psilac_full_seed6']['intercept'] = 0.112893660532015
    models['selbach_psilac_targetscan_seed6'] = {}
    models['selbach_psilac_targetscan_seed6']['slope_tgs_au'] = -0.261795564775942
    models['selbach_psilac_targetscan_seed6']['slope_tgs_position'] = 5.19785408156208e-05
    models['selbach_psilac_targetscan_seed6']['slope_tgs_pairing3p'] = -0.00411247449019235
    models['selbach_psilac_targetscan_seed6']['intercept'] = 0.109985998136486
    models['selbach_psilac_full_seed7'] = {}
    models['selbach_psilac_full_seed7']['slope_tgs_au'] = -0.507289469701232
    models['selbach_psilac_full_seed7']['slope_tgs_position'] = 2.40586733435545e-05
    models['selbach_psilac_full_seed7']['slope_tgs_pairing3p'] = -0.00425164703474552
    models['selbach_psilac_full_seed7']['slope_dg_duplex'] = 0.00600667567620741
    models['selbach_psilac_full_seed7']['slope_dg_binding'] = -0.00473504309814115
    models['selbach_psilac_full_seed7']['slope_dg_open'] = -0.000846719874537811
    models['selbach_psilac_full_seed7']['slope_prob_exact'] = 1.0750233580184
    models['selbach_psilac_full_seed7']['slope_prob_binomial'] = -0.607139311985669
    models['selbach_psilac_full_seed7']['slope_cons_bls'] = -0.00507735949058395
    models['selbach_psilac_full_seed7']['slope_selec_phylop'] = -0.0237087202093476
    models['selbach_psilac_full_seed7']['intercept'] = 0.212588485299256
    models['selbach_psilac_targetscan_seed7'] = {}
    models['selbach_psilac_targetscan_seed7']['slope_tgs_au'] = -0.466088596387115
    models['selbach_psilac_targetscan_seed7']['slope_tgs_position'] = 8.62629085137633e-05
    models['selbach_psilac_targetscan_seed7']['slope_tgs_pairing3p'] = -0.0132448810500069
    models['selbach_psilac_targetscan_seed7']['intercept'] = 0.165809518491651
    # -----------------------------------------------
    # Linsley et al dataset
    models['linsley_full_seed6'] = {}
    models['linsley_full_seed6']['slope_tgs_au'] = -0.222033577974352
    models['linsley_full_seed6']['slope_tgs_position'] = -2.69355145532793e-05
    models['linsley_full_seed6']['slope_tgs_pairing3p'] = -0.0139390679585005
    models['linsley_full_seed6']['slope_dg_duplex'] = -0.0108615051569127
    models['linsley_full_seed6']['slope_dg_binding'] = 0.0134336781495319
    models['linsley_full_seed6']['slope_dg_open'] = 0.000918130748643471
    models['linsley_full_seed6']['slope_prob_exact'] = -0.0335800676662317
    models['linsley_full_seed6']['slope_prob_binomial'] = 0.118998395011815
    models['linsley_full_seed6']['slope_cons_bls'] = -0.0133325474153919
    models['linsley_full_seed6']['slope_selec_phylop'] = -0.0133195616118572
    models['linsley_full_seed6']['intercept'] = 0.159366083623074
    models['linsley_targetscan_seed6'] = {}
    models['linsley_targetscan_seed6']['slope_tgs_au'] = -0.306697463530104
    models['linsley_targetscan_seed6']['slope_tgs_position'] = 1.0013008423451e-05
    models['linsley_targetscan_seed6']['slope_tgs_pairing3p'] = -0.0180002676530701
    models['linsley_targetscan_seed6']['intercept'] = 0.202891234374988
    models['linsley_full_seed7'] = {}
    models['linsley_full_seed7']['slope_tgs_au'] = -0.513793947628812
    models['linsley_full_seed7']['slope_tgs_position'] = 6.57610350835659e-05
    models['linsley_full_seed7']['slope_tgs_pairing3p'] = -0.00675376706692837
    models['linsley_full_seed7']['slope_dg_duplex'] = 0.00586267113783137
    models['linsley_full_seed7']['slope_dg_binding'] = -0.00786316319734974
    models['linsley_full_seed7']['slope_dg_open'] = -0.000389358839639988
    models['linsley_full_seed7']['slope_prob_exact'] = 0.381113950651266
    models['linsley_full_seed7']['slope_prob_binomial'] = -0.104067963956296
    models['linsley_full_seed7']['slope_cons_bls'] = -0.0151259194447267
    models['linsley_full_seed7']['slope_selec_phylop'] = -0.0171549119870622
    models['linsley_full_seed7']['intercept'] = 0.185446406358291
    models['linsley_targetscan_seed7'] = {}
    models['linsley_targetscan_seed7']['slope_tgs_au'] = -0.604323653510091
    models['linsley_targetscan_seed7']['slope_tgs_position'] = 0.000117553009700977
    models['linsley_targetscan_seed7']['slope_tgs_pairing3p'] = -0.00780595149004686
    models['linsley_targetscan_seed7']['intercept'] = 0.214564741916882
    # -----------------------------------------------
    # Chi et al dataset
    models['chi_clip_full_seed6'] = {}
    models['chi_clip_full_seed6']['slope_tgs_au'] = -3.4105589347059
    models['chi_clip_full_seed6']['slope_tgs_position'] = 2.3909956702737e-05
    models['chi_clip_full_seed6']['slope_tgs_pairing3p'] = -0.0251505378437593
    models['chi_clip_full_seed6']['slope_dg_duplex'] = 0.239813705783868
    models['chi_clip_full_seed6']['slope_dg_binding'] = -0.205672331187781
    models['chi_clip_full_seed6']['slope_dg_open'] = -0.0967848697878725
    models['chi_clip_full_seed6']['slope_prob_exact'] = 1.74213943569794
    models['chi_clip_full_seed6']['slope_prob_binomial'] = -0.993371035971244
    models['chi_clip_full_seed6']['slope_cons_bls'] = -0.498200485822226
    models['chi_clip_full_seed6']['slope_selec_phylop'] = 0.813353783900727
    models['chi_clip_full_seed6']['intercept'] = 12.3677944796927
    models['chi_clip_targetscan_seed6'] = {}
    models['chi_clip_targetscan_seed6']['slope_tgs_au'] = 1.63280062843948
    models['chi_clip_targetscan_seed6']['slope_tgs_position'] = -0.000497606846540528
    models['chi_clip_targetscan_seed6']['slope_tgs_pairing3p'] = -0.0362829278355419
    models['chi_clip_targetscan_seed6']['intercept'] = 7.60562144053847
    models['chi_clip_full_seed7'] = {}
    models['chi_clip_full_seed7']['slope_tgs_au'] = 6.9144196209403
    models['chi_clip_full_seed7']['slope_tgs_position'] = -5.70091981035941e-05
    models['chi_clip_full_seed7']['slope_tgs_pairing3p'] = 0.145039793084461
    models['chi_clip_full_seed7']['slope_dg_duplex'] = -0.248746645857738
    models['chi_clip_full_seed7']['slope_dg_binding'] = 0.164278158660717
    models['chi_clip_full_seed7']['slope_dg_open'] = -0.0271318452054673
    models['chi_clip_full_seed7']['slope_prob_exact'] = -3.77818992359216
    models['chi_clip_full_seed7']['slope_prob_binomial'] = 5.61636010175995
    models['chi_clip_full_seed7']['slope_cons_bls'] = -2.26531747643579
    models['chi_clip_full_seed7']['slope_selec_phylop'] = 1.48737297043884
    models['chi_clip_full_seed7']['intercept'] = 2.59788485703054
    models['chi_clip_targetscan_seed7'] = {}
    models['chi_clip_targetscan_seed7']['slope_tgs_au'] = 9.46540145235486
    models['chi_clip_targetscan_seed7']['slope_tgs_position'] = -0.000244104504394474
    models['chi_clip_targetscan_seed7']['slope_tgs_pairing3p'] = 0.345715860080735
    models['chi_clip_targetscan_seed7']['intercept'] = 2.2281423648033
    # -----------------------------------------------
    # Hendrickson et al IP dataset
    models['hendrickson_ip_full_seed6'] = {}
    models['hendrickson_ip_full_seed6']['slope_tgs_au'] = 0.20142997744822
    models['hendrickson_ip_full_seed6']['slope_tgs_position'] = -0.000103767159367576
    models['hendrickson_ip_full_seed6']['slope_tgs_pairing3p'] = 0.0241230613287024
    models['hendrickson_ip_full_seed6']['slope_dg_duplex'] = -0.00322660542544499
    models['hendrickson_ip_full_seed6']['slope_dg_binding'] = -0.00216572752811318
    models['hendrickson_ip_full_seed6']['slope_dg_open'] = 0.00291022506021806
    models['hendrickson_ip_full_seed6']['slope_prob_exact'] = -0.45056523437125
    models['hendrickson_ip_full_seed6']['slope_prob_binomial'] = -0.045738202429398
    models['hendrickson_ip_full_seed6']['slope_cons_bls'] = 0.0270165071011696
    models['hendrickson_ip_full_seed6']['slope_selec_phylop'] = 0.0196340342784427
    models['hendrickson_ip_full_seed6']['intercept'] = -0.170835163801025
    models['hendrickson_ip_targetscan_seed6'] = {}
    models['hendrickson_ip_targetscan_seed6']['slope_tgs_au'] = 0.209076677454877
    models['hendrickson_ip_targetscan_seed6']['slope_tgs_position'] = -0.000292286204635626
    models['hendrickson_ip_targetscan_seed6']['slope_tgs_pairing3p'] = 0.0260770655117346
    models['hendrickson_ip_targetscan_seed6']['intercept'] = -0.133871650030793
    models['hendrickson_ip_full_seed7'] = {}
    models['hendrickson_ip_full_seed7']['slope_tgs_au'] = 2.19435914165139
    models['hendrickson_ip_full_seed7']['slope_tgs_position'] = -0.000837366475583975
    models['hendrickson_ip_full_seed7']['slope_tgs_pairing3p'] = 0.120612949543222
    models['hendrickson_ip_full_seed7']['slope_dg_duplex'] = -0.138497097949444
    models['hendrickson_ip_full_seed7']['slope_dg_binding'] = 0.0885664453963726
    models['hendrickson_ip_full_seed7']['slope_dg_open'] = -0.0202562982328382
    models['hendrickson_ip_full_seed7']['slope_prob_exact'] = 1.62637249873114
    models['hendrickson_ip_full_seed7']['slope_prob_binomial'] = -1.32734202342928
    models['hendrickson_ip_full_seed7']['slope_cons_bls'] = 0.211172639817922
    models['hendrickson_ip_full_seed7']['slope_selec_phylop'] = -0.0591441881245265
    models['hendrickson_ip_full_seed7']['intercept'] = -1.04433898193709
    models['hendrickson_ip_targetscan_seed7'] = {}
    models['hendrickson_ip_targetscan_seed7']['slope_tgs_au'] = 2.23901000340425
    models['hendrickson_ip_targetscan_seed7']['slope_tgs_position'] = -0.0010201780415256
    models['hendrickson_ip_targetscan_seed7']['slope_tgs_pairing3p'] = 0.185947127889188
    models['hendrickson_ip_targetscan_seed7']['intercept'] = -0.454248075828581
    # -----------------------------------------------
    # Hendrickson et al Trans dataset
    models['hendrickson_exp_full_seed6'] = {}
    models['hendrickson_exp_full_seed6']['slope_tgs_au'] = -0.146527188452677
    models['hendrickson_exp_full_seed6']['slope_tgs_position'] = -7.46770344531951e-06
    models['hendrickson_exp_full_seed6']['slope_tgs_pairing3p'] = 0.00413183140941697
    models['hendrickson_exp_full_seed6']['slope_dg_duplex'] = 0.0124967832463282
    models['hendrickson_exp_full_seed6']['slope_dg_binding'] = -0.014107354410741
    models['hendrickson_exp_full_seed6']['slope_dg_open'] = -5.56136530226989e-05
    models['hendrickson_exp_full_seed6']['slope_prob_exact'] = 0.335281874513716
    models['hendrickson_exp_full_seed6']['slope_prob_binomial'] = -0.142677682234688
    models['hendrickson_exp_full_seed6']['slope_cons_bls'] = -0.0159576036595294
    models['hendrickson_exp_full_seed6']['slope_selec_phylop'] = 0.00487211279593042
    models['hendrickson_exp_full_seed6']['intercept'] = 0.036685361637679
    models['hendrickson_exp_targetscan_seed6'] = {}
    models['hendrickson_exp_targetscan_seed6']['slope_tgs_au'] = -0.10791361352
    models['hendrickson_exp_targetscan_seed6']['slope_tgs_position'] = 4.14933621215086e-05
    models['hendrickson_exp_targetscan_seed6']['slope_tgs_pairing3p'] = 0.00535263842670805
    models['hendrickson_exp_targetscan_seed6']['intercept'] = 0.0467334501087324
    models['hendrickson_exp_full_seed7'] = {}
    models['hendrickson_exp_full_seed7']['slope_tgs_au'] = -0.92198378294816
    models['hendrickson_exp_full_seed7']['slope_tgs_position'] = 0.000155188062135395
    models['hendrickson_exp_full_seed7']['slope_tgs_pairing3p'] = 0.00557373741571437
    models['hendrickson_exp_full_seed7']['slope_dg_duplex'] = 0.00279444311822454
    models['hendrickson_exp_full_seed7']['slope_dg_binding'] = 0.0224053629408548
    models['hendrickson_exp_full_seed7']['slope_dg_open'] = 0.00752138496416965
    models['hendrickson_exp_full_seed7']['slope_prob_exact'] = 0.045431931208893
    models['hendrickson_exp_full_seed7']['slope_prob_binomial'] = 0.365587675834314
    models['hendrickson_exp_full_seed7']['slope_cons_bls'] = -0.0573313788746663
    models['hendrickson_exp_full_seed7']['slope_selec_phylop'] = -0.00832908997583119
    models['hendrickson_exp_full_seed7']['intercept'] = 0.421951810550972
    models['hendrickson_exp_targetscan_seed7'] = {}
    models['hendrickson_exp_targetscan_seed7']['slope_tgs_au'] = -1.09961736571824
    models['hendrickson_exp_targetscan_seed7']['slope_tgs_position'] = 0.000309662811320741
    models['hendrickson_exp_targetscan_seed7']['slope_tgs_pairing3p'] = -0.0212241347166091
    models['hendrickson_exp_targetscan_seed7']['intercept'] = 0.212102786925411
    # -----------------------------------------------
    # Hendrickson et al RibN dataset
    models['hendrickson_ribnumber_full_seed6'] = {}
    models['hendrickson_ribnumber_full_seed6']['slope_tgs_au'] = -0.0205882262335487
    models['hendrickson_ribnumber_full_seed6']['slope_tgs_position'] = 1.33520117545911e-05
    models['hendrickson_ribnumber_full_seed6']['slope_tgs_pairing3p'] = -0.00396484131464407
    models['hendrickson_ribnumber_full_seed6']['slope_dg_duplex'] = -0.00185430990700682
    models['hendrickson_ribnumber_full_seed6']['slope_dg_binding'] = 0.00147790739987637
    models['hendrickson_ribnumber_full_seed6']['slope_dg_open'] = -0.000265994343441685
    models['hendrickson_ribnumber_full_seed6']['slope_prob_exact'] = -0.0455151299895328
    models['hendrickson_ribnumber_full_seed6']['slope_prob_binomial'] = 0.000136339066631797
    models['hendrickson_ribnumber_full_seed6']['slope_cons_bls'] = 0.00224954147926021
    models['hendrickson_ribnumber_full_seed6']['slope_selec_phylop'] = 5.88168410592548e-05
    models['hendrickson_ribnumber_full_seed6']['intercept'] = 1.03878009557721
    models['hendrickson_ribnumber_targetscan_seed6'] = {}
    models['hendrickson_ribnumber_targetscan_seed6']['slope_tgs_au'] = -0.0120444536493016
    models['hendrickson_ribnumber_targetscan_seed6']['slope_tgs_position'] = -3.23817167505417e-06
    models['hendrickson_ribnumber_targetscan_seed6']['slope_tgs_pairing3p'] = -0.00330849968664921
    models['hendrickson_ribnumber_targetscan_seed6']['intercept'] = 1.02293915111259
    models['hendrickson_ribnumber_full_seed7'] = {}
    models['hendrickson_ribnumber_full_seed7']['slope_tgs_au'] = -0.14808234945397
    models['hendrickson_ribnumber_full_seed7']['slope_tgs_position'] = 2.16298429155369e-05
    models['hendrickson_ribnumber_full_seed7']['slope_tgs_pairing3p'] = -0.00175155795162633
    models['hendrickson_ribnumber_full_seed7']['slope_dg_duplex'] = 0.000363000080124185
    models['hendrickson_ribnumber_full_seed7']['slope_dg_binding'] = 0.0025135196682497
    models['hendrickson_ribnumber_full_seed7']['slope_dg_open'] = 0.0011411356320409
    models['hendrickson_ribnumber_full_seed7']['slope_prob_exact'] = -0.493586356431284
    models['hendrickson_ribnumber_full_seed7']['slope_prob_binomial'] = 0.27767210962563
    models['hendrickson_ribnumber_full_seed7']['slope_cons_bls'] = -0.00410440907423636
    models['hendrickson_ribnumber_full_seed7']['slope_selec_phylop'] = -0.00164290928497366
    models['hendrickson_ribnumber_full_seed7']['intercept'] = 1.10731525548203
    models['hendrickson_ribnumber_targetscan_seed7'] = {}
    models['hendrickson_ribnumber_targetscan_seed7']['slope_tgs_au'] = -0.135909820409783
    models['hendrickson_ribnumber_targetscan_seed7']['slope_tgs_position'] = 2.98856306844393e-05
    models['hendrickson_ribnumber_targetscan_seed7']['slope_tgs_pairing3p'] = -0.00504340734361336
    models['hendrickson_ribnumber_targetscan_seed7']['intercept'] = 1.04945425805266
    # ---------------------------------------------------------------
    # Default model
    model = models['full_seed7']

    @staticmethod
    def model_select(mim, i):
        if mim.seed_lengths[i] == 6:
            return Defaults.models['full_seed6']
        elif mim.seed_lengths[i] >= 7:
            return Defaults.models['full_seed7']
