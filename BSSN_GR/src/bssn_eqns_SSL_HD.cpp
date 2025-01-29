// CODEGEN: SSL was enabled, adding term to gauge condition!
// CODEGEN: CAHD was enabled, adding damping term to chi!
// Codgen: generating unstage version
// Codgen: using standard gauge
// Codgen: using eta const damping
// Dendro: {{{
// Dendro: original ops: 623324
// Dendro: printing temp variables
const double DENDRO_0 = 2 * alpha[pp];
const double DENDRO_1 = sqrt(chi[pp]);
const double DENDRO_2 =
    (3.0 / 4.0) * alpha[pp] * lambda_f[1] + (3.0 / 4.0) * lambda_f[0];
const double DENDRO_3  = (4.0 / 3.0) * grad_0_beta0[pp];
const double DENDRO_4  = 2 * gt1[pp];
const double DENDRO_5  = 2 * gt2[pp];
const double DENDRO_6  = (2.0 / 3.0) * gt0[pp];
const double DENDRO_7  = (1.0 / 3.0) * gt1[pp];
const double DENDRO_8  = (2.0 / 3.0) * grad_2_beta2[pp];
const double DENDRO_9  = (1.0 / 3.0) * gt2[pp];
const double DENDRO_10 = (2.0 / 3.0) * grad_1_beta1[pp];
const double DENDRO_11 = (2.0 / 3.0) * grad_0_beta0[pp];
const double DENDRO_12 = (4.0 / 3.0) * grad_1_beta1[pp];
const double DENDRO_13 = 2 * gt4[pp];
const double DENDRO_14 = (1.0 / 3.0) * gt4[pp];
const double DENDRO_15 = (4.0 / 3.0) * grad_2_beta2[pp];
const double DENDRO_16 = (2.0 / 3.0) * chi[pp];
const double DENDRO_17 = grad_0_beta0[pp] + grad_1_beta1[pp] + grad_2_beta2[pp];
const double DENDRO_18 = 2 * At1[pp];
const double DENDRO_19 = 2 * At2[pp];
const double DENDRO_20 = gt2[pp] * gt4[pp];
const double DENDRO_21 = -DENDRO_20 + gt1[pp] * gt5[pp];
const double DENDRO_22 = At0[pp] * DENDRO_21;
const double DENDRO_23 = gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp];
const double DENDRO_24 = At2[pp] * DENDRO_23;
const double DENDRO_25 = pow(gt2[pp], 2);
const double DENDRO_26 = -DENDRO_25 + gt0[pp] * gt5[pp];
const double DENDRO_27 = pow(gt4[pp], 2);
const double DENDRO_28 = pow(gt1[pp], 2);
const double DENDRO_29 = gt3[pp] * gt5[pp];
const double DENDRO_30 = -DENDRO_20 * DENDRO_4 + DENDRO_25 * gt3[pp] +
                         DENDRO_27 * gt0[pp] + DENDRO_28 * gt5[pp] -
                         DENDRO_29 * gt0[pp];
const double DENDRO_31 = 1.0 / DENDRO_30;
const double DENDRO_32 = DENDRO_18 * DENDRO_31;
const double DENDRO_33 = gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp];
const double DENDRO_34 = At2[pp] * DENDRO_33;
const double DENDRO_35 = -DENDRO_27 + DENDRO_29;
const double DENDRO_36 = At0[pp] * DENDRO_35;
const double DENDRO_37 = 2 * DENDRO_31;
const double DENDRO_38 = At0[pp] * DENDRO_37;
const double DENDRO_39 = -DENDRO_28 + gt0[pp] * gt3[pp];
const double DENDRO_40 = At2[pp] * DENDRO_39;
const double DENDRO_41 = DENDRO_19 * DENDRO_31;
const double DENDRO_42 =
    DENDRO_21 * grad_0_chi[pp] + DENDRO_23 * grad_2_chi[pp];
const double DENDRO_43 = -DENDRO_26 * grad_1_chi[pp] + DENDRO_42;
const double DENDRO_44 = 0.5 * DENDRO_43;
const double DENDRO_45 = 1.0 / chi[pp];
const double DENDRO_46 = DENDRO_45 * gt0[pp];
const double DENDRO_47 = 1.0 * grad_0_gt1[pp];
const double DENDRO_48 = 0.5 * grad_1_gt0[pp];
const double DENDRO_49 = DENDRO_47 - DENDRO_48;
const double DENDRO_50 = 1.0 * grad_0_gt2[pp];
const double DENDRO_51 = 0.5 * grad_2_gt0[pp];
const double DENDRO_52 = DENDRO_50 - DENDRO_51;
const double DENDRO_53 = 0.5 * grad_0_gt0[pp];
const double DENDRO_54 = DENDRO_21 * DENDRO_53 + DENDRO_23 * DENDRO_52;
const double DENDRO_55 = -DENDRO_26 * DENDRO_49 + DENDRO_54;
const double DENDRO_56 = DENDRO_44 * DENDRO_46 + DENDRO_55;
const double DENDRO_57 = 12 * grad_1_alpha[pp];
const double DENDRO_58 = DENDRO_31 * DENDRO_57;
const double DENDRO_59 = -DENDRO_23 * grad_1_chi[pp] +
                         DENDRO_33 * grad_0_chi[pp] +
                         DENDRO_39 * grad_2_chi[pp];
const double DENDRO_60 = -DENDRO_59;
const double DENDRO_61 =
    -DENDRO_23 * DENDRO_49 + DENDRO_33 * DENDRO_53 + DENDRO_39 * DENDRO_52;
const double DENDRO_62 = -0.5 * DENDRO_45 * DENDRO_60 * gt0[pp] + DENDRO_61;
const double DENDRO_63 = 12 * grad_2_alpha[pp];
const double DENDRO_64 = DENDRO_31 * DENDRO_63;
const double DENDRO_65 = DENDRO_35 * DENDRO_53;
const double DENDRO_66 = DENDRO_33 * DENDRO_52;
const double DENDRO_67 = DENDRO_21 * DENDRO_49;
const double DENDRO_68 = DENDRO_31 * DENDRO_67;
const double DENDRO_69 = 1.0 * grad_0_chi[pp];
const double DENDRO_70 = DENDRO_31 * gt0[pp];
const double DENDRO_71 = -DENDRO_21 * grad_1_chi[pp] +
                         DENDRO_33 * grad_2_chi[pp] +
                         DENDRO_35 * grad_0_chi[pp];
const double DENDRO_72 = -DENDRO_71;
const double DENDRO_73 = 0.5 * DENDRO_72;
const double DENDRO_74 = DENDRO_31 * DENDRO_65 + DENDRO_31 * DENDRO_66 +
                         DENDRO_45 * (DENDRO_69 - DENDRO_70 * DENDRO_73) -
                         DENDRO_68;
const double DENDRO_75  = 12 * grad_0_alpha[pp];
const double DENDRO_76  = pow(chi[pp], -2);
const double DENDRO_77  = pow(grad_0_chi[pp], 2);
const double DENDRO_78  = DENDRO_76 * DENDRO_77;
const double DENDRO_79  = 4.0 * DENDRO_31;
const double DENDRO_80  = DENDRO_21 * DENDRO_79;
const double DENDRO_81  = DENDRO_80 * grad2_0_1_gt0[pp];
const double DENDRO_82  = DENDRO_23 * DENDRO_79;
const double DENDRO_83  = DENDRO_82 * grad2_1_2_gt0[pp];
const double DENDRO_84  = pow(DENDRO_30, -2);
const double DENDRO_85  = DENDRO_21 * grad_0_gt3[pp];
const double DENDRO_86  = DENDRO_35 * grad_1_gt0[pp];
const double DENDRO_87  = grad_0_gt4[pp] + grad_1_gt2[pp] - grad_2_gt1[pp];
const double DENDRO_88  = DENDRO_33 * DENDRO_87;
const double DENDRO_89  = -DENDRO_85 + DENDRO_86 + DENDRO_88;
const double DENDRO_90  = DENDRO_33 * grad_0_gt5[pp];
const double DENDRO_91  = DENDRO_35 * grad_2_gt0[pp];
const double DENDRO_92  = grad_0_gt4[pp] - grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_93  = DENDRO_21 * DENDRO_92;
const double DENDRO_94  = DENDRO_90 + DENDRO_91 - DENDRO_93;
const double DENDRO_95  = -DENDRO_21 * DENDRO_49 + DENDRO_65 + DENDRO_66;
const double DENDRO_96  = DENDRO_26 * grad_0_gt3[pp];
const double DENDRO_97  = DENDRO_21 * grad_1_gt0[pp];
const double DENDRO_98  = DENDRO_23 * DENDRO_87;
const double DENDRO_99  = DENDRO_97 + DENDRO_98;
const double DENDRO_100 = -DENDRO_96 + DENDRO_99;
const double DENDRO_101 = 0.25 * grad_0_gt3[pp];
const double DENDRO_102 = 1.0 * grad_1_gt1[pp];
const double DENDRO_103 = -DENDRO_102;
const double DENDRO_104 = DENDRO_26 * DENDRO_84;
const double DENDRO_105 = 4 * DENDRO_104;
const double DENDRO_106 = DENDRO_100 * DENDRO_105 * (-DENDRO_101 - DENDRO_103);
const double DENDRO_107 = DENDRO_33 * grad_2_gt0[pp];
const double DENDRO_108 = DENDRO_39 * grad_0_gt5[pp];
const double DENDRO_109 = DENDRO_107 + DENDRO_108 - DENDRO_23 * DENDRO_92;
const double DENDRO_110 = 0.25 * grad_0_gt5[pp];
const double DENDRO_111 = 1.0 * grad_2_gt2[pp];
const double DENDRO_112 = -DENDRO_111;
const double DENDRO_113 = DENDRO_110 + DENDRO_112;
const double DENDRO_114 = DENDRO_39 * DENDRO_84;
const double DENDRO_115 = 4 * DENDRO_114;
const double DENDRO_116 =
    DENDRO_21 * grad_2_gt0[pp] + DENDRO_23 * grad_0_gt5[pp];
const double DENDRO_117 = DENDRO_116 - DENDRO_26 * DENDRO_92;
const double DENDRO_118 = 0.25 * grad_0_gt4[pp];
const double DENDRO_119 = -DENDRO_118;
const double DENDRO_120 = 0.25 * grad_1_gt2[pp];
const double DENDRO_121 = 0.75 * grad_2_gt1[pp];
const double DENDRO_122 =
    DENDRO_115 * DENDRO_117 * (DENDRO_119 + DENDRO_120 + DENDRO_121);
const double DENDRO_123 = 0.75 * grad_1_gt2[pp];
const double DENDRO_124 = 0.25 * grad_2_gt1[pp];
const double DENDRO_125 = DENDRO_119 + DENDRO_123 + DENDRO_124;
const double DENDRO_126 = DENDRO_33 * grad_1_gt0[pp];
const double DENDRO_127 = DENDRO_39 * DENDRO_87;
const double DENDRO_128 = DENDRO_126 + DENDRO_127 - DENDRO_23 * grad_0_gt3[pp];
const double DENDRO_129 = DENDRO_35 * DENDRO_55;
const double DENDRO_130 = 4 * DENDRO_84;
const double DENDRO_131 = DENDRO_129 * DENDRO_130 * (DENDRO_47 + DENDRO_48);
const double DENDRO_132 = DENDRO_50 + DENDRO_51;
const double DENDRO_133 = 0.25 * grad_1_gt0[pp];
const double DENDRO_134 = DENDRO_133 * DENDRO_94;
const double DENDRO_135 = DENDRO_23 * DENDRO_84;
const double DENDRO_136 = 4 * DENDRO_135;
const double DENDRO_137 = 0.25 * grad_2_gt0[pp];
const double DENDRO_138 = DENDRO_137 * DENDRO_89;
const double DENDRO_139 = DENDRO_89 * grad_0_gt0[pp];
const double DENDRO_140 = DENDRO_95 * grad_1_gt0[pp];
const double DENDRO_141 = DENDRO_21 * DENDRO_84;
const double DENDRO_142 = 2.0 * DENDRO_141;
const double DENDRO_143 = DENDRO_94 * grad_0_gt0[pp];
const double DENDRO_144 = DENDRO_95 * grad_2_gt0[pp];
const double DENDRO_145 = DENDRO_100 * grad_1_gt0[pp];
const double DENDRO_146 = DENDRO_55 * grad_0_gt3[pp];
const double DENDRO_147 = DENDRO_145 + DENDRO_146;
const double DENDRO_148 = DENDRO_109 * grad_2_gt0[pp];
const double DENDRO_149 = DENDRO_61 * grad_0_gt5[pp];
const double DENDRO_150 = DENDRO_101 * DENDRO_117;
const double DENDRO_151 = -grad_0_gt4[pp] + grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_152 = 0.5 * DENDRO_151;
const double DENDRO_153 = DENDRO_100 * DENDRO_152 + DENDRO_150;
const double DENDRO_154 = DENDRO_110 * DENDRO_128;
const double DENDRO_155 = 0.5 * DENDRO_109;
const double DENDRO_156 = 0.25 * DENDRO_139;
const double DENDRO_157 = 4 * DENDRO_141;
const double DENDRO_158 = 0.25 * DENDRO_143;
const double DENDRO_159 = DENDRO_117 * grad_1_gt0[pp];
const double DENDRO_160 = DENDRO_55 * DENDRO_92;
const double DENDRO_161 = DENDRO_33 * DENDRO_84;
const double DENDRO_162 = 2.0 * DENDRO_161;
const double DENDRO_163 = DENDRO_162 * (DENDRO_159 + DENDRO_160);
const double DENDRO_164 = DENDRO_128 * grad_2_gt0[pp];
const double DENDRO_165 = DENDRO_61 * DENDRO_87;
const double DENDRO_166 = 0.5 * grad_0_gt3[pp];
const double DENDRO_167 = DENDRO_103 + DENDRO_166;
const double DENDRO_168 = DENDRO_117 * DENDRO_167;
const double DENDRO_169 = DENDRO_100 * DENDRO_92;
const double DENDRO_170 = 0.25 * DENDRO_169;
const double DENDRO_171 = -DENDRO_170;
const double DENDRO_172 = 0.5 * grad_0_gt5[pp];
const double DENDRO_173 = DENDRO_112 + DENDRO_172;
const double DENDRO_174 = -DENDRO_128;
const double DENDRO_175 = DENDRO_173 * DENDRO_174;
const double DENDRO_176 = -DENDRO_109;
const double DENDRO_177 = DENDRO_176 * DENDRO_87;
const double DENDRO_178 = -0.25 * DENDRO_177;
const double DENDRO_179 = 0.5 * DENDRO_49;
const double DENDRO_180 = DENDRO_117 * DENDRO_179;
const double DENDRO_181 = 4 * DENDRO_161;
const double DENDRO_182 = DENDRO_181 * (DENDRO_151 * DENDRO_55 + DENDRO_180);
const double DENDRO_183 = 0.5 * DENDRO_52;
const double DENDRO_184 = DENDRO_100 * DENDRO_179;
const double DENDRO_185 = -2 * DENDRO_167 * DENDRO_55 + DENDRO_184;
const double DENDRO_186 = -grad2_0_0_chi[pp];
const double DENDRO_187 = -DENDRO_95;
const double DENDRO_188 = DENDRO_31 * grad_0_chi[pp];
const double DENDRO_189 = DENDRO_31 * grad_1_chi[pp];
const double DENDRO_190 = -DENDRO_61;
const double DENDRO_191 = DENDRO_31 * grad_2_chi[pp];
const double DENDRO_192 = 2 * DENDRO_45;
const double DENDRO_193 = DENDRO_117 * DENDRO_33;
const double DENDRO_194 = DENDRO_26 * grad_2_gt3[pp];
const double DENDRO_195 = DENDRO_23 * grad_1_gt5[pp];
const double DENDRO_196 = DENDRO_151 * DENDRO_21;
const double DENDRO_197 = DENDRO_195 + DENDRO_196;
const double DENDRO_198 = -DENDRO_194 + DENDRO_197;
const double DENDRO_199 = 0.5 * grad_2_gt5[pp];
const double DENDRO_200 = DENDRO_199 * DENDRO_23;
const double DENDRO_201 = 0.5 * grad_1_gt5[pp];
const double DENDRO_202 = 1.0 * grad_2_gt4[pp];
const double DENDRO_203 = -DENDRO_202;
const double DENDRO_204 = DENDRO_201 + DENDRO_203;
const double DENDRO_205 =
    -DENDRO_173 * DENDRO_21 + DENDRO_200 + DENDRO_204 * DENDRO_26;
const double DENDRO_206 = DENDRO_205 * DENDRO_39;
const double DENDRO_207 = 0.5 * grad_1_gt3[pp];
const double DENDRO_208 = DENDRO_207 * DENDRO_26;
const double DENDRO_209 = DENDRO_167 * DENDRO_21;
const double DENDRO_210 = 1.0 * grad_1_gt4[pp];
const double DENDRO_211 = 0.5 * grad_2_gt3[pp];
const double DENDRO_212 = DENDRO_210 - DENDRO_211;
const double DENDRO_213 = DENDRO_208 + DENDRO_209 - DENDRO_212 * DENDRO_23;
const double DENDRO_214 = -DENDRO_213;
const double DENDRO_215 = DENDRO_214 * DENDRO_26;
const double DENDRO_216 = DENDRO_129 + DENDRO_206 + DENDRO_215;
const double DENDRO_217 =
    DENDRO_84 * (-1.0 * DENDRO_100 * DENDRO_21 + DENDRO_193 -
                 1.0 * DENDRO_198 * DENDRO_23 + DENDRO_216);
const double DENDRO_218 = 2.0 * grad_1_gt0[pp];
const double DENDRO_219 = -DENDRO_94;
const double DENDRO_220 = DENDRO_219 * DENDRO_33;
const double DENDRO_221 = DENDRO_151 * DENDRO_35 - DENDRO_21 * grad_2_gt3[pp] +
                          DENDRO_33 * grad_1_gt5[pp];
const double DENDRO_222 = -DENDRO_221;
const double DENDRO_223 = -DENDRO_89;
const double DENDRO_224 =
    -DENDRO_173 * DENDRO_35 + DENDRO_199 * DENDRO_33 + DENDRO_204 * DENDRO_21;
const double DENDRO_225 = -DENDRO_224;
const double DENDRO_226 = DENDRO_225 * DENDRO_39;
const double DENDRO_227 = DENDRO_207 * DENDRO_21;
const double DENDRO_228 =
    DENDRO_167 * DENDRO_35 - DENDRO_212 * DENDRO_33 + DENDRO_227;
const double DENDRO_229 = DENDRO_228 * DENDRO_26;
const double DENDRO_230 = DENDRO_187 * DENDRO_35;
const double DENDRO_231 = DENDRO_226 + DENDRO_229 + DENDRO_230;
const double DENDRO_232 =
    DENDRO_84 * (-1.0 * DENDRO_21 * DENDRO_223 + DENDRO_220 -
                 1.0 * DENDRO_222 * DENDRO_23 + DENDRO_231);
const double DENDRO_233 = 2.0 * grad_0_gt0[pp];
const double DENDRO_234 = DENDRO_176 * DENDRO_33;
const double DENDRO_235 = DENDRO_39 * grad_1_gt5[pp];
const double DENDRO_236 = DENDRO_151 * DENDRO_33;
const double DENDRO_237 = -DENDRO_23 * grad_2_gt3[pp] + DENDRO_235 + DENDRO_236;
const double DENDRO_238 = -DENDRO_237;
const double DENDRO_239 = DENDRO_199 * DENDRO_39;
const double DENDRO_240 = DENDRO_204 * DENDRO_23;
const double DENDRO_241 = -DENDRO_173 * DENDRO_33 + DENDRO_239 + DENDRO_240;
const double DENDRO_242 = -DENDRO_241;
const double DENDRO_243 = DENDRO_242 * DENDRO_39;
const double DENDRO_244 = DENDRO_207 * DENDRO_23;
const double DENDRO_245 =
    DENDRO_167 * DENDRO_33 - DENDRO_212 * DENDRO_39 + DENDRO_244;
const double DENDRO_246 = DENDRO_245 * DENDRO_26;
const double DENDRO_247 = DENDRO_190 * DENDRO_35;
const double DENDRO_248 = DENDRO_243 + DENDRO_246 + DENDRO_247;
const double DENDRO_249 =
    DENDRO_84 * (-1.0 * DENDRO_174 * DENDRO_21 - 1.0 * DENDRO_23 * DENDRO_238 +
                 DENDRO_234 + DENDRO_248);
const double DENDRO_250 = 2.0 * grad_2_gt0[pp];
const double DENDRO_251 = 3 * DENDRO_45;
const double DENDRO_252 = grad_0_chi[pp] * grad_1_chi[pp];
const double DENDRO_253 = 2 * DENDRO_21;
const double DENDRO_254 = DENDRO_251 * grad_2_chi[pp];
const double DENDRO_255 = 2 * DENDRO_33;
const double DENDRO_256 = 2 * DENDRO_23;
const double DENDRO_257 = pow(grad_1_chi[pp], 2);
const double DENDRO_258 = pow(grad_2_chi[pp], 2);
const double DENDRO_259 = DENDRO_100 * DENDRO_21 - 1.0 * DENDRO_193 +
                          DENDRO_198 * DENDRO_23 - DENDRO_216;
const double DENDRO_260 = DENDRO_21 * DENDRO_223 - 1.0 * DENDRO_220 +
                          DENDRO_222 * DENDRO_23 - DENDRO_231;
const double DENDRO_261 = DENDRO_174 * DENDRO_21 + DENDRO_23 * DENDRO_238 -
                          1.0 * DENDRO_234 - DENDRO_248;
const double DENDRO_262 =
    2 * DENDRO_188 * DENDRO_260 + 2 * DENDRO_189 * DENDRO_259 +
    2 * DENDRO_191 * DENDRO_261 -
    DENDRO_253 * (-DENDRO_251 * DENDRO_252 + 2 * grad2_0_1_chi[pp]) +
    DENDRO_255 * (-DENDRO_254 * grad_0_chi[pp] + 2 * grad2_0_2_chi[pp]) -
    DENDRO_256 * (-DENDRO_254 * grad_1_chi[pp] + 2 * grad2_1_2_chi[pp]) +
    DENDRO_26 * (-DENDRO_251 * DENDRO_257 + 2 * grad2_1_1_chi[pp]) +
    DENDRO_35 * (-DENDRO_251 * DENDRO_77 + 2 * grad2_0_0_chi[pp]) +
    DENDRO_39 * (-DENDRO_251 * DENDRO_258 + 2 * grad2_2_2_chi[pp]);
const double DENDRO_263 = DENDRO_45 * DENDRO_70;
const double DENDRO_264 = 3 * alpha[pp];
const double DENDRO_265 = DENDRO_45 * gt5[pp];
const double DENDRO_266 = DENDRO_205 + DENDRO_265 * DENDRO_44;
const double DENDRO_267 = 4 * grad_1_alpha[pp];
const double DENDRO_268 = DENDRO_267 * DENDRO_31;
const double DENDRO_269 = DENDRO_224 - 0.5 * DENDRO_45 * DENDRO_72 * gt5[pp];
const double DENDRO_270 = 4 * grad_0_alpha[pp];
const double DENDRO_271 = DENDRO_270 * DENDRO_31;
const double DENDRO_272 = 1.0 * grad_2_chi[pp];
const double DENDRO_273 = DENDRO_31 * gt5[pp];
const double DENDRO_274 = 0.5 * DENDRO_60;
const double DENDRO_275 = -DENDRO_173 * DENDRO_31 * DENDRO_33 +
                          DENDRO_239 * DENDRO_31 + DENDRO_240 * DENDRO_31 +
                          DENDRO_45 * (DENDRO_272 - DENDRO_273 * DENDRO_274);
const double DENDRO_276 = 4 * grad_2_alpha[pp];
const double DENDRO_277 = 4 * gt2[pp];
const double DENDRO_278 = 4 * gt4[pp];
const double DENDRO_279 = DENDRO_258 * DENDRO_76;
const double DENDRO_280 = DENDRO_80 * grad2_0_1_gt5[pp];
const double DENDRO_281 = DENDRO_82 * grad2_1_2_gt5[pp];
const double DENDRO_282 = DENDRO_31 * DENDRO_33;
const double DENDRO_283 = 4 * DENDRO_282;
const double DENDRO_284 = DENDRO_31 * DENDRO_35;
const double DENDRO_285 = 2.0 * DENDRO_284;
const double DENDRO_286 = DENDRO_26 * DENDRO_31;
const double DENDRO_287 = 2.0 * DENDRO_286;
const double DENDRO_288 = DENDRO_31 * DENDRO_39;
const double DENDRO_289 = 2.0 * DENDRO_288;
const double DENDRO_290 = DENDRO_176 * grad_0_gt5[pp];
const double DENDRO_291 = DENDRO_35 * DENDRO_84;
const double DENDRO_292 = 3.0 * DENDRO_291;
const double DENDRO_293 = DENDRO_238 * grad_1_gt5[pp];
const double DENDRO_294 = 3.0 * DENDRO_104;
const double DENDRO_295 = 6.0 * DENDRO_84;
const double DENDRO_296 = 0.25 * grad_2_gt3[pp];
const double DENDRO_297 = DENDRO_105 * DENDRO_198 * (DENDRO_210 - DENDRO_296);
const double DENDRO_298 = -DENDRO_137 + DENDRO_50;
const double DENDRO_299 = 4 * DENDRO_291;
const double DENDRO_300 = 0.75 * grad_0_gt4[pp];
const double DENDRO_301 = -DENDRO_124;
const double DENDRO_302 =
    DENDRO_117 * DENDRO_299 * (DENDRO_120 + DENDRO_300 + DENDRO_301);
const double DENDRO_303 = DENDRO_118 + DENDRO_123 + DENDRO_301;
const double DENDRO_304 = DENDRO_111 + DENDRO_172;
const double DENDRO_305 = DENDRO_130 * DENDRO_206 * (DENDRO_201 + DENDRO_202);
const double DENDRO_306 = DENDRO_110 * DENDRO_238;
const double DENDRO_307 = 0.25 * grad_1_gt5[pp];
const double DENDRO_308 = DENDRO_176 * DENDRO_307;
const double DENDRO_309 = DENDRO_198 * grad_1_gt5[pp];
const double DENDRO_310 = DENDRO_205 * grad_2_gt3[pp];
const double DENDRO_311 = DENDRO_309 + DENDRO_310;
const double DENDRO_312 = 2.0 * DENDRO_135;
const double DENDRO_313 = DENDRO_225 * grad_2_gt0[pp];
const double DENDRO_314 = DENDRO_176 * grad_2_gt5[pp];
const double DENDRO_315 = DENDRO_242 * grad_0_gt5[pp];
const double DENDRO_316 = DENDRO_238 * grad_2_gt5[pp];
const double DENDRO_317 = DENDRO_242 * grad_1_gt5[pp];
const double DENDRO_318 = DENDRO_117 * DENDRO_296;
const double DENDRO_319 = 0.5 * DENDRO_87;
const double DENDRO_320 = DENDRO_198 * DENDRO_319 + DENDRO_318;
const double DENDRO_321 = 0.25 * DENDRO_314;
const double DENDRO_322 = 0.25 * DENDRO_316;
const double DENDRO_323 = DENDRO_137 * DENDRO_222;
const double DENDRO_324 = DENDRO_117 * grad_1_gt5[pp];
const double DENDRO_325 = DENDRO_205 * DENDRO_92;
const double DENDRO_326 = DENDRO_162 * (DENDRO_324 + DENDRO_325);
const double DENDRO_327 = DENDRO_117 * DENDRO_212;
const double DENDRO_328 = DENDRO_198 * DENDRO_92;
const double DENDRO_329 = 0.25 * DENDRO_328;
const double DENDRO_330 = DENDRO_327 + DENDRO_329;
const double DENDRO_331 = DENDRO_151 * DENDRO_225;
const double DENDRO_332 = DENDRO_222 * DENDRO_52;
const double DENDRO_333 = DENDRO_151 * DENDRO_219;
const double DENDRO_334 = 0.25 * DENDRO_333;
const double DENDRO_335 = 0.5 * DENDRO_204;
const double DENDRO_336 = DENDRO_117 * DENDRO_335;
const double DENDRO_337 = 0.5 * DENDRO_173;
const double DENDRO_338 = DENDRO_222 * DENDRO_337;
const double DENDRO_339 = -DENDRO_198 * DENDRO_335;
const double DENDRO_340 = 2 * DENDRO_205 * DENDRO_212 + DENDRO_339;
const double DENDRO_341 = -DENDRO_219 * DENDRO_337;
const double DENDRO_342 = 2 * DENDRO_52;
const double DENDRO_343 = -grad2_2_2_chi[pp];
const double DENDRO_344 = -DENDRO_33;
const double DENDRO_345 = -DENDRO_173;
const double DENDRO_346 = -DENDRO_35;
const double DENDRO_347 = -DENDRO_204;
const double DENDRO_348 =
    DENDRO_199 * DENDRO_344 + DENDRO_21 * DENDRO_347 + DENDRO_345 * DENDRO_346;
const double DENDRO_349 = -DENDRO_26;
const double DENDRO_350 =
    DENDRO_200 + DENDRO_21 * DENDRO_345 + DENDRO_347 * DENDRO_349;
const double DENDRO_351 = -DENDRO_39;
const double DENDRO_352 = DENDRO_199 * DENDRO_351;
const double DENDRO_353 = DENDRO_344 * DENDRO_345;
const double DENDRO_354 = DENDRO_23 * DENDRO_347;
const double DENDRO_355 = DENDRO_259 * DENDRO_84;
const double DENDRO_356 = 2.0 * DENDRO_355;
const double DENDRO_357 = DENDRO_260 * DENDRO_84;
const double DENDRO_358 = 2.0 * DENDRO_357;
const double DENDRO_359 = DENDRO_261 * DENDRO_84;
const double DENDRO_360 = 2.0 * DENDRO_359;
const double DENDRO_361 = -DENDRO_262;
const double DENDRO_362 = DENDRO_361 * DENDRO_45;
const double DENDRO_363 = DENDRO_45 * gt3[pp];
const double DENDRO_364 = DENDRO_276 * DENDRO_31;
const double DENDRO_365 = 1.0 * grad_1_chi[pp];
const double DENDRO_366 = DENDRO_31 * gt3[pp];
const double DENDRO_367 = DENDRO_208 * DENDRO_31 + DENDRO_209 * DENDRO_31 -
                          DENDRO_212 * DENDRO_23 * DENDRO_31 +
                          DENDRO_45 * (DENDRO_365 - DENDRO_366 * DENDRO_44);
const double DENDRO_368 = -grad2_1_1_chi[pp];
const double DENDRO_369 = -DENDRO_167;
const double DENDRO_370 =
    DENDRO_212 * DENDRO_344 + DENDRO_227 + DENDRO_346 * DENDRO_369;
const double DENDRO_371 = DENDRO_207 * DENDRO_349;
const double DENDRO_372 = DENDRO_21 * DENDRO_369;
const double DENDRO_373 = DENDRO_212 * DENDRO_23;
const double DENDRO_374 =
    DENDRO_212 * DENDRO_351 + DENDRO_244 + DENDRO_344 * DENDRO_369;
const double DENDRO_375 = DENDRO_222 * DENDRO_49;
const double DENDRO_376 = DENDRO_151 * DENDRO_223;
const double DENDRO_377 = 0.25 * DENDRO_376;
const double DENDRO_378 = DENDRO_133 * DENDRO_222;
const double DENDRO_379 = 0.5 * DENDRO_92;
const double DENDRO_380 = DENDRO_174 * DENDRO_307;
const double DENDRO_381 = DENDRO_174 * DENDRO_204;
const double DENDRO_382 = DENDRO_238 * DENDRO_87;
const double DENDRO_383 = -0.25 * DENDRO_382;
const double DENDRO_384 = 0.5 * DENDRO_167;
const double DENDRO_385 = DENDRO_222 * DENDRO_384;
const double DENDRO_386 = 0.5 * DENDRO_212;
const double DENDRO_387 = DENDRO_238 * DENDRO_386;
const double DENDRO_388 = 2 * DENDRO_204 * DENDRO_245;
const double DENDRO_389 = DENDRO_198 * grad_1_gt3[pp];
const double DENDRO_390 = 0.25 * DENDRO_389;
const double DENDRO_391 = DENDRO_214 * grad_2_gt3[pp];
const double DENDRO_392 = DENDRO_174 * DENDRO_386;
const double DENDRO_393 = -DENDRO_223 * DENDRO_384;
const double DENDRO_394 = 2 * DENDRO_228 * DENDRO_49;
const double DENDRO_395 = DENDRO_100 * grad_1_gt3[pp];
const double DENDRO_396 = 0.25 * DENDRO_395;
const double DENDRO_397 = DENDRO_214 * grad_0_gt3[pp];
const double DENDRO_398 = DENDRO_228 * grad_1_gt0[pp];
const double DENDRO_399 = DENDRO_151 * DENDRO_228;
const double DENDRO_400 = DENDRO_245 * grad_1_gt5[pp];
const double DENDRO_401 = DENDRO_245 * DENDRO_87;
const double DENDRO_402 = DENDRO_203 + DENDRO_307;
const double DENDRO_403 = -DENDRO_120;
const double DENDRO_404 = DENDRO_115 * (DENDRO_118 + DENDRO_121 + DENDRO_403);
const double DENDRO_405 = DENDRO_299 * (-DENDRO_133 + DENDRO_47);
const double DENDRO_406 = DENDRO_299 * (DENDRO_124 + DENDRO_300 + DENDRO_403);
const double DENDRO_407 = 4 * gt1[pp];
const double DENDRO_408 = DENDRO_101 * DENDRO_198;
const double DENDRO_409 = DENDRO_100 * DENDRO_296;
const double DENDRO_410 = DENDRO_100 * grad_0_gt3[pp];
const double DENDRO_411 = DENDRO_198 * grad_2_gt3[pp];
const double DENDRO_412 = 3.0 * DENDRO_114;
const double DENDRO_413 =
    -DENDRO_130 * DENDRO_229 * (DENDRO_102 + DENDRO_166) -
    DENDRO_130 * DENDRO_246 * (DENDRO_210 + DENDRO_211) -
    DENDRO_181 * (DENDRO_100 * DENDRO_211 + DENDRO_408) -
    DENDRO_181 * (DENDRO_166 * DENDRO_198 + DENDRO_409) -
    DENDRO_257 * DENDRO_76 + DENDRO_278 * grad_1_Gt2[pp] +
    DENDRO_283 * grad2_0_2_gt3[pp] + DENDRO_285 * grad2_0_0_gt3[pp] +
    DENDRO_287 * grad2_1_1_gt3[pp] + DENDRO_289 * grad2_2_2_gt3[pp] -
    DENDRO_292 * DENDRO_410 + DENDRO_407 * grad_1_Gt0[pp] -
    DENDRO_411 * DENDRO_412 - DENDRO_80 * grad2_0_1_gt3[pp] -
    DENDRO_82 * grad2_1_2_gt3[pp] + 4 * grad_1_Gt1[pp] * gt3[pp];
const double DENDRO_414 = DENDRO_223 * grad_1_gt0[pp];
const double DENDRO_415 = DENDRO_219 * grad_2_gt0[pp];
const double DENDRO_416 = DENDRO_133 * DENDRO_219;
const double DENDRO_417 = DENDRO_137 * DENDRO_223;
const double DENDRO_418 = DENDRO_223 * grad_0_gt0[pp];
const double DENDRO_419 = DENDRO_187 * grad_1_gt0[pp];
const double DENDRO_420 = DENDRO_219 * grad_0_gt0[pp];
const double DENDRO_421 = DENDRO_187 * grad_2_gt0[pp];
const double DENDRO_422 = DENDRO_190 * grad_0_gt5[pp];
const double DENDRO_423 = 0.25 * DENDRO_418;
const double DENDRO_424 = 0.25 * DENDRO_420;
const double DENDRO_425 = DENDRO_110 * DENDRO_174;
const double DENDRO_426 = DENDRO_190 * DENDRO_87;
const double DENDRO_427 = DENDRO_174 * DENDRO_183;
const double DENDRO_428 = DENDRO_176 * DENDRO_183;
const double DENDRO_429 = DENDRO_346 * DENDRO_53;
const double DENDRO_430 = DENDRO_344 * DENDRO_52;
const double DENDRO_431 = DENDRO_349 * DENDRO_49 + DENDRO_54;
const double DENDRO_432 =
    DENDRO_23 * DENDRO_49 + DENDRO_344 * DENDRO_53 + DENDRO_351 * DENDRO_52;
const double DENDRO_433 = DENDRO_198 * grad_1_gt0[pp];
const double DENDRO_434 = DENDRO_117 * grad_0_gt3[pp];
const double DENDRO_435 = DENDRO_100 * DENDRO_87;
const double DENDRO_436 = DENDRO_434 + DENDRO_435;
const double DENDRO_437 = DENDRO_238 * grad_2_gt0[pp];
const double DENDRO_438 = DENDRO_174 * grad_0_gt5[pp];
const double DENDRO_439 = DENDRO_177 + DENDRO_438;
const double DENDRO_440 = DENDRO_110 * DENDRO_219;
const double DENDRO_441 = -DENDRO_173 * DENDRO_242;
const double DENDRO_442 = DENDRO_152 * DENDRO_205 + 0.25 * DENDRO_324;
const double DENDRO_443 = DENDRO_223 * DENDRO_87;
const double DENDRO_444 = 0.25 * DENDRO_443;
const double DENDRO_445 = DENDRO_100 * DENDRO_386;
const double DENDRO_446 = -DENDRO_198 * DENDRO_384;
const double DENDRO_447 = DENDRO_445 + DENDRO_446;
const double DENDRO_448 = DENDRO_187 * DENDRO_52;
const double DENDRO_449 = 0.25 * DENDRO_159 + DENDRO_319 * DENDRO_55;
const double DENDRO_450 = DENDRO_137 * DENDRO_176;
const double DENDRO_451 = 0.25 * DENDRO_117;
const double DENDRO_452 = DENDRO_151 * DENDRO_451 + DENDRO_201 * DENDRO_55;
const double DENDRO_453 = DENDRO_190 * DENDRO_199;
const double DENDRO_454 = -DENDRO_176 * DENDRO_337;
const double DENDRO_455 = DENDRO_137 * DENDRO_219;
const double DENDRO_456 = DENDRO_225 * DENDRO_53;
const double DENDRO_457 = DENDRO_183 * DENDRO_219 + DENDRO_456;
const double DENDRO_458 = DENDRO_451 * DENDRO_92;
const double DENDRO_459 = DENDRO_205 * DENDRO_48 + DENDRO_451 * DENDRO_87;
const double DENDRO_460 = DENDRO_167 * DENDRO_205;
const double DENDRO_461 = DENDRO_225 * DENDRO_48;
const double DENDRO_462 = DENDRO_110 * DENDRO_223 + DENDRO_461;
const double DENDRO_463 = DENDRO_152 * DENDRO_242;
const double DENDRO_464 = DENDRO_306 + DENDRO_308;
const double DENDRO_465 = DENDRO_151 * DENDRO_198;
const double DENDRO_466 = 0.25 * DENDRO_465;
const double DENDRO_467 = DENDRO_166 * DENDRO_205;
const double DENDRO_468 = DENDRO_100 * DENDRO_307 + DENDRO_467;
const double DENDRO_469 = 0.25 * DENDRO_87;
const double DENDRO_470 = DENDRO_219 * DENDRO_469 + DENDRO_461;
const double DENDRO_471 = DENDRO_238 * DENDRO_337;
const double DENDRO_472 = -DENDRO_471;
const double DENDRO_473 = 0.25 * DENDRO_174;
const double DENDRO_474 = DENDRO_473 * grad_2_gt5[pp];
const double DENDRO_475 = DENDRO_242 * DENDRO_319 + DENDRO_474;
const double DENDRO_476 = DENDRO_179 * DENDRO_198;
const double DENDRO_477 = -0.5 * DENDRO_168 + DENDRO_212 * DENDRO_55;
const double DENDRO_478 = 0.25 * DENDRO_222;
const double DENDRO_479 = DENDRO_478 * grad_0_gt0[pp];
const double DENDRO_480 = DENDRO_183 * DENDRO_223;
const double DENDRO_481 = DENDRO_479 + DENDRO_480;
const double DENDRO_482 = DENDRO_187 * DENDRO_319 + DENDRO_416;
const double DENDRO_483 = DENDRO_190 * DENDRO_201;
const double DENDRO_484 = 0.25 * DENDRO_151;
const double DENDRO_485 = DENDRO_176 * DENDRO_484;
const double DENDRO_486 = DENDRO_151 * DENDRO_238;
const double DENDRO_487 = DENDRO_174 * grad_1_gt5[pp];
const double DENDRO_488 = DENDRO_382 + DENDRO_487;
const double DENDRO_489 = 1.0 * DENDRO_104;
const double DENDRO_490 = -grad2_0_2_chi[pp];
const double DENDRO_491 = DENDRO_344 * grad_0_gt5[pp];
const double DENDRO_492 = DENDRO_346 * grad_2_gt0[pp];
const double DENDRO_493 = 0.5 * DENDRO_188;
const double DENDRO_494 = DENDRO_116 + DENDRO_349 * DENDRO_92;
const double DENDRO_495 = 0.5 * DENDRO_189;
const double DENDRO_496 = DENDRO_351 * grad_0_gt5[pp];
const double DENDRO_497 = DENDRO_344 * grad_2_gt0[pp];
const double DENDRO_498 = DENDRO_23 * DENDRO_92;
const double DENDRO_499 = 0.5 * DENDRO_191;
const double DENDRO_500 = 2.0 * gt2[pp];
const double DENDRO_501 = 2.0 * gt4[pp];
const double DENDRO_502 = 2.0 * gt5[pp];
const double DENDRO_503 = 2.0 * grad_2_Gt0[pp];
const double DENDRO_504 = 2.0 * grad_2_Gt1[pp];
const double DENDRO_505 = DENDRO_76 * grad_2_chi[pp];
const double DENDRO_506 = DENDRO_505 * grad_0_chi[pp];
const double DENDRO_507 = DENDRO_80 * grad2_0_1_gt2[pp];
const double DENDRO_508 = DENDRO_82 * grad2_1_2_gt2[pp];
const double DENDRO_509 = DENDRO_31 * gt2[pp];
const double DENDRO_510 =
    -DENDRO_192 *
        (DENDRO_490 + DENDRO_493 * (DENDRO_491 + DENDRO_492 + DENDRO_93) +
         DENDRO_494 * DENDRO_495 +
         DENDRO_499 * (DENDRO_496 + DENDRO_497 + DENDRO_498)) +
    DENDRO_283 * grad2_0_2_gt2[pp] + DENDRO_285 * grad2_0_0_gt2[pp] +
    DENDRO_287 * grad2_1_1_gt2[pp] + DENDRO_289 * grad2_2_2_gt2[pp] +
    DENDRO_356 * grad_1_gt2[pp] + DENDRO_358 * grad_0_gt2[pp] +
    DENDRO_360 * grad_2_gt2[pp] + DENDRO_362 * DENDRO_509 +
    DENDRO_500 * grad_0_Gt0[pp] + DENDRO_500 * grad_2_Gt2[pp] +
    DENDRO_501 * grad_0_Gt1[pp] + DENDRO_502 * grad_0_Gt2[pp] +
    DENDRO_503 * gt0[pp] + DENDRO_504 * gt1[pp] - DENDRO_506 - DENDRO_507 -
    DENDRO_508;
const double DENDRO_511 =
    -DENDRO_21 * DENDRO_31 * DENDRO_92 + DENDRO_31 * DENDRO_90 +
    DENDRO_31 * DENDRO_91 +
    DENDRO_45 * (-DENDRO_509 * DENDRO_72 + grad_2_chi[pp]);
const double DENDRO_512 = 2.0 * grad_0_alpha[pp];
const double DENDRO_513 =
    DENDRO_107 * DENDRO_31 + DENDRO_108 * DENDRO_31 -
    DENDRO_23 * DENDRO_31 * DENDRO_92 +
    DENDRO_45 * (-DENDRO_509 * DENDRO_60 + grad_0_chi[pp]);
const double DENDRO_514 = 2.0 * grad_2_alpha[pp];
const double DENDRO_515 = 2.0 * grad_1_alpha[pp];
const double DENDRO_516 = DENDRO_45 * gt2[pp];
const double DENDRO_517 = DENDRO_31 * (DENDRO_117 + DENDRO_43 * DENDRO_516);
const double DENDRO_518 = -DENDRO_511 * DENDRO_512 - DENDRO_513 * DENDRO_514 +
                          DENDRO_515 * DENDRO_517 - 4 * grad2_0_2_alpha[pp];
const double DENDRO_519 = DENDRO_222 * grad_2_gt0[pp];
const double DENDRO_520 = DENDRO_223 * grad_0_gt5[pp] + DENDRO_519;
const double DENDRO_521 = DENDRO_117 * grad_2_gt3[pp];
const double DENDRO_522 = DENDRO_100 * grad_1_gt5[pp] + DENDRO_521;
const double DENDRO_523 = DENDRO_465 + DENDRO_522;
const double DENDRO_524 = DENDRO_115 * (-DENDRO_336 + DENDRO_442);
const double DENDRO_525 = 0.25 * DENDRO_486;
const double DENDRO_526 = DENDRO_105 * (DENDRO_409 + DENDRO_447);
const double DENDRO_527 = DENDRO_299 * (0.5 * DENDRO_160 + DENDRO_449);
const double DENDRO_528 = DENDRO_181 * (-DENDRO_204 * DENDRO_55 + DENDRO_459);
const double DENDRO_529 = DENDRO_110 * DENDRO_176 + DENDRO_453;
const double DENDRO_530 = DENDRO_181 * (DENDRO_452 + DENDRO_458);
const double DENDRO_531 = DENDRO_100 * DENDRO_335;
const double DENDRO_532 =
    -0.5 * DENDRO_117 * DENDRO_212 + DENDRO_460 + DENDRO_531;
const double DENDRO_533 = DENDRO_308 + DENDRO_474;
const double DENDRO_534 = -DENDRO_223 * DENDRO_337;
const double DENDRO_535 = DENDRO_170 + DENDRO_477;
const double DENDRO_536 = DENDRO_152 * DENDRO_187;
const double DENDRO_537 = DENDRO_137 * DENDRO_238;
const double DENDRO_538 = DENDRO_425 + DENDRO_483;
const double DENDRO_539 = 0.25 * DENDRO_435;
const double DENDRO_540 = DENDRO_211 * DENDRO_55;
const double DENDRO_541 = DENDRO_133 * DENDRO_198 + DENDRO_540;
const double DENDRO_542 = DENDRO_539 + DENDRO_541;
const double DENDRO_543 = DENDRO_222 * grad_1_gt0[pp];
const double DENDRO_544 = DENDRO_443 + DENDRO_543;
const double DENDRO_545 = DENDRO_219 * grad_0_gt3[pp];
const double DENDRO_546 = DENDRO_176 * grad_2_gt3[pp];
const double DENDRO_547 = 0.25 * DENDRO_309;
const double DENDRO_548 = -DENDRO_204 * DENDRO_242;
const double DENDRO_549 = DENDRO_110 * DENDRO_222 + DENDRO_225 * DENDRO_379;
const double DENDRO_550 = DENDRO_212 * DENDRO_214;
const double DENDRO_551 = DENDRO_228 * DENDRO_319;
const double DENDRO_552 = DENDRO_101 * DENDRO_222 + DENDRO_551;
const double DENDRO_553 = DENDRO_238 * DENDRO_296;
const double DENDRO_554 = DENDRO_179 * DENDRO_219;
const double DENDRO_555 = DENDRO_480 + DENDRO_554;
const double DENDRO_556 = DENDRO_225 * DENDRO_49 + 0.5 * DENDRO_332;
const double DENDRO_557 = DENDRO_219 * DENDRO_92;
const double DENDRO_558 = 0.25 * DENDRO_557;
const double DENDRO_559 = DENDRO_242 * DENDRO_379;
const double DENDRO_560 = DENDRO_198 * DENDRO_469 + DENDRO_467;
const double DENDRO_561 = DENDRO_176 * DENDRO_335;
const double DENDRO_562 = -DENDRO_561;
const double DENDRO_563 = DENDRO_199 * DENDRO_245;
const double DENDRO_564 = -DENDRO_238 * DENDRO_335 + DENDRO_563;
const double DENDRO_565 = DENDRO_172 * DENDRO_228 + DENDRO_478 * DENDRO_92;
const double DENDRO_566 = DENDRO_205 * DENDRO_207;
const double DENDRO_567 = DENDRO_198 * DENDRO_296 + DENDRO_566;
const double DENDRO_568 = DENDRO_198 * DENDRO_386;
const double DENDRO_569 = DENDRO_151 * DENDRO_478;
const double DENDRO_570 = DENDRO_166 * DENDRO_225 + DENDRO_222 * DENDRO_469;
const double DENDRO_571 = -DENDRO_219 * DENDRO_384;
const double DENDRO_572 = DENDRO_228 * DENDRO_52 + 0.5 * DENDRO_375;
const double DENDRO_573 = DENDRO_451 * grad_1_gt3[pp];
const double DENDRO_574 = DENDRO_408 + DENDRO_573;
const double DENDRO_575 = DENDRO_214 * DENDRO_319;
const double DENDRO_576 = 0.25 * DENDRO_92;
const double DENDRO_577 = DENDRO_172 * DENDRO_245;
const double DENDRO_578 = DENDRO_238 * DENDRO_576 + DENDRO_577;
const double DENDRO_579 = DENDRO_176 * DENDRO_92;
const double DENDRO_580 = 1.0 * DENDRO_291;
const double DENDRO_581 = -grad2_1_2_chi[pp];
const double DENDRO_582 = DENDRO_151 * DENDRO_346 + DENDRO_21 * grad_2_gt3[pp] +
                          DENDRO_344 * grad_1_gt5[pp];
const double DENDRO_583 = DENDRO_349 * grad_2_gt3[pp];
const double DENDRO_584 = DENDRO_351 * grad_1_gt5[pp];
const double DENDRO_585 = DENDRO_23 * grad_2_gt3[pp];
const double DENDRO_586 = DENDRO_151 * DENDRO_344;
const double DENDRO_587 = DENDRO_31 * gt4[pp];
const double DENDRO_588 =
    DENDRO_283 * grad2_0_2_gt4[pp] + DENDRO_285 * grad2_0_0_gt4[pp] +
    DENDRO_287 * grad2_1_1_gt4[pp] + DENDRO_289 * grad2_2_2_gt4[pp] +
    DENDRO_500 * grad_1_Gt0[pp] + DENDRO_501 * grad_1_Gt1[pp] +
    DENDRO_501 * grad_2_Gt2[pp] + DENDRO_502 * grad_1_Gt2[pp] +
    DENDRO_503 * gt1[pp] + DENDRO_504 * gt3[pp] - DENDRO_505 * grad_1_chi[pp] -
    DENDRO_80 * grad2_0_1_gt4[pp] - DENDRO_82 * grad2_1_2_gt4[pp];
const double DENDRO_589 =
    -DENDRO_192 *
        (DENDRO_493 * DENDRO_582 + DENDRO_495 * (DENDRO_197 + DENDRO_583) +
         DENDRO_499 * (DENDRO_584 + DENDRO_585 + DENDRO_586) + DENDRO_581) +
    DENDRO_356 * grad_1_gt4[pp] + DENDRO_358 * grad_0_gt4[pp] +
    DENDRO_360 * grad_2_gt4[pp] + DENDRO_362 * DENDRO_587 + DENDRO_588;
const double DENDRO_590 = DENDRO_195 * DENDRO_31 + DENDRO_196 * DENDRO_31;
const double DENDRO_591 =
    -DENDRO_194 * DENDRO_31 -
    DENDRO_45 * (-DENDRO_43 * DENDRO_587 + grad_2_chi[pp]) + DENDRO_590;
const double DENDRO_592 =
    -DENDRO_23 * DENDRO_31 * grad_2_gt3[pp] + DENDRO_235 * DENDRO_31 +
    DENDRO_236 * DENDRO_31 +
    DENDRO_45 * (-DENDRO_587 * DENDRO_60 + grad_1_chi[pp]);
const double DENDRO_593 = DENDRO_221 - DENDRO_45 * DENDRO_72 * gt4[pp];
const double DENDRO_594 = -DENDRO_31 * DENDRO_512 * DENDRO_593 -
                          DENDRO_514 * DENDRO_592 + DENDRO_515 * DENDRO_591 -
                          4 * grad2_1_2_alpha[pp];
const double DENDRO_595 = 1.0 * DENDRO_400;
const double DENDRO_596 = 0.5 * DENDRO_399;
const double DENDRO_597 = 0.25 * DENDRO_579;
const double DENDRO_598 = DENDRO_306 + DENDRO_474;
const double DENDRO_599 = DENDRO_173 * DENDRO_228;
const double DENDRO_600 = DENDRO_566 + DENDRO_568;
const double DENDRO_601 = DENDRO_238 * DENDRO_307;
const double DENDRO_602 = DENDRO_408 + DENDRO_409;
const double DENDRO_603 = DENDRO_228 * DENDRO_51;
const double DENDRO_604 = DENDRO_101 * DENDRO_219 + DENDRO_603;
const double DENDRO_605 = DENDRO_214 * DENDRO_379 + DENDRO_573;
const double DENDRO_606 = DENDRO_176 * DENDRO_296 + DENDRO_577;
const double DENDRO_607 = 1.0 * DENDRO_161;
const double DENDRO_608 =
    -DENDRO_115 * (DENDRO_205 * DENDRO_211 + DENDRO_339 + DENDRO_547) -
    DENDRO_181 * (-DENDRO_531 + DENDRO_560) -
    DENDRO_580 * (DENDRO_169 + DENDRO_436) -
    DENDRO_607 * (DENDRO_328 + DENDRO_522);
const double DENDRO_609 = DENDRO_471 + DENDRO_561;
const double DENDRO_610 = DENDRO_101 * DENDRO_223;
const double DENDRO_611 = -DENDRO_167 * DENDRO_214;
const double DENDRO_612 = DENDRO_152 * DENDRO_245 + DENDRO_174 * DENDRO_296;
const double DENDRO_613 = DENDRO_187 * DENDRO_49;
const double DENDRO_614 = 0.25 * DENDRO_145;
const double DENDRO_615 = DENDRO_137 * DENDRO_174 + DENDRO_190 * DENDRO_379;
const double DENDRO_616 = 0.5 * DENDRO_175 + DENDRO_190 * DENDRO_204;
const double DENDRO_617 = DENDRO_479 + DENDRO_554;
const double DENDRO_618 = DENDRO_187 * DENDRO_379 + DENDRO_417;
const double DENDRO_619 = DENDRO_100 * DENDRO_484 + DENDRO_540;
const double DENDRO_620 = DENDRO_173 * DENDRO_245 + 0.5 * DENDRO_381;
const double DENDRO_621 = DENDRO_152 * DENDRO_214;
const double DENDRO_622 = DENDRO_223 * DENDRO_576 + DENDRO_603;
const double DENDRO_623 = -DENDRO_100 * DENDRO_384;
const double DENDRO_624 = DENDRO_207 * DENDRO_55;
const double DENDRO_625 = DENDRO_151 * DENDRO_473 + DENDRO_190 * DENDRO_211;
const double DENDRO_626 = DENDRO_133 * DENDRO_223;
const double DENDRO_627 = DENDRO_228 * DENDRO_53;
const double DENDRO_628 = DENDRO_179 * DENDRO_223 + DENDRO_627;
const double DENDRO_629 = DENDRO_473 * DENDRO_87;
const double DENDRO_630 = DENDRO_245 * DENDRO_51 + DENDRO_473 * DENDRO_92;
const double DENDRO_631 = 1.0 * DENDRO_114;
const double DENDRO_632 = -grad2_0_1_chi[pp];
const double DENDRO_633 = DENDRO_346 * grad_1_gt0[pp];
const double DENDRO_634 = DENDRO_344 * DENDRO_87;
const double DENDRO_635 = DENDRO_349 * grad_0_gt3[pp];
const double DENDRO_636 = DENDRO_23 * grad_0_gt3[pp];
const double DENDRO_637 =
    DENDRO_344 * grad_1_gt0[pp] + DENDRO_351 * DENDRO_87 + DENDRO_636;
const double DENDRO_638 = 2.0 * gt1[pp];
const double DENDRO_639 = DENDRO_638 * grad_0_Gt0[pp];
const double DENDRO_640 = 2.0 * grad_0_Gt1[pp] * gt3[pp];
const double DENDRO_641 = DENDRO_501 * grad_0_Gt2[pp];
const double DENDRO_642 = 2.0 * grad_1_Gt0[pp] * gt0[pp];
const double DENDRO_643 = DENDRO_638 * grad_1_Gt1[pp];
const double DENDRO_644 = DENDRO_500 * grad_1_Gt2[pp];
const double DENDRO_645 = -DENDRO_252 * DENDRO_76;
const double DENDRO_646 = DENDRO_283 * grad2_0_2_gt1[pp];
const double DENDRO_647 = DENDRO_285 * grad2_0_0_gt1[pp];
const double DENDRO_648 = DENDRO_287 * grad2_1_1_gt1[pp];
const double DENDRO_649 = DENDRO_289 * grad2_2_2_gt1[pp];
const double DENDRO_650 = -DENDRO_80 * grad2_0_1_gt1[pp];
const double DENDRO_651 = -DENDRO_82 * grad2_1_2_gt1[pp];
const double DENDRO_652 = DENDRO_31 * gt1[pp];
const double DENDRO_653 =
    -DENDRO_192 * (DENDRO_493 * (DENDRO_633 + DENDRO_634 + DENDRO_85) +
                   DENDRO_495 * (DENDRO_635 + DENDRO_99) +
                   DENDRO_499 * DENDRO_637 + DENDRO_632) +
    DENDRO_356 * grad_1_gt1[pp] + DENDRO_358 * grad_0_gt1[pp] +
    DENDRO_360 * grad_2_gt1[pp] + DENDRO_362 * DENDRO_652 + DENDRO_639 +
    DENDRO_640 + DENDRO_641 + DENDRO_642 + DENDRO_643 + DENDRO_644 +
    DENDRO_645 + DENDRO_646 + DENDRO_647 + DENDRO_648 + DENDRO_649 +
    DENDRO_650 + DENDRO_651;
const double DENDRO_654 =
    -DENDRO_21 * DENDRO_31 * grad_0_gt3[pp] + DENDRO_31 * DENDRO_86 +
    DENDRO_31 * DENDRO_88 +
    DENDRO_45 * (-DENDRO_652 * DENDRO_72 + grad_1_chi[pp]);
const double DENDRO_655 =
    -DENDRO_21 * DENDRO_31 * grad_1_gt0[pp] -
    DENDRO_23 * DENDRO_31 * DENDRO_87 + DENDRO_31 * DENDRO_96 +
    DENDRO_45 * (-DENDRO_43 * DENDRO_652 + grad_0_chi[pp]);
const double DENDRO_656 = DENDRO_45 * gt1[pp];
const double DENDRO_657 =
    DENDRO_31 * DENDRO_514 *
        (-DENDRO_126 - DENDRO_127 + DENDRO_60 * DENDRO_656 + DENDRO_636) -
    DENDRO_512 * DENDRO_654 - DENDRO_515 * DENDRO_655 - 4 * grad2_0_1_alpha[pp];
const double DENDRO_658 = -0.25 * DENDRO_176 * grad_1_gt5[pp] + DENDRO_609;
const double DENDRO_659 = 0.5 * DENDRO_395;
const double DENDRO_660 = DENDRO_228 * DENDRO_48;
const double DENDRO_661 = DENDRO_178 + DENDRO_616;
const double DENDRO_662 = -0.5 * DENDRO_176 * DENDRO_212 + DENDRO_620;
const double DENDRO_663 = DENDRO_409 + DENDRO_573;
const double DENDRO_664 = DENDRO_100 * DENDRO_101 + DENDRO_624;
const double DENDRO_665 = -DENDRO_115 * (DENDRO_117 * DENDRO_211 + DENDRO_466) +
                          DENDRO_136 * (DENDRO_446 + DENDRO_663) +
                          DENDRO_157 * (DENDRO_623 + DENDRO_664) -
                          DENDRO_181 * (DENDRO_150 + DENDRO_541) -
                          DENDRO_181 * (DENDRO_150 + DENDRO_619) -
                          DENDRO_299 * (1.0 * DENDRO_146 + DENDRO_614);
const double DENDRO_666 =
    -DENDRO_21 *
        (DENDRO_657 +
         alpha[pp] *
             (-DENDRO_105 * (DENDRO_392 + DENDRO_612) -
              DENDRO_105 * (DENDRO_611 + DENDRO_659) -
              DENDRO_105 * (DENDRO_393 + DENDRO_610 + DENDRO_660) +
              DENDRO_115 * DENDRO_658 +
              DENDRO_135 * (DENDRO_376 + DENDRO_543 + DENDRO_545) +
              DENDRO_135 * (DENDRO_486 + DENDRO_487 + DENDRO_546) -
              DENDRO_136 * DENDRO_662 + DENDRO_136 * (DENDRO_571 + DENDRO_622) +
              DENDRO_136 * (DENDRO_621 + DENDRO_663) +
              DENDRO_142 * (DENDRO_187 * grad_0_gt3[pp] + DENDRO_414) +
              DENDRO_157 * (DENDRO_625 + DENDRO_629) +
              DENDRO_157 * (-DENDRO_167 * DENDRO_187 + DENDRO_628) +
              DENDRO_157 * (DENDRO_190 * DENDRO_212 + DENDRO_630) +
              DENDRO_157 * (DENDRO_214 * DENDRO_48 + DENDRO_664) +
              DENDRO_181 * DENDRO_661 - DENDRO_181 * (DENDRO_416 + DENDRO_618) -
              DENDRO_181 * (DENDRO_536 + DENDRO_617) -
              DENDRO_181 * (DENDRO_483 + DENDRO_537 + DENDRO_597) -
              DENDRO_299 * (0.5 * DENDRO_426 + DENDRO_615) -
              DENDRO_299 * (DENDRO_187 * DENDRO_48 + DENDRO_423 + DENDRO_613) -
              DENDRO_631 * (DENDRO_333 + DENDRO_519 + DENDRO_557) + DENDRO_653 +
              DENDRO_665)) -
    DENDRO_21 *
        (DENDRO_657 +
         alpha[pp] *
             (-DENDRO_105 * (1.0 * DENDRO_398 + DENDRO_610) -
              DENDRO_105 * (0.5 * DENDRO_401 + DENDRO_612) -
              DENDRO_105 * (DENDRO_166 * DENDRO_214 + DENDRO_396 + DENDRO_611) -
              DENDRO_115 * (DENDRO_222 * DENDRO_51 + DENDRO_558) -
              DENDRO_115 * (0.25 * DENDRO_238 * grad_0_gt5[pp] - DENDRO_609) +
              DENDRO_136 * (DENDRO_378 + DENDRO_604) +
              DENDRO_136 * (DENDRO_378 + DENDRO_622) +
              DENDRO_136 * (-DENDRO_383 - DENDRO_620) +
              DENDRO_136 * (DENDRO_446 + DENDRO_605) +
              DENDRO_136 * (DENDRO_525 + DENDRO_606) +
              DENDRO_136 * (DENDRO_602 + DENDRO_621) +
              DENDRO_142 * (DENDRO_214 * grad_1_gt0[pp] + DENDRO_410) +
              DENDRO_157 * (DENDRO_626 + DENDRO_628) +
              DENDRO_157 * (DENDRO_629 + DENDRO_630) +
              DENDRO_157 * (DENDRO_245 * DENDRO_52 + DENDRO_625) +
              DENDRO_157 * (DENDRO_166 * DENDRO_187 + DENDRO_626 + DENDRO_627) +
              DENDRO_157 * (DENDRO_214 * DENDRO_49 + DENDRO_623 + DENDRO_624) -
              DENDRO_181 * (DENDRO_417 + DENDRO_617) -
              DENDRO_181 * (DENDRO_476 + DENDRO_619) -
              DENDRO_181 * (DENDRO_479 + DENDRO_618) -
              DENDRO_181 * (0.5 * DENDRO_238 * DENDRO_52 - DENDRO_616) -
              DENDRO_299 * (0.5 * DENDRO_418 + DENDRO_613) -
              DENDRO_299 * (DENDRO_427 + DENDRO_615) -
              DENDRO_299 * (DENDRO_166 * DENDRO_55 + DENDRO_184 + DENDRO_614) -
              DENDRO_607 * (DENDRO_169 + DENDRO_433 + DENDRO_434) -
              DENDRO_607 * (DENDRO_437 + DENDRO_438 + DENDRO_579) -
              DENDRO_631 * (DENDRO_328 + DENDRO_465 + DENDRO_521) +
              DENDRO_653)) -
    DENDRO_23 *
        (DENDRO_594 +
         alpha[pp] *
             (-DENDRO_105 * (DENDRO_552 + DENDRO_596) -
              DENDRO_105 * (DENDRO_553 + DENDRO_595) -
              DENDRO_105 * (DENDRO_211 * DENDRO_214 + DENDRO_390 + DENDRO_550) -
              DENDRO_115 * (0.5 * DENDRO_316 + DENDRO_548) -
              DENDRO_115 * (-DENDRO_338 + DENDRO_549) +
              DENDRO_136 * (DENDRO_564 + DENDRO_601) +
              DENDRO_136 * (DENDRO_565 + DENDRO_569) +
              DENDRO_136 * (DENDRO_570 - DENDRO_599) +
              DENDRO_136 * (-DENDRO_204 * DENDRO_214 + DENDRO_600) +
              DENDRO_136 * (DENDRO_211 * DENDRO_242 + DENDRO_563 + DENDRO_601) +
              DENDRO_157 * (DENDRO_377 + DENDRO_572) +
              DENDRO_157 * (DENDRO_380 + DENDRO_578) +
              DENDRO_157 * (DENDRO_380 + DENDRO_606) +
              DENDRO_157 * (DENDRO_444 + DENDRO_604) +
              DENDRO_157 * (DENDRO_445 + DENDRO_605) +
              DENDRO_157 * (DENDRO_575 + DENDRO_602) -
              DENDRO_181 * (DENDRO_534 + DENDRO_556) -
              DENDRO_181 * (DENDRO_559 + DENDRO_598) -
              DENDRO_181 * (DENDRO_562 + DENDRO_598) -
              DENDRO_299 * (DENDRO_417 + DENDRO_555) -
              DENDRO_299 * (DENDRO_172 * DENDRO_174 + DENDRO_597) +
              DENDRO_312 * (DENDRO_214 * grad_1_gt5[pp] + DENDRO_411) +
              DENDRO_589 - DENDRO_607 * (DENDRO_520 + DENDRO_557) +
              DENDRO_608)) -
    DENDRO_23 *
        (DENDRO_594 +
         alpha[pp] *
             (-DENDRO_105 * (-DENDRO_385 + DENDRO_552) -
              DENDRO_105 * (0.5 * DENDRO_389 + DENDRO_550) -
              DENDRO_105 * (DENDRO_201 * DENDRO_245 + DENDRO_387 + DENDRO_553) -
              DENDRO_115 * (1.0 * DENDRO_310 + DENDRO_547) -
              DENDRO_115 * (0.5 * DENDRO_331 + DENDRO_549) -
              DENDRO_115 * (DENDRO_201 * DENDRO_242 + DENDRO_322 + DENDRO_548) +
              DENDRO_136 * (DENDRO_567 + DENDRO_568) +
              DENDRO_136 * (DENDRO_569 + DENDRO_570) +
              DENDRO_136 * (-DENDRO_167 * DENDRO_225 + DENDRO_565) +
              DENDRO_136 * (DENDRO_201 * DENDRO_214 + DENDRO_567) +
              DENDRO_136 * (DENDRO_212 * DENDRO_242 + DENDRO_564) +
              DENDRO_141 * (DENDRO_488 + DENDRO_546) +
              DENDRO_141 * (DENDRO_544 + DENDRO_545) +
              DENDRO_157 * (DENDRO_445 + DENDRO_574) +
              DENDRO_157 * (DENDRO_571 + DENDRO_572) +
              DENDRO_157 * (DENDRO_574 + DENDRO_575) +
              DENDRO_157 * (DENDRO_176 * DENDRO_386 + DENDRO_578) -
              DENDRO_181 * (DENDRO_318 + DENDRO_468) -
              DENDRO_181 * (DENDRO_318 + DENDRO_560) -
              DENDRO_181 * (DENDRO_334 + DENDRO_556) -
              DENDRO_181 * (DENDRO_462 + DENDRO_558) -
              DENDRO_181 * (DENDRO_464 + DENDRO_559) -
              DENDRO_181 * (DENDRO_475 + DENDRO_562) -
              DENDRO_299 * (DENDRO_416 + DENDRO_555) -
              DENDRO_299 * (DENDRO_117 * DENDRO_166 + DENDRO_539) +
              DENDRO_312 * (DENDRO_242 * grad_2_gt3[pp] + DENDRO_293) -
              DENDRO_580 * (DENDRO_439 + DENDRO_579) + DENDRO_589)) +
    DENDRO_26 *
        (-DENDRO_267 * DENDRO_367 +
         DENDRO_271 * (DENDRO_228 + DENDRO_363 * DENDRO_73) +
         DENDRO_364 * (DENDRO_245 + DENDRO_274 * DENDRO_363) +
         alpha[pp] *
             (DENDRO_115 * DENDRO_238 * DENDRO_402 +
              DENDRO_136 * (DENDRO_387 - DENDRO_388) +
              DENDRO_136 * (DENDRO_390 + 1.0 * DENDRO_391) +
              DENDRO_136 * (DENDRO_228 * DENDRO_92 - DENDRO_385) +
              DENDRO_142 * (DENDRO_395 + DENDRO_397) +
              DENDRO_142 * (DENDRO_174 * grad_2_gt3[pp] + DENDRO_401) +
              DENDRO_142 * (DENDRO_223 * grad_0_gt3[pp] + DENDRO_398) +
              DENDRO_157 * (DENDRO_393 + DENDRO_394) +
              DENDRO_157 * (DENDRO_396 + 1.0 * DENDRO_397) +
              DENDRO_157 * (DENDRO_245 * DENDRO_92 + DENDRO_392) -
              DENDRO_174 * DENDRO_406 - DENDRO_181 * (DENDRO_375 + DENDRO_377) -
              DENDRO_181 * (-1.0 * DENDRO_381 - DENDRO_383) -
              DENDRO_181 * (DENDRO_223 * DENDRO_379 + DENDRO_378) -
              DENDRO_181 * (DENDRO_238 * DENDRO_379 + DENDRO_380) -
              DENDRO_192 *
                  (DENDRO_188 * DENDRO_370 +
                   DENDRO_189 * (DENDRO_371 + DENDRO_372 + DENDRO_373) +
                   DENDRO_191 * DENDRO_374 + DENDRO_368) -
              DENDRO_215 * DENDRO_295 * grad_1_gt3[pp] -
              DENDRO_222 * DENDRO_404 - DENDRO_223 * DENDRO_405 +
              DENDRO_312 * (DENDRO_389 + DENDRO_391) +
              DENDRO_312 * (DENDRO_222 * grad_0_gt3[pp] + DENDRO_399) +
              DENDRO_312 * (DENDRO_238 * grad_2_gt3[pp] + DENDRO_400) +
              DENDRO_356 * grad_1_gt3[pp] + DENDRO_358 * grad_0_gt3[pp] +
              DENDRO_360 * grad_2_gt3[pp] + DENDRO_362 * DENDRO_366 +
              DENDRO_413) -
         4 * grad2_1_1_alpha[pp]) +
    DENDRO_33 *
        (DENDRO_518 +
         alpha[pp] *
             (-DENDRO_105 * (DENDRO_408 + DENDRO_447) -
              DENDRO_105 * (DENDRO_222 * DENDRO_48 + DENDRO_444) -
              DENDRO_115 * (1.0 * DENDRO_313 + DENDRO_440) -
              DENDRO_115 * (0.5 * DENDRO_325 + DENDRO_442) -
              DENDRO_115 * (DENDRO_172 * DENDRO_242 + DENDRO_321 + DENDRO_441) +
              DENDRO_136 * (DENDRO_323 + DENDRO_462) +
              DENDRO_136 * (DENDRO_323 + DENDRO_470) +
              DENDRO_136 * (DENDRO_463 + DENDRO_464) +
              DENDRO_136 * (DENDRO_466 + DENDRO_468) +
              DENDRO_136 * (DENDRO_472 + DENDRO_475) +
              DENDRO_136 * (0.5 * DENDRO_327 + DENDRO_329 - DENDRO_460) +
              DENDRO_141 * (DENDRO_433 + DENDRO_436) +
              DENDRO_141 * (DENDRO_437 + DENDRO_439) +
              DENDRO_157 * (DENDRO_416 + DENDRO_481) +
              DENDRO_157 * (DENDRO_476 + DENDRO_477) +
              DENDRO_157 * (DENDRO_479 + DENDRO_482) +
              DENDRO_157 * (DENDRO_183 * DENDRO_238 + DENDRO_483 + DENDRO_485) -
              DENDRO_162 * (DENDRO_242 * grad_2_gt0[pp] + DENDRO_290) -
              DENDRO_181 * (DENDRO_455 + DENDRO_457) -
              DENDRO_181 * (DENDRO_458 + DENDRO_459) -
              DENDRO_181 * (DENDRO_205 * DENDRO_49 + DENDRO_452) -
              DENDRO_181 * (DENDRO_172 * DENDRO_187 + DENDRO_455 + DENDRO_456) -
              DENDRO_181 * (DENDRO_242 * DENDRO_52 + DENDRO_453 + DENDRO_454) -
              DENDRO_299 * (DENDRO_180 + DENDRO_449) -
              DENDRO_299 * (0.5 * DENDRO_420 + DENDRO_448) -
              DENDRO_299 * (DENDRO_172 * DENDRO_190 + DENDRO_428 + DENDRO_450) -
              DENDRO_489 * (DENDRO_486 + DENDRO_488) + DENDRO_510)) +
    DENDRO_33 *
        (DENDRO_518 +
         alpha[pp] *
             (-DENDRO_105 * (DENDRO_174 * DENDRO_201 + DENDRO_525) -
              DENDRO_115 * (0.5 * DENDRO_314 + DENDRO_441) -
              DENDRO_115 * (DENDRO_225 * DENDRO_51 + DENDRO_341 + DENDRO_440) +
              DENDRO_135 * DENDRO_523 + DENDRO_135 * (DENDRO_333 + DENDRO_520) -
              DENDRO_136 * DENDRO_532 + DENDRO_136 * (DENDRO_463 + DENDRO_533) +
              DENDRO_136 * (DENDRO_470 + DENDRO_534) +
              DENDRO_136 * (DENDRO_472 + DENDRO_533) + DENDRO_157 * DENDRO_535 +
              DENDRO_157 * DENDRO_542 + DENDRO_157 * (DENDRO_417 + DENDRO_482) +
              DENDRO_157 * (DENDRO_481 + DENDRO_536) +
              DENDRO_157 * (DENDRO_485 + DENDRO_538) +
              DENDRO_157 * (DENDRO_537 + DENDRO_538) -
              DENDRO_162 * (DENDRO_187 * grad_0_gt5[pp] + DENDRO_415) -
              DENDRO_181 * (DENDRO_454 + DENDRO_529) -
              DENDRO_181 * (-DENDRO_173 * DENDRO_187 + DENDRO_457) -
              DENDRO_181 * (DENDRO_242 * DENDRO_51 + DENDRO_529) -
              DENDRO_299 * (1.0 * DENDRO_422 + DENDRO_450) -
              DENDRO_299 * (DENDRO_187 * DENDRO_51 + DENDRO_424 + DENDRO_448) -
              DENDRO_489 * (DENDRO_376 + DENDRO_544) + DENDRO_510 - DENDRO_524 -
              DENDRO_526 - DENDRO_527 - DENDRO_528 - DENDRO_530)) +
    DENDRO_35 *
        (DENDRO_268 * DENDRO_56 - DENDRO_270 * DENDRO_74 -
         DENDRO_364 * DENDRO_62 +
         alpha[pp] *
             (-DENDRO_105 * DENDRO_125 * DENDRO_174 - DENDRO_106 +
              DENDRO_113 * DENDRO_115 * DENDRO_176 - DENDRO_122 -
              DENDRO_130 * DENDRO_132 * DENDRO_247 - DENDRO_131 +
              DENDRO_136 * DENDRO_153 +
              DENDRO_136 * (-1.0 * DENDRO_168 - DENDRO_171) +
              DENDRO_136 * (-1.0 * DENDRO_175 - DENDRO_178) +
              DENDRO_136 * (DENDRO_152 * DENDRO_176 + DENDRO_425) +
              DENDRO_136 * (DENDRO_219 * DENDRO_48 + DENDRO_417) +
              DENDRO_136 * (DENDRO_223 * DENDRO_51 + DENDRO_416) +
              DENDRO_142 * DENDRO_147 + DENDRO_142 * (DENDRO_418 + DENDRO_419) +
              DENDRO_142 * (DENDRO_174 * grad_2_gt0[pp] + DENDRO_426) +
              DENDRO_157 * DENDRO_185 +
              DENDRO_157 * (1.0 * DENDRO_419 + DENDRO_423) +
              DENDRO_157 * (DENDRO_151 * DENDRO_190 + DENDRO_427) -
              DENDRO_162 * (DENDRO_420 + DENDRO_421) -
              DENDRO_162 * (DENDRO_176 * grad_2_gt0[pp] + DENDRO_422) -
              DENDRO_163 - DENDRO_181 * (1.0 * DENDRO_421 + DENDRO_424) -
              DENDRO_181 * (-2 * DENDRO_173 * DENDRO_190 + DENDRO_428) -
              DENDRO_182 -
              DENDRO_192 * (DENDRO_186 +
                            DENDRO_188 * (DENDRO_429 + DENDRO_430 + DENDRO_67) +
                            DENDRO_189 * DENDRO_431 + DENDRO_191 * DENDRO_432) +
              DENDRO_218 * DENDRO_355 -
              6.0 * DENDRO_230 * DENDRO_84 * grad_0_gt0[pp] +
              DENDRO_233 * DENDRO_357 + DENDRO_250 * DENDRO_359 +
              DENDRO_263 * DENDRO_361 + DENDRO_277 * grad_0_Gt2[pp] +
              DENDRO_283 * grad2_0_2_gt0[pp] + DENDRO_285 * grad2_0_0_gt0[pp] +
              DENDRO_287 * grad2_1_1_gt0[pp] + DENDRO_289 * grad2_2_2_gt0[pp] -
              DENDRO_294 * DENDRO_414 + DENDRO_407 * grad_0_Gt1[pp] -
              DENDRO_412 * DENDRO_415 - DENDRO_78 - DENDRO_81 - DENDRO_83 +
              4 * grad_0_Gt0[pp] * gt0[pp]) -
         4 * grad2_0_0_alpha[pp]) +
    DENDRO_39 *
        (DENDRO_266 * DENDRO_268 - DENDRO_269 * DENDRO_271 -
         DENDRO_275 * DENDRO_276 +
         alpha[pp] *
             (-DENDRO_105 * DENDRO_222 * DENDRO_303 -
              DENDRO_130 * DENDRO_226 * DENDRO_304 + DENDRO_136 * DENDRO_340 +
              DENDRO_136 * (1.0 * DENDRO_317 + DENDRO_322) +
              DENDRO_136 * (DENDRO_225 * DENDRO_87 - DENDRO_338) +
              DENDRO_157 * DENDRO_320 + DENDRO_157 * DENDRO_330 +
              DENDRO_157 * (DENDRO_332 + DENDRO_334) +
              DENDRO_157 * (DENDRO_172 * DENDRO_238 + DENDRO_308) +
              DENDRO_157 * (DENDRO_176 * DENDRO_201 + DENDRO_306) +
              DENDRO_157 * (DENDRO_219 * DENDRO_319 + DENDRO_323) -
              DENDRO_162 * (DENDRO_314 + DENDRO_315) -
              DENDRO_162 * (DENDRO_219 * grad_0_gt5[pp] + DENDRO_313) -
              DENDRO_181 * (1.0 * DENDRO_315 + DENDRO_321) -
              DENDRO_181 * (DENDRO_205 * DENDRO_87 - DENDRO_336) -
              DENDRO_181 * (DENDRO_225 * DENDRO_342 + DENDRO_341) -
              DENDRO_192 *
                  (DENDRO_188 * DENDRO_348 + DENDRO_189 * DENDRO_350 +
                   DENDRO_191 * (DENDRO_352 + DENDRO_353 + DENDRO_354) +
                   DENDRO_343) -
              DENDRO_219 * DENDRO_298 * DENDRO_299 -
              DENDRO_243 * DENDRO_295 * grad_2_gt5[pp] +
              DENDRO_273 * DENDRO_362 + DENDRO_277 * grad_2_Gt0[pp] +
              DENDRO_278 * grad_2_Gt1[pp] - DENDRO_279 - DENDRO_280 -
              DENDRO_281 + DENDRO_283 * grad2_0_2_gt5[pp] +
              DENDRO_285 * grad2_0_0_gt5[pp] + DENDRO_287 * grad2_1_1_gt5[pp] +
              DENDRO_289 * grad2_2_2_gt5[pp] - DENDRO_290 * DENDRO_292 -
              DENDRO_293 * DENDRO_294 - DENDRO_297 - DENDRO_302 - DENDRO_305 +
              DENDRO_311 * DENDRO_312 + DENDRO_312 * (DENDRO_316 + DENDRO_317) +
              DENDRO_312 * (DENDRO_222 * grad_0_gt5[pp] + DENDRO_331) -
              DENDRO_326 + DENDRO_356 * grad_1_gt5[pp] +
              DENDRO_358 * grad_0_gt5[pp] + DENDRO_360 * grad_2_gt5[pp] +
              4 * grad_2_Gt2[pp] * gt5[pp]) -
         4 * grad2_2_2_alpha[pp]);
const double DENDRO_667 = DENDRO_31 * DENDRO_666;
const double DENDRO_668 = (1.0 / 12.0) * chi[pp];
const double DENDRO_669 = (1.0 / 3.0) * At1[pp];
const double DENDRO_670 = At1[pp] * DENDRO_21;
const double DENDRO_671 = At4[pp] * DENDRO_23;
const double DENDRO_672 = At3[pp] * DENDRO_26;
const double DENDRO_673 = DENDRO_670 + DENDRO_671 - DENDRO_672;
const double DENDRO_674 = At4[pp] * DENDRO_33;
const double DENDRO_675 =
    -At1[pp] * DENDRO_35 + At3[pp] * DENDRO_21 - DENDRO_674;
const double DENDRO_676 = At4[pp] * DENDRO_39;
const double DENDRO_677 =
    -At1[pp] * DENDRO_33 + At3[pp] * DENDRO_23 - DENDRO_676;
const double DENDRO_678 = 6.0 * grad_2_alpha[pp];
const double DENDRO_679 = 6.0 * grad_0_alpha[pp];
const double DENDRO_680 = 6.0 * grad_1_alpha[pp];
const double DENDRO_681 = DENDRO_92 * DENDRO_94;
const double DENDRO_682 = DENDRO_221 * grad_2_gt0[pp];
const double DENDRO_683 = DENDRO_151 * DENDRO_94;
const double DENDRO_684 = DENDRO_682 + DENDRO_683;
const double DENDRO_685 = DENDRO_384 * DENDRO_89;
const double DENDRO_686 = DENDRO_128 * DENDRO_386;
const double DENDRO_687 = DENDRO_179 * DENDRO_94;
const double DENDRO_688 =
    DENDRO_152 * DENDRO_95 + 0.25 * DENDRO_221 * grad_0_gt0[pp];
const double DENDRO_689 = DENDRO_134 + DENDRO_138;
const double DENDRO_690 = DENDRO_109 * DENDRO_576;
const double DENDRO_691 = DENDRO_201 * DENDRO_61;
const double DENDRO_692 = DENDRO_137 * DENDRO_237 + DENDRO_691;
const double DENDRO_693 = DENDRO_151 * DENDRO_89;
const double DENDRO_694 = DENDRO_221 * grad_1_gt0[pp] + DENDRO_693;
const double DENDRO_695 = 1.0 * DENDRO_135;
const double DENDRO_696 = DENDRO_151 * DENDRO_237;
const double DENDRO_697 = 2.0 * DENDRO_232;
const double DENDRO_698 = 2.0 * DENDRO_217;
const double DENDRO_699 = 2.0 * DENDRO_249;
const double DENDRO_700 = DENDRO_262 * DENDRO_45;
const double DENDRO_701 = (1.0 / 3.0) * At2[pp];
const double DENDRO_702 = At5[pp] * DENDRO_23;
const double DENDRO_703 =
    At2[pp] * DENDRO_21 - At4[pp] * DENDRO_26 + DENDRO_702;
const double DENDRO_704 = At5[pp] * DENDRO_33;
const double DENDRO_705 =
    -At2[pp] * DENDRO_35 + At4[pp] * DENDRO_21 - DENDRO_704;
const double DENDRO_706 = At4[pp] * DENDRO_23 - At5[pp] * DENDRO_39 - DENDRO_34;
const double DENDRO_707 = DENDRO_109 * grad_2_gt5[pp];
const double DENDRO_708 = DENDRO_89 * grad_0_gt5[pp];
const double DENDRO_709 = DENDRO_109 * DENDRO_307;
const double DENDRO_710 = 0.25 * DENDRO_128 * grad_2_gt5[pp];
const double DENDRO_711 = DENDRO_709 + DENDRO_710;
const double DENDRO_712 = DENDRO_87 * DENDRO_89;
const double DENDRO_713 = DENDRO_109 * DENDRO_110 + DENDRO_199 * DENDRO_61;
const double DENDRO_714 = -0.5 * DENDRO_173 * DENDRO_94;
const double DENDRO_715 = -0.5 * DENDRO_173 * DENDRO_89;
const double DENDRO_716 = DENDRO_183 * DENDRO_89;
const double DENDRO_717 = 2 * At4[pp];
const double DENDRO_718 = At3[pp] * DENDRO_37;
const double DENDRO_719 = DENDRO_31 * DENDRO_717;
const double DENDRO_720 = 0.5 * DENDRO_363;
const double DENDRO_721 = DENDRO_31 * DENDRO_75;
const double DENDRO_722 = DENDRO_213 * grad_2_gt3[pp];
const double DENDRO_723 = DENDRO_128 * DENDRO_307;
const double DENDRO_724 = 1.0 * DENDRO_221;
const double DENDRO_725 = 0.25 * DENDRO_693;
const double DENDRO_726 = DENDRO_213 * grad_0_gt3[pp];
const double DENDRO_727 = (1.0 / 3.0) * At4[pp];
const double DENDRO_728 = DENDRO_237 * grad_2_gt5[pp];
const double DENDRO_729 = DENDRO_110 * DENDRO_237;
const double DENDRO_730 = DENDRO_710 + DENDRO_729;
const double DENDRO_731 = DENDRO_237 * DENDRO_307;
const double DENDRO_732 = -0.5 * DENDRO_245 * grad_0_gt5[pp] + DENDRO_723;
const double DENDRO_733 = DENDRO_241 * grad_0_gt5[pp];
const double DENDRO_734 = DENDRO_241 * grad_1_gt5[pp];
const double DENDRO_735 = DENDRO_21 * grad_1_chi[pp] +
                          DENDRO_344 * grad_2_chi[pp] +
                          DENDRO_346 * grad_0_chi[pp];
const double DENDRO_736 = 0.5 * DENDRO_265;
const double DENDRO_737 = DENDRO_31 * grad_0_alpha[pp];
const double DENDRO_738 = DENDRO_349 * grad_1_chi[pp] + DENDRO_42;
const double DENDRO_739 = DENDRO_31 * grad_1_alpha[pp];
const double DENDRO_740 = DENDRO_23 * grad_1_chi[pp] +
                          DENDRO_344 * grad_0_chi[pp] +
                          DENDRO_351 * grad_2_chi[pp];
const double DENDRO_741 = 0.5 * DENDRO_740;
const double DENDRO_742 = DENDRO_31 * grad_2_alpha[pp];
const double DENDRO_743 = 0.5 * DENDRO_738;
const double DENDRO_744 = 0.5 * grad_1_alpha[pp];
const double DENDRO_745 = 0.5 * grad_2_alpha[pp];
const double DENDRO_746 = 0.5 * grad_0_alpha[pp];
const double DENDRO_747 = DENDRO_255 * DENDRO_31;
const double DENDRO_748 = (DENDRO_21 * DENDRO_21);
const double DENDRO_749 = (DENDRO_33 * DENDRO_33);
const double DENDRO_750 = 2 * DENDRO_35;
const double DENDRO_751 = At0[pp] * (DENDRO_35 * DENDRO_35) +
                          At3[pp] * DENDRO_748 + At5[pp] * DENDRO_749 -
                          DENDRO_253 * DENDRO_674 + DENDRO_34 * DENDRO_750 -
                          DENDRO_670 * DENDRO_750;
const double DENDRO_752 = 3 * DENDRO_84;
const double DENDRO_753 = (DENDRO_23 * DENDRO_23);
const double DENDRO_754 = 2 * DENDRO_26;
const double DENDRO_755 = At0[pp] * DENDRO_748 +
                          At3[pp] * (DENDRO_26 * DENDRO_26) +
                          At5[pp] * DENDRO_753 + DENDRO_24 * DENDRO_253 -
                          DENDRO_670 * DENDRO_754 - DENDRO_671 * DENDRO_754;
const double DENDRO_756 = At1[pp] * DENDRO_23;
const double DENDRO_757 = 2 * DENDRO_39;
const double DENDRO_758 = At0[pp] * DENDRO_749 + At3[pp] * DENDRO_753 +
                          At5[pp] * (DENDRO_39 * DENDRO_39) -
                          DENDRO_255 * DENDRO_756 + DENDRO_34 * DENDRO_757 -
                          DENDRO_671 * DENDRO_757;
const double DENDRO_759 = At3[pp] * DENDRO_21;
const double DENDRO_760 =
    At2[pp] * DENDRO_749 - DENDRO_21 * DENDRO_676 - DENDRO_23 * DENDRO_674 +
    DENDRO_23 * DENDRO_759 + DENDRO_33 * DENDRO_36 - DENDRO_33 * DENDRO_670 +
    DENDRO_35 * DENDRO_40 - DENDRO_35 * DENDRO_756 + DENDRO_39 * DENDRO_704;
const double DENDRO_761 = 6 * DENDRO_84;
const double DENDRO_762 =
    -At1[pp] * DENDRO_26 * DENDRO_35 - At1[pp] * DENDRO_748 -
    At4[pp] * DENDRO_21 * DENDRO_23 - At4[pp] * DENDRO_26 * DENDRO_33 +
    DENDRO_21 * DENDRO_34 + DENDRO_21 * DENDRO_36 + DENDRO_23 * DENDRO_704 +
    DENDRO_24 * DENDRO_35 + DENDRO_26 * DENDRO_759;
const double DENDRO_763 = -DENDRO_762;
const double DENDRO_764 =
    -At1[pp] * DENDRO_21 * DENDRO_23 - At1[pp] * DENDRO_26 * DENDRO_33 -
    At4[pp] * DENDRO_26 * DENDRO_39 - At4[pp] * DENDRO_753 +
    DENDRO_21 * DENDRO_40 + DENDRO_22 * DENDRO_33 + DENDRO_23 * DENDRO_34 +
    DENDRO_23 * DENDRO_672 + DENDRO_39 * DENDRO_702;
const double DENDRO_765 = -DENDRO_764;
const double DENDRO_766 = (1.0 / 3.0) * alpha[pp];
const double DENDRO_767 = (7.0 / 3.0) * DENDRO_282;
const double DENDRO_768 = (1.0 / 3.0) * DENDRO_282;
const double DENDRO_769 = (1.0 / 3.0) * DENDRO_284;
const double DENDRO_770 = 2 * DENDRO_84;
const double DENDRO_771 = DENDRO_770 * grad_0_alpha[pp];
const double DENDRO_772 = pow(DENDRO_30, -3);
const double DENDRO_773 = 4 * grad_0_K[pp];
const double DENDRO_774 = DENDRO_31 * DENDRO_766;
const double DENDRO_775 = DENDRO_770 * grad_2_alpha[pp];
const double DENDRO_776 = DENDRO_770 * grad_1_alpha[pp];
const double DENDRO_777 = 4 * grad_2_K[pp];
const double DENDRO_778 = 4 * grad_1_K[pp];
const double DENDRO_779 = 9 * DENDRO_45;
const double DENDRO_780 = DENDRO_189 * DENDRO_779;
const double DENDRO_781 =
    -2.0 / 3.0 * DENDRO_17 * DENDRO_260 * DENDRO_84 -
    2 * DENDRO_187 * DENDRO_751 * DENDRO_772 * alpha[pp] -
    7.0 / 3.0 * DENDRO_21 * DENDRO_31 * grad2_0_1_beta0[pp] -
    1.0 / 3.0 * DENDRO_21 * DENDRO_31 * grad2_1_1_beta1[pp] -
    1.0 / 3.0 * DENDRO_21 * DENDRO_31 * grad2_1_2_beta2[pp] -
    2.0 * DENDRO_219 * DENDRO_760 * DENDRO_772 * alpha[pp] -
    2.0 * DENDRO_222 * DENDRO_765 * DENDRO_772 * alpha[pp] -
    2.0 * DENDRO_223 * DENDRO_763 * DENDRO_772 * alpha[pp] -
    2 * DENDRO_225 * DENDRO_758 * DENDRO_772 * alpha[pp] -
    2 * DENDRO_228 * DENDRO_755 * DENDRO_772 * alpha[pp] -
    2 * DENDRO_23 * DENDRO_31 * grad2_1_2_beta0[pp] +
    (4.0 / 3.0) * DENDRO_284 * grad2_0_0_beta0[pp] +
    DENDRO_286 * grad2_1_1_beta0[pp] + DENDRO_288 * grad2_2_2_beta0[pp] +
    DENDRO_355 * grad_1_beta0[pp] + DENDRO_357 * grad_0_beta0[pp] +
    DENDRO_359 * grad_2_beta0[pp] + DENDRO_751 * DENDRO_771 +
    DENDRO_760 * DENDRO_775 + DENDRO_763 * DENDRO_776 +
    DENDRO_767 * grad2_0_2_beta0[pp] + DENDRO_768 * grad2_1_2_beta1[pp] +
    DENDRO_768 * grad2_2_2_beta2[pp] + DENDRO_769 * grad2_0_1_beta1[pp] +
    DENDRO_769 * grad2_0_2_beta2[pp] +
    DENDRO_774 * (DENDRO_21 * DENDRO_778 + DENDRO_763 * DENDRO_780) +
    DENDRO_774 * (9 * DENDRO_31 * DENDRO_45 * DENDRO_751 * grad_0_chi[pp] -
                  DENDRO_35 * DENDRO_773) +
    DENDRO_774 * (9 * DENDRO_31 * DENDRO_45 * DENDRO_760 * grad_2_chi[pp] -
                  DENDRO_33 * DENDRO_777) -
    beta0[pp] * grad_0_Gt0[pp] - beta1[pp] * grad_1_Gt0[pp] -
    beta2[pp] * grad_2_Gt0[pp];
const double DENDRO_782 = DENDRO_284 * grad2_0_0_beta1[pp];
const double DENDRO_783 = DENDRO_288 * grad2_2_2_beta1[pp];
const double DENDRO_784 = DENDRO_747 * grad2_0_2_beta1[pp];
const double DENDRO_785 = DENDRO_755 * DENDRO_776;
const double DENDRO_786 = (4.0 / 3.0) * DENDRO_286 * grad2_1_1_beta1[pp];
const double DENDRO_787 = DENDRO_21 * DENDRO_773;
const double DENDRO_788 = DENDRO_188 * DENDRO_779;
const double DENDRO_789 = DENDRO_23 * DENDRO_777;
const double DENDRO_790 = DENDRO_191 * DENDRO_779;
const double DENDRO_791 = (1.0 / 3.0) * DENDRO_286;
const double DENDRO_792 = DENDRO_791 * grad2_0_1_beta0[pp];
const double DENDRO_793 = DENDRO_791 * grad2_1_2_beta2[pp];
const double DENDRO_794 = DENDRO_26 * DENDRO_778 - DENDRO_755 * DENDRO_780;
const double DENDRO_795 = DENDRO_21 * DENDRO_31;
const double DENDRO_796 = (1.0 / 3.0) * DENDRO_795;
const double DENDRO_797 = DENDRO_23 * DENDRO_31;
const double DENDRO_798 = (1.0 / 3.0) * grad2_0_2_beta0[pp];
const double DENDRO_799 = DENDRO_0 * DENDRO_772;
const double DENDRO_800 = 2.0 * DENDRO_772 * alpha[pp];
const double DENDRO_801 = beta0[pp] * grad_0_Gt1[pp] +
                          beta1[pp] * grad_1_Gt1[pp] +
                          beta2[pp] * grad_2_Gt1[pp];
const double DENDRO_802 =
    -2.0 / 3.0 * DENDRO_17 * DENDRO_261 * DENDRO_84 -
    2.0 * DENDRO_174 * DENDRO_763 * DENDRO_772 * alpha[pp] -
    2.0 * DENDRO_176 * DENDRO_760 * DENDRO_772 * alpha[pp] -
    2 * DENDRO_190 * DENDRO_751 * DENDRO_772 * alpha[pp] -
    2 * DENDRO_21 * DENDRO_31 * grad2_0_1_beta2[pp] -
    1.0 / 3.0 * DENDRO_23 * DENDRO_31 * grad2_0_1_beta0[pp] -
    1.0 / 3.0 * DENDRO_23 * DENDRO_31 * grad2_1_1_beta1[pp] -
    7.0 / 3.0 * DENDRO_23 * DENDRO_31 * grad2_1_2_beta2[pp] -
    2.0 * DENDRO_238 * DENDRO_765 * DENDRO_772 * alpha[pp] -
    2 * DENDRO_242 * DENDRO_758 * DENDRO_772 * alpha[pp] -
    2 * DENDRO_245 * DENDRO_755 * DENDRO_772 * alpha[pp] +
    DENDRO_284 * grad2_0_0_beta2[pp] + DENDRO_286 * grad2_1_1_beta2[pp] +
    DENDRO_288 * DENDRO_798 + (1.0 / 3.0) * DENDRO_288 * grad2_1_2_beta1[pp] +
    (4.0 / 3.0) * DENDRO_288 * grad2_2_2_beta2[pp] +
    DENDRO_355 * grad_1_beta2[pp] + DENDRO_357 * grad_0_beta2[pp] +
    DENDRO_359 * grad_2_beta2[pp] + DENDRO_758 * DENDRO_775 +
    DENDRO_760 * DENDRO_771 + DENDRO_765 * DENDRO_776 +
    DENDRO_767 * grad2_0_2_beta2[pp] + DENDRO_768 * grad2_0_0_beta0[pp] +
    DENDRO_768 * grad2_0_1_beta1[pp] +
    DENDRO_774 * (DENDRO_23 * DENDRO_778 + DENDRO_765 * DENDRO_780) +
    DENDRO_774 * (9 * DENDRO_31 * DENDRO_45 * DENDRO_758 * grad_2_chi[pp] -
                  DENDRO_39 * DENDRO_777) +
    DENDRO_774 * (9 * DENDRO_31 * DENDRO_45 * DENDRO_760 * grad_0_chi[pp] -
                  DENDRO_33 * DENDRO_773) -
    beta0[pp] * grad_0_Gt2[pp] - beta1[pp] * grad_1_Gt2[pp] -
    beta2[pp] * grad_2_Gt2[pp];

// Dendro: printing variables
//--
a_rhs[pp] =
    -DENDRO_0 * K[pp] -
    0.59999999999999998 * DENDRO_1 * (-DENDRO_1 + alpha[pp]) *
        exp(-1.0 / 800.0 * (t * t)) +
    lambda[0] * (beta0[pp] * grad_0_alpha[pp] + beta1[pp] * grad_1_alpha[pp] +
                 beta2[pp] * grad_2_alpha[pp]);
//--
b_rhs0[pp]   = B0[pp] * DENDRO_2 + lambda[1] * (beta0[pp] * grad_0_beta0[pp] +
                                              beta1[pp] * grad_1_beta0[pp] +
                                              beta2[pp] * grad_2_beta0[pp]);
//--
b_rhs1[pp]   = B1[pp] * DENDRO_2 + lambda[1] * (beta0[pp] * grad_0_beta1[pp] +
                                              beta1[pp] * grad_1_beta1[pp] +
                                              beta2[pp] * grad_2_beta1[pp]);
//--
b_rhs2[pp]   = B2[pp] * DENDRO_2 + lambda[1] * (beta0[pp] * grad_0_beta2[pp] +
                                              beta1[pp] * grad_1_beta2[pp] +
                                              beta2[pp] * grad_2_beta2[pp]);
//--
gt_rhs00[pp] = -At0[pp] * DENDRO_0 + DENDRO_3 * gt0[pp] +
               DENDRO_4 * grad_0_beta1[pp] + DENDRO_5 * grad_0_beta2[pp] -
               DENDRO_6 * grad_1_beta1[pp] - DENDRO_6 * grad_2_beta2[pp] +
               beta0[pp] * grad_0_gt0[pp] + beta1[pp] * grad_1_gt0[pp] +
               beta2[pp] * grad_2_gt0[pp];
//--
gt_rhs01[pp] = -At1[pp] * DENDRO_0 + DENDRO_7 * grad_0_beta0[pp] +
               DENDRO_7 * grad_1_beta1[pp] - DENDRO_8 * gt1[pp] +
               beta0[pp] * grad_0_gt1[pp] + beta1[pp] * grad_1_gt1[pp] +
               beta2[pp] * grad_2_gt1[pp] + grad_0_beta1[pp] * gt3[pp] +
               grad_0_beta2[pp] * gt4[pp] + grad_1_beta0[pp] * gt0[pp] +
               grad_1_beta2[pp] * gt2[pp];
//--
gt_rhs02[pp] = -At2[pp] * DENDRO_0 - DENDRO_10 * gt2[pp] +
               DENDRO_9 * grad_0_beta0[pp] + DENDRO_9 * grad_2_beta2[pp] +
               beta0[pp] * grad_0_gt2[pp] + beta1[pp] * grad_1_gt2[pp] +
               beta2[pp] * grad_2_gt2[pp] + grad_0_beta1[pp] * gt4[pp] +
               grad_0_beta2[pp] * gt5[pp] + grad_2_beta0[pp] * gt0[pp] +
               grad_2_beta1[pp] * gt1[pp];
//--
gt_rhs11[pp] = -At3[pp] * DENDRO_0 - DENDRO_11 * gt3[pp] + DENDRO_12 * gt3[pp] +
               DENDRO_13 * grad_1_beta2[pp] + DENDRO_4 * grad_1_beta0[pp] -
               DENDRO_8 * gt3[pp] + beta0[pp] * grad_0_gt3[pp] +
               beta1[pp] * grad_1_gt3[pp] + beta2[pp] * grad_2_gt3[pp];
//--
gt_rhs12[pp] = -At4[pp] * DENDRO_0 - DENDRO_11 * gt4[pp] +
               DENDRO_14 * grad_1_beta1[pp] + DENDRO_14 * grad_2_beta2[pp] +
               beta0[pp] * grad_0_gt4[pp] + beta1[pp] * grad_1_gt4[pp] +
               beta2[pp] * grad_2_gt4[pp] + grad_1_beta0[pp] * gt2[pp] +
               grad_1_beta2[pp] * gt5[pp] + grad_2_beta0[pp] * gt1[pp] +
               grad_2_beta1[pp] * gt3[pp];
//--
gt_rhs22[pp] = -At5[pp] * DENDRO_0 - DENDRO_10 * gt5[pp] - DENDRO_11 * gt5[pp] +
               DENDRO_13 * grad_2_beta1[pp] + DENDRO_15 * gt5[pp] +
               DENDRO_5 * grad_2_beta0[pp] + beta0[pp] * grad_0_gt5[pp] +
               beta1[pp] * grad_1_gt5[pp] + beta2[pp] * grad_2_gt5[pp];
//--
chi_rhs[pp] = BSSN_CAHD_C * chi[pp] * (dx_i * dx_i) * ham[pp] / dt -
              DENDRO_16 * DENDRO_17 + DENDRO_16 * K[pp] * alpha[pp] +
              beta0[pp] * grad_0_chi[pp] + beta1[pp] * grad_1_chi[pp] +
              beta2[pp] * grad_2_chi[pp];
//--
At_rhs00[pp] =
    -At0[pp] * DENDRO_10 + At0[pp] * DENDRO_3 - At0[pp] * DENDRO_8 +
    DENDRO_18 * grad_0_beta1[pp] + DENDRO_19 * grad_0_beta2[pp] +
    DENDRO_668 *
        (DENDRO_264 *
             (-DENDRO_106 - DENDRO_109 * DENDRO_113 * DENDRO_115 - DENDRO_122 +
              4 * DENDRO_125 * DENDRO_128 * DENDRO_26 * DENDRO_84 - DENDRO_131 +
              4 * DENDRO_132 * DENDRO_35 * DENDRO_61 * DENDRO_84 -
              DENDRO_136 * (DENDRO_134 + DENDRO_51 * DENDRO_89) -
              DENDRO_136 * (DENDRO_138 + DENDRO_48 * DENDRO_94) -
              DENDRO_136 * (DENDRO_168 + DENDRO_171) -
              DENDRO_136 * (DENDRO_175 + DENDRO_178) -
              DENDRO_136 * (DENDRO_151 * DENDRO_155 + DENDRO_154) -
              DENDRO_142 * (DENDRO_139 + DENDRO_140) -
              DENDRO_142 * (DENDRO_164 + DENDRO_165) +
              2.0 * DENDRO_147 * DENDRO_21 * DENDRO_84 +
              4 * DENDRO_153 * DENDRO_23 * DENDRO_84 -
              DENDRO_157 * (1.0 * DENDRO_140 + DENDRO_156) -
              DENDRO_157 *
                  (DENDRO_128 * DENDRO_183 + 1.0 * DENDRO_151 * DENDRO_61) -
              DENDRO_163 -
              DENDRO_181 *
                  (-DENDRO_155 * DENDRO_52 + 2 * DENDRO_173 * DENDRO_61) -
              DENDRO_182 + 4 * DENDRO_185 * DENDRO_21 * DENDRO_84 -
              DENDRO_192 * (DENDRO_186 + DENDRO_187 * DENDRO_188 +
                            DENDRO_189 * DENDRO_55 + DENDRO_190 * DENDRO_191) -
              DENDRO_217 * DENDRO_218 - DENDRO_232 * DENDRO_233 -
              DENDRO_249 * DENDRO_250 +
              2.0 * DENDRO_26 * DENDRO_31 * grad2_1_1_gt0[pp] +
              3.0 * DENDRO_26 * DENDRO_84 * DENDRO_89 * grad_1_gt0[pp] -
              DENDRO_262 * DENDRO_263 +
              4 * DENDRO_31 * DENDRO_33 * grad2_0_2_gt0[pp] +
              2.0 * DENDRO_31 * DENDRO_35 * grad2_0_0_gt0[pp] +
              2.0 * DENDRO_31 * DENDRO_39 * grad2_2_2_gt0[pp] +
              2.0 * DENDRO_33 * DENDRO_84 * (DENDRO_143 + DENDRO_144) +
              4 * DENDRO_33 * DENDRO_84 * (1.0 * DENDRO_144 + DENDRO_158) +
              2.0 * DENDRO_33 * DENDRO_84 * (DENDRO_148 + DENDRO_149) +
              6.0 * DENDRO_35 * DENDRO_84 * DENDRO_95 * grad_0_gt0[pp] +
              3.0 * DENDRO_39 * DENDRO_84 * DENDRO_94 * grad_2_gt0[pp] -
              DENDRO_78 - DENDRO_81 - DENDRO_83 + 4 * grad_0_Gt0[pp] * gt0[pp] +
              4 * grad_0_Gt1[pp] * gt1[pp] + 4 * grad_0_Gt2[pp] * gt2[pp]) +
         DENDRO_56 * DENDRO_58 - DENDRO_62 * DENDRO_64 + DENDRO_667 * gt0[pp] -
         DENDRO_74 * DENDRO_75 - 12 * grad2_0_0_alpha[pp]) -
    alpha[pp] *
        (-At0[pp] * K[pp] +
         DENDRO_32 * (-At1[pp] * DENDRO_26 + DENDRO_22 + DENDRO_24) +
         DENDRO_38 * (At1[pp] * DENDRO_21 - DENDRO_34 - DENDRO_36) +
         DENDRO_41 * (-At0[pp] * DENDRO_33 + At1[pp] * DENDRO_23 - DENDRO_40)) +
    beta0[pp] * grad_0_At0[pp] + beta1[pp] * grad_1_At0[pp] +
    beta2[pp] * grad_2_At0[pp];
//--
At_rhs01[pp] =
    At0[pp] * grad_1_beta0[pp] - At1[pp] * DENDRO_8 +
    At2[pp] * grad_1_beta2[pp] + At3[pp] * grad_0_beta1[pp] +
    At4[pp] * grad_0_beta2[pp] +
    DENDRO_668 *
        (DENDRO_264 *
             (-DENDRO_105 * (DENDRO_167 * DENDRO_213 + DENDRO_659) -
              DENDRO_105 * (-DENDRO_101 * DENDRO_89 + DENDRO_660 + DENDRO_685) -
              DENDRO_105 * (-DENDRO_128 * DENDRO_296 +
                            0.5 * DENDRO_151 * DENDRO_245 - DENDRO_686) +
              DENDRO_114 * (DENDRO_681 + DENDRO_684) + DENDRO_115 * DENDRO_658 -
              DENDRO_136 * DENDRO_662 +
              DENDRO_136 * (-DENDRO_152 * DENDRO_213 + DENDRO_663) +
              DENDRO_136 * (DENDRO_384 * DENDRO_94 - DENDRO_576 * DENDRO_89 +
                            DENDRO_603) -
              DENDRO_142 *
                  (DENDRO_89 * grad_1_gt0[pp] + DENDRO_95 * grad_0_gt3[pp]) +
              DENDRO_157 * (-DENDRO_213 * DENDRO_48 + DENDRO_664) -
              DENDRO_157 * (DENDRO_128 * DENDRO_469 + DENDRO_128 * DENDRO_484 +
                            DENDRO_211 * DENDRO_61) +
              DENDRO_157 * (-DENDRO_128 * DENDRO_576 - DENDRO_212 * DENDRO_61 +
                            0.5 * DENDRO_245 * grad_2_gt0[pp]) +
              DENDRO_157 * (DENDRO_167 * DENDRO_95 - DENDRO_179 * DENDRO_89 +
                            DENDRO_627) +
              DENDRO_181 * DENDRO_661 + DENDRO_181 * (DENDRO_687 + DENDRO_688) +
              DENDRO_181 * (DENDRO_690 + DENDRO_692) +
              DENDRO_181 * (DENDRO_379 * DENDRO_95 + DENDRO_689) -
              DENDRO_192 * (DENDRO_100 * DENDRO_495 + DENDRO_174 * DENDRO_499 +
                            DENDRO_223 * DENDRO_493 + DENDRO_632) +
              DENDRO_299 *
                  (DENDRO_156 + DENDRO_48 * DENDRO_95 + DENDRO_49 * DENDRO_95) +
              DENDRO_299 * (0.25 * DENDRO_164 + 0.5 * DENDRO_165 +
                            DENDRO_379 * DENDRO_61) +
              DENDRO_639 + DENDRO_640 + DENDRO_641 + DENDRO_642 + DENDRO_643 +
              DENDRO_644 + DENDRO_645 + DENDRO_646 + DENDRO_647 + DENDRO_648 +
              DENDRO_649 + DENDRO_650 + DENDRO_651 - DENDRO_652 * DENDRO_700 +
              DENDRO_665 -
              DENDRO_695 * (DENDRO_694 + DENDRO_94 * grad_0_gt3[pp]) -
              DENDRO_695 * (DENDRO_109 * grad_2_gt3[pp] +
                            DENDRO_128 * grad_1_gt5[pp] + DENDRO_696) -
              DENDRO_697 * grad_0_gt1[pp] - DENDRO_698 * grad_1_gt1[pp] -
              DENDRO_699 * grad_2_gt1[pp]) +
         DENDRO_31 * DENDRO_678 * (-DENDRO_128 - DENDRO_59 * DENDRO_656) +
         DENDRO_652 * DENDRO_666 - DENDRO_654 * DENDRO_679 -
         DENDRO_655 * DENDRO_680 - 12 * grad2_0_1_alpha[pp]) +
    DENDRO_669 * grad_0_beta0[pp] + DENDRO_669 * grad_1_beta1[pp] -
    alpha[pp] * (-At1[pp] * K[pp] + DENDRO_32 * DENDRO_673 +
                 DENDRO_38 * DENDRO_675 + DENDRO_41 * DENDRO_677) +
    beta0[pp] * grad_0_At1[pp] + beta1[pp] * grad_1_At1[pp] +
    beta2[pp] * grad_2_At1[pp];
//--
At_rhs02[pp] =
    At0[pp] * grad_2_beta0[pp] + At1[pp] * grad_2_beta1[pp] -
    At2[pp] * DENDRO_10 + At4[pp] * grad_0_beta1[pp] +
    At5[pp] * grad_0_beta2[pp] +
    DENDRO_668 *
        (DENDRO_264 *
             (-DENDRO_115 * (DENDRO_173 * DENDRO_241 - 0.5 * DENDRO_707) -
              DENDRO_115 * (-DENDRO_110 * DENDRO_94 - DENDRO_224 * DENDRO_51 -
                            DENDRO_714) -
              DENDRO_136 * DENDRO_532 -
              DENDRO_136 * (DENDRO_152 * DENDRO_241 + DENDRO_711) -
              DENDRO_157 * (DENDRO_154 + DENDRO_692) -
              DENDRO_157 * (DENDRO_688 + DENDRO_716) -
              DENDRO_157 * (DENDRO_319 * DENDRO_95 + DENDRO_689) -
              DENDRO_157 * (DENDRO_109 * DENDRO_484 + DENDRO_154 + DENDRO_691) -
              DENDRO_181 * (0.5 * DENDRO_109 * DENDRO_173 - DENDRO_713) -
              DENDRO_181 * (DENDRO_173 * DENDRO_95 - DENDRO_183 * DENDRO_94 -
                            DENDRO_224 * DENDRO_53) -
              DENDRO_192 * (DENDRO_117 * DENDRO_495 + DENDRO_176 * DENDRO_499 +
                            DENDRO_219 * DENDRO_493 + DENDRO_490) +
              4 * DENDRO_21 * DENDRO_535 * DENDRO_84 +
              4 * DENDRO_21 * DENDRO_542 * DENDRO_84 +
              DENDRO_23 * DENDRO_523 * DENDRO_84 +
              4 * DENDRO_23 * DENDRO_84 *
                  (0.5 * DENDRO_173 * DENDRO_237 - DENDRO_711) +
              4 * DENDRO_23 * DENDRO_84 *
                  (-DENDRO_224 * DENDRO_48 - DENDRO_469 * DENDRO_94 -
                   DENDRO_715) +
              2.0 * DENDRO_26 * DENDRO_31 * grad2_1_1_gt2[pp] +
              DENDRO_26 * DENDRO_84 * (DENDRO_694 + DENDRO_712) +
              4 * DENDRO_26 * DENDRO_84 *
                  (DENDRO_128 * DENDRO_201 + 0.25 * DENDRO_696) +
              4 * DENDRO_31 * DENDRO_33 * grad2_0_2_gt2[pp] +
              2.0 * DENDRO_31 * DENDRO_35 * grad2_0_0_gt2[pp] +
              2.0 * DENDRO_31 * DENDRO_39 * grad2_2_2_gt2[pp] +
              4 * DENDRO_33 * DENDRO_84 *
                  (DENDRO_241 * DENDRO_51 + DENDRO_713) +
              2.0 * DENDRO_33 * DENDRO_84 *
                  (DENDRO_94 * grad_2_gt0[pp] + DENDRO_95 * grad_0_gt5[pp]) +
              4 * DENDRO_35 * DENDRO_84 *
                  (0.25 * DENDRO_148 + 1.0 * DENDRO_149) +
              4 * DENDRO_35 * DENDRO_84 *
                  (DENDRO_158 + DENDRO_51 * DENDRO_95 + DENDRO_52 * DENDRO_95) -
              DENDRO_506 - DENDRO_507 - DENDRO_508 - DENDRO_509 * DENDRO_700 -
              DENDRO_524 - DENDRO_526 - DENDRO_527 - DENDRO_528 - DENDRO_530 -
              DENDRO_695 * (DENDRO_684 + DENDRO_708) -
              DENDRO_697 * grad_0_gt2[pp] - DENDRO_698 * grad_1_gt2[pp] -
              DENDRO_699 * grad_2_gt2[pp] + 2.0 * grad_0_Gt0[pp] * gt2[pp] +
              2.0 * grad_0_Gt1[pp] * gt4[pp] + 2.0 * grad_0_Gt2[pp] * gt5[pp] +
              2.0 * grad_2_Gt0[pp] * gt0[pp] + 2.0 * grad_2_Gt1[pp] * gt1[pp] +
              2.0 * grad_2_Gt2[pp] * gt2[pp]) +
         DENDRO_509 * DENDRO_666 - DENDRO_511 * DENDRO_679 -
         DENDRO_513 * DENDRO_678 + DENDRO_517 * DENDRO_680 -
         12 * grad2_0_2_alpha[pp]) +
    DENDRO_701 * grad_0_beta0[pp] + DENDRO_701 * grad_2_beta2[pp] -
    alpha[pp] * (-At2[pp] * K[pp] + DENDRO_32 * DENDRO_703 +
                 DENDRO_38 * DENDRO_705 + DENDRO_41 * DENDRO_706) +
    beta0[pp] * grad_0_At2[pp] + beta1[pp] * grad_1_At2[pp] +
    beta2[pp] * grad_2_At2[pp];
//--
At_rhs11[pp] =
    -At3[pp] * DENDRO_11 + At3[pp] * DENDRO_12 - At3[pp] * DENDRO_8 +
    DENDRO_18 * grad_1_beta0[pp] +
    DENDRO_668 *
        (DENDRO_264 *
             (6.0 * DENDRO_104 * DENDRO_213 * grad_1_gt3[pp] -
              DENDRO_115 * DENDRO_237 * DENDRO_402 + DENDRO_128 * DENDRO_406 +
              DENDRO_136 * (DENDRO_390 - 1.0 * DENDRO_722) -
              DENDRO_136 * (-1.0 * DENDRO_228 * DENDRO_92 + DENDRO_385) -
              DENDRO_136 * (DENDRO_237 * DENDRO_386 + DENDRO_388) +
              DENDRO_142 * (DENDRO_100 * grad_1_gt3[pp] - DENDRO_726) +
              DENDRO_142 *
                  (-DENDRO_128 * grad_2_gt3[pp] + DENDRO_245 * DENDRO_87) +
              DENDRO_142 *
                  (DENDRO_228 * grad_1_gt0[pp] - DENDRO_89 * grad_0_gt3[pp]) +
              DENDRO_157 * (DENDRO_394 + DENDRO_685) +
              DENDRO_157 *
                  (0.25 * DENDRO_100 * grad_1_gt3[pp] - 1.0 * DENDRO_726) +
              DENDRO_157 * (DENDRO_245 * DENDRO_92 - DENDRO_686) +
              DENDRO_181 * (DENDRO_381 + DENDRO_383) +
              DENDRO_181 * (DENDRO_133 * DENDRO_221 + DENDRO_379 * DENDRO_89) +
              DENDRO_181 * (DENDRO_237 * DENDRO_379 + DENDRO_723) +
              DENDRO_181 * (DENDRO_49 * DENDRO_724 + DENDRO_725) -
              DENDRO_192 * (DENDRO_188 * DENDRO_228 + DENDRO_189 * DENDRO_214 +
                            DENDRO_191 * DENDRO_245 + DENDRO_368) +
              DENDRO_221 * DENDRO_404 + DENDRO_312 * (DENDRO_389 - DENDRO_722) +
              DENDRO_312 *
                  (DENDRO_151 * DENDRO_228 - DENDRO_221 * grad_0_gt3[pp]) +
              DENDRO_312 * (-DENDRO_237 * grad_2_gt3[pp] + DENDRO_400) -
              DENDRO_366 * DENDRO_700 + DENDRO_405 * DENDRO_89 + DENDRO_413 -
              DENDRO_697 * grad_0_gt3[pp] - DENDRO_698 * grad_1_gt3[pp] -
              DENDRO_699 * grad_2_gt3[pp]) -
         DENDRO_367 * DENDRO_57 +
         DENDRO_64 * (DENDRO_245 - DENDRO_59 * DENDRO_720) +
         DENDRO_667 * gt3[pp] +
         DENDRO_721 * (DENDRO_228 - DENDRO_71 * DENDRO_720) -
         12 * grad2_1_1_alpha[pp]) +
    DENDRO_717 * grad_1_beta2[pp] -
    alpha[pp] * (-At3[pp] * K[pp] + DENDRO_32 * DENDRO_675 +
                 DENDRO_673 * DENDRO_718 + DENDRO_677 * DENDRO_719) +
    beta0[pp] * grad_0_At3[pp] + beta1[pp] * grad_1_At3[pp] +
    beta2[pp] * grad_2_At3[pp];
//--
At_rhs12[pp] =
    At1[pp] * grad_2_beta0[pp] + At2[pp] * grad_1_beta0[pp] +
    At3[pp] * grad_2_beta1[pp] - At4[pp] * DENDRO_11 +
    At5[pp] * grad_1_beta2[pp] +
    DENDRO_668 *
        (DENDRO_264 *
             (-DENDRO_105 * (-DENDRO_237 * DENDRO_296 + DENDRO_595) -
              DENDRO_105 *
                  (-DENDRO_101 * DENDRO_221 + DENDRO_551 + DENDRO_596) -
              DENDRO_105 * (0.25 * DENDRO_198 * grad_1_gt3[pp] -
                            DENDRO_211 * DENDRO_213 - DENDRO_212 * DENDRO_213) -
              DENDRO_115 * (DENDRO_204 * DENDRO_241 - 0.5 * DENDRO_728) -
              DENDRO_115 *
                  (-DENDRO_110 * DENDRO_221 + 0.5 * DENDRO_173 * DENDRO_221 -
                   DENDRO_224 * DENDRO_379) +
              DENDRO_136 * (DENDRO_204 * DENDRO_213 + DENDRO_600) -
              DENDRO_136 * (DENDRO_166 * DENDRO_224 + DENDRO_221 * DENDRO_469 +
                            DENDRO_599) +
              DENDRO_136 * (-DENDRO_211 * DENDRO_241 +
                            0.5 * DENDRO_245 * grad_2_gt5[pp] - DENDRO_731) +
              DENDRO_136 * (-DENDRO_221 * DENDRO_484 - DENDRO_221 * DENDRO_576 +
                            0.5 * DENDRO_228 * grad_0_gt5[pp]) +
              DENDRO_136 * (DENDRO_237 * DENDRO_335 + DENDRO_563 - DENDRO_731) +
              DENDRO_157 * (-DENDRO_109 * DENDRO_296 - DENDRO_732) +
              DENDRO_157 * (-DENDRO_213 * DENDRO_319 + DENDRO_602) +
              DENDRO_157 * (-DENDRO_237 * DENDRO_576 - DENDRO_732) +
              DENDRO_157 *
                  (-DENDRO_101 * DENDRO_94 + 0.5 * DENDRO_228 * grad_2_gt0[pp] -
                   0.25 * DENDRO_712) +
              DENDRO_157 * (-DENDRO_179 * DENDRO_221 + DENDRO_228 * DENDRO_52 -
                            DENDRO_725) +
              DENDRO_157 *
                  (-DENDRO_213 * DENDRO_379 + DENDRO_445 + DENDRO_573) +
              DENDRO_161 * (DENDRO_681 + DENDRO_682 + DENDRO_708) -
              DENDRO_181 * (0.5 * DENDRO_109 * DENDRO_204 - DENDRO_730) +
              DENDRO_181 * (DENDRO_241 * DENDRO_379 + DENDRO_730) -
              DENDRO_181 * (-DENDRO_183 * DENDRO_221 - DENDRO_224 * DENDRO_49 -
                            DENDRO_715) -
              DENDRO_192 * (DENDRO_198 * DENDRO_495 + DENDRO_222 * DENDRO_493 +
                            DENDRO_238 * DENDRO_499 + DENDRO_581) +
              DENDRO_299 * (DENDRO_128 * DENDRO_172 + DENDRO_690) +
              DENDRO_299 * (DENDRO_138 + DENDRO_687 + DENDRO_716) +
              DENDRO_312 *
                  (DENDRO_198 * grad_2_gt3[pp] - DENDRO_213 * grad_1_gt5[pp]) -
              DENDRO_587 * DENDRO_700 + DENDRO_588 + DENDRO_608 -
              DENDRO_697 * grad_0_gt4[pp] - DENDRO_698 * grad_1_gt4[pp] -
              DENDRO_699 * grad_2_gt4[pp]) -
         DENDRO_31 * DENDRO_593 * DENDRO_679 + DENDRO_587 * DENDRO_666 +
         DENDRO_591 * DENDRO_680 - DENDRO_592 * DENDRO_678 -
         12 * grad2_1_2_alpha[pp]) +
    DENDRO_727 * grad_1_beta1[pp] + DENDRO_727 * grad_2_beta2[pp] -
    alpha[pp] * (-At4[pp] * K[pp] + DENDRO_32 * DENDRO_705 +
                 DENDRO_703 * DENDRO_718 + DENDRO_706 * DENDRO_719) +
    beta0[pp] * grad_0_At4[pp] + beta1[pp] * grad_1_At4[pp] +
    beta2[pp] * grad_2_At4[pp];
//--
At_rhs22[pp] =
    -At5[pp] * DENDRO_10 - At5[pp] * DENDRO_11 + At5[pp] * DENDRO_15 +
    DENDRO_19 * grad_2_beta0[pp] +
    DENDRO_668 *
        (DENDRO_264 *
             (3.0 * DENDRO_109 * DENDRO_35 * DENDRO_84 * grad_0_gt5[pp] +
              DENDRO_115 * DENDRO_224 * DENDRO_304 -
              DENDRO_136 * (0.25 * DENDRO_728 + 1.0 * DENDRO_734) -
              DENDRO_136 * (-1.0 * DENDRO_225 * DENDRO_87 + DENDRO_338) -
              DENDRO_157 * (DENDRO_109 * DENDRO_201 + DENDRO_729) -
              DENDRO_157 * (DENDRO_137 * DENDRO_221 + DENDRO_319 * DENDRO_94) -
              DENDRO_157 * (DENDRO_172 * DENDRO_237 + DENDRO_709) -
              DENDRO_157 * (DENDRO_52 * DENDRO_724 + 0.25 * DENDRO_683) -
              DENDRO_181 * (-DENDRO_224 * DENDRO_342 - DENDRO_714) -
              DENDRO_192 * (DENDRO_188 * DENDRO_225 + DENDRO_189 * DENDRO_205 +
                            DENDRO_191 * DENDRO_242 + DENDRO_343) +
              4 * DENDRO_21 * DENDRO_320 * DENDRO_84 +
              4 * DENDRO_21 * DENDRO_330 * DENDRO_84 +
              4 * DENDRO_221 * DENDRO_26 * DENDRO_303 * DENDRO_84 +
              2.0 * DENDRO_23 * DENDRO_311 * DENDRO_84 +
              4 * DENDRO_23 * DENDRO_340 * DENDRO_84 +
              3.0 * DENDRO_237 * DENDRO_26 * DENDRO_84 * grad_1_gt5[pp] +
              6.0 * DENDRO_241 * DENDRO_39 * DENDRO_84 * grad_2_gt5[pp] +
              2.0 * DENDRO_26 * DENDRO_31 * grad2_1_1_gt5[pp] -
              DENDRO_273 * DENDRO_700 - DENDRO_279 - DENDRO_280 - DENDRO_281 -
              DENDRO_297 + 4 * DENDRO_298 * DENDRO_35 * DENDRO_84 * DENDRO_94 -
              DENDRO_302 - DENDRO_305 +
              4 * DENDRO_31 * DENDRO_33 * grad2_0_2_gt5[pp] +
              2.0 * DENDRO_31 * DENDRO_35 * grad2_0_0_gt5[pp] +
              2.0 * DENDRO_31 * DENDRO_39 * grad2_2_2_gt5[pp] -
              DENDRO_312 * (DENDRO_728 + DENDRO_734) -
              DENDRO_312 *
                  (DENDRO_151 * DENDRO_224 + DENDRO_221 * grad_0_gt5[pp]) -
              DENDRO_326 +
              4 * DENDRO_33 * DENDRO_84 *
                  (0.25 * DENDRO_707 + 1.0 * DENDRO_733) +
              2.0 * DENDRO_33 * DENDRO_84 * (DENDRO_707 + DENDRO_733) +
              4 * DENDRO_33 * DENDRO_84 *
                  (-1.0 * DENDRO_205 * DENDRO_87 + DENDRO_336) +
              2.0 * DENDRO_33 * DENDRO_84 *
                  (DENDRO_224 * grad_2_gt0[pp] + DENDRO_94 * grad_0_gt5[pp]) -
              DENDRO_697 * grad_0_gt5[pp] - DENDRO_698 * grad_1_gt5[pp] -
              DENDRO_699 * grad_2_gt5[pp] + 4 * grad_2_Gt0[pp] * gt2[pp] +
              4 * grad_2_Gt1[pp] * gt4[pp] + 4 * grad_2_Gt2[pp] * gt5[pp]) +
         DENDRO_266 * DENDRO_58 - DENDRO_269 * DENDRO_721 -
         DENDRO_275 * DENDRO_63 + DENDRO_667 * gt5[pp] -
         12 * grad2_2_2_alpha[pp]) +
    DENDRO_717 * grad_2_beta1[pp] -
    alpha[pp] * (At5[pp] * DENDRO_37 * DENDRO_706 - At5[pp] * K[pp] +
                 DENDRO_41 * DENDRO_705 + DENDRO_703 * DENDRO_719) +
    beta0[pp] * grad_0_At5[pp] + beta1[pp] * grad_1_At5[pp] +
    beta2[pp] * grad_2_At5[pp];
//--
K_rhs[pp] =
    DENDRO_253 * DENDRO_31 * chi[pp] *
        (0.5 * DENDRO_742 * (DENDRO_637 + DENDRO_656 * DENDRO_740) +
         DENDRO_744 *
             (DENDRO_31 * DENDRO_635 + DENDRO_31 * DENDRO_97 +
              DENDRO_31 * DENDRO_98 -
              DENDRO_45 * (-DENDRO_652 * DENDRO_738 + grad_0_chi[pp])) +
         DENDRO_746 *
             (DENDRO_31 * DENDRO_633 + DENDRO_31 * DENDRO_634 +
              DENDRO_31 * DENDRO_85 -
              DENDRO_45 * (-DENDRO_652 * DENDRO_735 + grad_1_chi[pp])) -
         grad2_0_1_alpha[pp]) +
    DENDRO_256 * DENDRO_31 * chi[pp] *
        (0.5 * DENDRO_737 * (DENDRO_45 * DENDRO_735 * gt4[pp] + DENDRO_582) +
         DENDRO_744 * (DENDRO_31 * DENDRO_583 -
                       DENDRO_45 * (-DENDRO_587 * DENDRO_738 + grad_2_chi[pp]) +
                       DENDRO_590) +
         DENDRO_745 *
             (DENDRO_31 * DENDRO_584 + DENDRO_31 * DENDRO_585 +
              DENDRO_31 * DENDRO_586 -
              DENDRO_45 * (-DENDRO_587 * DENDRO_740 + grad_1_chi[pp])) -
         grad2_1_2_alpha[pp]) -
    DENDRO_284 * chi[pp] *
        (DENDRO_739 * (DENDRO_431 + DENDRO_46 * DENDRO_743) +
         DENDRO_742 * (DENDRO_432 + DENDRO_46 * DENDRO_741) -
         grad2_0_0_alpha[pp] +
         grad_0_alpha[pp] *
             (DENDRO_31 * DENDRO_429 + DENDRO_31 * DENDRO_430 -
              DENDRO_45 * (DENDRO_69 - 0.5 * DENDRO_70 * DENDRO_735) +
              DENDRO_68)) -
    DENDRO_286 * chi[pp] *
        (DENDRO_737 * (DENDRO_370 + DENDRO_720 * DENDRO_735) +
         DENDRO_742 * (DENDRO_374 + DENDRO_720 * DENDRO_740) -
         grad2_1_1_alpha[pp] +
         grad_1_alpha[pp] *
             (DENDRO_31 * DENDRO_371 + DENDRO_31 * DENDRO_372 +
              DENDRO_31 * DENDRO_373 -
              DENDRO_45 * (DENDRO_365 - DENDRO_366 * DENDRO_743))) -
    DENDRO_288 * chi[pp] *
        (DENDRO_737 * (DENDRO_348 + DENDRO_735 * DENDRO_736) +
         DENDRO_739 * (DENDRO_350 + DENDRO_736 * DENDRO_738) -
         grad2_2_2_alpha[pp] +
         grad_2_alpha[pp] *
             (DENDRO_31 * DENDRO_352 + DENDRO_31 * DENDRO_353 +
              DENDRO_31 * DENDRO_354 -
              DENDRO_45 * (DENDRO_272 - DENDRO_273 * DENDRO_741))) -
    DENDRO_747 * chi[pp] *
        (0.5 * DENDRO_739 * (DENDRO_494 + DENDRO_516 * DENDRO_738) +
         DENDRO_745 *
             (DENDRO_31 * DENDRO_496 + DENDRO_31 * DENDRO_497 +
              DENDRO_31 * DENDRO_498 -
              DENDRO_45 * (-DENDRO_509 * DENDRO_740 + grad_0_chi[pp])) +
         DENDRO_746 *
             (DENDRO_31 * DENDRO_491 + DENDRO_31 * DENDRO_492 +
              DENDRO_31 * DENDRO_93 -
              DENDRO_45 * (-DENDRO_509 * DENDRO_735 + grad_2_chi[pp])) -
         grad2_0_2_alpha[pp]) +
    DENDRO_766 *
        (At0[pp] * DENDRO_751 * DENDRO_752 + At1[pp] * DENDRO_761 * DENDRO_763 +
         At2[pp] * DENDRO_760 * DENDRO_761 + At3[pp] * DENDRO_752 * DENDRO_755 +
         At4[pp] * DENDRO_761 * DENDRO_765 + At5[pp] * DENDRO_752 * DENDRO_758 +
         pow(K[pp], 2)) +
    beta0[pp] * grad_0_K[pp] + beta1[pp] * grad_1_K[pp] +
    beta2[pp] * grad_2_K[pp];
//--
Gt_rhs0[pp] = -DENDRO_781;
//--
Gt_rhs1[pp] =
    -DENDRO_100 * DENDRO_762 * DENDRO_800 +
    DENDRO_117 * DENDRO_760 * DENDRO_800 - 2.0 / 3.0 * DENDRO_17 * DENDRO_217 -
    DENDRO_198 * DENDRO_764 * DENDRO_800 +
    DENDRO_205 * DENDRO_758 * DENDRO_799 -
    DENDRO_213 * DENDRO_755 * DENDRO_799 + DENDRO_217 * grad_1_beta1[pp] +
    DENDRO_232 * grad_0_beta1[pp] + DENDRO_249 * grad_2_beta1[pp] +
    DENDRO_55 * DENDRO_751 * DENDRO_799 + DENDRO_762 * DENDRO_771 +
    DENDRO_764 * DENDRO_775 + DENDRO_774 * DENDRO_794 -
    DENDRO_774 * (-DENDRO_762 * DENDRO_788 + DENDRO_787) -
    DENDRO_774 * (-DENDRO_764 * DENDRO_790 + DENDRO_789) - DENDRO_782 -
    DENDRO_783 - DENDRO_784 - DENDRO_785 - DENDRO_786 - DENDRO_792 -
    DENDRO_793 + (7.0 / 3.0) * DENDRO_795 * grad2_0_1_beta1[pp] +
    DENDRO_796 * grad2_0_0_beta0[pp] + DENDRO_796 * grad2_0_2_beta2[pp] +
    DENDRO_797 * DENDRO_798 + (7.0 / 3.0) * DENDRO_797 * grad2_1_2_beta1[pp] +
    (1.0 / 3.0) * DENDRO_797 * grad2_2_2_beta2[pp] + DENDRO_801;
//--
Gt_rhs2[pp] = -DENDRO_802;
//--
B_rhs0[pp] =
    -B0[pp] * eta - DENDRO_781 +
    lambda[2] * (beta0[pp] * grad_0_B0[pp] + beta1[pp] * grad_1_B0[pp] +
                 beta2[pp] * grad_2_B0[pp]) -
    lambda[3] * (beta0[pp] * grad_0_Gt0[pp] + beta1[pp] * grad_1_Gt0[pp] +
                 beta2[pp] * grad_2_Gt0[pp]);
//--
B_rhs1[pp] =
    -B1[pp] * eta + 2.0 * DENDRO_100 * DENDRO_763 * DENDRO_772 * alpha[pp] +
    2.0 * DENDRO_117 * DENDRO_760 * DENDRO_772 * alpha[pp] +
    (2.0 / 3.0) * DENDRO_17 * DENDRO_259 * DENDRO_84 +
    2.0 * DENDRO_198 * DENDRO_765 * DENDRO_772 * alpha[pp] +
    2 * DENDRO_205 * DENDRO_758 * DENDRO_772 * alpha[pp] +
    (1.0 / 3.0) * DENDRO_21 * DENDRO_31 * grad2_0_0_beta0[pp] +
    (7.0 / 3.0) * DENDRO_21 * DENDRO_31 * grad2_0_1_beta1[pp] +
    (1.0 / 3.0) * DENDRO_21 * DENDRO_31 * grad2_0_2_beta2[pp] +
    2 * DENDRO_214 * DENDRO_755 * DENDRO_772 * alpha[pp] +
    (1.0 / 3.0) * DENDRO_23 * DENDRO_31 * grad2_0_2_beta0[pp] +
    (7.0 / 3.0) * DENDRO_23 * DENDRO_31 * grad2_1_2_beta1[pp] +
    (1.0 / 3.0) * DENDRO_23 * DENDRO_31 * grad2_2_2_beta2[pp] -
    DENDRO_355 * grad_1_beta1[pp] - DENDRO_357 * grad_0_beta1[pp] -
    DENDRO_359 * grad_2_beta1[pp] +
    2 * DENDRO_55 * DENDRO_751 * DENDRO_772 * alpha[pp] -
    DENDRO_763 * DENDRO_771 - DENDRO_765 * DENDRO_775 +
    DENDRO_774 * DENDRO_794 -
    DENDRO_774 * (DENDRO_763 * DENDRO_788 + DENDRO_787) -
    DENDRO_774 * (DENDRO_765 * DENDRO_790 + DENDRO_789) - DENDRO_782 -
    DENDRO_783 - DENDRO_784 - DENDRO_785 - DENDRO_786 - DENDRO_792 -
    DENDRO_793 - DENDRO_801 * lambda[3] + beta0[pp] * grad_0_Gt1[pp] +
    beta1[pp] * grad_1_Gt1[pp] + beta2[pp] * grad_2_Gt1[pp] +
    lambda[2] * (beta0[pp] * grad_0_B1[pp] + beta1[pp] * grad_1_B1[pp] +
                 beta2[pp] * grad_2_B1[pp]);
//--
B_rhs2[pp] =
    -B2[pp] * eta - DENDRO_802 +
    lambda[2] * (beta0[pp] * grad_0_B2[pp] + beta1[pp] * grad_1_B2[pp] +
                 beta2[pp] * grad_2_B2[pp]) -
    lambda[3] * (beta0[pp] * grad_0_Gt2[pp] + beta1[pp] * grad_1_Gt2[pp] +
                 beta2[pp] * grad_2_Gt2[pp]);
// Dendro: reduced ops: 4312
// Dendro: }}}
