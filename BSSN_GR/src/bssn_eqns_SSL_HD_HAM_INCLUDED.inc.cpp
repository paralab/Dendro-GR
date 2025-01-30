// CODEGEN: SSL was enabled, adding term to gauge condition!
// CODEGEN: CAHD was enabled, adding damping term to chi!
// Codgen: generating unstage version
// Codgen: using standard gauge
// Codgen: using eta const damping
// Dendro: {{{
// Dendro: original ops: 708883
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
const double DENDRO_18 = pow(K[pp], 2);
const double DENDRO_19 = pow(gt4[pp], 2);
const double DENDRO_20 = pow(gt1[pp], 2);
const double DENDRO_21 = pow(gt2[pp], 2);
const double DENDRO_22 = gt3[pp] * gt5[pp];
const double DENDRO_23 = gt2[pp] * gt4[pp];
const double DENDRO_24 = DENDRO_19 * gt0[pp] + DENDRO_20 * gt5[pp] +
                         DENDRO_21 * gt3[pp] - DENDRO_22 * gt0[pp] -
                         DENDRO_23 * DENDRO_4;
const double DENDRO_25 = pow(DENDRO_24, -2);
const double DENDRO_26 = 12 * DENDRO_25;
const double DENDRO_27 = -DENDRO_23 + gt1[pp] * gt5[pp];
const double DENDRO_28 = (DENDRO_27 * DENDRO_27);
const double DENDRO_29 = gt1[pp] * gt4[pp] - gt2[pp] * gt3[pp];
const double DENDRO_30 = (DENDRO_29 * DENDRO_29);
const double DENDRO_31 = -DENDRO_19 + DENDRO_22;
const double DENDRO_32 = At4[pp] * DENDRO_29;
const double DENDRO_33 = 2 * DENDRO_27;
const double DENDRO_34 = At1[pp] * DENDRO_27;
const double DENDRO_35 = 2 * DENDRO_31;
const double DENDRO_36 = At2[pp] * DENDRO_29;
const double DENDRO_37 = At0[pp] * (DENDRO_31 * DENDRO_31) +
                         At3[pp] * DENDRO_28 + At5[pp] * DENDRO_30 -
                         DENDRO_32 * DENDRO_33 - DENDRO_34 * DENDRO_35 +
                         DENDRO_35 * DENDRO_36;
const double DENDRO_38 = At0[pp] * DENDRO_37;
const double DENDRO_39 = gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp];
const double DENDRO_40 = (DENDRO_39 * DENDRO_39);
const double DENDRO_41 = -DENDRO_21 + gt0[pp] * gt5[pp];
const double DENDRO_42 = At2[pp] * DENDRO_39;
const double DENDRO_43 = 2 * DENDRO_41;
const double DENDRO_44 = At4[pp] * DENDRO_39;
const double DENDRO_45 = At0[pp] * DENDRO_28 +
                         At3[pp] * (DENDRO_41 * DENDRO_41) +
                         At5[pp] * DENDRO_40 + DENDRO_33 * DENDRO_42 -
                         DENDRO_34 * DENDRO_43 - DENDRO_43 * DENDRO_44;
const double DENDRO_46 = At3[pp] * DENDRO_45;
const double DENDRO_47 = -DENDRO_20 + gt0[pp] * gt3[pp];
const double DENDRO_48 = At1[pp] * DENDRO_39;
const double DENDRO_49 = 2 * DENDRO_29;
const double DENDRO_50 = 2 * DENDRO_47;
const double DENDRO_51 = At0[pp] * DENDRO_30 + At3[pp] * DENDRO_40 +
                         At5[pp] * (DENDRO_47 * DENDRO_47) +
                         DENDRO_36 * DENDRO_50 - DENDRO_44 * DENDRO_50 -
                         DENDRO_48 * DENDRO_49;
const double DENDRO_52 = At5[pp] * DENDRO_51;
const double DENDRO_53 = 24 * DENDRO_25;
const double DENDRO_54 = At3[pp] * DENDRO_27;
const double DENDRO_55 = At0[pp] * DENDRO_31;
const double DENDRO_56 = At5[pp] * DENDRO_29;
const double DENDRO_57 = At2[pp] * DENDRO_47;
const double DENDRO_58 = At4[pp] * DENDRO_47;
const double DENDRO_59 =
    At2[pp] * DENDRO_30 - DENDRO_27 * DENDRO_58 - DENDRO_29 * DENDRO_34 +
    DENDRO_29 * DENDRO_55 - DENDRO_31 * DENDRO_48 + DENDRO_31 * DENDRO_57 -
    DENDRO_32 * DENDRO_39 + DENDRO_39 * DENDRO_54 + DENDRO_47 * DENDRO_56;
const double DENDRO_60 = At2[pp] * DENDRO_59;
const double DENDRO_61 =
    -At1[pp] * DENDRO_28 - At1[pp] * DENDRO_31 * DENDRO_41 -
    At4[pp] * DENDRO_27 * DENDRO_39 - At4[pp] * DENDRO_29 * DENDRO_41 +
    DENDRO_27 * DENDRO_36 + DENDRO_27 * DENDRO_55 + DENDRO_31 * DENDRO_42 +
    DENDRO_39 * DENDRO_56 + DENDRO_41 * DENDRO_54;
const double DENDRO_62 = -DENDRO_61;
const double DENDRO_63 = At1[pp] * DENDRO_62;
const double DENDRO_64 = At0[pp] * DENDRO_27;
const double DENDRO_65 = At3[pp] * DENDRO_41;
const double DENDRO_66 = At5[pp] * DENDRO_39;
const double DENDRO_67 =
    -At1[pp] * DENDRO_27 * DENDRO_39 - At1[pp] * DENDRO_29 * DENDRO_41 -
    At4[pp] * DENDRO_40 - At4[pp] * DENDRO_41 * DENDRO_47 +
    DENDRO_27 * DENDRO_57 + DENDRO_29 * DENDRO_64 + DENDRO_36 * DENDRO_39 +
    DENDRO_39 * DENDRO_65 + DENDRO_47 * DENDRO_66;
const double DENDRO_68  = -DENDRO_67;
const double DENDRO_69  = At4[pp] * DENDRO_68;
const double DENDRO_70  = 4.0 * gt2[pp];
const double DENDRO_71  = 4.0 * gt4[pp];
const double DENDRO_72  = pow(chi[pp], -2);
const double DENDRO_73  = pow(grad_2_chi[pp], 2);
const double DENDRO_74  = DENDRO_72 * DENDRO_73;
const double DENDRO_75  = 1.0 / DENDRO_24;
const double DENDRO_76  = 4.0 * DENDRO_75;
const double DENDRO_77  = DENDRO_27 * DENDRO_76;
const double DENDRO_78  = DENDRO_77 * grad2_0_1_gt5[pp];
const double DENDRO_79  = DENDRO_29 * DENDRO_76;
const double DENDRO_80  = DENDRO_39 * DENDRO_76;
const double DENDRO_81  = DENDRO_80 * grad2_1_2_gt5[pp];
const double DENDRO_82  = DENDRO_31 * DENDRO_75;
const double DENDRO_83  = 2.0 * DENDRO_82;
const double DENDRO_84  = DENDRO_41 * DENDRO_75;
const double DENDRO_85  = 2.0 * DENDRO_84;
const double DENDRO_86  = DENDRO_47 * DENDRO_75;
const double DENDRO_87  = 2.0 * DENDRO_86;
const double DENDRO_88  = DENDRO_29 * grad_2_gt0[pp];
const double DENDRO_89  = DENDRO_47 * grad_0_gt5[pp];
const double DENDRO_90  = grad_0_gt4[pp] - grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_91  = -DENDRO_39 * DENDRO_90 + DENDRO_88 + DENDRO_89;
const double DENDRO_92  = -DENDRO_91;
const double DENDRO_93  = DENDRO_92 * grad_0_gt5[pp];
const double DENDRO_94  = 3.0 * DENDRO_25;
const double DENDRO_95  = DENDRO_31 * DENDRO_94;
const double DENDRO_96  = DENDRO_47 * grad_1_gt5[pp];
const double DENDRO_97  = -grad_0_gt4[pp] + grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_98  = DENDRO_29 * DENDRO_97;
const double DENDRO_99  = -DENDRO_39 * grad_2_gt3[pp] + DENDRO_96 + DENDRO_98;
const double DENDRO_100 = -DENDRO_99;
const double DENDRO_101 = DENDRO_100 * grad_1_gt5[pp];
const double DENDRO_102 = DENDRO_41 * DENDRO_94;
const double DENDRO_103 = 0.5 * grad_2_gt5[pp];
const double DENDRO_104 = DENDRO_103 * DENDRO_47;
const double DENDRO_105 = 0.5 * grad_1_gt5[pp];
const double DENDRO_106 = 1.0 * grad_2_gt4[pp];
const double DENDRO_107 = -DENDRO_106;
const double DENDRO_108 = DENDRO_105 + DENDRO_107;
const double DENDRO_109 = DENDRO_108 * DENDRO_39;
const double DENDRO_110 = 0.5 * grad_0_gt5[pp];
const double DENDRO_111 = 1.0 * grad_2_gt2[pp];
const double DENDRO_112 = -DENDRO_111;
const double DENDRO_113 = DENDRO_110 + DENDRO_112;
const double DENDRO_114 = DENDRO_104 + DENDRO_109 - DENDRO_113 * DENDRO_29;
const double DENDRO_115 = -DENDRO_114;
const double DENDRO_116 = DENDRO_115 * DENDRO_47;
const double DENDRO_117 = 6.0 * DENDRO_25;
const double DENDRO_118 = DENDRO_41 * grad_2_gt3[pp];
const double DENDRO_119 = DENDRO_39 * grad_1_gt5[pp];
const double DENDRO_120 = DENDRO_27 * DENDRO_97;
const double DENDRO_121 = DENDRO_119 + DENDRO_120;
const double DENDRO_122 = -DENDRO_118 + DENDRO_121;
const double DENDRO_123 = 1.0 * grad_1_gt4[pp];
const double DENDRO_124 = 0.25 * grad_2_gt3[pp];
const double DENDRO_125 = 4 * DENDRO_25;
const double DENDRO_126 = DENDRO_125 * DENDRO_41;
const double DENDRO_127 = DENDRO_122 * DENDRO_126 * (DENDRO_123 - DENDRO_124);
const double DENDRO_128 = DENDRO_29 * grad_0_gt5[pp];
const double DENDRO_129 = DENDRO_31 * grad_2_gt0[pp];
const double DENDRO_130 = DENDRO_128 + DENDRO_129 - DENDRO_27 * DENDRO_90;
const double DENDRO_131 = -DENDRO_130;
const double DENDRO_132 = 1.0 * grad_0_gt2[pp];
const double DENDRO_133 = 0.25 * grad_2_gt0[pp];
const double DENDRO_134 = DENDRO_132 - DENDRO_133;
const double DENDRO_135 = DENDRO_125 * DENDRO_31;
const double DENDRO_136 =
    DENDRO_27 * grad_2_gt0[pp] + DENDRO_39 * grad_0_gt5[pp];
const double DENDRO_137 = DENDRO_136 - DENDRO_41 * DENDRO_90;
const double DENDRO_138 = 0.75 * grad_0_gt4[pp];
const double DENDRO_139 = 0.25 * grad_1_gt2[pp];
const double DENDRO_140 = 0.25 * grad_2_gt1[pp];
const double DENDRO_141 = -DENDRO_140;
const double DENDRO_142 =
    DENDRO_135 * DENDRO_137 * (DENDRO_138 + DENDRO_139 + DENDRO_141);
const double DENDRO_143 = -DENDRO_27 * grad_2_gt3[pp] +
                          DENDRO_29 * grad_1_gt5[pp] + DENDRO_31 * DENDRO_97;
const double DENDRO_144 = -DENDRO_143;
const double DENDRO_145 = 0.25 * grad_0_gt4[pp];
const double DENDRO_146 = 0.75 * grad_1_gt2[pp];
const double DENDRO_147 = DENDRO_141 + DENDRO_145 + DENDRO_146;
const double DENDRO_148 = DENDRO_110 + DENDRO_111;
const double DENDRO_149 =
    DENDRO_103 * DENDRO_29 + DENDRO_108 * DENDRO_27 - DENDRO_113 * DENDRO_31;
const double DENDRO_150 = -DENDRO_149;
const double DENDRO_151 = DENDRO_150 * DENDRO_47;
const double DENDRO_152 = DENDRO_103 * DENDRO_39;
const double DENDRO_153 =
    DENDRO_108 * DENDRO_41 - DENDRO_113 * DENDRO_27 + DENDRO_152;
const double DENDRO_154 = DENDRO_153 * DENDRO_47;
const double DENDRO_155 = DENDRO_125 * DENDRO_154 * (DENDRO_105 + DENDRO_106);
const double DENDRO_156 = 0.25 * grad_0_gt5[pp];
const double DENDRO_157 = DENDRO_100 * DENDRO_156;
const double DENDRO_158 = DENDRO_125 * DENDRO_27;
const double DENDRO_159 = 0.25 * grad_1_gt5[pp];
const double DENDRO_160 = DENDRO_159 * DENDRO_92;
const double DENDRO_161 = DENDRO_122 * grad_1_gt5[pp];
const double DENDRO_162 = DENDRO_153 * grad_2_gt3[pp];
const double DENDRO_163 = DENDRO_161 + DENDRO_162;
const double DENDRO_164 = 2.0 * DENDRO_25;
const double DENDRO_165 = DENDRO_164 * DENDRO_39;
const double DENDRO_166 = DENDRO_131 * grad_0_gt5[pp];
const double DENDRO_167 = DENDRO_150 * grad_2_gt0[pp];
const double DENDRO_168 = DENDRO_164 * DENDRO_29;
const double DENDRO_169 = DENDRO_92 * grad_2_gt5[pp];
const double DENDRO_170 = DENDRO_115 * grad_0_gt5[pp];
const double DENDRO_171 = DENDRO_100 * grad_2_gt5[pp];
const double DENDRO_172 = DENDRO_115 * grad_1_gt5[pp];
const double DENDRO_173 = DENDRO_137 * grad_2_gt3[pp];
const double DENDRO_174 = 0.25 * DENDRO_173;
const double DENDRO_175 = grad_0_gt4[pp] + grad_1_gt2[pp] - grad_2_gt1[pp];
const double DENDRO_176 = 0.5 * DENDRO_175;
const double DENDRO_177 = DENDRO_122 * DENDRO_176 + DENDRO_174;
const double DENDRO_178 = 0.25 * DENDRO_169;
const double DENDRO_179 = DENDRO_125 * DENDRO_29;
const double DENDRO_180 = 0.25 * DENDRO_171;
const double DENDRO_181 = DENDRO_125 * DENDRO_39;
const double DENDRO_182 = DENDRO_144 * grad_2_gt0[pp];
const double DENDRO_183 = 0.25 * DENDRO_182;
const double DENDRO_184 = DENDRO_137 * grad_1_gt5[pp];
const double DENDRO_185 = DENDRO_153 * DENDRO_90;
const double DENDRO_186 = DENDRO_168 * (DENDRO_184 + DENDRO_185);
const double DENDRO_187 = DENDRO_144 * grad_0_gt5[pp];
const double DENDRO_188 = DENDRO_150 * DENDRO_97;
const double DENDRO_189 = 0.5 * grad_2_gt3[pp];
const double DENDRO_190 = DENDRO_123 - DENDRO_189;
const double DENDRO_191 = 1.0 * DENDRO_137;
const double DENDRO_192 = DENDRO_122 * DENDRO_90;
const double DENDRO_193 = 0.25 * DENDRO_192;
const double DENDRO_194 = DENDRO_190 * DENDRO_191 + DENDRO_193;
const double DENDRO_195 = 0.5 * grad_2_gt0[pp];
const double DENDRO_196 = DENDRO_132 - DENDRO_195;
const double DENDRO_197 = 1.0 * DENDRO_144;
const double DENDRO_198 = DENDRO_131 * DENDRO_97;
const double DENDRO_199 = 0.25 * DENDRO_198;
const double DENDRO_200 = 0.5 * DENDRO_108;
const double DENDRO_201 = DENDRO_137 * DENDRO_200;
const double DENDRO_202 = -1.0 * DENDRO_153 * DENDRO_175 + DENDRO_201;
const double DENDRO_203 = -DENDRO_122 * DENDRO_200;
const double DENDRO_204 = 2 * DENDRO_153 * DENDRO_190 + DENDRO_203;
const double DENDRO_205 = 0.5 * DENDRO_144;
const double DENDRO_206 = DENDRO_113 * DENDRO_205;
const double DENDRO_207 = -1.0 * DENDRO_150 * DENDRO_175 + DENDRO_206;
const double DENDRO_208 = 0.5 * DENDRO_131;
const double DENDRO_209 = -DENDRO_113 * DENDRO_208;
const double DENDRO_210 = 2 * DENDRO_196;
const double DENDRO_211 = -grad2_2_2_chi[pp];
const double DENDRO_212 = -DENDRO_29;
const double DENDRO_213 = -DENDRO_113;
const double DENDRO_214 = -DENDRO_31;
const double DENDRO_215 = -DENDRO_108;
const double DENDRO_216 =
    DENDRO_103 * DENDRO_212 + DENDRO_213 * DENDRO_214 + DENDRO_215 * DENDRO_27;
const double DENDRO_217 = DENDRO_75 * grad_0_chi[pp];
const double DENDRO_218 = -DENDRO_41;
const double DENDRO_219 =
    DENDRO_152 + DENDRO_213 * DENDRO_27 + DENDRO_215 * DENDRO_218;
const double DENDRO_220 = DENDRO_75 * grad_1_chi[pp];
const double DENDRO_221 = -DENDRO_47;
const double DENDRO_222 = DENDRO_103 * DENDRO_221;
const double DENDRO_223 = DENDRO_212 * DENDRO_213;
const double DENDRO_224 = DENDRO_215 * DENDRO_39;
const double DENDRO_225 = DENDRO_75 * grad_2_chi[pp];
const double DENDRO_226 = 1.0 / chi[pp];
const double DENDRO_227 = 2 * DENDRO_226;
const double DENDRO_228 = 1.0 * DENDRO_29;
const double DENDRO_229 = DENDRO_41 * grad_0_gt3[pp];
const double DENDRO_230 = DENDRO_27 * grad_1_gt0[pp];
const double DENDRO_231 = DENDRO_175 * DENDRO_39;
const double DENDRO_232 = DENDRO_230 + DENDRO_231;
const double DENDRO_233 = -DENDRO_229 + DENDRO_232;
const double DENDRO_234 = 0.5 * grad_1_gt3[pp];
const double DENDRO_235 = DENDRO_234 * DENDRO_41;
const double DENDRO_236 = 0.5 * grad_0_gt3[pp];
const double DENDRO_237 = 1.0 * grad_1_gt1[pp];
const double DENDRO_238 = -DENDRO_237;
const double DENDRO_239 = DENDRO_236 + DENDRO_238;
const double DENDRO_240 = DENDRO_239 * DENDRO_27;
const double DENDRO_241 = -DENDRO_190 * DENDRO_39 + DENDRO_235 + DENDRO_240;
const double DENDRO_242 = -DENDRO_241;
const double DENDRO_243 = DENDRO_242 * DENDRO_41;
const double DENDRO_244 = 1.0 * grad_0_gt1[pp];
const double DENDRO_245 = 0.5 * grad_1_gt0[pp];
const double DENDRO_246 = DENDRO_244 - DENDRO_245;
const double DENDRO_247 = 0.5 * grad_0_gt0[pp];
const double DENDRO_248 = DENDRO_196 * DENDRO_39 + DENDRO_247 * DENDRO_27;
const double DENDRO_249 = -DENDRO_246 * DENDRO_41 + DENDRO_248;
const double DENDRO_250 = DENDRO_249 * DENDRO_31;
const double DENDRO_251 =
    -1.0 * DENDRO_122 * DENDRO_39 + DENDRO_137 * DENDRO_228 + DENDRO_154 -
    1.0 * DENDRO_233 * DENDRO_27 + DENDRO_243 + DENDRO_250;
const double DENDRO_252 = -DENDRO_251;
const double DENDRO_253 = DENDRO_25 * DENDRO_252;
const double DENDRO_254 = 2.0 * DENDRO_253;
const double DENDRO_255 = DENDRO_31 * grad_1_gt0[pp];
const double DENDRO_256 = DENDRO_175 * DENDRO_29;
const double DENDRO_257 = DENDRO_255 + DENDRO_256 - DENDRO_27 * grad_0_gt3[pp];
const double DENDRO_258 = -DENDRO_257;
const double DENDRO_259 = DENDRO_234 * DENDRO_27;
const double DENDRO_260 =
    -DENDRO_190 * DENDRO_29 + DENDRO_239 * DENDRO_31 + DENDRO_259;
const double DENDRO_261 = DENDRO_260 * DENDRO_41;
const double DENDRO_262 = DENDRO_247 * DENDRO_31;
const double DENDRO_263 = DENDRO_196 * DENDRO_29;
const double DENDRO_264 = -DENDRO_246 * DENDRO_27 + DENDRO_262 + DENDRO_263;
const double DENDRO_265 = -DENDRO_264;
const double DENDRO_266 = DENDRO_265 * DENDRO_31;
const double DENDRO_267 =
    DENDRO_131 * DENDRO_228 - 1.0 * DENDRO_144 * DENDRO_39 + DENDRO_151 -
    1.0 * DENDRO_258 * DENDRO_27 + DENDRO_261 + DENDRO_266;
const double DENDRO_268 = -DENDRO_267;
const double DENDRO_269 = DENDRO_25 * DENDRO_268;
const double DENDRO_270 = 2.0 * DENDRO_269;
const double DENDRO_271 = DENDRO_29 * grad_1_gt0[pp];
const double DENDRO_272 = DENDRO_175 * DENDRO_47;
const double DENDRO_273 = DENDRO_271 + DENDRO_272 - DENDRO_39 * grad_0_gt3[pp];
const double DENDRO_274 = -DENDRO_273;
const double DENDRO_275 = DENDRO_234 * DENDRO_39;
const double DENDRO_276 =
    -DENDRO_190 * DENDRO_47 + DENDRO_239 * DENDRO_29 + DENDRO_275;
const double DENDRO_277 = DENDRO_276 * DENDRO_41;
const double DENDRO_278 =
    DENDRO_196 * DENDRO_47 - DENDRO_246 * DENDRO_39 + DENDRO_247 * DENDRO_29;
const double DENDRO_279 = -DENDRO_278;
const double DENDRO_280 = DENDRO_279 * DENDRO_31;
const double DENDRO_281 =
    -1.0 * DENDRO_100 * DENDRO_39 + DENDRO_116 + DENDRO_228 * DENDRO_92 -
    1.0 * DENDRO_27 * DENDRO_274 + DENDRO_277 + DENDRO_280;
const double DENDRO_282 = -DENDRO_281;
const double DENDRO_283 = DENDRO_25 * DENDRO_282;
const double DENDRO_284 = 2.0 * DENDRO_283;
const double DENDRO_285 = 3 * DENDRO_226;
const double DENDRO_286 = DENDRO_285 * grad_1_chi[pp];
const double DENDRO_287 = grad_0_chi[pp] * grad_2_chi[pp];
const double DENDRO_288 = pow(grad_0_chi[pp], 2);
const double DENDRO_289 = pow(grad_1_chi[pp], 2);
const double DENDRO_290 =
    2 * DENDRO_217 * DENDRO_268 + 2 * DENDRO_220 * DENDRO_252 +
    2 * DENDRO_225 * DENDRO_282 -
    2 * DENDRO_27 * (-DENDRO_286 * grad_0_chi[pp] + 2 * grad2_0_1_chi[pp]) +
    DENDRO_31 * (-DENDRO_285 * DENDRO_288 + 2 * grad2_0_0_chi[pp]) -
    2 * DENDRO_39 * (-DENDRO_286 * grad_2_chi[pp] + 2 * grad2_1_2_chi[pp]) +
    DENDRO_41 * (-DENDRO_285 * DENDRO_289 + 2 * grad2_1_1_chi[pp]) +
    DENDRO_47 * (-DENDRO_285 * DENDRO_73 + 2 * grad2_2_2_chi[pp]) +
    DENDRO_49 * (-DENDRO_285 * DENDRO_287 + 2 * grad2_0_2_chi[pp]);
const double DENDRO_291 = -DENDRO_226 * DENDRO_290;
const double DENDRO_292 = DENDRO_291 * DENDRO_75;
const double DENDRO_293 =
    -DENDRO_101 * DENDRO_102 - DENDRO_116 * DENDRO_117 * grad_2_gt5[pp] -
    DENDRO_125 * DENDRO_148 * DENDRO_151 -
    DENDRO_126 * DENDRO_144 * DENDRO_147 - DENDRO_127 -
    DENDRO_131 * DENDRO_134 * DENDRO_135 - DENDRO_142 - DENDRO_155 +
    DENDRO_158 * DENDRO_177 + DENDRO_158 * DENDRO_194 +
    DENDRO_158 * (DENDRO_100 * DENDRO_110 + DENDRO_160) +
    DENDRO_158 * (DENDRO_105 * DENDRO_92 + DENDRO_157) +
    DENDRO_158 * (DENDRO_131 * DENDRO_176 + DENDRO_183) +
    DENDRO_158 * (DENDRO_196 * DENDRO_197 + DENDRO_199) +
    DENDRO_163 * DENDRO_165 + DENDRO_165 * (DENDRO_171 + DENDRO_172) +
    DENDRO_165 * (DENDRO_187 + DENDRO_188) -
    DENDRO_168 * (DENDRO_166 + DENDRO_167) -
    DENDRO_168 * (DENDRO_169 + DENDRO_170) + DENDRO_179 * DENDRO_202 -
    DENDRO_179 * (1.0 * DENDRO_170 + DENDRO_178) -
    DENDRO_179 * (DENDRO_150 * DENDRO_210 + DENDRO_209) +
    DENDRO_181 * DENDRO_204 - DENDRO_181 * DENDRO_207 +
    DENDRO_181 * (1.0 * DENDRO_172 + DENDRO_180) - DENDRO_186 -
    DENDRO_227 *
        (DENDRO_211 + DENDRO_216 * DENDRO_217 + DENDRO_219 * DENDRO_220 +
         DENDRO_225 * (DENDRO_222 + DENDRO_223 + DENDRO_224)) +
    DENDRO_254 * grad_1_gt5[pp] + DENDRO_270 * grad_0_gt5[pp] +
    DENDRO_284 * grad_2_gt5[pp] + DENDRO_292 * gt5[pp] +
    DENDRO_70 * grad_2_Gt0[pp] + DENDRO_71 * grad_2_Gt1[pp] - DENDRO_74 -
    DENDRO_78 + DENDRO_79 * grad2_0_2_gt5[pp] - DENDRO_81 +
    DENDRO_83 * grad2_0_0_gt5[pp] + DENDRO_85 * grad2_1_1_gt5[pp] +
    DENDRO_87 * grad2_2_2_gt5[pp] - DENDRO_93 * DENDRO_95 +
    4.0 * grad_2_Gt2[pp] * gt5[pp];
const double DENDRO_294 = DENDRO_86 * chi[pp];
const double DENDRO_295 = -grad2_1_1_chi[pp];
const double DENDRO_296 = -DENDRO_239;
const double DENDRO_297 =
    DENDRO_190 * DENDRO_212 + DENDRO_214 * DENDRO_296 + DENDRO_259;
const double DENDRO_298 = DENDRO_218 * DENDRO_234;
const double DENDRO_299 = DENDRO_27 * DENDRO_296;
const double DENDRO_300 = DENDRO_190 * DENDRO_39;
const double DENDRO_301 =
    DENDRO_190 * DENDRO_221 + DENDRO_212 * DENDRO_296 + DENDRO_275;
const double DENDRO_302 = DENDRO_144 * grad_1_gt0[pp];
const double DENDRO_303 = 0.25 * DENDRO_302;
const double DENDRO_304 = 0.5 * DENDRO_90;
const double DENDRO_305 = DENDRO_274 * grad_1_gt5[pp];
const double DENDRO_306 = 0.25 * DENDRO_305;
const double DENDRO_307 = 1.0 * DENDRO_274;
const double DENDRO_308 = -0.25 * DENDRO_100 * DENDRO_175;
const double DENDRO_309 = DENDRO_108 * DENDRO_307 + DENDRO_308;
const double DENDRO_310 = DENDRO_258 * DENDRO_97;
const double DENDRO_311 = 0.25 * DENDRO_310;
const double DENDRO_312 = 0.5 * DENDRO_190;
const double DENDRO_313 = DENDRO_100 * DENDRO_312;
const double DENDRO_314 = 2 * DENDRO_108 * DENDRO_276;
const double DENDRO_315 = DENDRO_122 * grad_1_gt3[pp];
const double DENDRO_316 = 0.25 * DENDRO_315;
const double DENDRO_317 = DENDRO_242 * grad_2_gt3[pp];
const double DENDRO_318 = DENDRO_205 * DENDRO_239;
const double DENDRO_319 = -1.0 * DENDRO_260 * DENDRO_90 + DENDRO_318;
const double DENDRO_320 = 0.5 * DENDRO_258;
const double DENDRO_321 = -DENDRO_239 * DENDRO_320;
const double DENDRO_322 = 2 * DENDRO_246 * DENDRO_260;
const double DENDRO_323 = DENDRO_233 * grad_1_gt3[pp];
const double DENDRO_324 = 0.25 * DENDRO_323;
const double DENDRO_325 = DENDRO_242 * grad_0_gt3[pp];
const double DENDRO_326 = DENDRO_274 * DENDRO_312;
const double DENDRO_327 = DENDRO_258 * grad_0_gt3[pp];
const double DENDRO_328 = DENDRO_260 * grad_1_gt0[pp];
const double DENDRO_329 = DENDRO_164 * DENDRO_27;
const double DENDRO_330 = DENDRO_144 * grad_0_gt3[pp];
const double DENDRO_331 = DENDRO_260 * DENDRO_97;
const double DENDRO_332 = DENDRO_276 * grad_1_gt5[pp];
const double DENDRO_333 = DENDRO_175 * DENDRO_276;
const double DENDRO_334 = DENDRO_107 + DENDRO_159;
const double DENDRO_335 = DENDRO_125 * DENDRO_47;
const double DENDRO_336 = -DENDRO_139;
const double DENDRO_337 = 0.75 * grad_2_gt1[pp];
const double DENDRO_338 = DENDRO_335 * (DENDRO_145 + DENDRO_336 + DENDRO_337);
const double DENDRO_339 = 0.25 * grad_1_gt0[pp];
const double DENDRO_340 = DENDRO_135 * (DENDRO_244 - DENDRO_339);
const double DENDRO_341 = DENDRO_135 * (DENDRO_138 + DENDRO_140 + DENDRO_336);
const double DENDRO_342 = DENDRO_117 * grad_1_gt3[pp];
const double DENDRO_343 = 4.0 * gt1[pp];
const double DENDRO_344 = 0.25 * grad_0_gt3[pp];
const double DENDRO_345 = DENDRO_122 * DENDRO_344;
const double DENDRO_346 = DENDRO_124 * DENDRO_233;
const double DENDRO_347 = DENDRO_233 * grad_0_gt3[pp];
const double DENDRO_348 = DENDRO_122 * grad_2_gt3[pp];
const double DENDRO_349 = DENDRO_47 * DENDRO_94;
const double DENDRO_350 =
    -DENDRO_125 * DENDRO_261 * (DENDRO_236 + DENDRO_237) -
    DENDRO_125 * DENDRO_277 * (DENDRO_123 + DENDRO_189) -
    DENDRO_179 * (DENDRO_122 * DENDRO_236 + DENDRO_346) -
    DENDRO_179 * (DENDRO_189 * DENDRO_233 + DENDRO_345) -
    DENDRO_289 * DENDRO_72 + DENDRO_343 * grad_1_Gt0[pp] -
    DENDRO_347 * DENDRO_95 - DENDRO_348 * DENDRO_349 +
    DENDRO_71 * grad_1_Gt2[pp] - DENDRO_77 * grad2_0_1_gt3[pp] +
    DENDRO_79 * grad2_0_2_gt3[pp] - DENDRO_80 * grad2_1_2_gt3[pp] +
    DENDRO_83 * grad2_0_0_gt3[pp] + DENDRO_85 * grad2_1_1_gt3[pp] +
    DENDRO_87 * grad2_2_2_gt3[pp] + 4.0 * grad_1_Gt1[pp] * gt3[pp];
const double DENDRO_351 =
    DENDRO_100 * DENDRO_334 * DENDRO_335 - DENDRO_144 * DENDRO_338 +
    DENDRO_158 * (DENDRO_321 + DENDRO_322) +
    DENDRO_158 * (DENDRO_324 + 1.0 * DENDRO_325) +
    DENDRO_158 * (1.0 * DENDRO_276 * DENDRO_90 + DENDRO_326) +
    DENDRO_165 * (DENDRO_315 + DENDRO_317) +
    DENDRO_165 * (DENDRO_330 + DENDRO_331) +
    DENDRO_165 * (DENDRO_100 * grad_2_gt3[pp] + DENDRO_332) +
    DENDRO_179 * DENDRO_309 -
    DENDRO_179 * (DENDRO_100 * DENDRO_304 + DENDRO_306) -
    DENDRO_179 * (DENDRO_197 * DENDRO_246 + DENDRO_311) -
    DENDRO_179 * (DENDRO_258 * DENDRO_304 + DENDRO_303) -
    DENDRO_181 * DENDRO_319 + DENDRO_181 * (DENDRO_313 - DENDRO_314) +
    DENDRO_181 * (DENDRO_316 + 1.0 * DENDRO_317) -
    DENDRO_227 * (DENDRO_217 * DENDRO_297 +
                  DENDRO_220 * (DENDRO_298 + DENDRO_299 + DENDRO_300) +
                  DENDRO_225 * DENDRO_301 + DENDRO_295) -
    DENDRO_243 * DENDRO_342 + DENDRO_254 * grad_1_gt3[pp] -
    DENDRO_258 * DENDRO_340 + DENDRO_270 * grad_0_gt3[pp] -
    DENDRO_274 * DENDRO_341 + DENDRO_284 * grad_2_gt3[pp] +
    DENDRO_292 * gt3[pp] + DENDRO_329 * (DENDRO_323 + DENDRO_325) +
    DENDRO_329 * (DENDRO_327 + DENDRO_328) +
    DENDRO_329 * (DENDRO_274 * grad_2_gt3[pp] + DENDRO_333) + DENDRO_350;
const double DENDRO_352 = DENDRO_84 * chi[pp];
const double DENDRO_353 = DENDRO_288 * DENDRO_72;
const double DENDRO_354 = DENDRO_77 * grad2_0_1_gt0[pp];
const double DENDRO_355 = DENDRO_80 * grad2_1_2_gt0[pp];
const double DENDRO_356 = DENDRO_258 * grad_1_gt0[pp];
const double DENDRO_357 = DENDRO_131 * grad_2_gt0[pp];
const double DENDRO_358 = DENDRO_126 * DENDRO_233 * (-DENDRO_238 - DENDRO_344);
const double DENDRO_359 = DENDRO_112 + DENDRO_156;
const double DENDRO_360 = -DENDRO_145;
const double DENDRO_361 =
    DENDRO_137 * DENDRO_335 * (DENDRO_139 + DENDRO_337 + DENDRO_360);
const double DENDRO_362 = DENDRO_140 + DENDRO_146 + DENDRO_360;
const double DENDRO_363 = DENDRO_125 * DENDRO_250 * (DENDRO_244 + DENDRO_245);
const double DENDRO_364 = DENDRO_132 + DENDRO_195;
const double DENDRO_365 = DENDRO_131 * DENDRO_339;
const double DENDRO_366 = DENDRO_133 * DENDRO_258;
const double DENDRO_367 = DENDRO_249 * grad_0_gt3[pp];
const double DENDRO_368 = DENDRO_233 * grad_1_gt0[pp] + DENDRO_367;
const double DENDRO_369 = DENDRO_258 * grad_0_gt0[pp];
const double DENDRO_370 = DENDRO_265 * grad_1_gt0[pp];
const double DENDRO_371 = DENDRO_131 * grad_0_gt0[pp];
const double DENDRO_372 = DENDRO_265 * grad_2_gt0[pp];
const double DENDRO_373 = DENDRO_137 * grad_0_gt3[pp];
const double DENDRO_374 = 0.25 * DENDRO_373;
const double DENDRO_375 = 0.5 * DENDRO_97;
const double DENDRO_376 = DENDRO_233 * DENDRO_375 + DENDRO_374;
const double DENDRO_377 = DENDRO_279 * grad_0_gt5[pp];
const double DENDRO_378 = 0.25 * DENDRO_369;
const double DENDRO_379 = 0.25 * DENDRO_371;
const double DENDRO_380 = DENDRO_274 * grad_0_gt5[pp];
const double DENDRO_381 = 0.25 * DENDRO_380;
const double DENDRO_382 = DENDRO_249 * DENDRO_90;
const double DENDRO_383 =
    DENDRO_168 * (DENDRO_137 * grad_1_gt0[pp] + DENDRO_382);
const double DENDRO_384 = DENDRO_175 * DENDRO_279;
const double DENDRO_385 =
    DENDRO_191 * DENDRO_239 - 0.25 * DENDRO_233 * DENDRO_90;
const double DENDRO_386 = -0.25 * DENDRO_175 * DENDRO_92;
const double DENDRO_387 = DENDRO_113 * DENDRO_307 + DENDRO_386;
const double DENDRO_388 = 0.5 * DENDRO_246;
const double DENDRO_389 = DENDRO_137 * DENDRO_388;
const double DENDRO_390 = 1.0 * DENDRO_97;
const double DENDRO_391 = DENDRO_179 * (DENDRO_249 * DENDRO_390 + DENDRO_389);
const double DENDRO_392 = DENDRO_233 * DENDRO_388;
const double DENDRO_393 = -2 * DENDRO_239 * DENDRO_249 + DENDRO_392;
const double DENDRO_394 = 0.5 * DENDRO_274;
const double DENDRO_395 = DENDRO_196 * DENDRO_394;
const double DENDRO_396 = 0.5 * DENDRO_92;
const double DENDRO_397 = DENDRO_196 * DENDRO_396;
const double DENDRO_398 = -grad2_0_0_chi[pp];
const double DENDRO_399 = DENDRO_214 * DENDRO_247;
const double DENDRO_400 = DENDRO_246 * DENDRO_27;
const double DENDRO_401 = DENDRO_196 * DENDRO_212;
const double DENDRO_402 = DENDRO_218 * DENDRO_246 + DENDRO_248;
const double DENDRO_403 =
    DENDRO_196 * DENDRO_221 + DENDRO_212 * DENDRO_247 + DENDRO_246 * DENDRO_39;
const double DENDRO_404 =
    -DENDRO_102 * DENDRO_356 - DENDRO_117 * DENDRO_266 * grad_0_gt0[pp] -
    DENDRO_125 * DENDRO_280 * DENDRO_364 -
    DENDRO_126 * DENDRO_274 * DENDRO_362 + DENDRO_158 * DENDRO_393 +
    DENDRO_158 * (1.0 * DENDRO_370 + DENDRO_378) +
    DENDRO_158 * (DENDRO_279 * DENDRO_390 + DENDRO_395) -
    DENDRO_168 * (DENDRO_371 + DENDRO_372) -
    DENDRO_168 * (DENDRO_377 + DENDRO_92 * grad_2_gt0[pp]) -
    DENDRO_179 * (1.0 * DENDRO_372 + DENDRO_379) -
    DENDRO_179 * (-2 * DENDRO_113 * DENDRO_279 + DENDRO_397) +
    DENDRO_181 * DENDRO_376 - DENDRO_181 * DENDRO_385 -
    DENDRO_181 * DENDRO_387 +
    DENDRO_181 * (DENDRO_131 * DENDRO_245 + DENDRO_366) +
    DENDRO_181 * (DENDRO_195 * DENDRO_258 + DENDRO_365) +
    DENDRO_181 * (DENDRO_375 * DENDRO_92 + DENDRO_381) -
    DENDRO_227 *
        (DENDRO_217 * (DENDRO_399 + DENDRO_400 + DENDRO_401) +
         DENDRO_220 * DENDRO_402 + DENDRO_225 * DENDRO_403 + DENDRO_398) +
    DENDRO_254 * grad_1_gt0[pp] + DENDRO_270 * grad_0_gt0[pp] +
    DENDRO_284 * grad_2_gt0[pp] + DENDRO_292 * gt0[pp] +
    DENDRO_329 * DENDRO_368 + DENDRO_329 * (DENDRO_369 + DENDRO_370) +
    DENDRO_329 * (DENDRO_274 * grad_2_gt0[pp] + DENDRO_384) +
    DENDRO_335 * DENDRO_359 * DENDRO_92 + DENDRO_343 * grad_0_Gt1[pp] -
    DENDRO_349 * DENDRO_357 - DENDRO_353 - DENDRO_354 - DENDRO_355 -
    DENDRO_358 - DENDRO_361 - DENDRO_363 - DENDRO_383 - DENDRO_391 +
    DENDRO_70 * grad_0_Gt2[pp] + DENDRO_79 * grad2_0_2_gt0[pp] +
    DENDRO_83 * grad2_0_0_gt0[pp] + DENDRO_85 * grad2_1_1_gt0[pp] +
    DENDRO_87 * grad2_2_2_gt0[pp] + 4.0 * grad_0_Gt0[pp] * gt0[pp];
const double DENDRO_405 = DENDRO_82 * chi[pp];
const double DENDRO_406 = 0.25 * DENDRO_166;
const double DENDRO_407 = -DENDRO_113 * DENDRO_115;
const double DENDRO_408 = DENDRO_153 * DENDRO_375 + 0.25 * DENDRO_184;
const double DENDRO_409 = DENDRO_175 * DENDRO_258;
const double DENDRO_410 = 0.25 * DENDRO_409;
const double DENDRO_411 = DENDRO_233 * DENDRO_312;
const double DENDRO_412 = 0.5 * DENDRO_239;
const double DENDRO_413 = -DENDRO_122 * DENDRO_412;
const double DENDRO_414 = DENDRO_411 + DENDRO_413;
const double DENDRO_415 = DENDRO_196 * DENDRO_265;
const double DENDRO_416 = DENDRO_137 * DENDRO_339 + DENDRO_176 * DENDRO_249;
const double DENDRO_417 = DENDRO_133 * DENDRO_92;
const double DENDRO_418 = 0.25 * DENDRO_137;
const double DENDRO_419 = DENDRO_105 * DENDRO_249 + DENDRO_418 * DENDRO_97;
const double DENDRO_420 = DENDRO_103 * DENDRO_279;
const double DENDRO_421 = -DENDRO_113 * DENDRO_396;
const double DENDRO_422 = DENDRO_131 * DENDRO_133;
const double DENDRO_423 = DENDRO_150 * DENDRO_247;
const double DENDRO_424 = DENDRO_196 * DENDRO_208 + DENDRO_423;
const double DENDRO_425 = DENDRO_418 * DENDRO_90;
const double DENDRO_426 = 0.25 * DENDRO_175;
const double DENDRO_427 = DENDRO_137 * DENDRO_426 + DENDRO_153 * DENDRO_245;
const double DENDRO_428 = DENDRO_153 * DENDRO_239;
const double DENDRO_429 = DENDRO_258 * grad_0_gt5[pp];
const double DENDRO_430 = DENDRO_150 * DENDRO_245;
const double DENDRO_431 = 0.25 * DENDRO_429 + DENDRO_430;
const double DENDRO_432 = DENDRO_115 * DENDRO_375;
const double DENDRO_433 = DENDRO_157 + DENDRO_160;
const double DENDRO_434 = DENDRO_122 * DENDRO_97;
const double DENDRO_435 = 0.25 * DENDRO_434;
const double DENDRO_436 = DENDRO_153 * DENDRO_236;
const double DENDRO_437 = DENDRO_233 * grad_1_gt5[pp];
const double DENDRO_438 = DENDRO_436 + 0.25 * DENDRO_437;
const double DENDRO_439 = DENDRO_131 * DENDRO_426 + DENDRO_430;
const double DENDRO_440 = 0.5 * DENDRO_100;
const double DENDRO_441 = DENDRO_113 * DENDRO_440;
const double DENDRO_442 = -DENDRO_441;
const double DENDRO_443 = 0.25 * DENDRO_274;
const double DENDRO_444 = DENDRO_443 * grad_2_gt5[pp];
const double DENDRO_445 = DENDRO_115 * DENDRO_176 + DENDRO_444;
const double DENDRO_446 = DENDRO_122 * DENDRO_388;
const double DENDRO_447 = -DENDRO_137 * DENDRO_412 + DENDRO_190 * DENDRO_249;
const double DENDRO_448 = DENDRO_196 * DENDRO_320;
const double DENDRO_449 = 0.25 * DENDRO_144;
const double DENDRO_450 = DENDRO_449 * grad_0_gt0[pp];
const double DENDRO_451 = DENDRO_365 + DENDRO_450;
const double DENDRO_452 = DENDRO_176 * DENDRO_265;
const double DENDRO_453 = DENDRO_105 * DENDRO_279;
const double DENDRO_454 = 0.25 * DENDRO_97;
const double DENDRO_455 = DENDRO_453 + DENDRO_454 * DENDRO_92;
const double DENDRO_456 = DENDRO_122 * grad_1_gt0[pp];
const double DENDRO_457 = DENDRO_175 * DENDRO_233;
const double DENDRO_458 = 1.0 * DENDRO_25 * DENDRO_27;
const double DENDRO_459 = DENDRO_100 * grad_2_gt0[pp];
const double DENDRO_460 = DENDRO_175 * DENDRO_92;
const double DENDRO_461 = DENDRO_100 * DENDRO_97;
const double DENDRO_462 = DENDRO_100 * DENDRO_175;
const double DENDRO_463 = 1.0 * DENDRO_25;
const double DENDRO_464 = DENDRO_41 * DENDRO_463;
const double DENDRO_465 = -grad2_0_2_chi[pp];
const double DENDRO_466 = DENDRO_212 * grad_0_gt5[pp];
const double DENDRO_467 = DENDRO_214 * grad_2_gt0[pp];
const double DENDRO_468 = DENDRO_27 * DENDRO_90;
const double DENDRO_469 = 0.5 * DENDRO_217;
const double DENDRO_470 = DENDRO_136 + DENDRO_218 * DENDRO_90;
const double DENDRO_471 = 0.5 * DENDRO_220;
const double DENDRO_472 = DENDRO_221 * grad_0_gt5[pp];
const double DENDRO_473 = DENDRO_212 * grad_2_gt0[pp];
const double DENDRO_474 = DENDRO_39 * DENDRO_90;
const double DENDRO_475 = 0.5 * DENDRO_225;
const double DENDRO_476 = 2.0 * gt2[pp];
const double DENDRO_477 = 2.0 * gt4[pp];
const double DENDRO_478 = 2.0 * gt1[pp];
const double DENDRO_479 = DENDRO_287 * DENDRO_72;
const double DENDRO_480 = DENDRO_77 * grad2_0_1_gt2[pp];
const double DENDRO_481 = DENDRO_80 * grad2_1_2_gt2[pp];
const double DENDRO_482 = DENDRO_75 * gt2[pp];
const double DENDRO_483 =
    -DENDRO_227 *
        (DENDRO_465 + DENDRO_469 * (DENDRO_466 + DENDRO_467 + DENDRO_468) +
         DENDRO_470 * DENDRO_471 +
         DENDRO_475 * (DENDRO_472 + DENDRO_473 + DENDRO_474)) +
    DENDRO_254 * grad_1_gt2[pp] + DENDRO_270 * grad_0_gt2[pp] +
    DENDRO_284 * grad_2_gt2[pp] + DENDRO_291 * DENDRO_482 +
    DENDRO_476 * grad_0_Gt0[pp] + DENDRO_476 * grad_2_Gt2[pp] +
    DENDRO_477 * grad_0_Gt1[pp] + DENDRO_478 * grad_2_Gt1[pp] - DENDRO_479 -
    DENDRO_480 - DENDRO_481 + DENDRO_79 * grad2_0_2_gt2[pp] +
    DENDRO_83 * grad2_0_0_gt2[pp] + DENDRO_85 * grad2_1_1_gt2[pp] +
    DENDRO_87 * grad2_2_2_gt2[pp] + 2.0 * grad_0_Gt2[pp] * gt5[pp] +
    2.0 * grad_2_Gt0[pp] * gt0[pp];
const double DENDRO_484 =
    -DENDRO_126 * (DENDRO_345 + DENDRO_414) -
    DENDRO_126 * (DENDRO_144 * DENDRO_245 + DENDRO_410) -
    DENDRO_135 * (0.5 * DENDRO_371 + DENDRO_415) -
    DENDRO_135 * (DENDRO_389 + DENDRO_416) -
    DENDRO_135 * (DENDRO_110 * DENDRO_279 + DENDRO_397 + DENDRO_417) +
    DENDRO_158 * (DENDRO_446 + DENDRO_447) +
    DENDRO_158 * (DENDRO_448 + DENDRO_451) +
    DENDRO_158 * (DENDRO_451 + DENDRO_452) +
    DENDRO_158 * (DENDRO_196 * DENDRO_440 + DENDRO_455) -
    DENDRO_168 * (DENDRO_115 * grad_2_gt0[pp] + DENDRO_93) -
    DENDRO_179 * (DENDRO_422 + DENDRO_424) -
    DENDRO_179 * (DENDRO_425 + DENDRO_427) -
    DENDRO_179 * (DENDRO_153 * DENDRO_246 + DENDRO_419) -
    DENDRO_179 * (DENDRO_110 * DENDRO_265 + DENDRO_422 + DENDRO_423) -
    DENDRO_179 * (DENDRO_115 * DENDRO_196 + DENDRO_420 + DENDRO_421) +
    DENDRO_181 * (DENDRO_183 + DENDRO_431) +
    DENDRO_181 * (DENDRO_183 + DENDRO_439) +
    DENDRO_181 * (DENDRO_432 + DENDRO_433) +
    DENDRO_181 * (DENDRO_435 + DENDRO_438) +
    DENDRO_181 * (DENDRO_442 + DENDRO_445) +
    DENDRO_181 * (DENDRO_137 * DENDRO_312 + DENDRO_193 - DENDRO_428) -
    DENDRO_335 * (1.0 * DENDRO_167 + DENDRO_406) -
    DENDRO_335 * (0.5 * DENDRO_185 + DENDRO_408) -
    DENDRO_335 * (DENDRO_110 * DENDRO_115 + DENDRO_178 + DENDRO_407) +
    DENDRO_458 * (DENDRO_373 + DENDRO_456 + DENDRO_457) +
    DENDRO_458 * (DENDRO_380 + DENDRO_459 + DENDRO_460) -
    DENDRO_464 * (DENDRO_305 + DENDRO_461 + DENDRO_462) + DENDRO_483;
const double DENDRO_485 = DENDRO_29 * DENDRO_75;
const double DENDRO_486 = 3 * chi[pp];
const double DENDRO_487 = DENDRO_485 * DENDRO_486;
const double DENDRO_488 = 2.0 * grad_1_Gt0[pp];
const double DENDRO_489 = DENDRO_488 * gt2[pp];
const double DENDRO_490 = DENDRO_477 * grad_1_Gt1[pp];
const double DENDRO_491 = 2.0 * grad_1_Gt2[pp];
const double DENDRO_492 = DENDRO_491 * gt5[pp];
const double DENDRO_493 = DENDRO_478 * grad_2_Gt0[pp];
const double DENDRO_494 = 2.0 * gt3[pp];
const double DENDRO_495 = DENDRO_494 * grad_2_Gt1[pp];
const double DENDRO_496 = DENDRO_477 * grad_2_Gt2[pp];
const double DENDRO_497 = DENDRO_72 * grad_1_chi[pp];
const double DENDRO_498 = -DENDRO_497 * grad_2_chi[pp];
const double DENDRO_499 = DENDRO_83 * grad2_0_0_gt4[pp];
const double DENDRO_500 = DENDRO_85 * grad2_1_1_gt4[pp];
const double DENDRO_501 = DENDRO_87 * grad2_2_2_gt4[pp];
const double DENDRO_502 = DENDRO_79 * grad2_0_2_gt4[pp];
const double DENDRO_503 = -DENDRO_77 * grad2_0_1_gt4[pp];
const double DENDRO_504 = -DENDRO_80 * grad2_1_2_gt4[pp];
const double DENDRO_505 = 0.25 * DENDRO_161;
const double DENDRO_506 = -DENDRO_108 * DENDRO_115;
const double DENDRO_507 = DENDRO_150 * DENDRO_304 + 0.25 * DENDRO_187;
const double DENDRO_508 = DENDRO_190 * DENDRO_242;
const double DENDRO_509 = DENDRO_176 * DENDRO_260;
const double DENDRO_510 = 0.25 * DENDRO_330 + DENDRO_509;
const double DENDRO_511 = DENDRO_100 * DENDRO_124;
const double DENDRO_512 = 0.25 * DENDRO_457;
const double DENDRO_513 = DENDRO_208 * DENDRO_246;
const double DENDRO_514 = DENDRO_448 + DENDRO_513;
const double DENDRO_515 = DENDRO_150 * DENDRO_246 + DENDRO_196 * DENDRO_205;
const double DENDRO_516 = DENDRO_131 * DENDRO_90;
const double DENDRO_517 = 0.25 * DENDRO_516;
const double DENDRO_518 = DENDRO_115 * DENDRO_304;
const double DENDRO_519 = DENDRO_122 * DENDRO_426 + DENDRO_436;
const double DENDRO_520 = DENDRO_200 * DENDRO_92;
const double DENDRO_521 = -DENDRO_520;
const double DENDRO_522 = DENDRO_103 * DENDRO_276;
const double DENDRO_523 = -DENDRO_100 * DENDRO_200;
const double DENDRO_524 = DENDRO_110 * DENDRO_260 + DENDRO_449 * DENDRO_90;
const double DENDRO_525 = 0.25 * DENDRO_348;
const double DENDRO_526 = DENDRO_153 * DENDRO_234;
const double DENDRO_527 = DENDRO_122 * DENDRO_312 + DENDRO_526;
const double DENDRO_528 = DENDRO_449 * DENDRO_97;
const double DENDRO_529 = DENDRO_144 * DENDRO_426 + DENDRO_150 * DENDRO_236;
const double DENDRO_530 = -DENDRO_208 * DENDRO_239;
const double DENDRO_531 = DENDRO_196 * DENDRO_260 + DENDRO_205 * DENDRO_246;
const double DENDRO_532 = DENDRO_418 * grad_1_gt3[pp];
const double DENDRO_533 = DENDRO_411 + DENDRO_532;
const double DENDRO_534 = DENDRO_176 * DENDRO_242 + DENDRO_345;
const double DENDRO_535 = DENDRO_110 * DENDRO_276;
const double DENDRO_536 = 0.25 * DENDRO_90;
const double DENDRO_537 = DENDRO_100 * DENDRO_536;
const double DENDRO_538 = DENDRO_131 * grad_0_gt3[pp];
const double DENDRO_539 = DENDRO_302 + DENDRO_538;
const double DENDRO_540 = DENDRO_92 * grad_2_gt3[pp];
const double DENDRO_541 = DENDRO_305 + DENDRO_540;
const double DENDRO_542 = DENDRO_90 * DENDRO_92;
const double DENDRO_543 = DENDRO_380 + DENDRO_542;
const double DENDRO_544 = DENDRO_31 * DENDRO_463;
const double DENDRO_545 = -grad2_1_2_chi[pp];
const double DENDRO_546 = DENDRO_212 * grad_1_gt5[pp] + DENDRO_214 * DENDRO_97 +
                          DENDRO_27 * grad_2_gt3[pp];
const double DENDRO_547 = DENDRO_218 * grad_2_gt3[pp];
const double DENDRO_548 = DENDRO_221 * grad_1_gt5[pp];
const double DENDRO_549 = DENDRO_39 * grad_2_gt3[pp];
const double DENDRO_550 = DENDRO_212 * DENDRO_97;
const double DENDRO_551 = DENDRO_75 * gt4[pp];
const double DENDRO_552 =
    -DENDRO_227 *
        (DENDRO_469 * DENDRO_546 + DENDRO_471 * (DENDRO_121 + DENDRO_547) +
         DENDRO_475 * (DENDRO_548 + DENDRO_549 + DENDRO_550) + DENDRO_545) +
    DENDRO_254 * grad_1_gt4[pp] + DENDRO_270 * grad_0_gt4[pp] +
    DENDRO_284 * grad_2_gt4[pp] + DENDRO_291 * DENDRO_551;
const double DENDRO_553 =
    -DENDRO_126 * (0.5 * DENDRO_315 + DENDRO_508) -
    DENDRO_126 * (-DENDRO_318 + DENDRO_510) -
    DENDRO_126 * (DENDRO_105 * DENDRO_276 + DENDRO_313 + DENDRO_511) -
    DENDRO_135 * (DENDRO_365 + DENDRO_514) -
    DENDRO_135 * (DENDRO_137 * DENDRO_236 + DENDRO_512) +
    DENDRO_158 * (DENDRO_345 + DENDRO_533) +
    DENDRO_158 * (DENDRO_530 + DENDRO_531) +
    DENDRO_158 * (DENDRO_532 + DENDRO_534) +
    DENDRO_158 * (DENDRO_312 * DENDRO_92 + DENDRO_535 + DENDRO_537) +
    DENDRO_165 * (DENDRO_101 + DENDRO_115 * grad_2_gt3[pp]) -
    DENDRO_179 * (DENDRO_174 + DENDRO_438) -
    DENDRO_179 * (DENDRO_174 + DENDRO_519) -
    DENDRO_179 * (DENDRO_199 + DENDRO_515) -
    DENDRO_179 * (DENDRO_431 + DENDRO_517) -
    DENDRO_179 * (DENDRO_433 + DENDRO_518) -
    DENDRO_179 * (DENDRO_445 + DENDRO_521) +
    DENDRO_181 * (DENDRO_525 + DENDRO_527) +
    DENDRO_181 * (DENDRO_528 + DENDRO_529) +
    DENDRO_181 * (-DENDRO_150 * DENDRO_239 + DENDRO_524) +
    DENDRO_181 * (DENDRO_105 * DENDRO_242 + DENDRO_525 + DENDRO_526) +
    DENDRO_181 * (DENDRO_115 * DENDRO_190 + DENDRO_522 + DENDRO_523) -
    DENDRO_335 * (1.0 * DENDRO_162 + DENDRO_505) -
    DENDRO_335 * (0.5 * DENDRO_188 + DENDRO_507) -
    DENDRO_335 * (DENDRO_105 * DENDRO_115 + DENDRO_180 + DENDRO_506) +
    DENDRO_458 * (DENDRO_409 + DENDRO_539) +
    DENDRO_458 * (DENDRO_462 + DENDRO_541) + DENDRO_489 + DENDRO_490 +
    DENDRO_492 + DENDRO_493 + DENDRO_495 + DENDRO_496 + DENDRO_498 +
    DENDRO_499 + DENDRO_500 + DENDRO_501 + DENDRO_502 + DENDRO_503 +
    DENDRO_504 - DENDRO_544 * (DENDRO_460 + DENDRO_543) + DENDRO_552;
const double DENDRO_554 = DENDRO_39 * DENDRO_75;
const double DENDRO_555 = DENDRO_486 * DENDRO_554;
const double DENDRO_556 = 1.0 * DENDRO_332;
const double DENDRO_557 = 0.5 * DENDRO_331;
const double DENDRO_558 = 0.25 * DENDRO_542;
const double DENDRO_559 = -DENDRO_113 * DENDRO_320;
const double DENDRO_560 = DENDRO_157 + DENDRO_444;
const double DENDRO_561 = DENDRO_113 * DENDRO_260;
const double DENDRO_562 = 0.25 * DENDRO_101 + DENDRO_522;
const double DENDRO_563 = DENDRO_195 * DENDRO_260;
const double DENDRO_564 = 0.25 * DENDRO_538 + DENDRO_563;
const double DENDRO_565 = DENDRO_242 * DENDRO_304;
const double DENDRO_566 = 0.25 * DENDRO_540;
const double DENDRO_567 = DENDRO_306 + DENDRO_535;
const double DENDRO_568 = DENDRO_182 + DENDRO_516;
const double DENDRO_569 = DENDRO_228 * DENDRO_25;
const double DENDRO_570 = DENDRO_200 * DENDRO_233;
const double DENDRO_571 = DENDRO_233 * DENDRO_90;
const double DENDRO_572 = DENDRO_373 + DENDRO_571;
const double DENDRO_573 = DENDRO_173 + DENDRO_192;
const double DENDRO_574 =
    -DENDRO_179 * (DENDRO_519 - DENDRO_570) -
    DENDRO_335 * (DENDRO_153 * DENDRO_189 + DENDRO_203 + DENDRO_505) +
    DENDRO_489 + DENDRO_490 + DENDRO_492 + DENDRO_493 + DENDRO_495 +
    DENDRO_496 + DENDRO_498 + DENDRO_499 + DENDRO_500 + DENDRO_501 +
    DENDRO_502 + DENDRO_503 + DENDRO_504 -
    DENDRO_544 * (DENDRO_457 + DENDRO_572) -
    DENDRO_569 * (DENDRO_437 + DENDRO_573);
const double DENDRO_575 =
    -DENDRO_126 * (DENDRO_510 + DENDRO_557) -
    DENDRO_126 * (DENDRO_511 + DENDRO_556) -
    DENDRO_126 * (DENDRO_189 * DENDRO_242 + DENDRO_316 + DENDRO_508) -
    DENDRO_135 * (DENDRO_366 + DENDRO_514) -
    DENDRO_135 * (DENDRO_110 * DENDRO_274 + DENDRO_558) +
    DENDRO_158 * (DENDRO_311 + DENDRO_531) +
    DENDRO_158 * (DENDRO_346 + DENDRO_534) +
    DENDRO_158 * (DENDRO_410 + DENDRO_564) +
    DENDRO_158 * (DENDRO_533 + DENDRO_565) +
    DENDRO_158 * (DENDRO_537 + DENDRO_567) +
    DENDRO_158 * (DENDRO_566 + DENDRO_567) +
    DENDRO_165 * (DENDRO_242 * grad_1_gt5[pp] + DENDRO_348) -
    DENDRO_179 * (DENDRO_515 + DENDRO_559) -
    DENDRO_179 * (DENDRO_518 + DENDRO_560) -
    DENDRO_179 * (DENDRO_521 + DENDRO_560) +
    DENDRO_181 * (DENDRO_523 + DENDRO_562) +
    DENDRO_181 * (DENDRO_524 + DENDRO_528) +
    DENDRO_181 * (DENDRO_529 - DENDRO_561) +
    DENDRO_181 * (-DENDRO_108 * DENDRO_242 + DENDRO_527) +
    DENDRO_181 * (DENDRO_115 * DENDRO_189 + DENDRO_562) -
    DENDRO_335 * (0.5 * DENDRO_171 + DENDRO_506) -
    DENDRO_335 * (-DENDRO_206 + DENDRO_507) + DENDRO_552 -
    DENDRO_569 * (DENDRO_429 + DENDRO_568) + DENDRO_574;
const double DENDRO_576 = DENDRO_335 * (-DENDRO_201 + DENDRO_408);
const double DENDRO_577 = 0.25 * DENDRO_461;
const double DENDRO_578 = DENDRO_126 * (DENDRO_346 + DENDRO_414);
const double DENDRO_579 = DENDRO_135 * (0.5 * DENDRO_382 + DENDRO_416);
const double DENDRO_580 = DENDRO_179 * (-DENDRO_108 * DENDRO_249 + DENDRO_427);
const double DENDRO_581 = DENDRO_156 * DENDRO_92 + DENDRO_420;
const double DENDRO_582 = DENDRO_179 * (DENDRO_419 + DENDRO_425);
const double DENDRO_583 =
    -0.5 * DENDRO_137 * DENDRO_190 + DENDRO_428 + DENDRO_570;
const double DENDRO_584 = DENDRO_160 + DENDRO_444;
const double DENDRO_585 = DENDRO_447 + 0.25 * DENDRO_571;
const double DENDRO_586 = DENDRO_265 * DENDRO_375 + DENDRO_450;
const double DENDRO_587 = DENDRO_453 + 0.25 * DENDRO_459;
const double DENDRO_588 = DENDRO_189 * DENDRO_249;
const double DENDRO_589 = 0.25 * DENDRO_456 + DENDRO_588;
const double DENDRO_590 = DENDRO_512 + DENDRO_589;
const double DENDRO_591 = DENDRO_365 + DENDRO_366;
const double DENDRO_592 = 1.0 * DENDRO_25 * DENDRO_39;
const double DENDRO_593 = DENDRO_173 + DENDRO_434 + DENDRO_437;
const double DENDRO_594 =
    -DENDRO_126 * (DENDRO_105 * DENDRO_274 + DENDRO_577) -
    DENDRO_135 * (1.0 * DENDRO_377 + DENDRO_417) -
    DENDRO_135 * (DENDRO_195 * DENDRO_265 + DENDRO_379 + DENDRO_415) +
    DENDRO_158 * DENDRO_585 + DENDRO_158 * DENDRO_590 +
    DENDRO_158 * (DENDRO_381 + DENDRO_455) +
    DENDRO_158 * (DENDRO_381 + DENDRO_587) +
    DENDRO_158 * (DENDRO_448 + DENDRO_586) +
    DENDRO_158 * (DENDRO_452 + DENDRO_591) -
    DENDRO_168 * (DENDRO_265 * grad_0_gt5[pp] + DENDRO_357) -
    DENDRO_179 * (DENDRO_421 + DENDRO_581) -
    DENDRO_179 * (-DENDRO_113 * DENDRO_265 + DENDRO_424) -
    DENDRO_179 * (DENDRO_115 * DENDRO_195 + DENDRO_581) -
    DENDRO_181 * DENDRO_583 + DENDRO_181 * (DENDRO_432 + DENDRO_584) +
    DENDRO_181 * (DENDRO_439 + DENDRO_559) +
    DENDRO_181 * (DENDRO_442 + DENDRO_584) -
    DENDRO_335 * (0.5 * DENDRO_169 + DENDRO_407) -
    DENDRO_335 * (DENDRO_150 * DENDRO_195 + DENDRO_209 + DENDRO_406) -
    DENDRO_464 * (DENDRO_302 + DENDRO_310 + DENDRO_409) + DENDRO_483 -
    DENDRO_576 - DENDRO_578 - DENDRO_579 - DENDRO_580 - DENDRO_582 +
    DENDRO_592 * DENDRO_593 +
    DENDRO_592 * (DENDRO_182 + DENDRO_198 + DENDRO_429);
const double DENDRO_595 = DENDRO_441 + DENDRO_520;
const double DENDRO_596 = 0.25 * DENDRO_327;
const double DENDRO_597 = -DENDRO_239 * DENDRO_242;
const double DENDRO_598 = DENDRO_124 * DENDRO_274 + DENDRO_276 * DENDRO_375;
const double DENDRO_599 = DENDRO_246 * DENDRO_265;
const double DENDRO_600 = DENDRO_233 * DENDRO_339;
const double DENDRO_601 = DENDRO_133 * DENDRO_274 + DENDRO_279 * DENDRO_304;
const double DENDRO_602 = DENDRO_108 * DENDRO_279 + DENDRO_113 * DENDRO_394;
const double DENDRO_603 = DENDRO_366 + DENDRO_450;
const double DENDRO_604 = DENDRO_265 * DENDRO_304;
const double DENDRO_605 = DENDRO_233 * DENDRO_454 + DENDRO_588;
const double DENDRO_606 = DENDRO_113 * DENDRO_276 + DENDRO_200 * DENDRO_274;
const double DENDRO_607 = DENDRO_242 * DENDRO_375 + DENDRO_346;
const double DENDRO_608 = DENDRO_258 * DENDRO_536 + DENDRO_563;
const double DENDRO_609 = DENDRO_413 + DENDRO_532;
const double DENDRO_610 = DENDRO_234 * DENDRO_249;
const double DENDRO_611 = -DENDRO_233 * DENDRO_412 + DENDRO_610;
const double DENDRO_612 = DENDRO_189 * DENDRO_279 + DENDRO_443 * DENDRO_97;
const double DENDRO_613 = DENDRO_247 * DENDRO_260;
const double DENDRO_614 = 0.25 * DENDRO_356 + DENDRO_613;
const double DENDRO_615 = DENDRO_246 * DENDRO_320;
const double DENDRO_616 = DENDRO_175 * DENDRO_443;
const double DENDRO_617 = DENDRO_195 * DENDRO_276 + DENDRO_443 * DENDRO_90;
const double DENDRO_618 = DENDRO_463 * DENDRO_47;
const double DENDRO_619 = -grad2_0_1_chi[pp];
const double DENDRO_620 = DENDRO_27 * grad_0_gt3[pp];
const double DENDRO_621 = DENDRO_214 * grad_1_gt0[pp];
const double DENDRO_622 = DENDRO_175 * DENDRO_212;
const double DENDRO_623 = DENDRO_218 * grad_0_gt3[pp];
const double DENDRO_624 = DENDRO_39 * grad_0_gt3[pp];
const double DENDRO_625 =
    DENDRO_175 * DENDRO_221 + DENDRO_212 * grad_1_gt0[pp] + DENDRO_624;
const double DENDRO_626 = DENDRO_75 * gt1[pp];
const double DENDRO_627 =
    DENDRO_477 * grad_0_Gt2[pp] + DENDRO_478 * grad_0_Gt0[pp] +
    DENDRO_478 * grad_1_Gt1[pp] + DENDRO_488 * gt0[pp] + DENDRO_491 * gt2[pp] +
    DENDRO_494 * grad_0_Gt1[pp] - DENDRO_497 * grad_0_chi[pp] -
    DENDRO_77 * grad2_0_1_gt1[pp] + DENDRO_79 * grad2_0_2_gt1[pp] -
    DENDRO_80 * grad2_1_2_gt1[pp] + DENDRO_83 * grad2_0_0_gt1[pp] +
    DENDRO_85 * grad2_1_1_gt1[pp] + DENDRO_87 * grad2_2_2_gt1[pp];
const double DENDRO_628 =
    -DENDRO_227 * (DENDRO_469 * (DENDRO_620 + DENDRO_621 + DENDRO_622) +
                   DENDRO_471 * (DENDRO_232 + DENDRO_623) +
                   DENDRO_475 * DENDRO_625 + DENDRO_619) +
    DENDRO_254 * grad_1_gt1[pp] + DENDRO_270 * grad_0_gt1[pp] +
    DENDRO_284 * grad_2_gt1[pp] + DENDRO_291 * DENDRO_626 + DENDRO_627;
const double DENDRO_629 =
    -DENDRO_126 * (1.0 * DENDRO_328 + DENDRO_596) -
    DENDRO_126 * (0.5 * DENDRO_333 + DENDRO_598) -
    DENDRO_126 * (DENDRO_236 * DENDRO_242 + DENDRO_324 + DENDRO_597) -
    DENDRO_135 * (0.5 * DENDRO_369 + DENDRO_599) -
    DENDRO_135 * (DENDRO_395 + DENDRO_601) -
    DENDRO_135 * (DENDRO_236 * DENDRO_249 + DENDRO_392 + DENDRO_600) +
    DENDRO_158 * (DENDRO_614 + DENDRO_615) +
    DENDRO_158 * (DENDRO_616 + DENDRO_617) +
    DENDRO_158 * (DENDRO_196 * DENDRO_276 + DENDRO_612) +
    DENDRO_158 * (DENDRO_236 * DENDRO_265 + DENDRO_614) +
    DENDRO_158 * (DENDRO_242 * DENDRO_246 + DENDRO_611) -
    DENDRO_179 * (DENDRO_446 + DENDRO_605) -
    DENDRO_179 * (DENDRO_513 + DENDRO_603) -
    DENDRO_179 * (DENDRO_603 + DENDRO_604) -
    DENDRO_179 * (0.5 * DENDRO_100 * DENDRO_196 - DENDRO_602) +
    DENDRO_181 * (DENDRO_303 + DENDRO_564) +
    DENDRO_181 * (DENDRO_303 + DENDRO_608) +
    DENDRO_181 * (-DENDRO_308 - DENDRO_606) +
    DENDRO_181 * (DENDRO_345 + DENDRO_607) +
    DENDRO_181 * (DENDRO_565 + DENDRO_609) +
    DENDRO_181 * (DENDRO_535 + DENDRO_566 + DENDRO_577) +
    DENDRO_329 * (DENDRO_242 * grad_1_gt0[pp] + DENDRO_347) -
    DENDRO_335 * (0.25 * DENDRO_100 * grad_0_gt5[pp] - DENDRO_595) -
    DENDRO_335 * (DENDRO_144 * DENDRO_195 + DENDRO_517) -
    DENDRO_569 * (DENDRO_456 + DENDRO_572) -
    DENDRO_569 * (DENDRO_459 + DENDRO_543) -
    DENDRO_618 * (DENDRO_434 + DENDRO_573) + DENDRO_628;
const double DENDRO_630 = DENDRO_27 * DENDRO_75;
const double DENDRO_631 = DENDRO_486 * DENDRO_630;
const double DENDRO_632 = DENDRO_595 - 0.25 * DENDRO_92 * grad_1_gt5[pp];
const double DENDRO_633 = 0.5 * DENDRO_323;
const double DENDRO_634 = DENDRO_245 * DENDRO_260;
const double DENDRO_635 = DENDRO_386 + DENDRO_602;
const double DENDRO_636 = -0.5 * DENDRO_190 * DENDRO_92 + DENDRO_606;
const double DENDRO_637 = 0.25 * DENDRO_347;
const double DENDRO_638 = DENDRO_610 + DENDRO_637;
const double DENDRO_639 = -DENDRO_135 * (1.0 * DENDRO_367 + DENDRO_600) +
                          DENDRO_158 * (DENDRO_611 + DENDRO_637) -
                          DENDRO_179 * (DENDRO_374 + DENDRO_589) -
                          DENDRO_179 * (DENDRO_374 + DENDRO_605) +
                          DENDRO_181 * (DENDRO_346 + DENDRO_609) -
                          DENDRO_335 * (DENDRO_137 * DENDRO_189 + DENDRO_435);
const double DENDRO_640 =
    -DENDRO_126 * (DENDRO_326 + DENDRO_598) -
    DENDRO_126 * (DENDRO_597 + DENDRO_633) -
    DENDRO_126 * (DENDRO_321 + DENDRO_596 + DENDRO_634) -
    DENDRO_135 * (0.5 * DENDRO_384 + DENDRO_601) -
    DENDRO_135 * (DENDRO_245 * DENDRO_265 + DENDRO_378 + DENDRO_599) +
    DENDRO_158 * (DENDRO_612 + DENDRO_616) +
    DENDRO_158 * (DENDRO_190 * DENDRO_279 + DENDRO_617) +
    DENDRO_158 * (DENDRO_242 * DENDRO_245 + DENDRO_638) +
    DENDRO_158 * (-DENDRO_239 * DENDRO_265 + DENDRO_613 + DENDRO_615) +
    DENDRO_179 * DENDRO_635 - DENDRO_179 * (DENDRO_513 + DENDRO_586) -
    DENDRO_179 * (DENDRO_558 + DENDRO_587) -
    DENDRO_179 * (DENDRO_591 + DENDRO_604) - DENDRO_181 * DENDRO_636 +
    DENDRO_181 * (DENDRO_530 + DENDRO_608) +
    DENDRO_181 * (DENDRO_532 + DENDRO_607) +
    DENDRO_329 * (DENDRO_265 * grad_0_gt3[pp] + DENDRO_356) +
    DENDRO_335 * DENDRO_632 + DENDRO_592 * (DENDRO_310 + DENDRO_539) +
    DENDRO_592 * (DENDRO_461 + DENDRO_541) -
    DENDRO_618 * (DENDRO_198 + DENDRO_568) + DENDRO_628 + DENDRO_639;
const double DENDRO_641 = (1.0 / 12.0) * chi[pp];
const double DENDRO_642 = 2 * At1[pp];
const double DENDRO_643 = 2 * At2[pp];
const double DENDRO_644 = DENDRO_642 * DENDRO_75;
const double DENDRO_645 = 2 * DENDRO_75;
const double DENDRO_646 = At0[pp] * DENDRO_645;
const double DENDRO_647 = DENDRO_643 * DENDRO_75;
const double DENDRO_648 =
    DENDRO_27 * grad_0_chi[pp] + DENDRO_39 * grad_2_chi[pp];
const double DENDRO_649 = -DENDRO_41 * grad_1_chi[pp] + DENDRO_648;
const double DENDRO_650 = DENDRO_226 * DENDRO_649;
const double DENDRO_651 = 0.5 * gt0[pp];
const double DENDRO_652 = DENDRO_249 + DENDRO_650 * DENDRO_651;
const double DENDRO_653 = 12 * grad_1_alpha[pp];
const double DENDRO_654 = DENDRO_653 * DENDRO_75;
const double DENDRO_655 = DENDRO_29 * grad_0_chi[pp] -
                          DENDRO_39 * grad_1_chi[pp] +
                          DENDRO_47 * grad_2_chi[pp];
const double DENDRO_656 = -DENDRO_655;
const double DENDRO_657 = -0.5 * DENDRO_226 * DENDRO_656 * gt0[pp] + DENDRO_278;
const double DENDRO_658 = 12 * grad_2_alpha[pp];
const double DENDRO_659 = DENDRO_658 * DENDRO_75;
const double DENDRO_660 = DENDRO_400 * DENDRO_75;
const double DENDRO_661 = 1.0 * grad_0_chi[pp];
const double DENDRO_662 = DENDRO_75 * gt0[pp];
const double DENDRO_663 = -DENDRO_27 * grad_1_chi[pp] +
                          DENDRO_29 * grad_2_chi[pp] +
                          DENDRO_31 * grad_0_chi[pp];
const double DENDRO_664 = -DENDRO_663;
const double DENDRO_665 = 0.5 * DENDRO_664;
const double DENDRO_666 = DENDRO_226 * (DENDRO_661 - DENDRO_662 * DENDRO_665) +
                          DENDRO_262 * DENDRO_75 + DENDRO_263 * DENDRO_75 -
                          DENDRO_660;
const double DENDRO_667 = 12 * grad_0_alpha[pp];
const double DENDRO_668 = DENDRO_130 * DENDRO_339;
const double DENDRO_669 = DENDRO_133 * DENDRO_257;
const double DENDRO_670 = DENDRO_257 * grad_0_gt0[pp];
const double DENDRO_671 = DENDRO_264 * grad_1_gt0[pp];
const double DENDRO_672 = DENDRO_130 * grad_0_gt0[pp];
const double DENDRO_673 = DENDRO_264 * grad_2_gt0[pp];
const double DENDRO_674 = DENDRO_278 * grad_0_gt5[pp];
const double DENDRO_675 = DENDRO_156 * DENDRO_273;
const double DENDRO_676 = 0.25 * DENDRO_670;
const double DENDRO_677 = 0.25 * DENDRO_672;
const double DENDRO_678 = DENDRO_175 * DENDRO_278;
const double DENDRO_679 = 0.5 * DENDRO_196;
const double DENDRO_680 = DENDRO_25 * DENDRO_251;
const double DENDRO_681 = 2.0 * DENDRO_680;
const double DENDRO_682 = DENDRO_25 * DENDRO_267;
const double DENDRO_683 = 2.0 * DENDRO_682;
const double DENDRO_684 = DENDRO_25 * DENDRO_281;
const double DENDRO_685 = 2.0 * DENDRO_684;
const double DENDRO_686 = DENDRO_226 * DENDRO_290;
const double DENDRO_687 = 3 * alpha[pp];
const double DENDRO_688 = 0.5 * gt5[pp];
const double DENDRO_689 = DENDRO_153 + DENDRO_650 * DENDRO_688;
const double DENDRO_690 = 4 * grad_1_alpha[pp];
const double DENDRO_691 = DENDRO_690 * DENDRO_75;
const double DENDRO_692 = DENDRO_149 - 0.5 * DENDRO_226 * DENDRO_664 * gt5[pp];
const double DENDRO_693 = 4 * grad_0_alpha[pp];
const double DENDRO_694 = DENDRO_693 * DENDRO_75;
const double DENDRO_695 = 1.0 * grad_2_chi[pp];
const double DENDRO_696 = DENDRO_688 * DENDRO_75;
const double DENDRO_697 = DENDRO_104 * DENDRO_75 + DENDRO_109 * DENDRO_75 -
                          DENDRO_113 * DENDRO_29 * DENDRO_75 +
                          DENDRO_226 * (-DENDRO_656 * DENDRO_696 + DENDRO_695);
const double DENDRO_698 = 4 * grad_2_alpha[pp];
const double DENDRO_699 = 0.5 * gt3[pp];
const double DENDRO_700 = DENDRO_698 * DENDRO_75;
const double DENDRO_701 = 1.0 * grad_1_chi[pp];
const double DENDRO_702 = DENDRO_699 * DENDRO_75;
const double DENDRO_703 = -DENDRO_190 * DENDRO_39 * DENDRO_75 +
                          DENDRO_226 * (-DENDRO_649 * DENDRO_702 + DENDRO_701) +
                          DENDRO_235 * DENDRO_75 + DENDRO_240 * DENDRO_75;
const double DENDRO_704 =
    DENDRO_128 * DENDRO_75 + DENDRO_129 * DENDRO_75 +
    DENDRO_226 * (-DENDRO_482 * DENDRO_664 + grad_2_chi[pp]) -
    DENDRO_27 * DENDRO_75 * DENDRO_90;
const double DENDRO_705 = 2.0 * grad_0_alpha[pp];
const double DENDRO_706 =
    DENDRO_226 * (-DENDRO_482 * DENDRO_656 + grad_0_chi[pp]) -
    DENDRO_39 * DENDRO_75 * DENDRO_90 + DENDRO_75 * DENDRO_88 +
    DENDRO_75 * DENDRO_89;
const double DENDRO_707 = 2.0 * grad_2_alpha[pp];
const double DENDRO_708 = 2.0 * grad_1_alpha[pp];
const double DENDRO_709 = DENDRO_226 * gt2[pp];
const double DENDRO_710 = DENDRO_75 * (DENDRO_137 + DENDRO_649 * DENDRO_709);
const double DENDRO_711 = -DENDRO_704 * DENDRO_705 - DENDRO_706 * DENDRO_707 +
                          DENDRO_708 * DENDRO_710 - 4 * grad2_0_2_alpha[pp];
const double DENDRO_712 = DENDRO_119 * DENDRO_75 + DENDRO_120 * DENDRO_75;
const double DENDRO_713 =
    -DENDRO_118 * DENDRO_75 -
    DENDRO_226 * (-DENDRO_551 * DENDRO_649 + grad_2_chi[pp]) + DENDRO_712;
const double DENDRO_714 =
    DENDRO_226 * (-DENDRO_551 * DENDRO_656 + grad_1_chi[pp]) -
    DENDRO_39 * DENDRO_75 * grad_2_gt3[pp] + DENDRO_75 * DENDRO_96 +
    DENDRO_75 * DENDRO_98;
const double DENDRO_715 = DENDRO_143 - DENDRO_226 * DENDRO_664 * gt4[pp];
const double DENDRO_716 = -DENDRO_705 * DENDRO_715 * DENDRO_75 -
                          DENDRO_707 * DENDRO_714 + DENDRO_708 * DENDRO_713 -
                          4 * grad2_1_2_alpha[pp];
const double DENDRO_717 =
    DENDRO_226 * (-DENDRO_626 * DENDRO_664 + grad_1_chi[pp]) +
    DENDRO_255 * DENDRO_75 + DENDRO_256 * DENDRO_75 -
    DENDRO_27 * DENDRO_75 * grad_0_gt3[pp];
const double DENDRO_718 =
    -DENDRO_175 * DENDRO_39 * DENDRO_75 +
    DENDRO_226 * (-DENDRO_626 * DENDRO_649 + grad_0_chi[pp]) +
    DENDRO_229 * DENDRO_75 - DENDRO_27 * DENDRO_75 * grad_1_gt0[pp];
const double DENDRO_719 = DENDRO_226 * gt1[pp];
const double DENDRO_720 =
    -DENDRO_705 * DENDRO_717 +
    DENDRO_707 * DENDRO_75 *
        (-DENDRO_271 - DENDRO_272 + DENDRO_624 + DENDRO_656 * DENDRO_719) -
    DENDRO_708 * DENDRO_718 - 4 * grad2_0_1_alpha[pp];
const double DENDRO_721 =
    -DENDRO_27 * (DENDRO_629 * alpha[pp] + DENDRO_720) -
    DENDRO_27 * (DENDRO_640 * alpha[pp] + DENDRO_720) +
    DENDRO_29 * (DENDRO_484 * alpha[pp] + DENDRO_711) +
    DENDRO_29 * (DENDRO_594 * alpha[pp] + DENDRO_711) +
    DENDRO_31 * (DENDRO_404 * alpha[pp] + DENDRO_652 * DENDRO_691 -
                 DENDRO_657 * DENDRO_700 - DENDRO_666 * DENDRO_693 -
                 4 * grad2_0_0_alpha[pp]) -
    DENDRO_39 * (DENDRO_553 * alpha[pp] + DENDRO_716) -
    DENDRO_39 * (DENDRO_575 * alpha[pp] + DENDRO_716) +
    DENDRO_41 *
        (DENDRO_351 * alpha[pp] - DENDRO_690 * DENDRO_703 +
         DENDRO_694 * (DENDRO_226 * DENDRO_665 * gt3[pp] + DENDRO_260) +
         DENDRO_700 * (DENDRO_226 * DENDRO_656 * DENDRO_699 + DENDRO_276) -
         4 * grad2_1_1_alpha[pp]) +
    DENDRO_47 * (DENDRO_293 * alpha[pp] + DENDRO_689 * DENDRO_691 -
                 DENDRO_692 * DENDRO_694 - DENDRO_697 * DENDRO_698 -
                 4 * grad2_2_2_alpha[pp]);
const double DENDRO_722 = DENDRO_721 * DENDRO_75;
const double DENDRO_723 = (1.0 / 3.0) * At1[pp];
const double DENDRO_724 = DENDRO_34 + DENDRO_44 - DENDRO_65;
const double DENDRO_725 =
    -At1[pp] * DENDRO_31 + At3[pp] * DENDRO_27 - DENDRO_32;
const double DENDRO_726 =
    -At1[pp] * DENDRO_29 + At3[pp] * DENDRO_39 - DENDRO_58;
const double DENDRO_727 = 6.0 * grad_2_alpha[pp];
const double DENDRO_728 = 6.0 * grad_0_alpha[pp];
const double DENDRO_729 = 6.0 * grad_1_alpha[pp];
const double DENDRO_730 = DENDRO_257 * DENDRO_412;
const double DENDRO_731 = DENDRO_273 * DENDRO_312;
const double DENDRO_732 = DENDRO_130 * DENDRO_388;
const double DENDRO_733 =
    0.25 * DENDRO_143 * grad_0_gt0[pp] + DENDRO_264 * DENDRO_375;
const double DENDRO_734 = DENDRO_668 + DENDRO_669;
const double DENDRO_735 = DENDRO_536 * DENDRO_91;
const double DENDRO_736 = DENDRO_105 * DENDRO_278;
const double DENDRO_737 = DENDRO_133 * DENDRO_99 + DENDRO_736;
const double DENDRO_738 = DENDRO_130 * DENDRO_90;
const double DENDRO_739 = DENDRO_143 * grad_2_gt0[pp];
const double DENDRO_740 = DENDRO_130 * DENDRO_97;
const double DENDRO_741 = DENDRO_739 + DENDRO_740;
const double DENDRO_742 = DENDRO_257 * DENDRO_97;
const double DENDRO_743 = DENDRO_143 * grad_1_gt0[pp] + DENDRO_742;
const double DENDRO_744 = DENDRO_97 * DENDRO_99;
const double DENDRO_745 = (1.0 / 3.0) * At2[pp];
const double DENDRO_746 = At2[pp] * DENDRO_27 - At4[pp] * DENDRO_41 + DENDRO_66;
const double DENDRO_747 =
    -At2[pp] * DENDRO_31 + At4[pp] * DENDRO_27 - DENDRO_56;
const double DENDRO_748 = At4[pp] * DENDRO_39 - At5[pp] * DENDRO_47 - DENDRO_36;
const double DENDRO_749 = DENDRO_91 * grad_2_gt5[pp];
const double DENDRO_750 = DENDRO_257 * grad_0_gt5[pp];
const double DENDRO_751 = DENDRO_159 * DENDRO_91;
const double DENDRO_752 = 0.25 * DENDRO_273 * grad_2_gt5[pp];
const double DENDRO_753 = DENDRO_751 + DENDRO_752;
const double DENDRO_754 = DENDRO_175 * DENDRO_257;
const double DENDRO_755 = DENDRO_103 * DENDRO_278 + DENDRO_156 * DENDRO_91;
const double DENDRO_756 = -0.5 * DENDRO_113 * DENDRO_130;
const double DENDRO_757 = -0.5 * DENDRO_113 * DENDRO_257;
const double DENDRO_758 = DENDRO_257 * DENDRO_679;
const double DENDRO_759 = 2 * At4[pp];
const double DENDRO_760 = At3[pp] * DENDRO_645;
const double DENDRO_761 = DENDRO_75 * DENDRO_759;
const double DENDRO_762 = DENDRO_226 * DENDRO_699;
const double DENDRO_763 = DENDRO_667 * DENDRO_75;
const double DENDRO_764 = DENDRO_241 * grad_2_gt3[pp];
const double DENDRO_765 = DENDRO_159 * DENDRO_273;
const double DENDRO_766 = 1.0 * DENDRO_143;
const double DENDRO_767 = 0.25 * DENDRO_742;
const double DENDRO_768 = DENDRO_241 * grad_0_gt3[pp];
const double DENDRO_769 = DENDRO_686 * DENDRO_75;
const double DENDRO_770 = (1.0 / 3.0) * At4[pp];
const double DENDRO_771 = DENDRO_99 * grad_2_gt5[pp];
const double DENDRO_772 = DENDRO_156 * DENDRO_99;
const double DENDRO_773 = DENDRO_752 + DENDRO_772;
const double DENDRO_774 = DENDRO_159 * DENDRO_99;
const double DENDRO_775 = -0.5 * DENDRO_276 * grad_0_gt5[pp] + DENDRO_765;
const double DENDRO_776 = DENDRO_114 * grad_0_gt5[pp];
const double DENDRO_777 = DENDRO_114 * grad_1_gt5[pp];
const double DENDRO_778 = DENDRO_212 * grad_2_chi[pp] +
                          DENDRO_214 * grad_0_chi[pp] +
                          DENDRO_27 * grad_1_chi[pp];
const double DENDRO_779 = DENDRO_226 * DENDRO_688;
const double DENDRO_780 = DENDRO_75 * grad_0_alpha[pp];
const double DENDRO_781 = DENDRO_218 * grad_1_chi[pp] + DENDRO_648;
const double DENDRO_782 = DENDRO_75 * grad_1_alpha[pp];
const double DENDRO_783 = DENDRO_212 * grad_0_chi[pp] +
                          DENDRO_221 * grad_2_chi[pp] +
                          DENDRO_39 * grad_1_chi[pp];
const double DENDRO_784 = DENDRO_75 * grad_2_alpha[pp];
const double DENDRO_785 = DENDRO_226 * DENDRO_651;
const double DENDRO_786 = 0.5 * grad_1_alpha[pp];
const double DENDRO_787 = 0.5 * grad_2_alpha[pp];
const double DENDRO_788 = 0.5 * grad_0_alpha[pp];
const double DENDRO_789 = DENDRO_49 * DENDRO_75;
const double DENDRO_790 = 3 * DENDRO_25;
const double DENDRO_791 = 6 * DENDRO_25;
const double DENDRO_792 = (1.0 / 3.0) * alpha[pp];
const double DENDRO_793 = (7.0 / 3.0) * DENDRO_485;
const double DENDRO_794 = (1.0 / 3.0) * DENDRO_485;
const double DENDRO_795 = (1.0 / 3.0) * DENDRO_82;
const double DENDRO_796 = 2 * DENDRO_25;
const double DENDRO_797 = DENDRO_796 * grad_0_alpha[pp];
const double DENDRO_798 = pow(DENDRO_24, -3);
const double DENDRO_799 = 4 * grad_0_K[pp];
const double DENDRO_800 = DENDRO_75 * DENDRO_792;
const double DENDRO_801 = DENDRO_796 * grad_2_alpha[pp];
const double DENDRO_802 = DENDRO_796 * grad_1_alpha[pp];
const double DENDRO_803 = 4 * grad_2_K[pp];
const double DENDRO_804 = 4 * grad_1_K[pp];
const double DENDRO_805 = 9 * DENDRO_226;
const double DENDRO_806 = DENDRO_220 * DENDRO_805;
const double DENDRO_807 =
    -2.0 * DENDRO_131 * DENDRO_59 * DENDRO_798 * alpha[pp] -
    2.0 * DENDRO_144 * DENDRO_68 * DENDRO_798 * alpha[pp] -
    2 * DENDRO_150 * DENDRO_51 * DENDRO_798 * alpha[pp] -
    2.0 / 3.0 * DENDRO_17 * DENDRO_25 * DENDRO_268 +
    DENDRO_253 * grad_1_beta0[pp] -
    2.0 * DENDRO_258 * DENDRO_62 * DENDRO_798 * alpha[pp] -
    2 * DENDRO_260 * DENDRO_45 * DENDRO_798 * alpha[pp] -
    2 * DENDRO_265 * DENDRO_37 * DENDRO_798 * alpha[pp] +
    DENDRO_269 * grad_0_beta0[pp] -
    7.0 / 3.0 * DENDRO_27 * DENDRO_75 * grad2_0_1_beta0[pp] -
    1.0 / 3.0 * DENDRO_27 * DENDRO_75 * grad2_1_1_beta1[pp] -
    1.0 / 3.0 * DENDRO_27 * DENDRO_75 * grad2_1_2_beta2[pp] +
    DENDRO_283 * grad_2_beta0[pp] + DENDRO_37 * DENDRO_797 -
    2 * DENDRO_39 * DENDRO_75 * grad2_1_2_beta0[pp] + DENDRO_59 * DENDRO_801 +
    DENDRO_62 * DENDRO_802 + DENDRO_793 * grad2_0_2_beta0[pp] +
    DENDRO_794 * grad2_1_2_beta1[pp] + DENDRO_794 * grad2_2_2_beta2[pp] +
    DENDRO_795 * grad2_0_1_beta1[pp] + DENDRO_795 * grad2_0_2_beta2[pp] +
    DENDRO_800 * (DENDRO_27 * DENDRO_804 + DENDRO_62 * DENDRO_806) +
    DENDRO_800 * (9 * DENDRO_226 * DENDRO_37 * DENDRO_75 * grad_0_chi[pp] -
                  DENDRO_31 * DENDRO_799) +
    DENDRO_800 * (9 * DENDRO_226 * DENDRO_59 * DENDRO_75 * grad_2_chi[pp] -
                  DENDRO_29 * DENDRO_803) +
    (4.0 / 3.0) * DENDRO_82 * grad2_0_0_beta0[pp] +
    DENDRO_84 * grad2_1_1_beta0[pp] + DENDRO_86 * grad2_2_2_beta0[pp] -
    beta0[pp] * grad_0_Gt0[pp] - beta1[pp] * grad_1_Gt0[pp] -
    beta2[pp] * grad_2_Gt0[pp];
const double DENDRO_808 = DENDRO_82 * grad2_0_0_beta1[pp];
const double DENDRO_809 = DENDRO_86 * grad2_2_2_beta1[pp];
const double DENDRO_810 = DENDRO_789 * grad2_0_2_beta1[pp];
const double DENDRO_811 = DENDRO_45 * DENDRO_802;
const double DENDRO_812 = (4.0 / 3.0) * DENDRO_84 * grad2_1_1_beta1[pp];
const double DENDRO_813 = DENDRO_27 * DENDRO_799;
const double DENDRO_814 = DENDRO_217 * DENDRO_805;
const double DENDRO_815 = DENDRO_39 * DENDRO_803;
const double DENDRO_816 = DENDRO_225 * DENDRO_805;
const double DENDRO_817 = (1.0 / 3.0) * DENDRO_84;
const double DENDRO_818 = DENDRO_817 * grad2_0_1_beta0[pp];
const double DENDRO_819 = DENDRO_817 * grad2_1_2_beta2[pp];
const double DENDRO_820 = DENDRO_41 * DENDRO_804 - DENDRO_45 * DENDRO_806;
const double DENDRO_821 = (1.0 / 3.0) * DENDRO_630;
const double DENDRO_822 = (1.0 / 3.0) * grad2_0_2_beta0[pp];
const double DENDRO_823 = DENDRO_0 * DENDRO_798;
const double DENDRO_824 = 2.0 * DENDRO_798 * alpha[pp];
const double DENDRO_825 = beta0[pp] * grad_0_Gt1[pp] +
                          beta1[pp] * grad_1_Gt1[pp] +
                          beta2[pp] * grad_2_Gt1[pp];
const double DENDRO_826 =
    -2.0 * DENDRO_100 * DENDRO_68 * DENDRO_798 * alpha[pp] -
    2 * DENDRO_115 * DENDRO_51 * DENDRO_798 * alpha[pp] -
    2.0 / 3.0 * DENDRO_17 * DENDRO_25 * DENDRO_282 +
    DENDRO_253 * grad_1_beta2[pp] + DENDRO_269 * grad_0_beta2[pp] -
    2 * DENDRO_27 * DENDRO_75 * grad2_0_1_beta2[pp] -
    2.0 * DENDRO_274 * DENDRO_62 * DENDRO_798 * alpha[pp] -
    2 * DENDRO_276 * DENDRO_45 * DENDRO_798 * alpha[pp] -
    2 * DENDRO_279 * DENDRO_37 * DENDRO_798 * alpha[pp] +
    DENDRO_283 * grad_2_beta2[pp] -
    1.0 / 3.0 * DENDRO_39 * DENDRO_75 * grad2_0_1_beta0[pp] -
    1.0 / 3.0 * DENDRO_39 * DENDRO_75 * grad2_1_1_beta1[pp] -
    7.0 / 3.0 * DENDRO_39 * DENDRO_75 * grad2_1_2_beta2[pp] +
    DENDRO_51 * DENDRO_801 + DENDRO_59 * DENDRO_797 -
    2.0 * DENDRO_59 * DENDRO_798 * DENDRO_92 * alpha[pp] +
    DENDRO_68 * DENDRO_802 + DENDRO_793 * grad2_0_2_beta2[pp] +
    DENDRO_794 * grad2_0_0_beta0[pp] + DENDRO_794 * grad2_0_1_beta1[pp] +
    DENDRO_800 * (DENDRO_39 * DENDRO_804 + DENDRO_68 * DENDRO_806) +
    DENDRO_800 * (9 * DENDRO_226 * DENDRO_51 * DENDRO_75 * grad_2_chi[pp] -
                  DENDRO_47 * DENDRO_803) +
    DENDRO_800 * (9 * DENDRO_226 * DENDRO_59 * DENDRO_75 * grad_0_chi[pp] -
                  DENDRO_29 * DENDRO_799) +
    DENDRO_82 * grad2_0_0_beta2[pp] + DENDRO_822 * DENDRO_86 +
    DENDRO_84 * grad2_1_1_beta2[pp] +
    (1.0 / 3.0) * DENDRO_86 * grad2_1_2_beta1[pp] +
    (4.0 / 3.0) * DENDRO_86 * grad2_2_2_beta2[pp] - beta0[pp] * grad_0_Gt2[pp] -
    beta1[pp] * grad_1_Gt2[pp] - beta2[pp] * grad_2_Gt2[pp];

// Dendro: printing variables
//--
a_rhs[pp] =
    -DENDRO_0 * K[pp] -
    DENDRO_1 * h_ssl * (-DENDRO_1 + alpha[pp]) *
        exp(-1.0 / 2.0 * (t * t) / (sig_ssl * sig_ssl)) +
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
chi_rhs[pp] =
    -BSSN_CAHD_C * DENDRO_641 * (dx_i * dx_i) *
        (-8 * DENDRO_18 + DENDRO_26 * DENDRO_38 + DENDRO_26 * DENDRO_46 +
         DENDRO_26 * DENDRO_52 + 3 * DENDRO_293 * DENDRO_294 +
         3 * DENDRO_351 * DENDRO_352 + 3 * DENDRO_404 * DENDRO_405 +
         DENDRO_484 * DENDRO_487 + DENDRO_487 * DENDRO_594 +
         DENDRO_53 * DENDRO_60 + DENDRO_53 * DENDRO_63 + DENDRO_53 * DENDRO_69 -
         DENDRO_553 * DENDRO_555 - DENDRO_555 * DENDRO_575 -
         DENDRO_629 * DENDRO_631 - DENDRO_631 * DENDRO_640) /
        dt -
    DENDRO_16 * DENDRO_17 + DENDRO_16 * K[pp] * alpha[pp] +
    beta0[pp] * grad_0_chi[pp] + beta1[pp] * grad_1_chi[pp] +
    beta2[pp] * grad_2_chi[pp];
//--
At_rhs00[pp] =
    -At0[pp] * DENDRO_10 + At0[pp] * DENDRO_3 - At0[pp] * DENDRO_8 +
    DENDRO_641 *
        (DENDRO_652 * DENDRO_654 - DENDRO_657 * DENDRO_659 -
         DENDRO_666 * DENDRO_667 +
         DENDRO_687 *
             (3.0 * DENDRO_130 * DENDRO_25 * DENDRO_47 * grad_2_gt0[pp] -
              DENDRO_158 * (1.0 * DENDRO_671 + DENDRO_676) -
              DENDRO_158 * (DENDRO_273 * DENDRO_679 + DENDRO_278 * DENDRO_390) -
              DENDRO_179 *
                  (2 * DENDRO_113 * DENDRO_278 - DENDRO_679 * DENDRO_91) -
              DENDRO_181 * DENDRO_385 - DENDRO_181 * DENDRO_387 -
              DENDRO_181 * (DENDRO_130 * DENDRO_245 + DENDRO_669) -
              DENDRO_181 * (DENDRO_195 * DENDRO_257 + DENDRO_668) -
              DENDRO_181 * (DENDRO_375 * DENDRO_91 + DENDRO_675) -
              DENDRO_227 * (DENDRO_217 * DENDRO_265 + DENDRO_220 * DENDRO_249 +
                            DENDRO_225 * DENDRO_279 + DENDRO_398) +
              3.0 * DENDRO_25 * DENDRO_257 * DENDRO_41 * grad_1_gt0[pp] +
              6.0 * DENDRO_25 * DENDRO_264 * DENDRO_31 * grad_0_gt0[pp] +
              2.0 * DENDRO_25 * DENDRO_27 * DENDRO_368 +
              4 * DENDRO_25 * DENDRO_27 * DENDRO_393 +
              4 * DENDRO_25 * DENDRO_273 * DENDRO_362 * DENDRO_41 +
              4 * DENDRO_25 * DENDRO_278 * DENDRO_31 * DENDRO_364 +
              2.0 * DENDRO_25 * DENDRO_29 * (DENDRO_672 + DENDRO_673) +
              4 * DENDRO_25 * DENDRO_29 * (1.0 * DENDRO_673 + DENDRO_677) +
              2.0 * DENDRO_25 * DENDRO_29 *
                  (DENDRO_674 + DENDRO_91 * grad_2_gt0[pp]) +
              4 * DENDRO_25 * DENDRO_376 * DENDRO_39 +
              4.0 * DENDRO_29 * DENDRO_75 * grad2_0_2_gt0[pp] +
              2.0 * DENDRO_31 * DENDRO_75 * grad2_0_0_gt0[pp] -
              DENDRO_329 * (DENDRO_670 + DENDRO_671) -
              DENDRO_329 * (DENDRO_273 * grad_2_gt0[pp] + DENDRO_678) -
              DENDRO_335 * DENDRO_359 * DENDRO_91 - DENDRO_353 - DENDRO_354 -
              DENDRO_355 - DENDRO_358 - DENDRO_361 - DENDRO_363 - DENDRO_383 -
              DENDRO_391 + 2.0 * DENDRO_41 * DENDRO_75 * grad2_1_1_gt0[pp] +
              2.0 * DENDRO_47 * DENDRO_75 * grad2_2_2_gt0[pp] -
              DENDRO_662 * DENDRO_686 - DENDRO_681 * grad_1_gt0[pp] -
              DENDRO_683 * grad_0_gt0[pp] - DENDRO_685 * grad_2_gt0[pp] +
              4.0 * grad_0_Gt0[pp] * gt0[pp] + 4.0 * grad_0_Gt1[pp] * gt1[pp] +
              4.0 * grad_0_Gt2[pp] * gt2[pp]) +
         DENDRO_722 * gt0[pp] - 12 * grad2_0_0_alpha[pp]) +
    DENDRO_642 * grad_0_beta1[pp] + DENDRO_643 * grad_0_beta2[pp] -
    alpha[pp] * (-At0[pp] * K[pp] +
                 DENDRO_644 * (-At1[pp] * DENDRO_41 + DENDRO_42 + DENDRO_64) +
                 DENDRO_646 * (At1[pp] * DENDRO_27 - DENDRO_36 - DENDRO_55) +
                 DENDRO_647 *
                     (-At0[pp] * DENDRO_29 + At1[pp] * DENDRO_39 - DENDRO_57)) +
    beta0[pp] * grad_0_At0[pp] + beta1[pp] * grad_1_At0[pp] +
    beta2[pp] * grad_2_At0[pp];
//--
At_rhs01[pp] =
    At0[pp] * grad_1_beta0[pp] - At1[pp] * DENDRO_8 +
    At2[pp] * grad_1_beta2[pp] + At3[pp] * grad_0_beta1[pp] +
    At4[pp] * grad_0_beta2[pp] +
    DENDRO_641 *
        (DENDRO_626 * DENDRO_721 +
         DENDRO_687 *
             (-DENDRO_126 * (DENDRO_239 * DENDRO_241 + DENDRO_633) -
              DENDRO_126 * (-DENDRO_124 * DENDRO_273 +
                            0.5 * DENDRO_276 * DENDRO_97 - DENDRO_731) -
              DENDRO_126 *
                  (-DENDRO_257 * DENDRO_344 + DENDRO_634 + DENDRO_730) +
              DENDRO_135 * (DENDRO_133 * DENDRO_273 + DENDRO_278 * DENDRO_304 +
                            0.5 * DENDRO_678) +
              DENDRO_135 * (DENDRO_245 * DENDRO_264 + DENDRO_246 * DENDRO_264 +
                            DENDRO_676) +
              DENDRO_158 * (-DENDRO_241 * DENDRO_245 + DENDRO_638) -
              DENDRO_158 * (DENDRO_189 * DENDRO_278 + DENDRO_273 * DENDRO_426 +
                            DENDRO_273 * DENDRO_454) +
              DENDRO_158 * (-DENDRO_190 * DENDRO_278 - DENDRO_273 * DENDRO_536 +
                            0.5 * DENDRO_276 * grad_2_gt0[pp]) +
              DENDRO_158 * (DENDRO_239 * DENDRO_264 - DENDRO_257 * DENDRO_388 +
                            DENDRO_613) +
              DENDRO_179 * DENDRO_635 + DENDRO_179 * (DENDRO_732 + DENDRO_733) +
              DENDRO_179 * (DENDRO_735 + DENDRO_737) +
              DENDRO_179 * (DENDRO_264 * DENDRO_304 + DENDRO_734) -
              DENDRO_181 * DENDRO_636 +
              DENDRO_181 * (DENDRO_130 * DENDRO_412 - DENDRO_257 * DENDRO_536 +
                            DENDRO_563) +
              DENDRO_181 *
                  (-DENDRO_241 * DENDRO_375 + DENDRO_346 + DENDRO_532) -
              DENDRO_227 * (DENDRO_217 * DENDRO_320 + DENDRO_225 * DENDRO_394 +
                            DENDRO_233 * DENDRO_471 + DENDRO_619) -
              DENDRO_329 *
                  (DENDRO_257 * grad_1_gt0[pp] + DENDRO_264 * grad_0_gt3[pp]) +
              DENDRO_335 * DENDRO_632 -
              DENDRO_592 * (DENDRO_130 * grad_0_gt3[pp] + DENDRO_743) -
              DENDRO_592 * (DENDRO_273 * grad_1_gt5[pp] + DENDRO_744 +
                            DENDRO_91 * grad_2_gt3[pp]) +
              DENDRO_618 * (DENDRO_738 + DENDRO_741) - DENDRO_626 * DENDRO_686 +
              DENDRO_627 + DENDRO_639 - DENDRO_681 * grad_1_gt1[pp] -
              DENDRO_683 * grad_0_gt1[pp] - DENDRO_685 * grad_2_gt1[pp]) -
         DENDRO_717 * DENDRO_728 - DENDRO_718 * DENDRO_729 +
         DENDRO_727 * DENDRO_75 * (-DENDRO_273 - DENDRO_655 * DENDRO_719) -
         12 * grad2_0_1_alpha[pp]) +
    DENDRO_723 * grad_0_beta0[pp] + DENDRO_723 * grad_1_beta1[pp] -
    alpha[pp] * (-At1[pp] * K[pp] + DENDRO_644 * DENDRO_724 +
                 DENDRO_646 * DENDRO_725 + DENDRO_647 * DENDRO_726) +
    beta0[pp] * grad_0_At1[pp] + beta1[pp] * grad_1_At1[pp] +
    beta2[pp] * grad_2_At1[pp];
//--
At_rhs02[pp] =
    At0[pp] * grad_2_beta0[pp] + At1[pp] * grad_2_beta1[pp] -
    At2[pp] * DENDRO_10 + At4[pp] * grad_0_beta1[pp] +
    At5[pp] * grad_0_beta2[pp] +
    DENDRO_641 *
        (DENDRO_482 * DENDRO_721 +
         DENDRO_687 *
             (-DENDRO_158 * (DENDRO_675 + DENDRO_737) -
              DENDRO_158 * (DENDRO_733 + DENDRO_758) -
              DENDRO_158 * (DENDRO_176 * DENDRO_264 + DENDRO_734) -
              DENDRO_158 * (DENDRO_454 * DENDRO_91 + DENDRO_675 + DENDRO_736) -
              DENDRO_179 * (0.5 * DENDRO_113 * DENDRO_91 - DENDRO_755) -
              DENDRO_179 * (DENDRO_113 * DENDRO_264 - DENDRO_130 * DENDRO_679 -
                            DENDRO_149 * DENDRO_247) -
              DENDRO_181 * DENDRO_583 -
              DENDRO_181 * (DENDRO_114 * DENDRO_375 + DENDRO_753) -
              DENDRO_227 * (DENDRO_137 * DENDRO_471 + DENDRO_208 * DENDRO_217 +
                            DENDRO_225 * DENDRO_396 + DENDRO_465) +
              4 * DENDRO_25 * DENDRO_27 * DENDRO_585 +
              4 * DENDRO_25 * DENDRO_27 * DENDRO_590 +
              4 * DENDRO_25 * DENDRO_29 *
                  (DENDRO_114 * DENDRO_195 + DENDRO_755) +
              2.0 * DENDRO_25 * DENDRO_29 *
                  (DENDRO_130 * grad_2_gt0[pp] + DENDRO_264 * grad_0_gt5[pp]) +
              4 * DENDRO_25 * DENDRO_31 *
                  (DENDRO_133 * DENDRO_91 + 1.0 * DENDRO_674) +
              4 * DENDRO_25 * DENDRO_31 *
                  (DENDRO_195 * DENDRO_264 + DENDRO_196 * DENDRO_264 +
                   DENDRO_677) +
              1.0 * DENDRO_25 * DENDRO_39 * DENDRO_593 +
              4 * DENDRO_25 * DENDRO_39 *
                  (0.5 * DENDRO_113 * DENDRO_99 - DENDRO_753) +
              4 * DENDRO_25 * DENDRO_39 *
                  (-DENDRO_130 * DENDRO_426 - DENDRO_149 * DENDRO_245 -
                   DENDRO_757) +
              1.0 * DENDRO_25 * DENDRO_41 * (DENDRO_743 + DENDRO_754) +
              4 * DENDRO_25 * DENDRO_41 *
                  (DENDRO_105 * DENDRO_273 + 0.25 * DENDRO_744) +
              4.0 * DENDRO_29 * DENDRO_75 * grad2_0_2_gt2[pp] +
              2.0 * DENDRO_31 * DENDRO_75 * grad2_0_0_gt2[pp] -
              DENDRO_335 * (DENDRO_113 * DENDRO_114 - 0.5 * DENDRO_749) -
              DENDRO_335 * (-DENDRO_130 * DENDRO_156 - DENDRO_149 * DENDRO_195 -
                            DENDRO_756) +
              2.0 * DENDRO_41 * DENDRO_75 * grad2_1_1_gt2[pp] +
              2.0 * DENDRO_47 * DENDRO_75 * grad2_2_2_gt2[pp] - DENDRO_479 -
              DENDRO_480 - DENDRO_481 - DENDRO_482 * DENDRO_686 - DENDRO_576 -
              DENDRO_578 - DENDRO_579 - DENDRO_580 - DENDRO_582 -
              DENDRO_592 * (DENDRO_741 + DENDRO_750) -
              DENDRO_681 * grad_1_gt2[pp] - DENDRO_683 * grad_0_gt2[pp] -
              DENDRO_685 * grad_2_gt2[pp] + 2.0 * grad_0_Gt0[pp] * gt2[pp] +
              2.0 * grad_0_Gt1[pp] * gt4[pp] + 2.0 * grad_0_Gt2[pp] * gt5[pp] +
              2.0 * grad_2_Gt0[pp] * gt0[pp] + 2.0 * grad_2_Gt1[pp] * gt1[pp] +
              2.0 * grad_2_Gt2[pp] * gt2[pp]) -
         DENDRO_704 * DENDRO_728 - DENDRO_706 * DENDRO_727 +
         DENDRO_710 * DENDRO_729 - 12 * grad2_0_2_alpha[pp]) +
    DENDRO_745 * grad_0_beta0[pp] + DENDRO_745 * grad_2_beta2[pp] -
    alpha[pp] * (-At2[pp] * K[pp] + DENDRO_644 * DENDRO_746 +
                 DENDRO_646 * DENDRO_747 + DENDRO_647 * DENDRO_748) +
    beta0[pp] * grad_0_At2[pp] + beta1[pp] * grad_1_At2[pp] +
    beta2[pp] * grad_2_At2[pp];
//--
At_rhs11[pp] =
    -At3[pp] * DENDRO_11 + At3[pp] * DENDRO_12 - At3[pp] * DENDRO_8 +
    DENDRO_641 *
        (-DENDRO_653 * DENDRO_703 +
         DENDRO_659 * (DENDRO_276 - DENDRO_655 * DENDRO_762) +
         DENDRO_687 *
             (DENDRO_143 * DENDRO_338 + DENDRO_158 * (DENDRO_322 + DENDRO_730) +
              DENDRO_158 *
                  (0.25 * DENDRO_233 * grad_1_gt3[pp] - 1.0 * DENDRO_768) +
              DENDRO_158 * (1.0 * DENDRO_276 * DENDRO_90 - DENDRO_731) +
              DENDRO_165 * (DENDRO_315 - DENDRO_764) +
              DENDRO_165 * (DENDRO_332 - DENDRO_99 * grad_2_gt3[pp]) +
              DENDRO_165 *
                  (-DENDRO_143 * grad_0_gt3[pp] + DENDRO_260 * DENDRO_97) +
              DENDRO_179 * DENDRO_309 +
              DENDRO_179 * (DENDRO_143 * DENDRO_339 + DENDRO_257 * DENDRO_304) +
              DENDRO_179 * (DENDRO_246 * DENDRO_766 + DENDRO_767) +
              DENDRO_179 * (DENDRO_304 * DENDRO_99 + DENDRO_765) -
              DENDRO_181 * DENDRO_319 +
              DENDRO_181 * (DENDRO_316 - 1.0 * DENDRO_764) -
              DENDRO_181 * (DENDRO_312 * DENDRO_99 + DENDRO_314) -
              DENDRO_227 * (DENDRO_217 * DENDRO_260 + DENDRO_220 * DENDRO_242 +
                            DENDRO_225 * DENDRO_276 + DENDRO_295) +
              DENDRO_241 * DENDRO_342 * DENDRO_41 + DENDRO_257 * DENDRO_340 +
              DENDRO_273 * DENDRO_341 +
              DENDRO_329 *
                  (DENDRO_175 * DENDRO_276 - DENDRO_273 * grad_2_gt3[pp]) +
              DENDRO_329 * (DENDRO_233 * grad_1_gt3[pp] - DENDRO_768) +
              DENDRO_329 *
                  (-DENDRO_257 * grad_0_gt3[pp] + DENDRO_260 * grad_1_gt0[pp]) -
              DENDRO_334 * DENDRO_335 * DENDRO_99 + DENDRO_350 -
              DENDRO_681 * grad_1_gt3[pp] - DENDRO_683 * grad_0_gt3[pp] -
              DENDRO_685 * grad_2_gt3[pp] - DENDRO_769 * gt3[pp]) +
         DENDRO_722 * gt3[pp] +
         DENDRO_763 * (DENDRO_260 - DENDRO_663 * DENDRO_762) -
         12 * grad2_1_1_alpha[pp]) +
    DENDRO_642 * grad_1_beta0[pp] + DENDRO_759 * grad_1_beta2[pp] -
    alpha[pp] * (-At3[pp] * K[pp] + DENDRO_644 * DENDRO_725 +
                 DENDRO_724 * DENDRO_760 + DENDRO_726 * DENDRO_761) +
    beta0[pp] * grad_0_At3[pp] + beta1[pp] * grad_1_At3[pp] +
    beta2[pp] * grad_2_At3[pp];
//--
At_rhs12[pp] =
    At1[pp] * grad_2_beta0[pp] + At2[pp] * grad_1_beta0[pp] +
    At3[pp] * grad_2_beta1[pp] - At4[pp] * DENDRO_11 +
    At5[pp] * grad_1_beta2[pp] +
    DENDRO_641 *
        (DENDRO_551 * DENDRO_721 +
         DENDRO_687 *
             (-DENDRO_126 * (-DENDRO_124 * DENDRO_99 + DENDRO_556) -
              DENDRO_126 * (0.25 * DENDRO_122 * grad_1_gt3[pp] -
                            DENDRO_189 * DENDRO_241 - DENDRO_190 * DENDRO_241) -
              DENDRO_126 *
                  (-DENDRO_143 * DENDRO_344 + DENDRO_509 + DENDRO_557) +
              DENDRO_135 * (DENDRO_110 * DENDRO_273 + DENDRO_735) +
              DENDRO_135 * (DENDRO_669 + DENDRO_732 + DENDRO_758) +
              DENDRO_158 * (-DENDRO_124 * DENDRO_91 - DENDRO_775) +
              DENDRO_158 * (-DENDRO_241 * DENDRO_304 + DENDRO_533) +
              DENDRO_158 * (-DENDRO_536 * DENDRO_99 - DENDRO_775) +
              DENDRO_158 *
                  (-DENDRO_130 * DENDRO_344 +
                   0.5 * DENDRO_260 * grad_2_gt0[pp] - 0.25 * DENDRO_754) +
              DENDRO_158 * (-DENDRO_143 * DENDRO_388 + DENDRO_196 * DENDRO_260 -
                            DENDRO_767) +
              DENDRO_158 *
                  (-DENDRO_176 * DENDRO_241 + DENDRO_345 + DENDRO_346) +
              DENDRO_165 *
                  (DENDRO_122 * grad_2_gt3[pp] - DENDRO_241 * grad_1_gt5[pp]) -
              DENDRO_179 * (0.5 * DENDRO_108 * DENDRO_91 - DENDRO_773) +
              DENDRO_179 * (DENDRO_114 * DENDRO_304 + DENDRO_773) -
              DENDRO_179 * (-DENDRO_143 * DENDRO_679 - DENDRO_149 * DENDRO_246 -
                            DENDRO_757) +
              DENDRO_181 * (DENDRO_108 * DENDRO_241 + DENDRO_527) +
              DENDRO_181 * (-DENDRO_114 * DENDRO_189 +
                            0.5 * DENDRO_276 * grad_2_gt5[pp] - DENDRO_774) -
              DENDRO_181 * (DENDRO_143 * DENDRO_426 + DENDRO_149 * DENDRO_236 +
                            DENDRO_561) +
              DENDRO_181 * (-DENDRO_143 * DENDRO_454 - DENDRO_143 * DENDRO_536 +
                            0.5 * DENDRO_260 * grad_0_gt5[pp]) +
              DENDRO_181 * (DENDRO_200 * DENDRO_99 + DENDRO_522 - DENDRO_774) -
              DENDRO_227 * (DENDRO_122 * DENDRO_471 + DENDRO_205 * DENDRO_217 +
                            DENDRO_225 * DENDRO_440 + DENDRO_545) -
              DENDRO_335 * (DENDRO_108 * DENDRO_114 - 0.5 * DENDRO_771) -
              DENDRO_335 * (0.5 * DENDRO_113 * DENDRO_143 -
                            DENDRO_143 * DENDRO_156 - DENDRO_149 * DENDRO_304) -
              DENDRO_551 * DENDRO_686 +
              DENDRO_569 * (DENDRO_738 + DENDRO_739 + DENDRO_750) + DENDRO_574 -
              DENDRO_681 * grad_1_gt4[pp] - DENDRO_683 * grad_0_gt4[pp] -
              DENDRO_685 * grad_2_gt4[pp]) +
         DENDRO_713 * DENDRO_729 - DENDRO_714 * DENDRO_727 -
         DENDRO_715 * DENDRO_728 * DENDRO_75 - 12 * grad2_1_2_alpha[pp]) +
    DENDRO_770 * grad_1_beta1[pp] + DENDRO_770 * grad_2_beta2[pp] -
    alpha[pp] * (-At4[pp] * K[pp] + DENDRO_644 * DENDRO_747 +
                 DENDRO_746 * DENDRO_760 + DENDRO_748 * DENDRO_761) +
    beta0[pp] * grad_0_At4[pp] + beta1[pp] * grad_1_At4[pp] +
    beta2[pp] * grad_2_At4[pp];
//--
At_rhs22[pp] =
    -At5[pp] * DENDRO_10 - At5[pp] * DENDRO_11 + At5[pp] * DENDRO_15 +
    DENDRO_641 *
        (DENDRO_654 * DENDRO_689 - DENDRO_658 * DENDRO_697 +
         DENDRO_687 *
             (6.0 * DENDRO_114 * DENDRO_25 * DENDRO_47 * grad_2_gt5[pp] -
              DENDRO_127 + 4 * DENDRO_130 * DENDRO_134 * DENDRO_25 * DENDRO_31 -
              DENDRO_142 + 4 * DENDRO_143 * DENDRO_147 * DENDRO_25 * DENDRO_41 +
              DENDRO_148 * DENDRO_149 * DENDRO_335 - DENDRO_155 -
              DENDRO_158 * (DENDRO_105 * DENDRO_91 + DENDRO_772) -
              DENDRO_158 * (DENDRO_110 * DENDRO_99 + DENDRO_751) -
              DENDRO_158 * (DENDRO_130 * DENDRO_176 + DENDRO_133 * DENDRO_143) -
              DENDRO_158 * (DENDRO_196 * DENDRO_766 + 0.25 * DENDRO_740) +
              2.0 * DENDRO_163 * DENDRO_25 * DENDRO_39 -
              DENDRO_165 * (DENDRO_771 + DENDRO_777) -
              DENDRO_165 *
                  (DENDRO_143 * grad_0_gt5[pp] + DENDRO_149 * DENDRO_97) +
              4 * DENDRO_177 * DENDRO_25 * DENDRO_27 -
              DENDRO_179 * (-DENDRO_149 * DENDRO_210 - DENDRO_756) -
              DENDRO_181 * DENDRO_207 -
              DENDRO_181 * (0.25 * DENDRO_771 + 1.0 * DENDRO_777) - DENDRO_186 +
              4 * DENDRO_194 * DENDRO_25 * DENDRO_27 +
              4 * DENDRO_202 * DENDRO_25 * DENDRO_29 +
              4 * DENDRO_204 * DENDRO_25 * DENDRO_39 -
              DENDRO_227 * (DENDRO_115 * DENDRO_225 + DENDRO_150 * DENDRO_217 +
                            DENDRO_153 * DENDRO_220 + DENDRO_211) +
              4 * DENDRO_25 * DENDRO_29 *
                  (0.25 * DENDRO_749 + 1.0 * DENDRO_776) +
              2.0 * DENDRO_25 * DENDRO_29 * (DENDRO_749 + DENDRO_776) +
              2.0 * DENDRO_25 * DENDRO_29 *
                  (DENDRO_130 * grad_0_gt5[pp] + DENDRO_149 * grad_2_gt0[pp]) +
              3.0 * DENDRO_25 * DENDRO_31 * DENDRO_91 * grad_0_gt5[pp] +
              3.0 * DENDRO_25 * DENDRO_41 * DENDRO_99 * grad_1_gt5[pp] +
              4.0 * DENDRO_29 * DENDRO_75 * grad2_0_2_gt5[pp] +
              2.0 * DENDRO_31 * DENDRO_75 * grad2_0_0_gt5[pp] +
              2.0 * DENDRO_41 * DENDRO_75 * grad2_1_1_gt5[pp] +
              2.0 * DENDRO_47 * DENDRO_75 * grad2_2_2_gt5[pp] -
              DENDRO_681 * grad_1_gt5[pp] - DENDRO_683 * grad_0_gt5[pp] -
              DENDRO_685 * grad_2_gt5[pp] - DENDRO_74 - DENDRO_769 * gt5[pp] -
              DENDRO_78 - DENDRO_81 + 4.0 * grad_2_Gt0[pp] * gt2[pp] +
              4.0 * grad_2_Gt1[pp] * gt4[pp] + 4.0 * grad_2_Gt2[pp] * gt5[pp]) -
         DENDRO_692 * DENDRO_763 + DENDRO_722 * gt5[pp] -
         12 * grad2_2_2_alpha[pp]) +
    DENDRO_643 * grad_2_beta0[pp] + DENDRO_759 * grad_2_beta1[pp] -
    alpha[pp] * (At5[pp] * DENDRO_645 * DENDRO_748 - At5[pp] * K[pp] +
                 DENDRO_647 * DENDRO_747 + DENDRO_746 * DENDRO_761) +
    beta0[pp] * grad_0_At5[pp] + beta1[pp] * grad_1_At5[pp] +
    beta2[pp] * grad_2_At5[pp];
//--
K_rhs[pp] =
    -DENDRO_294 * (DENDRO_780 * (DENDRO_216 + DENDRO_778 * DENDRO_779) +
                   DENDRO_782 * (DENDRO_219 + DENDRO_779 * DENDRO_781) -
                   grad2_2_2_alpha[pp] +
                   grad_2_alpha[pp] *
                       (DENDRO_222 * DENDRO_75 + DENDRO_223 * DENDRO_75 +
                        DENDRO_224 * DENDRO_75 -
                        DENDRO_226 * (DENDRO_695 - DENDRO_696 * DENDRO_783))) +
    DENDRO_33 * DENDRO_75 * chi[pp] *
        (0.5 * DENDRO_784 * (DENDRO_625 + DENDRO_719 * DENDRO_783) +
         DENDRO_786 *
             (-DENDRO_226 * (-DENDRO_626 * DENDRO_781 + grad_0_chi[pp]) +
              DENDRO_230 * DENDRO_75 + DENDRO_231 * DENDRO_75 +
              DENDRO_623 * DENDRO_75) +
         DENDRO_788 *
             (-DENDRO_226 * (-DENDRO_626 * DENDRO_778 + grad_1_chi[pp]) +
              DENDRO_620 * DENDRO_75 + DENDRO_621 * DENDRO_75 +
              DENDRO_622 * DENDRO_75) -
         grad2_0_1_alpha[pp]) -
    DENDRO_352 * (DENDRO_780 * (DENDRO_297 + DENDRO_762 * DENDRO_778) +
                  DENDRO_784 * (DENDRO_301 + DENDRO_762 * DENDRO_783) -
                  grad2_1_1_alpha[pp] +
                  grad_1_alpha[pp] *
                      (-DENDRO_226 * (DENDRO_701 - DENDRO_702 * DENDRO_781) +
                       DENDRO_298 * DENDRO_75 + DENDRO_299 * DENDRO_75 +
                       DENDRO_300 * DENDRO_75)) +
    2 * DENDRO_39 * DENDRO_75 * chi[pp] *
        (0.5 * DENDRO_780 * (DENDRO_226 * DENDRO_778 * gt4[pp] + DENDRO_546) +
         DENDRO_786 *
             (-DENDRO_226 * (-DENDRO_551 * DENDRO_781 + grad_2_chi[pp]) +
              DENDRO_547 * DENDRO_75 + DENDRO_712) +
         DENDRO_787 *
             (-DENDRO_226 * (-DENDRO_551 * DENDRO_783 + grad_1_chi[pp]) +
              DENDRO_548 * DENDRO_75 + DENDRO_549 * DENDRO_75 +
              DENDRO_550 * DENDRO_75) -
         grad2_1_2_alpha[pp]) -
    DENDRO_405 *
        (DENDRO_782 * (DENDRO_402 + DENDRO_781 * DENDRO_785) +
         DENDRO_784 * (DENDRO_403 + DENDRO_783 * DENDRO_785) -
         grad2_0_0_alpha[pp] +
         grad_0_alpha[pp] *
             (-DENDRO_226 * (DENDRO_661 - 0.5 * DENDRO_662 * DENDRO_778) +
              DENDRO_399 * DENDRO_75 + DENDRO_401 * DENDRO_75 + DENDRO_660)) -
    DENDRO_789 * chi[pp] *
        (0.5 * DENDRO_782 * (DENDRO_470 + DENDRO_709 * DENDRO_781) +
         DENDRO_787 *
             (-DENDRO_226 * (-DENDRO_482 * DENDRO_783 + grad_0_chi[pp]) +
              DENDRO_472 * DENDRO_75 + DENDRO_473 * DENDRO_75 +
              DENDRO_474 * DENDRO_75) +
         DENDRO_788 *
             (-DENDRO_226 * (-DENDRO_482 * DENDRO_778 + grad_2_chi[pp]) +
              DENDRO_466 * DENDRO_75 + DENDRO_467 * DENDRO_75 +
              DENDRO_468 * DENDRO_75) -
         grad2_0_2_alpha[pp]) +
    DENDRO_792 * (DENDRO_18 + DENDRO_38 * DENDRO_790 + DENDRO_46 * DENDRO_790 +
                  DENDRO_52 * DENDRO_790 + DENDRO_60 * DENDRO_791 +
                  DENDRO_63 * DENDRO_791 + DENDRO_69 * DENDRO_791) +
    beta0[pp] * grad_0_K[pp] + beta1[pp] * grad_1_K[pp] +
    beta2[pp] * grad_2_K[pp];
//--
Gt_rhs0[pp] = -DENDRO_807;
//--
Gt_rhs1[pp] =
    -DENDRO_122 * DENDRO_67 * DENDRO_824 + DENDRO_137 * DENDRO_59 * DENDRO_824 +
    DENDRO_153 * DENDRO_51 * DENDRO_823 - 2.0 / 3.0 * DENDRO_17 * DENDRO_680 -
    DENDRO_233 * DENDRO_61 * DENDRO_824 - DENDRO_241 * DENDRO_45 * DENDRO_823 +
    DENDRO_249 * DENDRO_37 * DENDRO_823 + DENDRO_554 * DENDRO_822 +
    (7.0 / 3.0) * DENDRO_554 * grad2_1_2_beta1[pp] +
    (1.0 / 3.0) * DENDRO_554 * grad2_2_2_beta2[pp] + DENDRO_61 * DENDRO_797 +
    (7.0 / 3.0) * DENDRO_630 * grad2_0_1_beta1[pp] + DENDRO_67 * DENDRO_801 +
    DENDRO_680 * grad_1_beta1[pp] + DENDRO_682 * grad_0_beta1[pp] +
    DENDRO_684 * grad_2_beta1[pp] + DENDRO_800 * DENDRO_820 -
    DENDRO_800 * (-DENDRO_61 * DENDRO_814 + DENDRO_813) -
    DENDRO_800 * (-DENDRO_67 * DENDRO_816 + DENDRO_815) - DENDRO_808 -
    DENDRO_809 - DENDRO_810 - DENDRO_811 - DENDRO_812 - DENDRO_818 -
    DENDRO_819 + DENDRO_821 * grad2_0_0_beta0[pp] +
    DENDRO_821 * grad2_0_2_beta2[pp] + DENDRO_825;
//--
Gt_rhs2[pp] = -DENDRO_826;
//--
B_rhs0[pp] =
    -B0[pp] * eta - DENDRO_807 +
    lambda[2] * (beta0[pp] * grad_0_B0[pp] + beta1[pp] * grad_1_B0[pp] +
                 beta2[pp] * grad_2_B0[pp]) -
    lambda[3] * (beta0[pp] * grad_0_Gt0[pp] + beta1[pp] * grad_1_Gt0[pp] +
                 beta2[pp] * grad_2_Gt0[pp]);
//--
B_rhs1[pp] =
    -B1[pp] * eta + 2.0 * DENDRO_122 * DENDRO_68 * DENDRO_798 * alpha[pp] +
    2.0 * DENDRO_137 * DENDRO_59 * DENDRO_798 * alpha[pp] +
    2 * DENDRO_153 * DENDRO_51 * DENDRO_798 * alpha[pp] +
    (2.0 / 3.0) * DENDRO_17 * DENDRO_25 * DENDRO_252 +
    2.0 * DENDRO_233 * DENDRO_62 * DENDRO_798 * alpha[pp] +
    2 * DENDRO_242 * DENDRO_45 * DENDRO_798 * alpha[pp] +
    2 * DENDRO_249 * DENDRO_37 * DENDRO_798 * alpha[pp] -
    DENDRO_253 * grad_1_beta1[pp] - DENDRO_269 * grad_0_beta1[pp] +
    (1.0 / 3.0) * DENDRO_27 * DENDRO_75 * grad2_0_0_beta0[pp] +
    (7.0 / 3.0) * DENDRO_27 * DENDRO_75 * grad2_0_1_beta1[pp] +
    (1.0 / 3.0) * DENDRO_27 * DENDRO_75 * grad2_0_2_beta2[pp] -
    DENDRO_283 * grad_2_beta1[pp] +
    (1.0 / 3.0) * DENDRO_39 * DENDRO_75 * grad2_0_2_beta0[pp] +
    (7.0 / 3.0) * DENDRO_39 * DENDRO_75 * grad2_1_2_beta1[pp] +
    (1.0 / 3.0) * DENDRO_39 * DENDRO_75 * grad2_2_2_beta2[pp] -
    DENDRO_62 * DENDRO_797 - DENDRO_68 * DENDRO_801 + DENDRO_800 * DENDRO_820 -
    DENDRO_800 * (DENDRO_62 * DENDRO_814 + DENDRO_813) -
    DENDRO_800 * (DENDRO_68 * DENDRO_816 + DENDRO_815) - DENDRO_808 -
    DENDRO_809 - DENDRO_810 - DENDRO_811 - DENDRO_812 - DENDRO_818 -
    DENDRO_819 - DENDRO_825 * lambda[3] + beta0[pp] * grad_0_Gt1[pp] +
    beta1[pp] * grad_1_Gt1[pp] + beta2[pp] * grad_2_Gt1[pp] +
    lambda[2] * (beta0[pp] * grad_0_B1[pp] + beta1[pp] * grad_1_B1[pp] +
                 beta2[pp] * grad_2_B1[pp]);
//--
B_rhs2[pp] =
    -B2[pp] * eta - DENDRO_826 +
    lambda[2] * (beta0[pp] * grad_0_B2[pp] + beta1[pp] * grad_1_B2[pp] +
                 beta2[pp] * grad_2_B2[pp]) -
    lambda[3] * (beta0[pp] * grad_0_Gt2[pp] + beta1[pp] * grad_1_Gt2[pp] +
                 beta2[pp] * grad_2_Gt2[pp]);
// Dendro: reduced ops: 4339
// Dendro: }}}
