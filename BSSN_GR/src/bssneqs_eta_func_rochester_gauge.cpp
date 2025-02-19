// Codgen: generating unstage version
// Codgen: using rochester gauge
// Codgen: using eta func damping
//  Dendro: {{{
//  Dendro: original ops: 607507
//  Dendro: printing temp variables
const double DENDRO_0 = 2 * alpha[pp];
const double DENDRO_1 = (3.0 / 4.0) * BSSN_XI[2];
const double DENDRO_2 = pow(gt4[pp], 2);
const double DENDRO_3 = pow(gt1[pp], 2);
const double DENDRO_4 = pow(gt2[pp], 2);
const double DENDRO_5 = gt3[pp] * gt5[pp];
const double DENDRO_6 = gt1[pp] * gt4[pp];
const double DENDRO_7 = 2 * gt2[pp];
const double DENDRO_8 = DENDRO_2 * gt0[pp] + DENDRO_3 * gt5[pp] +
                        DENDRO_4 * gt3[pp] - DENDRO_5 * gt0[pp] -
                        DENDRO_6 * DENDRO_7;
const double DENDRO_9  = 1.0 / DENDRO_8;
const double DENDRO_10 = gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp];
const double DENDRO_11 = DENDRO_10 * grad_1_chi[pp];
const double DENDRO_12 = 2 * grad_0_chi[pp];
const double DENDRO_13 = gt2[pp] * gt3[pp];
const double DENDRO_14 = -DENDRO_13 + DENDRO_6;
const double DENDRO_15 = DENDRO_14 * grad_2_chi[pp];
const double DENDRO_16 = gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp];
const double DENDRO_17 = DENDRO_16 * grad_1_chi[pp];
const double DENDRO_18 = pow(grad_0_chi[pp], 2);
const double DENDRO_19 = -DENDRO_2 + DENDRO_5;
const double DENDRO_20 = pow(grad_1_chi[pp], 2);
const double DENDRO_21 = gt0[pp] * gt5[pp];
const double DENDRO_22 = DENDRO_21 - DENDRO_4;
const double DENDRO_23 = pow(grad_2_chi[pp], 2);
const double DENDRO_24 = gt0[pp] * gt3[pp];
const double DENDRO_25 = DENDRO_24 - DENDRO_3;
const double DENDRO_26 =
    BSSN_ETA_R0 *
    sqrt(DENDRO_9 * (DENDRO_11 * DENDRO_12 - DENDRO_12 * DENDRO_15 +
                     2 * DENDRO_17 * grad_2_chi[pp] - DENDRO_18 * DENDRO_19 -
                     DENDRO_20 * DENDRO_22 - DENDRO_23 * DENDRO_25)) *
    pow(1 - pow(chi[pp], BSSN_ETA_POWER[0]), -BSSN_ETA_POWER[1]);
const double DENDRO_27 = (4.0 / 3.0) * grad_0_beta0[pp];
const double DENDRO_28 = 2 * gt1[pp];
const double DENDRO_29 = (2.0 / 3.0) * gt0[pp];
const double DENDRO_30 = (1.0 / 3.0) * gt1[pp];
const double DENDRO_31 = (2.0 / 3.0) * grad_2_beta2[pp];
const double DENDRO_32 = (1.0 / 3.0) * gt2[pp];
const double DENDRO_33 = (2.0 / 3.0) * grad_1_beta1[pp];
const double DENDRO_34 = (2.0 / 3.0) * grad_0_beta0[pp];
const double DENDRO_35 = (4.0 / 3.0) * grad_1_beta1[pp];
const double DENDRO_36 = 2 * gt4[pp];
const double DENDRO_37 = (1.0 / 3.0) * gt4[pp];
const double DENDRO_38 = (4.0 / 3.0) * grad_2_beta2[pp];
const double DENDRO_39 = (2.0 / 3.0) * chi[pp];
const double DENDRO_40 = grad_0_beta0[pp] + grad_1_beta1[pp] + grad_2_beta2[pp];
const double DENDRO_41 = 2 * At1[pp];
const double DENDRO_42 = 2 * At2[pp];
const double DENDRO_43 = At0[pp] * DENDRO_10;
const double DENDRO_44 = At2[pp] * DENDRO_16;
const double DENDRO_45 = At1[pp] * DENDRO_22;
const double DENDRO_46 = DENDRO_41 * DENDRO_9;
const double DENDRO_47 = At1[pp] * DENDRO_10;
const double DENDRO_48 = At2[pp] * DENDRO_14;
const double DENDRO_49 = -DENDRO_48;
const double DENDRO_50 = At0[pp] * DENDRO_19;
const double DENDRO_51 = 2 * DENDRO_9;
const double DENDRO_52 = At0[pp] * DENDRO_51;
const double DENDRO_53 = At1[pp] * DENDRO_16;
const double DENDRO_54 = At2[pp] * DENDRO_25;
const double DENDRO_55 = DENDRO_42 * DENDRO_9;
const double DENDRO_56 =
    DENDRO_10 * grad_0_chi[pp] + DENDRO_16 * grad_2_chi[pp];
const double DENDRO_57  = -DENDRO_22 * grad_1_chi[pp] + DENDRO_56;
const double DENDRO_58  = 0.5 * DENDRO_57;
const double DENDRO_59  = 1.0 / chi[pp];
const double DENDRO_60  = DENDRO_59 * gt0[pp];
const double DENDRO_61  = 1.0 * grad_0_gt1[pp];
const double DENDRO_62  = 0.5 * grad_1_gt0[pp];
const double DENDRO_63  = DENDRO_61 - DENDRO_62;
const double DENDRO_64  = 1.0 * grad_0_gt2[pp];
const double DENDRO_65  = 0.5 * grad_2_gt0[pp];
const double DENDRO_66  = DENDRO_64 - DENDRO_65;
const double DENDRO_67  = 0.5 * grad_0_gt0[pp];
const double DENDRO_68  = DENDRO_10 * DENDRO_67 + DENDRO_16 * DENDRO_66;
const double DENDRO_69  = -DENDRO_22 * DENDRO_63 + DENDRO_68;
const double DENDRO_70  = DENDRO_58 * DENDRO_60 + DENDRO_69;
const double DENDRO_71  = 12 * grad_1_alpha[pp];
const double DENDRO_72  = DENDRO_71 * DENDRO_9;
const double DENDRO_73  = DENDRO_14 * grad_0_chi[pp];
const double DENDRO_74  = DENDRO_25 * grad_2_chi[pp];
const double DENDRO_75  = DENDRO_17 - DENDRO_73 - DENDRO_74;
const double DENDRO_76  = 0.5 * DENDRO_75;
const double DENDRO_77  = DENDRO_60 * DENDRO_76;
const double DENDRO_78  = DENDRO_14 * DENDRO_67;
const double DENDRO_79  = DENDRO_16 * DENDRO_63;
const double DENDRO_80  = DENDRO_25 * DENDRO_66;
const double DENDRO_81  = DENDRO_78 - DENDRO_79 + DENDRO_80;
const double DENDRO_82  = 12 * grad_2_alpha[pp];
const double DENDRO_83  = DENDRO_82 * DENDRO_9;
const double DENDRO_84  = DENDRO_19 * DENDRO_67;
const double DENDRO_85  = DENDRO_84 * DENDRO_9;
const double DENDRO_86  = DENDRO_14 * DENDRO_66;
const double DENDRO_87  = DENDRO_86 * DENDRO_9;
const double DENDRO_88  = DENDRO_10 * DENDRO_63;
const double DENDRO_89  = DENDRO_88 * DENDRO_9;
const double DENDRO_90  = 1.0 * grad_0_chi[pp];
const double DENDRO_91  = DENDRO_9 * gt0[pp];
const double DENDRO_92  = DENDRO_19 * grad_0_chi[pp];
const double DENDRO_93  = DENDRO_11 - DENDRO_15 - DENDRO_92;
const double DENDRO_94  = 0.5 * DENDRO_93;
const double DENDRO_95  = DENDRO_59 * (DENDRO_90 - DENDRO_91 * DENDRO_94);
const double DENDRO_96  = 12 * grad_0_alpha[pp];
const double DENDRO_97  = -grad2_0_0_chi[pp];
const double DENDRO_98  = -DENDRO_84 - DENDRO_86 + DENDRO_88;
const double DENDRO_99  = DENDRO_9 * grad_0_chi[pp];
const double DENDRO_100 = DENDRO_9 * grad_1_chi[pp];
const double DENDRO_101 = -DENDRO_78 + DENDRO_79 - DENDRO_80;
const double DENDRO_102 = DENDRO_9 * grad_2_chi[pp];
const double DENDRO_103 = 2 * DENDRO_59;
const double DENDRO_104 = 0.5 * grad_0_gt3[pp];
const double DENDRO_105 = 1.0 * grad_1_gt1[pp];
const double DENDRO_106 = DENDRO_104 - DENDRO_105;
const double DENDRO_107 = grad_0_gt4[pp] - grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_108 =
    DENDRO_10 * grad_2_gt0[pp] + DENDRO_16 * grad_0_gt5[pp];
const double DENDRO_109 = -DENDRO_107 * DENDRO_22 + DENDRO_108;
const double DENDRO_110 = DENDRO_106 * DENDRO_109;
const double DENDRO_111 = DENDRO_22 * grad_0_gt3[pp];
const double DENDRO_112 = DENDRO_10 * grad_1_gt0[pp];
const double DENDRO_113 = grad_0_gt4[pp] + grad_1_gt2[pp] - grad_2_gt1[pp];
const double DENDRO_114 = DENDRO_113 * DENDRO_16;
const double DENDRO_115 = DENDRO_112 + DENDRO_114;
const double DENDRO_116 = -DENDRO_111 + DENDRO_115;
const double DENDRO_117 = DENDRO_107 * DENDRO_116;
const double DENDRO_118 = 0.25 * DENDRO_117;
const double DENDRO_119 = pow(DENDRO_8, -2);
const double DENDRO_120 = DENDRO_119 * DENDRO_16;
const double DENDRO_121 = 4 * DENDRO_120;
const double DENDRO_122 = 0.5 * grad_0_gt5[pp];
const double DENDRO_123 = 1.0 * grad_2_gt2[pp];
const double DENDRO_124 = -DENDRO_123;
const double DENDRO_125 = DENDRO_122 + DENDRO_124;
const double DENDRO_126 = DENDRO_16 * grad_0_gt3[pp];
const double DENDRO_127 = DENDRO_14 * grad_1_gt0[pp];
const double DENDRO_128 = DENDRO_113 * DENDRO_25;
const double DENDRO_129 = DENDRO_126 - DENDRO_127 - DENDRO_128;
const double DENDRO_130 = DENDRO_125 * DENDRO_129;
const double DENDRO_131 = DENDRO_14 * grad_2_gt0[pp];
const double DENDRO_132 = DENDRO_25 * grad_0_gt5[pp];
const double DENDRO_133 = DENDRO_107 * DENDRO_16;
const double DENDRO_134 = -DENDRO_131 - DENDRO_132 + DENDRO_133;
const double DENDRO_135 = DENDRO_113 * DENDRO_134;
const double DENDRO_136 = 0.25 * DENDRO_135;
const double DENDRO_137 = -DENDRO_136;
const double DENDRO_138 = -DENDRO_126 + DENDRO_127 + DENDRO_128;
const double DENDRO_139 = 0.25 * grad_0_gt5[pp];
const double DENDRO_140 = DENDRO_138 * DENDRO_139;
const double DENDRO_141 = -grad_0_gt4[pp] + grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_142 = DENDRO_131 + DENDRO_132 - DENDRO_133;
const double DENDRO_143 = 0.5 * DENDRO_142;
const double DENDRO_144 = DENDRO_14 * grad_0_gt5[pp];
const double DENDRO_145 = DENDRO_19 * grad_2_gt0[pp];
const double DENDRO_146 = DENDRO_10 * DENDRO_107;
const double DENDRO_147 = DENDRO_144 + DENDRO_145 - DENDRO_146;
const double DENDRO_148 = 0.25 * grad_1_gt0[pp];
const double DENDRO_149 = DENDRO_147 * DENDRO_148;
const double DENDRO_150 = DENDRO_10 * grad_0_gt3[pp];
const double DENDRO_151 = DENDRO_19 * grad_1_gt0[pp];
const double DENDRO_152 = DENDRO_113 * DENDRO_14;
const double DENDRO_153 = -DENDRO_150 + DENDRO_151 + DENDRO_152;
const double DENDRO_154 = 0.25 * grad_2_gt0[pp];
const double DENDRO_155 = DENDRO_153 * DENDRO_154;
const double DENDRO_156 = 2 * DENDRO_125;
const double DENDRO_157 = DENDRO_119 * DENDRO_14;
const double DENDRO_158 = 4 * DENDRO_157;
const double DENDRO_159 = DENDRO_153 * grad_0_gt0[pp];
const double DENDRO_160 = 0.25 * DENDRO_159;
const double DENDRO_161 = DENDRO_84 + DENDRO_86 - DENDRO_88;
const double DENDRO_162 = DENDRO_161 * grad_1_gt0[pp];
const double DENDRO_163 = DENDRO_10 * DENDRO_119;
const double DENDRO_164 = 4 * DENDRO_163;
const double DENDRO_165 = 0.5 * DENDRO_66;
const double DENDRO_166 = DENDRO_147 * grad_0_gt0[pp];
const double DENDRO_167 = 0.25 * DENDRO_166;
const double DENDRO_168 = DENDRO_161 * grad_2_gt0[pp];
const double DENDRO_169 = 2.0 * DENDRO_157;
const double DENDRO_170 = DENDRO_142 * grad_2_gt0[pp];
const double DENDRO_171 = DENDRO_81 * grad_0_gt5[pp];
const double DENDRO_172 = -DENDRO_144 - DENDRO_145 + DENDRO_146;
const double DENDRO_173 = DENDRO_14 * DENDRO_172;
const double DENDRO_174 = DENDRO_10 * grad_2_gt3[pp];
const double DENDRO_175 = DENDRO_14 * grad_1_gt5[pp];
const double DENDRO_176 = DENDRO_141 * DENDRO_19;
const double DENDRO_177 = DENDRO_174 - DENDRO_175 - DENDRO_176;
const double DENDRO_178 = DENDRO_16 * DENDRO_177;
const double DENDRO_179 = DENDRO_150 - DENDRO_151 - DENDRO_152;
const double DENDRO_180 = DENDRO_10 * DENDRO_179;
const double DENDRO_181 = 0.5 * grad_2_gt5[pp];
const double DENDRO_182 = DENDRO_14 * DENDRO_181;
const double DENDRO_183 = DENDRO_125 * DENDRO_19;
const double DENDRO_184 = 0.5 * grad_1_gt5[pp];
const double DENDRO_185 = 1.0 * grad_2_gt4[pp];
const double DENDRO_186 = -DENDRO_185;
const double DENDRO_187 = DENDRO_184 + DENDRO_186;
const double DENDRO_188 = DENDRO_10 * DENDRO_187;
const double DENDRO_189 = -DENDRO_182 + DENDRO_183 - DENDRO_188;
const double DENDRO_190 = DENDRO_189 * DENDRO_25;
const double DENDRO_191 = 0.5 * grad_1_gt3[pp];
const double DENDRO_192 = DENDRO_10 * DENDRO_191;
const double DENDRO_193 = 1.0 * grad_1_gt4[pp];
const double DENDRO_194 = 0.5 * grad_2_gt3[pp];
const double DENDRO_195 = DENDRO_193 - DENDRO_194;
const double DENDRO_196 =
    DENDRO_106 * DENDRO_19 - DENDRO_14 * DENDRO_195 + DENDRO_192;
const double DENDRO_197 = DENDRO_196 * DENDRO_22;
const double DENDRO_198 = DENDRO_19 * DENDRO_98;
const double DENDRO_199 =
    DENDRO_119 * (DENDRO_173 - 1.0 * DENDRO_178 - 1.0 * DENDRO_180 +
                  DENDRO_190 + DENDRO_197 + DENDRO_198);
const double DENDRO_200 = 2.0 * grad_0_gt0[pp];
const double DENDRO_201 = DENDRO_109 * DENDRO_14;
const double DENDRO_202 = DENDRO_22 * grad_2_gt3[pp];
const double DENDRO_203 = DENDRO_16 * grad_1_gt5[pp];
const double DENDRO_204 = DENDRO_10 * DENDRO_141;
const double DENDRO_205 = DENDRO_203 + DENDRO_204;
const double DENDRO_206 = -DENDRO_202 + DENDRO_205;
const double DENDRO_207 = DENDRO_16 * DENDRO_206;
const double DENDRO_208 = DENDRO_10 * DENDRO_116;
const double DENDRO_209 = DENDRO_16 * DENDRO_181;
const double DENDRO_210 =
    -DENDRO_10 * DENDRO_125 + DENDRO_187 * DENDRO_22 + DENDRO_209;
const double DENDRO_211 = DENDRO_210 * DENDRO_25;
const double DENDRO_212 = DENDRO_191 * DENDRO_22;
const double DENDRO_213 = DENDRO_16 * DENDRO_195;
const double DENDRO_214 = DENDRO_10 * DENDRO_106;
const double DENDRO_215 = -DENDRO_212 + DENDRO_213 - DENDRO_214;
const double DENDRO_216 = DENDRO_215 * DENDRO_22;
const double DENDRO_217 = DENDRO_19 * DENDRO_69;
const double DENDRO_218 =
    DENDRO_119 * (DENDRO_201 - 1.0 * DENDRO_207 - 1.0 * DENDRO_208 +
                  DENDRO_211 + DENDRO_216 + DENDRO_217);
const double DENDRO_219 = 2.0 * grad_1_gt0[pp];
const double DENDRO_220 = DENDRO_134 * DENDRO_14;
const double DENDRO_221 = DENDRO_16 * grad_2_gt3[pp];
const double DENDRO_222 = DENDRO_25 * grad_1_gt5[pp];
const double DENDRO_223 = DENDRO_14 * DENDRO_141;
const double DENDRO_224 = DENDRO_221 - DENDRO_222 - DENDRO_223;
const double DENDRO_225 = DENDRO_16 * DENDRO_224;
const double DENDRO_226 = DENDRO_10 * DENDRO_129;
const double DENDRO_227 = DENDRO_181 * DENDRO_25;
const double DENDRO_228 = DENDRO_125 * DENDRO_14;
const double DENDRO_229 = DENDRO_16 * DENDRO_187;
const double DENDRO_230 = -DENDRO_227 + DENDRO_228 - DENDRO_229;
const double DENDRO_231 = DENDRO_230 * DENDRO_25;
const double DENDRO_232 = DENDRO_16 * DENDRO_191;
const double DENDRO_233 =
    DENDRO_106 * DENDRO_14 - DENDRO_195 * DENDRO_25 + DENDRO_232;
const double DENDRO_234 = DENDRO_22 * DENDRO_233;
const double DENDRO_235 = DENDRO_101 * DENDRO_19;
const double DENDRO_236 =
    DENDRO_119 * (DENDRO_220 - 1.0 * DENDRO_225 - 1.0 * DENDRO_226 +
                  DENDRO_231 + DENDRO_234 + DENDRO_235);
const double DENDRO_237 = 2.0 * grad_2_gt0[pp];
const double DENDRO_238 = 2.0 * DENDRO_163;
const double DENDRO_239 = DENDRO_138 * grad_2_gt0[pp];
const double DENDRO_240 = DENDRO_113 * DENDRO_81;
const double DENDRO_241 = 3 * DENDRO_59;
const double DENDRO_242 = grad_0_chi[pp] * grad_1_chi[pp];
const double DENDRO_243 = 2 * DENDRO_10;
const double DENDRO_244 =
    DENDRO_243 * (-DENDRO_241 * DENDRO_242 + 2 * grad2_0_1_chi[pp]);
const double DENDRO_245 = DENDRO_241 * grad_2_chi[pp];
const double DENDRO_246 = 2 * DENDRO_14;
const double DENDRO_247 =
    DENDRO_246 * (-DENDRO_245 * grad_0_chi[pp] + 2 * grad2_0_2_chi[pp]);
const double DENDRO_248 = 2 * DENDRO_16;
const double DENDRO_249 =
    DENDRO_248 * (-DENDRO_245 * grad_1_chi[pp] + 2 * grad2_1_2_chi[pp]);
const double DENDRO_250 =
    DENDRO_19 * (-DENDRO_18 * DENDRO_241 + 2 * grad2_0_0_chi[pp]);
const double DENDRO_251 =
    DENDRO_22 * (-DENDRO_20 * DENDRO_241 + 2 * grad2_1_1_chi[pp]);
const double DENDRO_252 =
    DENDRO_25 * (-DENDRO_23 * DENDRO_241 + 2 * grad2_2_2_chi[pp]);
const double DENDRO_253 = -1.0 * DENDRO_201 + DENDRO_207 + DENDRO_208 -
                          DENDRO_211 - DENDRO_216 - DENDRO_217;
const double DENDRO_254 = 2 * DENDRO_100 * DENDRO_253;
const double DENDRO_255 = -1.0 * DENDRO_173 + DENDRO_178 + DENDRO_180 -
                          DENDRO_190 - DENDRO_197 - DENDRO_198;
const double DENDRO_256 = 2 * DENDRO_255 * DENDRO_99;
const double DENDRO_257 = -1.0 * DENDRO_220 + DENDRO_225 + DENDRO_226 -
                          DENDRO_231 - DENDRO_234 - DENDRO_235;
const double DENDRO_258 = 2 * DENDRO_102 * DENDRO_257;
const double DENDRO_259 = -DENDRO_244 + DENDRO_247 - DENDRO_249 + DENDRO_250 +
                          DENDRO_251 + DENDRO_252 + DENDRO_254 + DENDRO_256 +
                          DENDRO_258;
const double DENDRO_260 = DENDRO_59 * DENDRO_91;
const double DENDRO_261 = DENDRO_119 * DENDRO_25;
const double DENDRO_262 = 4 * DENDRO_261;
const double DENDRO_263 = 0.25 * grad_0_gt4[pp];
const double DENDRO_264 = -DENDRO_263;
const double DENDRO_265 = 0.75 * grad_1_gt2[pp];
const double DENDRO_266 = 0.25 * grad_2_gt1[pp];
const double DENDRO_267 = DENDRO_119 * DENDRO_22;
const double DENDRO_268 = 4 * DENDRO_267;
const double DENDRO_269 = DENDRO_268 * (DENDRO_264 + DENDRO_265 + DENDRO_266);
const double DENDRO_270 = DENDRO_64 + DENDRO_65;
const double DENDRO_271 = DENDRO_119 * DENDRO_19;
const double DENDRO_272 = 4 * DENDRO_271;
const double DENDRO_273 = DENDRO_153 * grad_1_gt0[pp];
const double DENDRO_274 = 3.0 * DENDRO_267;
const double DENDRO_275 = DENDRO_147 * grad_2_gt0[pp];
const double DENDRO_276 = 3.0 * DENDRO_261;
const double DENDRO_277 = 6.0 * grad_0_gt0[pp];
const double DENDRO_278 = pow(chi[pp], -2);
const double DENDRO_279 = 4 * gt1[pp];
const double DENDRO_280 = 4 * gt2[pp];
const double DENDRO_281 = 0.5 * DENDRO_63;
const double DENDRO_282 = DENDRO_109 * DENDRO_281;
const double DENDRO_283 = DENDRO_14 * DENDRO_9;
const double DENDRO_284 = 4 * DENDRO_283;
const double DENDRO_285 = 0.25 * grad_0_gt3[pp];
const double DENDRO_286 = DENDRO_109 * DENDRO_285;
const double DENDRO_287 = 0.5 * DENDRO_141;
const double DENDRO_288 = DENDRO_116 * DENDRO_281;
const double DENDRO_289 = DENDRO_19 * DENDRO_9;
const double DENDRO_290 = 2.0 * DENDRO_289;
const double DENDRO_291 = DENDRO_22 * DENDRO_9;
const double DENDRO_292 = 2.0 * DENDRO_291;
const double DENDRO_293 = DENDRO_25 * DENDRO_9;
const double DENDRO_294 = 2.0 * DENDRO_293;
const double DENDRO_295 = DENDRO_116 * grad_1_gt0[pp];
const double DENDRO_296 = DENDRO_69 * grad_0_gt3[pp];
const double DENDRO_297 = DENDRO_109 * grad_1_gt0[pp];
const double DENDRO_298 = DENDRO_107 * DENDRO_69;
const double DENDRO_299 = 4.0 * DENDRO_9;
const double DENDRO_300 = DENDRO_10 * DENDRO_299;
const double DENDRO_301 = DENDRO_16 * DENDRO_299;
const double DENDRO_302 = 0.25 * grad_1_gt2[pp];
const double DENDRO_303 = 0.75 * grad_2_gt1[pp];
const double DENDRO_304 = 4 * DENDRO_119;
const double DENDRO_305 =
    -DENDRO_109 * DENDRO_262 * (DENDRO_264 + DENDRO_302 + DENDRO_303) -
    DENDRO_116 * DENDRO_268 * (DENDRO_105 - DENDRO_285) +
    DENDRO_121 * (DENDRO_116 * DENDRO_287 + DENDRO_286) -
    DENDRO_158 * (DENDRO_141 * DENDRO_69 + DENDRO_282) +
    DENDRO_164 * (-2 * DENDRO_106 * DENDRO_69 + DENDRO_288) -
    DENDRO_169 * (DENDRO_297 + DENDRO_298) - DENDRO_18 * DENDRO_278 -
    DENDRO_217 * DENDRO_304 * (DENDRO_61 + DENDRO_62) +
    DENDRO_238 * (DENDRO_295 + DENDRO_296) + DENDRO_279 * grad_0_Gt1[pp] +
    DENDRO_280 * grad_0_Gt2[pp] + DENDRO_284 * grad2_0_2_gt0[pp] +
    DENDRO_290 * grad2_0_0_gt0[pp] + DENDRO_292 * grad2_1_1_gt0[pp] +
    DENDRO_294 * grad2_2_2_gt0[pp] - DENDRO_300 * grad2_0_1_gt0[pp] -
    DENDRO_301 * grad2_1_2_gt0[pp] + 4 * grad_0_Gt0[pp] * gt0[pp];
const double DENDRO_306 = 3 * alpha[pp];
const double DENDRO_307 = DENDRO_59 * gt5[pp];
const double DENDRO_308 = DENDRO_210 + DENDRO_307 * DENDRO_58;
const double DENDRO_309 = 4 * grad_1_alpha[pp];
const double DENDRO_310 = DENDRO_309 * DENDRO_9;
const double DENDRO_311 = DENDRO_307 * DENDRO_94;
const double DENDRO_312 = 4 * grad_0_alpha[pp];
const double DENDRO_313 = DENDRO_312 * DENDRO_9;
const double DENDRO_314 = DENDRO_227 * DENDRO_9;
const double DENDRO_315 = DENDRO_228 * DENDRO_9;
const double DENDRO_316 = DENDRO_229 * DENDRO_9;
const double DENDRO_317 = 1.0 * grad_2_chi[pp];
const double DENDRO_318 = DENDRO_9 * gt5[pp];
const double DENDRO_319 = DENDRO_59 * (DENDRO_317 - DENDRO_318 * DENDRO_76);
const double DENDRO_320 = 4 * grad_2_alpha[pp];
const double DENDRO_321 = -grad2_2_2_chi[pp];
const double DENDRO_322 = DENDRO_13 - DENDRO_6;
const double DENDRO_323 = -DENDRO_122;
const double DENDRO_324 = DENDRO_123 + DENDRO_323;
const double DENDRO_325 = DENDRO_2 - DENDRO_5;
const double DENDRO_326 = -DENDRO_184 + DENDRO_185;
const double DENDRO_327 =
    DENDRO_10 * DENDRO_326 + DENDRO_181 * DENDRO_322 + DENDRO_324 * DENDRO_325;
const double DENDRO_328 = -DENDRO_21 + DENDRO_4;
const double DENDRO_329 =
    DENDRO_10 * DENDRO_324 + DENDRO_209 + DENDRO_326 * DENDRO_328;
const double DENDRO_330 = -DENDRO_24 + DENDRO_3;
const double DENDRO_331 = DENDRO_181 * DENDRO_330;
const double DENDRO_332 = DENDRO_322 * DENDRO_324;
const double DENDRO_333 = DENDRO_16 * DENDRO_326;
const double DENDRO_334 = 0.5 * DENDRO_187;
const double DENDRO_335 = DENDRO_109 * DENDRO_334;
const double DENDRO_336 = -DENDRO_335;
const double DENDRO_337 = DENDRO_113 * DENDRO_210;
const double DENDRO_338 = 0.5 * DENDRO_125;
const double DENDRO_339 = -DENDRO_172 * DENDRO_338;
const double DENDRO_340 = 2 * DENDRO_66;
const double DENDRO_341 = DENDRO_134 * grad_2_gt5[pp];
const double DENDRO_342 = 0.25 * DENDRO_341;
const double DENDRO_343 = DENDRO_230 * grad_0_gt5[pp];
const double DENDRO_344 = DENDRO_177 * DENDRO_338;
const double DENDRO_345 = -DENDRO_344;
const double DENDRO_346 = DENDRO_113 * DENDRO_189;
const double DENDRO_347 = DENDRO_224 * grad_2_gt5[pp];
const double DENDRO_348 = 0.25 * DENDRO_347;
const double DENDRO_349 = DENDRO_230 * grad_1_gt5[pp];
const double DENDRO_350 = DENDRO_177 * DENDRO_66;
const double DENDRO_351 = DENDRO_141 * DENDRO_172;
const double DENDRO_352 = 0.25 * DENDRO_351;
const double DENDRO_353 = DENDRO_139 * DENDRO_224;
const double DENDRO_354 = 0.25 * grad_1_gt5[pp];
const double DENDRO_355 = DENDRO_134 * DENDRO_354;
const double DENDRO_356 = DENDRO_154 * DENDRO_177;
const double DENDRO_357 = 0.5 * DENDRO_113;
const double DENDRO_358 = DENDRO_119 * DENDRO_255;
const double DENDRO_359 = 2.0 * DENDRO_358;
const double DENDRO_360 = DENDRO_119 * DENDRO_253;
const double DENDRO_361 = 2.0 * DENDRO_360;
const double DENDRO_362 = DENDRO_119 * DENDRO_257;
const double DENDRO_363 = 2.0 * DENDRO_362;
const double DENDRO_364 = DENDRO_141 * DENDRO_189;
const double DENDRO_365 = 2.0 * DENDRO_120;
const double DENDRO_366 = DENDRO_189 * grad_2_gt0[pp];
const double DENDRO_367 = DENDRO_244 - DENDRO_247 + DENDRO_249 - DENDRO_250 -
                          DENDRO_251 - DENDRO_252 - DENDRO_254 - DENDRO_256 -
                          DENDRO_258;
const double DENDRO_368 = DENDRO_367 * DENDRO_59;
const double DENDRO_369 = -DENDRO_266;
const double DENDRO_370 = DENDRO_268 * (DENDRO_263 + DENDRO_265 + DENDRO_369);
const double DENDRO_371 = DENDRO_272 * (-DENDRO_154 + DENDRO_64);
const double DENDRO_372 = DENDRO_134 * grad_0_gt5[pp];
const double DENDRO_373 = 3.0 * DENDRO_271;
const double DENDRO_374 = DENDRO_224 * grad_1_gt5[pp];
const double DENDRO_375 = 6.0 * DENDRO_119;
const double DENDRO_376 = 4 * gt4[pp];
const double DENDRO_377 = -DENDRO_206 * DENDRO_334;
const double DENDRO_378 = DENDRO_109 * DENDRO_195;
const double DENDRO_379 = DENDRO_107 * DENDRO_206;
const double DENDRO_380 = 0.25 * DENDRO_379;
const double DENDRO_381 = 0.25 * grad_2_gt3[pp];
const double DENDRO_382 = DENDRO_109 * DENDRO_381;
const double DENDRO_383 = DENDRO_206 * grad_1_gt5[pp];
const double DENDRO_384 = DENDRO_210 * grad_2_gt3[pp];
const double DENDRO_385 = DENDRO_109 * grad_1_gt5[pp];
const double DENDRO_386 = DENDRO_107 * DENDRO_210;
const double DENDRO_387 = 0.75 * grad_0_gt4[pp];
const double DENDRO_388 =
    -DENDRO_109 * DENDRO_272 * (DENDRO_302 + DENDRO_369 + DENDRO_387) +
    DENDRO_121 * (2 * DENDRO_195 * DENDRO_210 + DENDRO_377) +
    DENDRO_164 * (DENDRO_378 + DENDRO_380) +
    DENDRO_164 * (DENDRO_206 * DENDRO_357 + DENDRO_382) -
    DENDRO_169 * (DENDRO_385 + DENDRO_386) -
    DENDRO_206 * DENDRO_268 * (DENDRO_193 - DENDRO_381) -
    DENDRO_211 * DENDRO_304 * (DENDRO_184 + DENDRO_185) -
    DENDRO_23 * DENDRO_278 + DENDRO_280 * grad_2_Gt0[pp] +
    DENDRO_284 * grad2_0_2_gt5[pp] + DENDRO_290 * grad2_0_0_gt5[pp] +
    DENDRO_292 * grad2_1_1_gt5[pp] + DENDRO_294 * grad2_2_2_gt5[pp] -
    DENDRO_300 * grad2_0_1_gt5[pp] - DENDRO_301 * grad2_1_2_gt5[pp] +
    DENDRO_365 * (DENDRO_383 + DENDRO_384) + DENDRO_376 * grad_2_Gt1[pp] +
    4 * grad_2_Gt2[pp] * gt5[pp];
const double DENDRO_389 = DENDRO_59 * gt3[pp];
const double DENDRO_390 = DENDRO_320 * DENDRO_9;
const double DENDRO_391 = DENDRO_212 * DENDRO_9;
const double DENDRO_392 = DENDRO_213 * DENDRO_9;
const double DENDRO_393 = DENDRO_214 * DENDRO_9;
const double DENDRO_394 = 1.0 * grad_1_chi[pp];
const double DENDRO_395 = DENDRO_9 * gt3[pp];
const double DENDRO_396 = DENDRO_59 * (DENDRO_394 - DENDRO_395 * DENDRO_58);
const double DENDRO_397 = -grad2_1_1_chi[pp];
const double DENDRO_398 = -DENDRO_104 + DENDRO_105;
const double DENDRO_399 =
    DENDRO_192 + DENDRO_195 * DENDRO_322 + DENDRO_325 * DENDRO_398;
const double DENDRO_400 = DENDRO_191 * DENDRO_328;
const double DENDRO_401 = DENDRO_10 * DENDRO_398;
const double DENDRO_402 =
    DENDRO_195 * DENDRO_330 + DENDRO_232 + DENDRO_322 * DENDRO_398;
const double DENDRO_403 = DENDRO_177 * DENDRO_63;
const double DENDRO_404 = DENDRO_141 * DENDRO_179;
const double DENDRO_405 = 0.25 * DENDRO_404;
const double DENDRO_406 = DENDRO_148 * DENDRO_177;
const double DENDRO_407 = 0.5 * DENDRO_107;
const double DENDRO_408 = DENDRO_129 * DENDRO_354;
const double DENDRO_409 = DENDRO_113 * DENDRO_224;
const double DENDRO_410 = 0.25 * DENDRO_409;
const double DENDRO_411 = DENDRO_129 * DENDRO_187;
const double DENDRO_412 = 0.5 * DENDRO_106;
const double DENDRO_413 = DENDRO_177 * DENDRO_412;
const double DENDRO_414 = -DENDRO_413;
const double DENDRO_415 = DENDRO_107 * DENDRO_196;
const double DENDRO_416 = 0.5 * DENDRO_195;
const double DENDRO_417 = DENDRO_224 * DENDRO_416;
const double DENDRO_418 = 2 * DENDRO_187 * DENDRO_233;
const double DENDRO_419 = DENDRO_206 * grad_1_gt3[pp];
const double DENDRO_420 = 0.25 * DENDRO_419;
const double DENDRO_421 = DENDRO_215 * grad_2_gt3[pp];
const double DENDRO_422 = DENDRO_129 * DENDRO_416;
const double DENDRO_423 = DENDRO_107 * DENDRO_233;
const double DENDRO_424 = -DENDRO_179 * DENDRO_412;
const double DENDRO_425 = 2 * DENDRO_196 * DENDRO_63;
const double DENDRO_426 = DENDRO_116 * grad_1_gt3[pp];
const double DENDRO_427 = 0.25 * DENDRO_426;
const double DENDRO_428 = DENDRO_215 * grad_0_gt3[pp];
const double DENDRO_429 = DENDRO_196 * grad_1_gt0[pp];
const double DENDRO_430 = DENDRO_141 * DENDRO_196;
const double DENDRO_431 = DENDRO_233 * grad_1_gt5[pp];
const double DENDRO_432 = DENDRO_113 * DENDRO_233;
const double DENDRO_433 = -DENDRO_302;
const double DENDRO_434 = DENDRO_262 * (DENDRO_263 + DENDRO_303 + DENDRO_433);
const double DENDRO_435 = DENDRO_272 * (-DENDRO_148 + DENDRO_61);
const double DENDRO_436 = DENDRO_272 * (DENDRO_266 + DENDRO_387 + DENDRO_433);
const double DENDRO_437 = DENDRO_206 * DENDRO_285;
const double DENDRO_438 = DENDRO_116 * DENDRO_381;
const double DENDRO_439 = DENDRO_116 * grad_0_gt3[pp];
const double DENDRO_440 = DENDRO_206 * grad_2_gt3[pp];
const double DENDRO_441 =
    -DENDRO_158 * (DENDRO_104 * DENDRO_206 + DENDRO_438) -
    DENDRO_158 * (DENDRO_116 * DENDRO_194 + DENDRO_437) -
    DENDRO_197 * DENDRO_304 * (DENDRO_104 + DENDRO_105) -
    DENDRO_20 * DENDRO_278 -
    DENDRO_234 * DENDRO_304 * (DENDRO_193 + DENDRO_194) -
    DENDRO_276 * DENDRO_440 + DENDRO_279 * grad_1_Gt0[pp] +
    DENDRO_284 * grad2_0_2_gt3[pp] + DENDRO_290 * grad2_0_0_gt3[pp] +
    DENDRO_292 * grad2_1_1_gt3[pp] + DENDRO_294 * grad2_2_2_gt3[pp] -
    DENDRO_300 * grad2_0_1_gt3[pp] - DENDRO_301 * grad2_1_2_gt3[pp] -
    DENDRO_373 * DENDRO_439 + DENDRO_376 * grad_1_Gt2[pp] +
    4 * grad_1_Gt1[pp] * gt3[pp];
const double DENDRO_442 = DENDRO_325 * DENDRO_67;
const double DENDRO_443 = DENDRO_322 * DENDRO_66;
const double DENDRO_444 = DENDRO_328 * DENDRO_63 + DENDRO_68;
const double DENDRO_445 =
    DENDRO_322 * DENDRO_67 + DENDRO_330 * DENDRO_66 + DENDRO_79;
const double DENDRO_446 = DENDRO_134 * DENDRO_165;
const double DENDRO_447 = DENDRO_172 * grad_0_gt0[pp];
const double DENDRO_448 = 0.25 * DENDRO_447;
const double DENDRO_449 = DENDRO_98 * grad_2_gt0[pp];
const double DENDRO_450 = DENDRO_129 * DENDRO_139;
const double DENDRO_451 = DENDRO_148 * DENDRO_172;
const double DENDRO_452 = DENDRO_154 * DENDRO_179;
const double DENDRO_453 = DENDRO_129 * DENDRO_165;
const double DENDRO_454 = DENDRO_179 * grad_0_gt0[pp];
const double DENDRO_455 = 0.25 * DENDRO_454;
const double DENDRO_456 = DENDRO_98 * grad_1_gt0[pp];
const double DENDRO_457 = DENDRO_101 * DENDRO_113;
const double DENDRO_458 = DENDRO_101 * grad_0_gt5[pp];
const double DENDRO_459 = DENDRO_179 * grad_1_gt0[pp];
const double DENDRO_460 = DENDRO_172 * grad_2_gt0[pp];
const double DENDRO_461 = DENDRO_206 * grad_1_gt0[pp];
const double DENDRO_462 = DENDRO_109 * grad_0_gt3[pp];
const double DENDRO_463 = DENDRO_113 * DENDRO_116;
const double DENDRO_464 = DENDRO_462 + DENDRO_463;
const double DENDRO_465 = DENDRO_224 * grad_2_gt0[pp];
const double DENDRO_466 = DENDRO_129 * grad_0_gt5[pp];
const double DENDRO_467 = DENDRO_135 + DENDRO_466;
const double DENDRO_468 = DENDRO_139 * DENDRO_172;
const double DENDRO_469 = -DENDRO_125 * DENDRO_230;
const double DENDRO_470 = DENDRO_210 * DENDRO_287 + 0.25 * DENDRO_385;
const double DENDRO_471 = DENDRO_113 * DENDRO_179;
const double DENDRO_472 = 0.25 * DENDRO_471;
const double DENDRO_473 = DENDRO_116 * DENDRO_416;
const double DENDRO_474 = -DENDRO_206 * DENDRO_412;
const double DENDRO_475 = DENDRO_473 + DENDRO_474;
const double DENDRO_476 = DENDRO_66 * DENDRO_98;
const double DENDRO_477 = 0.25 * DENDRO_297 + DENDRO_357 * DENDRO_69;
const double DENDRO_478 = DENDRO_134 * DENDRO_154;
const double DENDRO_479 = 0.25 * DENDRO_109;
const double DENDRO_480 = DENDRO_141 * DENDRO_479 + DENDRO_184 * DENDRO_69;
const double DENDRO_481 = -DENDRO_134 * DENDRO_338;
const double DENDRO_482 = DENDRO_101 * DENDRO_181;
const double DENDRO_483 = DENDRO_154 * DENDRO_172;
const double DENDRO_484 = DENDRO_189 * DENDRO_67;
const double DENDRO_485 = DENDRO_165 * DENDRO_172 + DENDRO_484;
const double DENDRO_486 = DENDRO_107 * DENDRO_479;
const double DENDRO_487 = DENDRO_113 * DENDRO_479 + DENDRO_210 * DENDRO_62;
const double DENDRO_488 = DENDRO_106 * DENDRO_210;
const double DENDRO_489 = 0.5 * DENDRO_378;
const double DENDRO_490 = -DENDRO_488 + DENDRO_489;
const double DENDRO_491 = DENDRO_189 * DENDRO_62;
const double DENDRO_492 = DENDRO_139 * DENDRO_179 + DENDRO_491;
const double DENDRO_493 = DENDRO_230 * DENDRO_287;
const double DENDRO_494 = DENDRO_353 + DENDRO_355;
const double DENDRO_495 = DENDRO_141 * DENDRO_206;
const double DENDRO_496 = 0.25 * DENDRO_495;
const double DENDRO_497 = DENDRO_104 * DENDRO_210;
const double DENDRO_498 = DENDRO_116 * DENDRO_354 + DENDRO_497;
const double DENDRO_499 = 0.25 * DENDRO_113;
const double DENDRO_500 = DENDRO_172 * DENDRO_499 + DENDRO_491;
const double DENDRO_501 = DENDRO_224 * DENDRO_338;
const double DENDRO_502 = -DENDRO_501;
const double DENDRO_503 = 0.25 * DENDRO_129;
const double DENDRO_504 = DENDRO_503 * grad_2_gt5[pp];
const double DENDRO_505 = DENDRO_230 * DENDRO_357 + DENDRO_504;
const double DENDRO_506 = DENDRO_206 * DENDRO_281;
const double DENDRO_507 = -0.5 * DENDRO_110 + DENDRO_195 * DENDRO_69;
const double DENDRO_508 = 0.25 * DENDRO_177;
const double DENDRO_509 = DENDRO_508 * grad_0_gt0[pp];
const double DENDRO_510 = DENDRO_165 * DENDRO_179;
const double DENDRO_511 = DENDRO_509 + DENDRO_510;
const double DENDRO_512 = DENDRO_357 * DENDRO_98 + DENDRO_451;
const double DENDRO_513 = 0.25 * DENDRO_141;
const double DENDRO_514 = DENDRO_134 * DENDRO_513;
const double DENDRO_515 = DENDRO_165 * DENDRO_224;
const double DENDRO_516 = DENDRO_101 * DENDRO_184;
const double DENDRO_517 = DENDRO_141 * DENDRO_224;
const double DENDRO_518 = DENDRO_129 * grad_1_gt5[pp];
const double DENDRO_519 = DENDRO_409 + DENDRO_518;
const double DENDRO_520 = 1.0 * DENDRO_267;
const double DENDRO_521 = -grad2_0_2_chi[pp];
const double DENDRO_522 = DENDRO_322 * grad_0_gt5[pp];
const double DENDRO_523 = DENDRO_325 * grad_2_gt0[pp];
const double DENDRO_524 = 0.5 * DENDRO_99;
const double DENDRO_525 = DENDRO_107 * DENDRO_328 + DENDRO_108;
const double DENDRO_526 = 0.5 * DENDRO_100;
const double DENDRO_527 = DENDRO_330 * grad_0_gt5[pp];
const double DENDRO_528 = DENDRO_322 * grad_2_gt0[pp];
const double DENDRO_529 = 0.5 * DENDRO_102;
const double DENDRO_530 = 2.0 * gt2[pp];
const double DENDRO_531 = DENDRO_530 * grad_0_Gt0[pp];
const double DENDRO_532 = 2.0 * gt4[pp];
const double DENDRO_533 = DENDRO_532 * grad_0_Gt1[pp];
const double DENDRO_534 = 2.0 * gt5[pp];
const double DENDRO_535 = DENDRO_534 * grad_0_Gt2[pp];
const double DENDRO_536 = 2.0 * grad_2_Gt0[pp];
const double DENDRO_537 = DENDRO_536 * gt0[pp];
const double DENDRO_538 = 2.0 * grad_2_Gt1[pp];
const double DENDRO_539 = DENDRO_538 * gt1[pp];
const double DENDRO_540 = DENDRO_530 * grad_2_Gt2[pp];
const double DENDRO_541 = DENDRO_278 * grad_2_chi[pp];
const double DENDRO_542 = -DENDRO_541 * grad_0_chi[pp];
const double DENDRO_543 = DENDRO_284 * grad2_0_2_gt2[pp];
const double DENDRO_544 = DENDRO_290 * grad2_0_0_gt2[pp];
const double DENDRO_545 = DENDRO_292 * grad2_1_1_gt2[pp];
const double DENDRO_546 = DENDRO_294 * grad2_2_2_gt2[pp];
const double DENDRO_547 = -DENDRO_300 * grad2_0_1_gt2[pp];
const double DENDRO_548 = -DENDRO_301 * grad2_1_2_gt2[pp];
const double DENDRO_549 = DENDRO_9 * gt2[pp];
const double DENDRO_550 =
    -DENDRO_103 *
        (DENDRO_521 + DENDRO_524 * (DENDRO_146 + DENDRO_522 + DENDRO_523) +
         DENDRO_525 * DENDRO_526 +
         DENDRO_529 * (DENDRO_133 + DENDRO_527 + DENDRO_528)) +
    DENDRO_359 * grad_0_gt2[pp] + DENDRO_361 * grad_1_gt2[pp] +
    DENDRO_363 * grad_2_gt2[pp] + DENDRO_368 * DENDRO_549 + DENDRO_531 +
    DENDRO_533 + DENDRO_535 + DENDRO_537 + DENDRO_539 + DENDRO_540 +
    DENDRO_542 + DENDRO_543 + DENDRO_544 + DENDRO_545 + DENDRO_546 +
    DENDRO_547 + DENDRO_548;
const double DENDRO_551 = DENDRO_144 * DENDRO_9;
const double DENDRO_552 = DENDRO_145 * DENDRO_9;
const double DENDRO_553 = DENDRO_146 * DENDRO_9;
const double DENDRO_554 =
    DENDRO_59 * (-DENDRO_549 * DENDRO_93 + grad_2_chi[pp]);
const double DENDRO_555 = 2.0 * grad_0_alpha[pp];
const double DENDRO_556 = DENDRO_131 * DENDRO_9;
const double DENDRO_557 = DENDRO_132 * DENDRO_9;
const double DENDRO_558 = DENDRO_133 * DENDRO_9;
const double DENDRO_559 =
    DENDRO_59 * (-DENDRO_549 * DENDRO_75 + grad_0_chi[pp]);
const double DENDRO_560 = 2.0 * grad_2_alpha[pp];
const double DENDRO_561 = 2.0 * grad_1_alpha[pp];
const double DENDRO_562 = DENDRO_59 * gt2[pp];
const double DENDRO_563 = DENDRO_9 * (DENDRO_109 + DENDRO_562 * DENDRO_57);
const double DENDRO_564 =
    DENDRO_555 * (-DENDRO_551 - DENDRO_552 + DENDRO_553 - DENDRO_554) +
    DENDRO_560 * (-DENDRO_556 - DENDRO_557 + DENDRO_558 - DENDRO_559) +
    DENDRO_561 * DENDRO_563 - 4 * grad2_0_2_alpha[pp];
const double DENDRO_565 = DENDRO_177 * grad_2_gt0[pp];
const double DENDRO_566 = DENDRO_179 * grad_0_gt5[pp] + DENDRO_565;
const double DENDRO_567 = 0.25 * DENDRO_517;
const double DENDRO_568 = DENDRO_134 * DENDRO_139 + DENDRO_482;
const double DENDRO_569 = DENDRO_116 * DENDRO_334;
const double DENDRO_570 = -DENDRO_569;
const double DENDRO_571 = DENDRO_355 + DENDRO_504;
const double DENDRO_572 = -DENDRO_179 * DENDRO_338;
const double DENDRO_573 = DENDRO_287 * DENDRO_98;
const double DENDRO_574 = DENDRO_154 * DENDRO_224;
const double DENDRO_575 = DENDRO_450 + DENDRO_516;
const double DENDRO_576 = DENDRO_177 * grad_1_gt0[pp];
const double DENDRO_577 = DENDRO_471 + DENDRO_576;
const double DENDRO_578 = DENDRO_109 * grad_2_gt3[pp];
const double DENDRO_579 = DENDRO_116 * grad_1_gt5[pp] + DENDRO_578;
const double DENDRO_580 = 0.25 * DENDRO_463;
const double DENDRO_581 = DENDRO_194 * DENDRO_69;
const double DENDRO_582 = DENDRO_148 * DENDRO_206 + DENDRO_581;
const double DENDRO_583 = DENDRO_120 * (DENDRO_495 + DENDRO_579) -
                          DENDRO_158 * (DENDRO_480 + DENDRO_486) -
                          DENDRO_158 * (-DENDRO_187 * DENDRO_69 + DENDRO_487) +
                          DENDRO_164 * (DENDRO_118 + DENDRO_507) +
                          DENDRO_164 * (DENDRO_580 + DENDRO_582) -
                          DENDRO_262 * (DENDRO_336 + DENDRO_470) -
                          DENDRO_268 * (DENDRO_438 + DENDRO_475) -
                          DENDRO_272 * (0.5 * DENDRO_298 + DENDRO_477);
const double DENDRO_584 = DENDRO_172 * grad_0_gt3[pp];
const double DENDRO_585 = DENDRO_134 * grad_2_gt3[pp];
const double DENDRO_586 = 0.25 * DENDRO_383;
const double DENDRO_587 = -DENDRO_187 * DENDRO_230;
const double DENDRO_588 = DENDRO_139 * DENDRO_177 + DENDRO_189 * DENDRO_407;
const double DENDRO_589 = DENDRO_195 * DENDRO_215;
const double DENDRO_590 = DENDRO_196 * DENDRO_357;
const double DENDRO_591 = DENDRO_177 * DENDRO_285 + DENDRO_590;
const double DENDRO_592 = DENDRO_224 * DENDRO_381;
const double DENDRO_593 = DENDRO_172 * DENDRO_281;
const double DENDRO_594 = DENDRO_510 + DENDRO_593;
const double DENDRO_595 = DENDRO_189 * DENDRO_63 + 0.5 * DENDRO_350;
const double DENDRO_596 = DENDRO_107 * DENDRO_172;
const double DENDRO_597 = 0.25 * DENDRO_596;
const double DENDRO_598 = DENDRO_230 * DENDRO_407;
const double DENDRO_599 = DENDRO_206 * DENDRO_499 + DENDRO_497;
const double DENDRO_600 = DENDRO_134 * DENDRO_334;
const double DENDRO_601 = -DENDRO_600;
const double DENDRO_602 = DENDRO_181 * DENDRO_233;
const double DENDRO_603 = -DENDRO_224 * DENDRO_334 + DENDRO_602;
const double DENDRO_604 = DENDRO_122 * DENDRO_196;
const double DENDRO_605 = DENDRO_107 * DENDRO_508 + DENDRO_604;
const double DENDRO_606 = DENDRO_191 * DENDRO_210;
const double DENDRO_607 = DENDRO_206 * DENDRO_381 + DENDRO_606;
const double DENDRO_608 = DENDRO_206 * DENDRO_416;
const double DENDRO_609 = DENDRO_141 * DENDRO_508;
const double DENDRO_610 = DENDRO_104 * DENDRO_189 + DENDRO_177 * DENDRO_499;
const double DENDRO_611 = -DENDRO_172 * DENDRO_412;
const double DENDRO_612 = DENDRO_196 * DENDRO_66;
const double DENDRO_613 = 0.5 * DENDRO_403 + DENDRO_612;
const double DENDRO_614 = DENDRO_479 * grad_1_gt3[pp];
const double DENDRO_615 = DENDRO_437 + DENDRO_614;
const double DENDRO_616 = DENDRO_215 * DENDRO_357;
const double DENDRO_617 = DENDRO_134 * DENDRO_416;
const double DENDRO_618 = 0.25 * DENDRO_107;
const double DENDRO_619 = DENDRO_122 * DENDRO_233;
const double DENDRO_620 = DENDRO_224 * DENDRO_618 + DENDRO_619;
const double DENDRO_621 = DENDRO_107 * DENDRO_134;
const double DENDRO_622 = 1.0 * DENDRO_271;
const double DENDRO_623 = -grad2_1_2_chi[pp];
const double DENDRO_624 =
    DENDRO_141 * DENDRO_325 + DENDRO_174 + DENDRO_322 * grad_1_gt5[pp];
const double DENDRO_625 = DENDRO_328 * grad_2_gt3[pp];
const double DENDRO_626 = DENDRO_330 * grad_1_gt5[pp];
const double DENDRO_627 = DENDRO_141 * DENDRO_322;
const double DENDRO_628 = DENDRO_9 * gt4[pp];
const double DENDRO_629 =
    DENDRO_284 * grad2_0_2_gt4[pp] + DENDRO_290 * grad2_0_0_gt4[pp] +
    DENDRO_292 * grad2_1_1_gt4[pp] + DENDRO_294 * grad2_2_2_gt4[pp] -
    DENDRO_300 * grad2_0_1_gt4[pp] - DENDRO_301 * grad2_1_2_gt4[pp] +
    DENDRO_530 * grad_1_Gt0[pp] + DENDRO_532 * grad_1_Gt1[pp] +
    DENDRO_532 * grad_2_Gt2[pp] + DENDRO_534 * grad_1_Gt2[pp] +
    DENDRO_536 * gt1[pp] + DENDRO_538 * gt3[pp] - DENDRO_541 * grad_1_chi[pp];
const double DENDRO_630 =
    -DENDRO_103 *
        (DENDRO_524 * DENDRO_624 + DENDRO_526 * (DENDRO_205 + DENDRO_625) +
         DENDRO_529 * (DENDRO_221 + DENDRO_626 + DENDRO_627) + DENDRO_623) +
    DENDRO_359 * grad_0_gt4[pp] + DENDRO_361 * grad_1_gt4[pp] +
    DENDRO_363 * grad_2_gt4[pp] + DENDRO_368 * DENDRO_628 + DENDRO_629;
const double DENDRO_631 = DENDRO_203 * DENDRO_9 + DENDRO_204 * DENDRO_9;
const double DENDRO_632 =
    -DENDRO_202 * DENDRO_9 -
    DENDRO_59 * (-DENDRO_57 * DENDRO_628 + grad_2_chi[pp]) + DENDRO_631;
const double DENDRO_633 = DENDRO_221 * DENDRO_9;
const double DENDRO_634 = DENDRO_222 * DENDRO_9;
const double DENDRO_635 = DENDRO_223 * DENDRO_9;
const double DENDRO_636 =
    DENDRO_59 * (-DENDRO_628 * DENDRO_75 + grad_1_chi[pp]);
const double DENDRO_637 = DENDRO_59 * gt4[pp];
const double DENDRO_638 = DENDRO_637 * DENDRO_93;
const double DENDRO_639 =
    DENDRO_555 * DENDRO_9 * (DENDRO_177 + DENDRO_638) +
    DENDRO_560 * (DENDRO_633 - DENDRO_634 - DENDRO_635 - DENDRO_636) +
    DENDRO_561 * DENDRO_632 - 4 * grad2_1_2_alpha[pp];
const double DENDRO_640 = 1.0 * DENDRO_431;
const double DENDRO_641 = 0.5 * DENDRO_430;
const double DENDRO_642 = 0.25 * DENDRO_621;
const double DENDRO_643 = DENDRO_353 + DENDRO_504;
const double DENDRO_644 = DENDRO_125 * DENDRO_196;
const double DENDRO_645 = DENDRO_606 + DENDRO_608;
const double DENDRO_646 = DENDRO_224 * DENDRO_354;
const double DENDRO_647 = DENDRO_437 + DENDRO_438;
const double DENDRO_648 = DENDRO_196 * DENDRO_65;
const double DENDRO_649 = DENDRO_172 * DENDRO_285 + DENDRO_648;
const double DENDRO_650 = DENDRO_215 * DENDRO_407 + DENDRO_614;
const double DENDRO_651 = DENDRO_134 * DENDRO_381 + DENDRO_619;
const double DENDRO_652 = 1.0 * DENDRO_157;
const double DENDRO_653 =
    -DENDRO_158 * (DENDRO_570 + DENDRO_599) -
    DENDRO_262 * (DENDRO_194 * DENDRO_210 + DENDRO_377 + DENDRO_586) -
    DENDRO_622 * (DENDRO_117 + DENDRO_464) -
    DENDRO_652 * (DENDRO_379 + DENDRO_579);
const double DENDRO_654 = DENDRO_502 + DENDRO_601;
const double DENDRO_655 = DENDRO_179 * DENDRO_285;
const double DENDRO_656 = -DENDRO_106 * DENDRO_215;
const double DENDRO_657 = DENDRO_233 * DENDRO_287;
const double DENDRO_658 = DENDRO_129 * DENDRO_381 + DENDRO_657;
const double DENDRO_659 = DENDRO_63 * DENDRO_98;
const double DENDRO_660 = 0.25 * DENDRO_295;
const double DENDRO_661 = DENDRO_101 * DENDRO_407 + DENDRO_129 * DENDRO_154;
const double DENDRO_662 = DENDRO_101 * DENDRO_187;
const double DENDRO_663 = 0.5 * DENDRO_130;
const double DENDRO_664 = -DENDRO_662 - DENDRO_663;
const double DENDRO_665 = DENDRO_509 + DENDRO_593;
const double DENDRO_666 = DENDRO_407 * DENDRO_98 + DENDRO_452;
const double DENDRO_667 = DENDRO_116 * DENDRO_513 + DENDRO_581;
const double DENDRO_668 = DENDRO_125 * DENDRO_233;
const double DENDRO_669 = 0.5 * DENDRO_411;
const double DENDRO_670 = -DENDRO_668 - DENDRO_669;
const double DENDRO_671 = DENDRO_215 * DENDRO_287;
const double DENDRO_672 = DENDRO_179 * DENDRO_618 + DENDRO_648;
const double DENDRO_673 = -DENDRO_116 * DENDRO_412;
const double DENDRO_674 = DENDRO_191 * DENDRO_69;
const double DENDRO_675 = DENDRO_101 * DENDRO_194 + DENDRO_141 * DENDRO_503;
const double DENDRO_676 = DENDRO_148 * DENDRO_179;
const double DENDRO_677 = DENDRO_196 * DENDRO_67;
const double DENDRO_678 = DENDRO_179 * DENDRO_281 + DENDRO_677;
const double DENDRO_679 = DENDRO_113 * DENDRO_503;
const double DENDRO_680 = DENDRO_233 * DENDRO_65;
const double DENDRO_681 = DENDRO_107 * DENDRO_503 + DENDRO_680;
const double DENDRO_682 = 1.0 * DENDRO_261;
const double DENDRO_683 = -grad2_0_1_chi[pp];
const double DENDRO_684 = DENDRO_325 * grad_1_gt0[pp];
const double DENDRO_685 = DENDRO_113 * DENDRO_322;
const double DENDRO_686 = DENDRO_328 * grad_0_gt3[pp];
const double DENDRO_687 =
    DENDRO_113 * DENDRO_330 + DENDRO_126 + DENDRO_322 * grad_1_gt0[pp];
const double DENDRO_688 = 2.0 * gt1[pp];
const double DENDRO_689 = DENDRO_688 * grad_0_Gt0[pp];
const double DENDRO_690 = 2.0 * grad_0_Gt1[pp] * gt3[pp];
const double DENDRO_691 = DENDRO_532 * grad_0_Gt2[pp];
const double DENDRO_692 = 2.0 * grad_1_Gt0[pp] * gt0[pp];
const double DENDRO_693 = DENDRO_688 * grad_1_Gt1[pp];
const double DENDRO_694 = DENDRO_530 * grad_1_Gt2[pp];
const double DENDRO_695 = -DENDRO_242 * DENDRO_278;
const double DENDRO_696 = DENDRO_284 * grad2_0_2_gt1[pp];
const double DENDRO_697 = DENDRO_290 * grad2_0_0_gt1[pp];
const double DENDRO_698 = DENDRO_292 * grad2_1_1_gt1[pp];
const double DENDRO_699 = DENDRO_294 * grad2_2_2_gt1[pp];
const double DENDRO_700 = -DENDRO_300 * grad2_0_1_gt1[pp];
const double DENDRO_701 = -DENDRO_301 * grad2_1_2_gt1[pp];
const double DENDRO_702 = DENDRO_9 * gt1[pp];
const double DENDRO_703 =
    -DENDRO_103 * (DENDRO_524 * (DENDRO_150 + DENDRO_684 + DENDRO_685) +
                   DENDRO_526 * (DENDRO_115 + DENDRO_686) +
                   DENDRO_529 * DENDRO_687 + DENDRO_683) +
    DENDRO_359 * grad_0_gt1[pp] + DENDRO_361 * grad_1_gt1[pp] +
    DENDRO_363 * grad_2_gt1[pp] + DENDRO_368 * DENDRO_702 + DENDRO_689 +
    DENDRO_690 + DENDRO_691 + DENDRO_692 + DENDRO_693 + DENDRO_694 +
    DENDRO_695 + DENDRO_696 + DENDRO_697 + DENDRO_698 + DENDRO_699 +
    DENDRO_700 + DENDRO_701;
const double DENDRO_704 = DENDRO_150 * DENDRO_9;
const double DENDRO_705 = DENDRO_151 * DENDRO_9;
const double DENDRO_706 = DENDRO_152 * DENDRO_9;
const double DENDRO_707 =
    DENDRO_59 * (-DENDRO_702 * DENDRO_93 + grad_1_chi[pp]);
const double DENDRO_708 =
    DENDRO_59 * (-DENDRO_57 * DENDRO_702 + grad_0_chi[pp]);
const double DENDRO_709 = DENDRO_111 * DENDRO_9;
const double DENDRO_710 = DENDRO_112 * DENDRO_9;
const double DENDRO_711 = DENDRO_114 * DENDRO_9;
const double DENDRO_712 = DENDRO_710 + DENDRO_711;
const double DENDRO_713 = DENDRO_59 * gt1[pp];
const double DENDRO_714 =
    DENDRO_555 * (DENDRO_704 - DENDRO_705 - DENDRO_706 - DENDRO_707) +
    DENDRO_560 * DENDRO_9 * (DENDRO_129 + DENDRO_713 * DENDRO_75) +
    DENDRO_561 * (-DENDRO_708 - DENDRO_709 + DENDRO_712) -
    4 * grad2_0_1_alpha[pp];
const double DENDRO_715 = 0.5 * DENDRO_426;
const double DENDRO_716 = DENDRO_196 * DENDRO_62;
const double DENDRO_717 = DENDRO_438 + DENDRO_614;
const double DENDRO_718 = DENDRO_116 * DENDRO_285 + DENDRO_674;
const double DENDRO_719 = DENDRO_121 * (DENDRO_474 + DENDRO_717) -
                          DENDRO_158 * (DENDRO_286 + DENDRO_582) -
                          DENDRO_158 * (DENDRO_286 + DENDRO_667) +
                          DENDRO_164 * (DENDRO_673 + DENDRO_718) -
                          DENDRO_262 * (DENDRO_109 * DENDRO_194 + DENDRO_496) -
                          DENDRO_272 * (1.0 * DENDRO_296 + DENDRO_660);
const double DENDRO_720 =
    -DENDRO_10 *
        (DENDRO_714 +
         alpha[pp] *
             (DENDRO_120 * (DENDRO_404 + DENDRO_576 + DENDRO_584) +
              DENDRO_120 * (DENDRO_517 + DENDRO_518 + DENDRO_585) +
              DENDRO_121 * (DENDRO_611 + DENDRO_672) +
              DENDRO_121 * (DENDRO_617 + DENDRO_670) +
              DENDRO_121 * (DENDRO_671 + DENDRO_717) -
              DENDRO_158 * (DENDRO_136 + DENDRO_664) -
              DENDRO_158 * (DENDRO_451 + DENDRO_666) -
              DENDRO_158 * (DENDRO_573 + DENDRO_665) -
              DENDRO_158 * (DENDRO_516 + DENDRO_574 + DENDRO_642) +
              DENDRO_164 * (DENDRO_675 + DENDRO_679) +
              DENDRO_164 * (DENDRO_101 * DENDRO_195 + DENDRO_681) +
              DENDRO_164 * (-DENDRO_106 * DENDRO_98 + DENDRO_678) +
              DENDRO_164 * (DENDRO_215 * DENDRO_62 + DENDRO_718) +
              DENDRO_238 * (DENDRO_459 + DENDRO_98 * grad_0_gt3[pp]) -
              DENDRO_262 * (DENDRO_355 + DENDRO_654) -
              DENDRO_268 * (DENDRO_422 + DENDRO_658) -
              DENDRO_268 * (DENDRO_656 + DENDRO_715) -
              DENDRO_268 * (DENDRO_424 + DENDRO_655 + DENDRO_716) -
              DENDRO_272 * (0.5 * DENDRO_457 + DENDRO_661) -
              DENDRO_272 * (DENDRO_455 + DENDRO_62 * DENDRO_98 + DENDRO_659) -
              DENDRO_682 * (DENDRO_351 + DENDRO_565 + DENDRO_596) + DENDRO_703 +
              DENDRO_719)) -
    DENDRO_10 *
        (DENDRO_714 +
         alpha[pp] *
             (DENDRO_121 * (DENDRO_406 + DENDRO_649) +
              DENDRO_121 * (DENDRO_406 + DENDRO_672) +
              DENDRO_121 * (DENDRO_410 + DENDRO_670) +
              DENDRO_121 * (DENDRO_474 + DENDRO_650) +
              DENDRO_121 * (DENDRO_567 + DENDRO_651) +
              DENDRO_121 * (DENDRO_647 + DENDRO_671) -
              DENDRO_158 * (DENDRO_452 + DENDRO_665) -
              DENDRO_158 * (DENDRO_506 + DENDRO_667) -
              DENDRO_158 * (DENDRO_509 + DENDRO_666) -
              DENDRO_158 * (DENDRO_515 + DENDRO_664) +
              DENDRO_164 * (DENDRO_676 + DENDRO_678) +
              DENDRO_164 * (DENDRO_679 + DENDRO_681) +
              DENDRO_164 * (DENDRO_233 * DENDRO_66 + DENDRO_675) +
              DENDRO_164 * (DENDRO_104 * DENDRO_98 + DENDRO_676 + DENDRO_677) +
              DENDRO_164 * (DENDRO_215 * DENDRO_63 + DENDRO_673 + DENDRO_674) +
              DENDRO_238 * (DENDRO_215 * grad_1_gt0[pp] + DENDRO_439) -
              DENDRO_262 * (DENDRO_353 + DENDRO_654) -
              DENDRO_262 * (DENDRO_177 * DENDRO_65 + DENDRO_597) -
              DENDRO_268 * (1.0 * DENDRO_429 + DENDRO_655) -
              DENDRO_268 * (0.5 * DENDRO_432 + DENDRO_658) -
              DENDRO_268 * (DENDRO_104 * DENDRO_215 + DENDRO_427 + DENDRO_656) -
              DENDRO_272 * (DENDRO_453 + DENDRO_661) -
              DENDRO_272 * (0.5 * DENDRO_454 + DENDRO_659) -
              DENDRO_272 * (DENDRO_104 * DENDRO_69 + DENDRO_288 + DENDRO_660) -
              DENDRO_652 * (DENDRO_117 + DENDRO_461 + DENDRO_462) -
              DENDRO_652 * (DENDRO_465 + DENDRO_466 + DENDRO_621) -
              DENDRO_682 * (DENDRO_379 + DENDRO_495 + DENDRO_578) +
              DENDRO_703)) +
    DENDRO_14 *
        (DENDRO_564 +
         alpha[pp] *
             (DENDRO_120 * (DENDRO_351 + DENDRO_566) +
              DENDRO_121 * (DENDRO_490 + DENDRO_570) +
              DENDRO_121 * (DENDRO_493 + DENDRO_571) +
              DENDRO_121 * (DENDRO_500 + DENDRO_572) +
              DENDRO_121 * (DENDRO_502 + DENDRO_571) -
              DENDRO_158 * (DENDRO_481 + DENDRO_568) -
              DENDRO_158 * (-DENDRO_125 * DENDRO_98 + DENDRO_485) -
              DENDRO_158 * (DENDRO_230 * DENDRO_65 + DENDRO_568) +
              DENDRO_164 * (DENDRO_452 + DENDRO_512) +
              DENDRO_164 * (DENDRO_511 + DENDRO_573) +
              DENDRO_164 * (DENDRO_514 + DENDRO_575) +
              DENDRO_164 * (DENDRO_574 + DENDRO_575) -
              DENDRO_169 * (DENDRO_460 + DENDRO_98 * grad_0_gt5[pp]) -
              DENDRO_262 * (0.5 * DENDRO_341 + DENDRO_469) -
              DENDRO_262 * (DENDRO_189 * DENDRO_65 + DENDRO_339 + DENDRO_468) -
              DENDRO_268 * (DENDRO_129 * DENDRO_184 + DENDRO_567) -
              DENDRO_272 * (1.0 * DENDRO_458 + DENDRO_478) -
              DENDRO_272 * (DENDRO_448 + DENDRO_476 + DENDRO_65 * DENDRO_98) -
              DENDRO_520 * (DENDRO_404 + DENDRO_577) + DENDRO_550 +
              DENDRO_583)) +
    DENDRO_14 *
        (DENDRO_564 +
         alpha[pp] *
             (DENDRO_121 * (DENDRO_356 + DENDRO_492) +
              DENDRO_121 * (DENDRO_356 + DENDRO_500) +
              DENDRO_121 * (DENDRO_380 + DENDRO_490) +
              DENDRO_121 * (DENDRO_493 + DENDRO_494) +
              DENDRO_121 * (DENDRO_496 + DENDRO_498) +
              DENDRO_121 * (DENDRO_502 + DENDRO_505) -
              DENDRO_158 * (DENDRO_483 + DENDRO_485) -
              DENDRO_158 * (DENDRO_486 + DENDRO_487) -
              DENDRO_158 * (DENDRO_210 * DENDRO_63 + DENDRO_480) -
              DENDRO_158 * (DENDRO_122 * DENDRO_98 + DENDRO_483 + DENDRO_484) -
              DENDRO_158 * (DENDRO_230 * DENDRO_66 + DENDRO_481 + DENDRO_482) +
              DENDRO_163 * (DENDRO_461 + DENDRO_464) +
              DENDRO_163 * (DENDRO_465 + DENDRO_467) +
              DENDRO_164 * (DENDRO_451 + DENDRO_511) +
              DENDRO_164 * (DENDRO_506 + DENDRO_507) +
              DENDRO_164 * (DENDRO_509 + DENDRO_512) +
              DENDRO_164 * (DENDRO_514 + DENDRO_515 + DENDRO_516) -
              DENDRO_169 * (DENDRO_230 * grad_2_gt0[pp] + DENDRO_372) -
              DENDRO_262 * (1.0 * DENDRO_366 + DENDRO_468) -
              DENDRO_262 * (0.5 * DENDRO_386 + DENDRO_470) -
              DENDRO_262 * (DENDRO_122 * DENDRO_230 + DENDRO_342 + DENDRO_469) -
              DENDRO_268 * (DENDRO_437 + DENDRO_475) -
              DENDRO_268 * (DENDRO_177 * DENDRO_62 + DENDRO_472) -
              DENDRO_272 * (DENDRO_282 + DENDRO_477) -
              DENDRO_272 * (0.5 * DENDRO_447 + DENDRO_476) -
              DENDRO_272 * (DENDRO_101 * DENDRO_122 + DENDRO_446 + DENDRO_478) -
              DENDRO_520 * (DENDRO_517 + DENDRO_519) + DENDRO_550)) -
    DENDRO_16 *
        (DENDRO_639 +
         alpha[pp] *
             (DENDRO_121 * (DENDRO_603 + DENDRO_646) +
              DENDRO_121 * (DENDRO_605 + DENDRO_609) +
              DENDRO_121 * (DENDRO_610 - DENDRO_644) +
              DENDRO_121 * (-DENDRO_187 * DENDRO_215 + DENDRO_645) +
              DENDRO_121 * (DENDRO_194 * DENDRO_230 + DENDRO_602 + DENDRO_646) -
              DENDRO_158 * (DENDRO_572 + DENDRO_595) -
              DENDRO_158 * (DENDRO_598 + DENDRO_643) -
              DENDRO_158 * (DENDRO_601 + DENDRO_643) +
              DENDRO_164 * (DENDRO_405 + DENDRO_613) +
              DENDRO_164 * (DENDRO_408 + DENDRO_620) +
              DENDRO_164 * (DENDRO_408 + DENDRO_651) +
              DENDRO_164 * (DENDRO_472 + DENDRO_649) +
              DENDRO_164 * (DENDRO_473 + DENDRO_650) +
              DENDRO_164 * (DENDRO_616 + DENDRO_647) -
              DENDRO_262 * (DENDRO_345 + DENDRO_588) -
              DENDRO_262 * (0.5 * DENDRO_347 + DENDRO_587) -
              DENDRO_268 * (DENDRO_591 + DENDRO_641) -
              DENDRO_268 * (DENDRO_592 + DENDRO_640) -
              DENDRO_268 * (DENDRO_194 * DENDRO_215 + DENDRO_420 + DENDRO_589) -
              DENDRO_272 * (DENDRO_452 + DENDRO_594) -
              DENDRO_272 * (DENDRO_122 * DENDRO_129 + DENDRO_642) +
              DENDRO_365 * (DENDRO_215 * grad_1_gt5[pp] + DENDRO_440) +
              DENDRO_630 - DENDRO_652 * (DENDRO_566 + DENDRO_596) +
              DENDRO_653)) -
    DENDRO_16 *
        (DENDRO_639 +
         alpha[pp] *
             (DENDRO_121 * (DENDRO_607 + DENDRO_608) +
              DENDRO_121 * (DENDRO_609 + DENDRO_610) +
              DENDRO_121 * (-DENDRO_106 * DENDRO_189 + DENDRO_605) +
              DENDRO_121 * (DENDRO_184 * DENDRO_215 + DENDRO_607) +
              DENDRO_121 * (DENDRO_195 * DENDRO_230 + DENDRO_603) -
              DENDRO_158 * (DENDRO_352 + DENDRO_595) -
              DENDRO_158 * (DENDRO_382 + DENDRO_498) -
              DENDRO_158 * (DENDRO_382 + DENDRO_599) -
              DENDRO_158 * (DENDRO_492 + DENDRO_597) -
              DENDRO_158 * (DENDRO_494 + DENDRO_598) -
              DENDRO_158 * (DENDRO_505 + DENDRO_601) +
              DENDRO_163 * (DENDRO_519 + DENDRO_585) +
              DENDRO_163 * (DENDRO_577 + DENDRO_584) +
              DENDRO_164 * (DENDRO_473 + DENDRO_615) +
              DENDRO_164 * (DENDRO_611 + DENDRO_613) +
              DENDRO_164 * (DENDRO_615 + DENDRO_616) +
              DENDRO_164 * (DENDRO_617 + DENDRO_620) -
              DENDRO_262 * (0.5 * DENDRO_364 + DENDRO_588) -
              DENDRO_262 * (1.0 * DENDRO_384 + DENDRO_586) -
              DENDRO_262 * (DENDRO_184 * DENDRO_230 + DENDRO_348 + DENDRO_587) -
              DENDRO_268 * (DENDRO_414 + DENDRO_591) -
              DENDRO_268 * (0.5 * DENDRO_419 + DENDRO_589) -
              DENDRO_268 * (DENDRO_184 * DENDRO_233 + DENDRO_417 + DENDRO_592) -
              DENDRO_272 * (DENDRO_451 + DENDRO_594) -
              DENDRO_272 * (DENDRO_104 * DENDRO_109 + DENDRO_580) +
              DENDRO_365 * (DENDRO_230 * grad_2_gt3[pp] + DENDRO_374) -
              DENDRO_622 * (DENDRO_467 + DENDRO_621) + DENDRO_630)) +
    DENDRO_19 *
        (DENDRO_310 * DENDRO_70 +
         DENDRO_312 * (-DENDRO_85 - DENDRO_87 + DENDRO_89 - DENDRO_95) +
         DENDRO_390 * (DENDRO_101 + DENDRO_77) +
         alpha[pp] *
             (-DENDRO_103 *
                  (DENDRO_100 * DENDRO_444 + DENDRO_102 * DENDRO_445 +
                   DENDRO_97 +
                   DENDRO_99 * (DENDRO_442 + DENDRO_443 + DENDRO_88)) -
              DENDRO_119 * DENDRO_198 * DENDRO_277 +
              DENDRO_121 * (-1.0 * DENDRO_110 + DENDRO_118) +
              DENDRO_121 * (-1.0 * DENDRO_130 + DENDRO_136) +
              DENDRO_121 * (DENDRO_134 * DENDRO_287 + DENDRO_450) +
              DENDRO_121 * (DENDRO_172 * DENDRO_62 + DENDRO_452) +
              DENDRO_121 * (DENDRO_179 * DENDRO_65 + DENDRO_451) -
              DENDRO_129 * DENDRO_269 -
              DENDRO_134 * DENDRO_262 * (DENDRO_123 - DENDRO_139) -
              DENDRO_158 * (DENDRO_448 + 1.0 * DENDRO_449) -
              DENDRO_158 * (-DENDRO_101 * DENDRO_156 + DENDRO_446) +
              DENDRO_164 * (DENDRO_455 + 1.0 * DENDRO_456) +
              DENDRO_164 * (DENDRO_101 * DENDRO_141 + DENDRO_453) -
              DENDRO_169 * (DENDRO_447 + DENDRO_449) -
              DENDRO_169 * (DENDRO_134 * grad_2_gt0[pp] + DENDRO_458) +
              DENDRO_200 * DENDRO_358 + DENDRO_219 * DENDRO_360 -
              DENDRO_235 * DENDRO_270 * DENDRO_304 + DENDRO_237 * DENDRO_362 +
              DENDRO_238 * (DENDRO_454 + DENDRO_456) +
              DENDRO_238 * (DENDRO_129 * grad_2_gt0[pp] + DENDRO_457) +
              DENDRO_260 * DENDRO_367 - DENDRO_274 * DENDRO_459 -
              DENDRO_276 * DENDRO_460 + DENDRO_305) -
         4 * grad2_0_0_alpha[pp]) +
    DENDRO_22 *
        (DENDRO_309 * (-DENDRO_391 + DENDRO_392 - DENDRO_393 - DENDRO_396) +
         DENDRO_313 * (DENDRO_196 + DENDRO_389 * DENDRO_94) +
         DENDRO_390 * (DENDRO_233 + DENDRO_389 * DENDRO_76) +
         alpha[pp] *
             (-DENDRO_103 *
                  (DENDRO_100 * (DENDRO_213 + DENDRO_400 + DENDRO_401) +
                   DENDRO_102 * DENDRO_402 + DENDRO_397 +
                   DENDRO_399 * DENDRO_99) +
              DENDRO_121 * (DENDRO_414 + DENDRO_415) +
              DENDRO_121 * (DENDRO_417 - DENDRO_418) +
              DENDRO_121 * (DENDRO_420 + 1.0 * DENDRO_421) -
              DENDRO_129 * DENDRO_436 - DENDRO_158 * (DENDRO_403 + DENDRO_405) -
              DENDRO_158 * (DENDRO_410 - 1.0 * DENDRO_411) -
              DENDRO_158 * (DENDRO_179 * DENDRO_407 + DENDRO_406) -
              DENDRO_158 * (DENDRO_224 * DENDRO_407 + DENDRO_408) +
              DENDRO_164 * (DENDRO_422 + DENDRO_423) +
              DENDRO_164 * (DENDRO_424 + DENDRO_425) +
              DENDRO_164 * (DENDRO_427 + 1.0 * DENDRO_428) -
              DENDRO_177 * DENDRO_434 - DENDRO_179 * DENDRO_435 -
              DENDRO_216 * DENDRO_375 * grad_1_gt3[pp] -
              DENDRO_224 * DENDRO_262 * (DENDRO_185 - DENDRO_354) +
              DENDRO_238 * (DENDRO_426 + DENDRO_428) +
              DENDRO_238 * (DENDRO_129 * grad_2_gt3[pp] + DENDRO_432) +
              DENDRO_238 * (DENDRO_179 * grad_0_gt3[pp] + DENDRO_429) +
              DENDRO_359 * grad_0_gt3[pp] + DENDRO_361 * grad_1_gt3[pp] +
              DENDRO_363 * grad_2_gt3[pp] +
              DENDRO_365 * (DENDRO_419 + DENDRO_421) +
              DENDRO_365 * (DENDRO_177 * grad_0_gt3[pp] + DENDRO_430) +
              DENDRO_365 * (DENDRO_224 * grad_2_gt3[pp] + DENDRO_431) +
              DENDRO_368 * DENDRO_395 + DENDRO_441) -
         4 * grad2_1_1_alpha[pp]) +
    DENDRO_25 *
        (DENDRO_308 * DENDRO_310 + DENDRO_313 * (DENDRO_189 + DENDRO_311) +
         DENDRO_320 * (-DENDRO_314 + DENDRO_315 - DENDRO_316 - DENDRO_319) +
         alpha[pp] *
             (-DENDRO_103 *
                  (DENDRO_100 * DENDRO_329 +
                   DENDRO_102 * (DENDRO_331 + DENDRO_332 + DENDRO_333) +
                   DENDRO_321 + DENDRO_327 * DENDRO_99) +
              DENDRO_121 * (DENDRO_345 + DENDRO_346) +
              DENDRO_121 * (DENDRO_348 + 1.0 * DENDRO_349) -
              DENDRO_158 * (DENDRO_336 + DENDRO_337) -
              DENDRO_158 * (DENDRO_342 + 1.0 * DENDRO_343) -
              DENDRO_158 * (DENDRO_189 * DENDRO_340 + DENDRO_339) +
              DENDRO_164 * (DENDRO_350 + DENDRO_352) +
              DENDRO_164 * (DENDRO_122 * DENDRO_224 + DENDRO_355) +
              DENDRO_164 * (DENDRO_134 * DENDRO_184 + DENDRO_353) +
              DENDRO_164 * (DENDRO_172 * DENDRO_357 + DENDRO_356) -
              DENDRO_169 * (DENDRO_341 + DENDRO_343) -
              DENDRO_169 * (DENDRO_172 * grad_0_gt5[pp] + DENDRO_366) -
              DENDRO_172 * DENDRO_371 - DENDRO_177 * DENDRO_370 -
              DENDRO_190 * DENDRO_304 * (DENDRO_122 + DENDRO_123) -
              DENDRO_231 * DENDRO_375 * grad_2_gt5[pp] -
              DENDRO_274 * DENDRO_374 + DENDRO_318 * DENDRO_368 +
              DENDRO_359 * grad_0_gt5[pp] + DENDRO_361 * grad_1_gt5[pp] +
              DENDRO_363 * grad_2_gt5[pp] +
              DENDRO_365 * (DENDRO_347 + DENDRO_349) +
              DENDRO_365 * (DENDRO_177 * grad_0_gt5[pp] + DENDRO_364) -
              DENDRO_372 * DENDRO_373 + DENDRO_388) -
         4 * grad2_2_2_alpha[pp]);
const double DENDRO_721 = DENDRO_720 * DENDRO_9;
const double DENDRO_722 = (1.0 / 12.0) * chi[pp];
const double DENDRO_723 = (1.0 / 3.0) * At1[pp];
const double DENDRO_724 = At4[pp] * DENDRO_16;
const double DENDRO_725 = At3[pp] * DENDRO_22;
const double DENDRO_726 = DENDRO_47 + DENDRO_724 - DENDRO_725;
const double DENDRO_727 = At3[pp] * DENDRO_10;
const double DENDRO_728 = At4[pp] * DENDRO_14;
const double DENDRO_729 = -At1[pp] * DENDRO_19 + DENDRO_727 - DENDRO_728;
const double DENDRO_730 = At4[pp] * DENDRO_25;
const double DENDRO_731 =
    -At1[pp] * DENDRO_14 + At3[pp] * DENDRO_16 - DENDRO_730;
const double DENDRO_732 = -DENDRO_17 + DENDRO_73 + DENDRO_74;
const double DENDRO_733 = 6.0 * grad_2_alpha[pp];
const double DENDRO_734 = 6.0 * grad_0_alpha[pp];
const double DENDRO_735 = 6.0 * grad_1_alpha[pp];
const double DENDRO_736 = DENDRO_107 * DENDRO_147;
const double DENDRO_737 = -DENDRO_174 + DENDRO_175 + DENDRO_176;
const double DENDRO_738 = DENDRO_737 * grad_2_gt0[pp];
const double DENDRO_739 = DENDRO_141 * DENDRO_147;
const double DENDRO_740 = DENDRO_738 + DENDRO_739;
const double DENDRO_741 = DENDRO_212 - DENDRO_213 + DENDRO_214;
const double DENDRO_742 = DENDRO_153 * DENDRO_412;
const double DENDRO_743 = -DENDRO_138 * DENDRO_416;
const double DENDRO_744 = DENDRO_147 * DENDRO_281;
const double DENDRO_745 =
    DENDRO_161 * DENDRO_287 + 0.25 * DENDRO_737 * grad_0_gt0[pp];
const double DENDRO_746 = DENDRO_149 + DENDRO_155;
const double DENDRO_747 = DENDRO_142 * DENDRO_618;
const double DENDRO_748 = -DENDRO_221 + DENDRO_222 + DENDRO_223;
const double DENDRO_749 = DENDRO_184 * DENDRO_81;
const double DENDRO_750 = DENDRO_154 * DENDRO_748 + DENDRO_749;
const double DENDRO_751 = DENDRO_141 * DENDRO_153;
const double DENDRO_752 = DENDRO_737 * grad_1_gt0[pp] + DENDRO_751;
const double DENDRO_753 = 1.0 * DENDRO_120;
const double DENDRO_754 = DENDRO_141 * DENDRO_748;
const double DENDRO_755 = 2.0 * DENDRO_199;
const double DENDRO_756 = 2.0 * DENDRO_218;
const double DENDRO_757 = 2.0 * DENDRO_236;
const double DENDRO_758 = DENDRO_259 * DENDRO_59;
const double DENDRO_759 = (1.0 / 3.0) * At2[pp];
const double DENDRO_760 = At5[pp] * DENDRO_16;
const double DENDRO_761 =
    At2[pp] * DENDRO_10 - At4[pp] * DENDRO_22 + DENDRO_760;
const double DENDRO_762 = At5[pp] * DENDRO_14;
const double DENDRO_763 =
    -At2[pp] * DENDRO_19 + At4[pp] * DENDRO_10 - DENDRO_762;
const double DENDRO_764 = -At5[pp] * DENDRO_25 + DENDRO_49 + DENDRO_724;
const double DENDRO_765 = DENDRO_113 * DENDRO_153;
const double DENDRO_766 = DENDRO_142 * grad_2_gt5[pp];
const double DENDRO_767 = DENDRO_227 - DENDRO_228 + DENDRO_229;
const double DENDRO_768 = DENDRO_147 * DENDRO_338;
const double DENDRO_769 = DENDRO_182 - DENDRO_183 + DENDRO_188;
const double DENDRO_770 = DENDRO_142 * DENDRO_354;
const double DENDRO_771 = 0.25 * DENDRO_138 * grad_2_gt5[pp];
const double DENDRO_772 = DENDRO_139 * DENDRO_142;
const double DENDRO_773 = DENDRO_181 * DENDRO_81;
const double DENDRO_774 = DENDRO_153 * DENDRO_165;
const double DENDRO_775 = DENDRO_153 * DENDRO_338;
const double DENDRO_776 = -DENDRO_771;
const double DENDRO_777 = DENDRO_153 * grad_0_gt5[pp];
const double DENDRO_778 = 2 * At4[pp];
const double DENDRO_779 = At3[pp] * DENDRO_51;
const double DENDRO_780 = DENDRO_778 * DENDRO_9;
const double DENDRO_781 = 0.5 * DENDRO_389;
const double DENDRO_782 = DENDRO_9 * DENDRO_96;
const double DENDRO_783 = DENDRO_741 * grad_2_gt3[pp];
const double DENDRO_784 = DENDRO_138 * DENDRO_354;
const double DENDRO_785 = 0.25 * DENDRO_751;
const double DENDRO_786 = 1.0 * DENDRO_737;
const double DENDRO_787 = DENDRO_741 * grad_0_gt3[pp];
const double DENDRO_788 = (1.0 / 3.0) * At4[pp];
const double DENDRO_789 = DENDRO_748 * grad_2_gt5[pp];
const double DENDRO_790 = DENDRO_139 * DENDRO_748;
const double DENDRO_791 = -DENDRO_354 * DENDRO_748 + DENDRO_602;
const double DENDRO_792 = DENDRO_619 - DENDRO_784;
const double DENDRO_793 = DENDRO_767 * grad_1_gt5[pp];
const double DENDRO_794 = DENDRO_767 * grad_0_gt5[pp];
const double DENDRO_795 =
    DENDRO_11 + DENDRO_322 * grad_2_chi[pp] + DENDRO_325 * grad_0_chi[pp];
const double DENDRO_796 = 0.5 * DENDRO_307;
const double DENDRO_797 = DENDRO_9 * grad_0_alpha[pp];
const double DENDRO_798 = DENDRO_328 * grad_1_chi[pp] + DENDRO_56;
const double DENDRO_799 = DENDRO_9 * grad_1_alpha[pp];
const double DENDRO_800 =
    DENDRO_17 + DENDRO_322 * grad_0_chi[pp] + DENDRO_330 * grad_2_chi[pp];
const double DENDRO_801 = 0.5 * DENDRO_800;
const double DENDRO_802 = DENDRO_9 * grad_2_alpha[pp];
const double DENDRO_803 = 0.5 * DENDRO_798;
const double DENDRO_804 = 0.5 * grad_1_alpha[pp];
const double DENDRO_805 = 0.5 * grad_2_alpha[pp];
const double DENDRO_806 = DENDRO_248 * DENDRO_9;
const double DENDRO_807 = 0.5 * grad_0_alpha[pp];
const double DENDRO_808 = DENDRO_246 * DENDRO_9;
const double DENDRO_809 = DENDRO_243 * DENDRO_9;
const double DENDRO_810 = (DENDRO_10 * DENDRO_10);
const double DENDRO_811 = (DENDRO_14 * DENDRO_14);
const double DENDRO_812 = 2 * DENDRO_19;
const double DENDRO_813 = At0[pp] * (DENDRO_19 * DENDRO_19) +
                          At3[pp] * DENDRO_810 + At5[pp] * DENDRO_811 -
                          DENDRO_243 * DENDRO_728 - DENDRO_47 * DENDRO_812 +
                          DENDRO_48 * DENDRO_812;
const double DENDRO_814 = 3 * DENDRO_119;
const double DENDRO_815 = (DENDRO_16 * DENDRO_16);
const double DENDRO_816 = 2 * DENDRO_22;
const double DENDRO_817 = At0[pp] * DENDRO_810 +
                          At3[pp] * (DENDRO_22 * DENDRO_22) +
                          At5[pp] * DENDRO_815 + DENDRO_243 * DENDRO_44 -
                          DENDRO_47 * DENDRO_816 - DENDRO_724 * DENDRO_816;
const double DENDRO_818 = 2 * DENDRO_25;
const double DENDRO_819 = At0[pp] * DENDRO_811 + At3[pp] * DENDRO_815 +
                          At5[pp] * (DENDRO_25 * DENDRO_25) -
                          DENDRO_246 * DENDRO_53 + DENDRO_48 * DENDRO_818 -
                          DENDRO_724 * DENDRO_818;
const double DENDRO_820 =
    At2[pp] * DENDRO_811 - DENDRO_10 * DENDRO_730 - DENDRO_14 * DENDRO_47 +
    DENDRO_14 * DENDRO_50 + DENDRO_16 * DENDRO_727 - DENDRO_16 * DENDRO_728 -
    DENDRO_19 * DENDRO_53 + DENDRO_19 * DENDRO_54 + DENDRO_25 * DENDRO_762;
const double DENDRO_821 = 6 * DENDRO_119;
const double DENDRO_822 = At1[pp] * DENDRO_810;
const double DENDRO_823 = DENDRO_10 * DENDRO_724;
const double DENDRO_824 = DENDRO_10 * DENDRO_48;
const double DENDRO_825 = DENDRO_22 * DENDRO_728;
const double DENDRO_826 = DENDRO_16 * DENDRO_762;
const double DENDRO_827 = DENDRO_10 * DENDRO_50;
const double DENDRO_828 = DENDRO_19 * DENDRO_45;
const double DENDRO_829 = DENDRO_19 * DENDRO_44;
const double DENDRO_830 = DENDRO_22 * DENDRO_727;
const double DENDRO_831 = DENDRO_822 + DENDRO_823 - DENDRO_824 + DENDRO_825 -
                          DENDRO_826 - DENDRO_827 + DENDRO_828 - DENDRO_829 -
                          DENDRO_830;
const double DENDRO_832 = At4[pp] * DENDRO_815;
const double DENDRO_833 = DENDRO_16 * DENDRO_47;
const double DENDRO_834 = DENDRO_14 * DENDRO_43;
const double DENDRO_835 = DENDRO_14 * DENDRO_45;
const double DENDRO_836 = DENDRO_16 * DENDRO_48;
const double DENDRO_837 = DENDRO_10 * DENDRO_54;
const double DENDRO_838 = DENDRO_16 * DENDRO_725;
const double DENDRO_839 = DENDRO_22 * DENDRO_730;
const double DENDRO_840 = DENDRO_25 * DENDRO_760;
const double DENDRO_841 = DENDRO_832 + DENDRO_833 - DENDRO_834 + DENDRO_835 -
                          DENDRO_836 - DENDRO_837 - DENDRO_838 + DENDRO_839 -
                          DENDRO_840;
const double DENDRO_842 = (1.0 / 3.0) * alpha[pp];
const double DENDRO_843 = DENDRO_10 * DENDRO_9;
const double DENDRO_844 = (7.0 / 3.0) * DENDRO_843;
const double DENDRO_845 = (7.0 / 3.0) * DENDRO_283;
const double DENDRO_846 = (1.0 / 3.0) * DENDRO_843;
const double DENDRO_847 = (1.0 / 3.0) * DENDRO_283;
const double DENDRO_848 = (1.0 / 3.0) * DENDRO_289;
const double DENDRO_849 = 2 * DENDRO_119;
const double DENDRO_850 = DENDRO_849 * grad_0_alpha[pp];
const double DENDRO_851 = pow(DENDRO_8, -3);
const double DENDRO_852 = DENDRO_0 * DENDRO_851;
const double DENDRO_853 = DENDRO_813 * DENDRO_852;
const double DENDRO_854 = DENDRO_817 * DENDRO_852;
const double DENDRO_855 = DENDRO_819 * DENDRO_852;
const double DENDRO_856 = 4 * grad_0_K[pp];
const double DENDRO_857 = 9 * DENDRO_59;
const double DENDRO_858 = DENDRO_857 * DENDRO_99;
const double DENDRO_859 = DENDRO_842 * DENDRO_9;
const double DENDRO_860 = DENDRO_849 * grad_2_alpha[pp];
const double DENDRO_861 = DENDRO_849 * grad_1_alpha[pp];
const double DENDRO_862 = 2.0 * DENDRO_851 * alpha[pp];
const double DENDRO_863 = DENDRO_820 * DENDRO_862;
const double DENDRO_864 = DENDRO_831 * DENDRO_862;
const double DENDRO_865 = DENDRO_841 * DENDRO_862;
const double DENDRO_866 = 4 * grad_2_K[pp];
const double DENDRO_867 = DENDRO_102 * DENDRO_857;
const double DENDRO_868 = 4 * grad_1_K[pp];
const double DENDRO_869 = DENDRO_100 * DENDRO_857;
const double DENDRO_870 = (2.0 / 3.0) * DENDRO_40;
const double DENDRO_871 = DENDRO_16 * DENDRO_9;
const double DENDRO_872 = (1.0 / 3.0) * grad2_0_2_beta0[pp];
const double DENDRO_873 = (7.0 / 3.0) * DENDRO_871;
const double DENDRO_874 = (1.0 / 3.0) * DENDRO_871;
const double DENDRO_875 = (1.0 / 3.0) * DENDRO_291;
const double DENDRO_876 = -DENDRO_822 - DENDRO_823 + DENDRO_824 - DENDRO_825 +
                          DENDRO_826 + DENDRO_827 - DENDRO_828 + DENDRO_829 +
                          DENDRO_830;
const double DENDRO_877 = -DENDRO_832 - DENDRO_833 + DENDRO_834 - DENDRO_835 +
                          DENDRO_836 + DENDRO_837 + DENDRO_838 - DENDRO_839 +
                          DENDRO_840;

// Dendro: printing variables
//--
a_rhs[pp] = -DENDRO_0 * K[pp] + lambda[0] * (beta0[pp] * grad_0_alpha[pp] +
                                             beta1[pp] * grad_1_alpha[pp] +
                                             beta2[pp] * grad_2_alpha[pp]);
//--
b_rhs0[pp] =
    BSSN_XI[1] * (beta0[pp] * grad_0_beta0[pp] + beta1[pp] * grad_1_beta0[pp] +
                  beta2[pp] * grad_2_beta0[pp]) +
    DENDRO_1 * Gt0[pp] - DENDRO_26 * beta0[pp];
//--
b_rhs1[pp] =
    BSSN_XI[1] * (beta0[pp] * grad_0_beta1[pp] + beta1[pp] * grad_1_beta1[pp] +
                  beta2[pp] * grad_2_beta1[pp]) +
    DENDRO_1 * Gt1[pp] - DENDRO_26 * beta1[pp];
//--
b_rhs2[pp] =
    BSSN_XI[1] * (beta0[pp] * grad_0_beta2[pp] + beta1[pp] * grad_1_beta2[pp] +
                  beta2[pp] * grad_2_beta2[pp]) +
    DENDRO_1 * Gt2[pp] - DENDRO_26 * beta2[pp];
//--
gt_rhs00[pp] = -At0[pp] * DENDRO_0 + DENDRO_27 * gt0[pp] +
               DENDRO_28 * grad_0_beta1[pp] - DENDRO_29 * grad_1_beta1[pp] -
               DENDRO_29 * grad_2_beta2[pp] + DENDRO_7 * grad_0_beta2[pp] +
               beta0[pp] * grad_0_gt0[pp] + beta1[pp] * grad_1_gt0[pp] +
               beta2[pp] * grad_2_gt0[pp];
//--
gt_rhs01[pp] = -At1[pp] * DENDRO_0 + DENDRO_30 * grad_0_beta0[pp] +
               DENDRO_30 * grad_1_beta1[pp] - DENDRO_31 * gt1[pp] +
               beta0[pp] * grad_0_gt1[pp] + beta1[pp] * grad_1_gt1[pp] +
               beta2[pp] * grad_2_gt1[pp] + grad_0_beta1[pp] * gt3[pp] +
               grad_0_beta2[pp] * gt4[pp] + grad_1_beta0[pp] * gt0[pp] +
               grad_1_beta2[pp] * gt2[pp];
//--
gt_rhs02[pp] = -At2[pp] * DENDRO_0 + DENDRO_32 * grad_0_beta0[pp] +
               DENDRO_32 * grad_2_beta2[pp] - DENDRO_33 * gt2[pp] +
               beta0[pp] * grad_0_gt2[pp] + beta1[pp] * grad_1_gt2[pp] +
               beta2[pp] * grad_2_gt2[pp] + grad_0_beta1[pp] * gt4[pp] +
               grad_0_beta2[pp] * gt5[pp] + grad_2_beta0[pp] * gt0[pp] +
               grad_2_beta1[pp] * gt1[pp];
//--
gt_rhs11[pp] = -At3[pp] * DENDRO_0 + DENDRO_28 * grad_1_beta0[pp] -
               DENDRO_31 * gt3[pp] - DENDRO_34 * gt3[pp] + DENDRO_35 * gt3[pp] +
               DENDRO_36 * grad_1_beta2[pp] + beta0[pp] * grad_0_gt3[pp] +
               beta1[pp] * grad_1_gt3[pp] + beta2[pp] * grad_2_gt3[pp];
//--
gt_rhs12[pp] = -At4[pp] * DENDRO_0 - DENDRO_34 * gt4[pp] +
               DENDRO_37 * grad_1_beta1[pp] + DENDRO_37 * grad_2_beta2[pp] +
               beta0[pp] * grad_0_gt4[pp] + beta1[pp] * grad_1_gt4[pp] +
               beta2[pp] * grad_2_gt4[pp] + grad_1_beta0[pp] * gt2[pp] +
               grad_1_beta2[pp] * gt5[pp] + grad_2_beta0[pp] * gt1[pp] +
               grad_2_beta1[pp] * gt3[pp];
//--
gt_rhs22[pp] = -At5[pp] * DENDRO_0 - DENDRO_33 * gt5[pp] - DENDRO_34 * gt5[pp] +
               DENDRO_36 * grad_2_beta1[pp] + DENDRO_38 * gt5[pp] +
               DENDRO_7 * grad_2_beta0[pp] + beta0[pp] * grad_0_gt5[pp] +
               beta1[pp] * grad_1_gt5[pp] + beta2[pp] * grad_2_gt5[pp];
//--
chi_rhs[pp] = -DENDRO_39 * DENDRO_40 + DENDRO_39 * K[pp] * alpha[pp] +
              beta0[pp] * grad_0_chi[pp] + beta1[pp] * grad_1_chi[pp] +
              beta2[pp] * grad_2_chi[pp];
//--
At_rhs00[pp] =
    At0[pp] * DENDRO_27 - At0[pp] * DENDRO_31 - At0[pp] * DENDRO_33 +
    DENDRO_41 * grad_0_beta1[pp] + DENDRO_42 * grad_0_beta2[pp] +
    DENDRO_722 *
        (DENDRO_306 *
             (-DENDRO_103 * (DENDRO_100 * DENDRO_69 + DENDRO_101 * DENDRO_102 +
                             DENDRO_97 + DENDRO_98 * DENDRO_99) -
              DENDRO_121 * (DENDRO_110 - DENDRO_118) -
              DENDRO_121 * (DENDRO_130 + DENDRO_137) -
              DENDRO_121 * (DENDRO_140 + DENDRO_141 * DENDRO_143) -
              DENDRO_121 * (DENDRO_149 + DENDRO_153 * DENDRO_65) -
              DENDRO_121 * (DENDRO_147 * DENDRO_62 + DENDRO_155) +
              DENDRO_138 * DENDRO_269 -
              DENDRO_142 * DENDRO_262 * (DENDRO_124 + DENDRO_139) +
              DENDRO_158 * (DENDRO_167 + 1.0 * DENDRO_168) -
              DENDRO_158 * (-DENDRO_143 * DENDRO_66 + DENDRO_156 * DENDRO_81) +
              DENDRO_161 * DENDRO_271 * DENDRO_277 -
              DENDRO_164 * (DENDRO_160 + 1.0 * DENDRO_162) -
              DENDRO_164 *
                  (DENDRO_138 * DENDRO_165 + 1.0 * DENDRO_141 * DENDRO_81) +
              DENDRO_169 * (DENDRO_166 + DENDRO_168) +
              DENDRO_169 * (DENDRO_170 + DENDRO_171) - DENDRO_199 * DENDRO_200 -
              DENDRO_218 * DENDRO_219 - DENDRO_236 * DENDRO_237 -
              DENDRO_238 * (DENDRO_159 + DENDRO_162) -
              DENDRO_238 * (DENDRO_239 + DENDRO_240) - DENDRO_259 * DENDRO_260 +
              DENDRO_270 * DENDRO_272 * DENDRO_81 + DENDRO_273 * DENDRO_274 +
              DENDRO_275 * DENDRO_276 + DENDRO_305) +
         DENDRO_70 * DENDRO_72 + DENDRO_721 * gt0[pp] -
         DENDRO_83 * (-DENDRO_77 + DENDRO_81) -
         DENDRO_96 * (DENDRO_85 + DENDRO_87 - DENDRO_89 + DENDRO_95) -
         12 * grad2_0_0_alpha[pp]) -
    alpha[pp] *
        (-At0[pp] * K[pp] + DENDRO_46 * (DENDRO_43 + DENDRO_44 - DENDRO_45) +
         DENDRO_52 * (DENDRO_47 + DENDRO_49 - DENDRO_50) +
         DENDRO_55 * (-At0[pp] * DENDRO_14 + DENDRO_53 - DENDRO_54)) +
    beta0[pp] * grad_0_At0[pp] + beta1[pp] * grad_1_At0[pp] +
    beta2[pp] * grad_2_At0[pp];
//--
At_rhs01[pp] =
    At0[pp] * grad_1_beta0[pp] - At1[pp] * DENDRO_31 +
    At2[pp] * grad_1_beta2[pp] + At3[pp] * grad_0_beta1[pp] +
    At4[pp] * grad_0_beta2[pp] +
    DENDRO_722 *
        (DENDRO_306 *
             (-DENDRO_103 * (DENDRO_116 * DENDRO_526 + DENDRO_129 * DENDRO_529 +
                             DENDRO_179 * DENDRO_524 + DENDRO_683) +
              DENDRO_121 * (-DENDRO_287 * DENDRO_741 + DENDRO_717) -
              DENDRO_121 * (-DENDRO_617 + DENDRO_668 + DENDRO_669) +
              DENDRO_121 * (DENDRO_147 * DENDRO_412 - DENDRO_153 * DENDRO_618 +
                            DENDRO_648) +
              DENDRO_158 * (DENDRO_744 + DENDRO_745) +
              DENDRO_158 * (DENDRO_747 + DENDRO_750) +
              DENDRO_158 * (DENDRO_161 * DENDRO_407 + DENDRO_746) +
              DENDRO_158 * (DENDRO_137 + DENDRO_662 + DENDRO_663) +
              DENDRO_164 * (-DENDRO_62 * DENDRO_741 + DENDRO_718) +
              DENDRO_164 * (DENDRO_106 * DENDRO_161 - DENDRO_153 * DENDRO_281 +
                            DENDRO_677) -
              DENDRO_164 * (DENDRO_138 * DENDRO_499 + DENDRO_138 * DENDRO_513 +
                            DENDRO_194 * DENDRO_81) +
              DENDRO_164 * (-DENDRO_138 * DENDRO_618 - DENDRO_195 * DENDRO_81 +
                            DENDRO_680) -
              DENDRO_238 * (DENDRO_161 * grad_0_gt3[pp] + DENDRO_273) +
              DENDRO_261 * (DENDRO_736 + DENDRO_740) +
              DENDRO_262 * (-DENDRO_355 + DENDRO_501 + DENDRO_600) -
              DENDRO_268 * (DENDRO_106 * DENDRO_741 + DENDRO_715) -
              DENDRO_268 *
                  (-DENDRO_138 * DENDRO_381 + DENDRO_657 + DENDRO_743) -
              DENDRO_268 *
                  (-DENDRO_153 * DENDRO_285 + DENDRO_716 + DENDRO_742) +
              DENDRO_272 * (DENDRO_160 + DENDRO_161 * DENDRO_62 +
                            DENDRO_161 * DENDRO_63) +
              DENDRO_272 * (0.25 * DENDRO_239 + 0.5 * DENDRO_240 +
                            DENDRO_407 * DENDRO_81) +
              DENDRO_689 + DENDRO_690 + DENDRO_691 + DENDRO_692 + DENDRO_693 +
              DENDRO_694 + DENDRO_695 + DENDRO_696 + DENDRO_697 + DENDRO_698 +
              DENDRO_699 + DENDRO_700 + DENDRO_701 - DENDRO_702 * DENDRO_758 +
              DENDRO_719 -
              DENDRO_753 * (DENDRO_147 * grad_0_gt3[pp] + DENDRO_752) -
              DENDRO_753 * (DENDRO_138 * grad_1_gt5[pp] +
                            DENDRO_142 * grad_2_gt3[pp] + DENDRO_754) -
              DENDRO_755 * grad_0_gt1[pp] - DENDRO_756 * grad_1_gt1[pp] -
              DENDRO_757 * grad_2_gt1[pp]) +
         DENDRO_702 * DENDRO_720 +
         DENDRO_733 * DENDRO_9 * (DENDRO_129 - DENDRO_713 * DENDRO_732) -
         DENDRO_734 * (-DENDRO_704 + DENDRO_705 + DENDRO_706 + DENDRO_707) -
         DENDRO_735 * (DENDRO_708 + DENDRO_709 - DENDRO_710 - DENDRO_711) -
         12 * grad2_0_1_alpha[pp]) +
    DENDRO_723 * grad_0_beta0[pp] + DENDRO_723 * grad_1_beta1[pp] -
    alpha[pp] * (-At1[pp] * K[pp] + DENDRO_46 * DENDRO_726 +
                 DENDRO_52 * DENDRO_729 + DENDRO_55 * DENDRO_731) +
    beta0[pp] * grad_0_At1[pp] + beta1[pp] * grad_1_At1[pp] +
    beta2[pp] * grad_2_At1[pp];
//--
At_rhs02[pp] =
    At0[pp] * grad_2_beta0[pp] + At1[pp] * grad_2_beta1[pp] -
    At2[pp] * DENDRO_33 + At4[pp] * grad_0_beta1[pp] +
    At5[pp] * grad_0_beta2[pp] +
    DENDRO_722 *
        (DENDRO_306 *
             (-DENDRO_103 * (DENDRO_109 * DENDRO_526 + DENDRO_134 * DENDRO_529 +
                             DENDRO_172 * DENDRO_524 + DENDRO_521) -
              DENDRO_121 * (DENDRO_488 - DENDRO_489 + DENDRO_569) +
              DENDRO_121 * (-DENDRO_147 * DENDRO_499 - DENDRO_62 * DENDRO_769 +
                            DENDRO_775) -
              DENDRO_121 * (DENDRO_287 * DENDRO_767 + DENDRO_770 + DENDRO_771) +
              DENDRO_121 * (DENDRO_338 * DENDRO_748 - DENDRO_770 + DENDRO_776) -
              DENDRO_158 * (DENDRO_125 * DENDRO_143 - DENDRO_772 - DENDRO_773) -
              DENDRO_158 * (DENDRO_125 * DENDRO_161 - DENDRO_147 * DENDRO_165 -
                            DENDRO_67 * DENDRO_769) +
              DENDRO_158 * (DENDRO_65 * DENDRO_767 + DENDRO_772 + DENDRO_773) -
              DENDRO_164 * (DENDRO_140 + DENDRO_750) -
              DENDRO_164 * (DENDRO_745 + DENDRO_774) -
              DENDRO_164 * (DENDRO_161 * DENDRO_357 + DENDRO_746) -
              DENDRO_164 * (DENDRO_140 + DENDRO_142 * DENDRO_513 + DENDRO_749) +
              DENDRO_169 * (DENDRO_161 * grad_0_gt5[pp] + DENDRO_275) -
              DENDRO_262 * (DENDRO_125 * DENDRO_767 - 0.5 * DENDRO_766) -
              DENDRO_262 * (-DENDRO_139 * DENDRO_147 - DENDRO_65 * DENDRO_769 +
                            DENDRO_768) +
              DENDRO_267 * (DENDRO_752 + DENDRO_765) +
              DENDRO_268 * (DENDRO_138 * DENDRO_184 + 0.25 * DENDRO_754) +
              DENDRO_272 * (0.25 * DENDRO_170 + 1.0 * DENDRO_171) +
              DENDRO_272 * (DENDRO_161 * DENDRO_65 + DENDRO_161 * DENDRO_66 +
                            DENDRO_167) +
              DENDRO_531 + DENDRO_533 + DENDRO_535 + DENDRO_537 + DENDRO_539 +
              DENDRO_540 + DENDRO_542 + DENDRO_543 + DENDRO_544 + DENDRO_545 +
              DENDRO_546 + DENDRO_547 + DENDRO_548 - DENDRO_549 * DENDRO_758 +
              DENDRO_583 - DENDRO_753 * (DENDRO_740 + DENDRO_777) -
              DENDRO_755 * grad_0_gt2[pp] - DENDRO_756 * grad_1_gt2[pp] -
              DENDRO_757 * grad_2_gt2[pp]) +
         DENDRO_549 * DENDRO_720 + DENDRO_563 * DENDRO_735 -
         DENDRO_733 * (DENDRO_556 + DENDRO_557 - DENDRO_558 + DENDRO_559) -
         DENDRO_734 * (DENDRO_551 + DENDRO_552 - DENDRO_553 + DENDRO_554) -
         12 * grad2_0_2_alpha[pp]) +
    DENDRO_759 * grad_0_beta0[pp] + DENDRO_759 * grad_2_beta2[pp] -
    alpha[pp] * (-At2[pp] * K[pp] + DENDRO_46 * DENDRO_761 +
                 DENDRO_52 * DENDRO_763 + DENDRO_55 * DENDRO_764) +
    beta0[pp] * grad_0_At2[pp] + beta1[pp] * grad_1_At2[pp] +
    beta2[pp] * grad_2_At2[pp];
//--
At_rhs11[pp] =
    -At3[pp] * DENDRO_31 - At3[pp] * DENDRO_34 + At3[pp] * DENDRO_35 +
    DENDRO_41 * grad_1_beta0[pp] +
    DENDRO_722 *
        (DENDRO_306 *
             (-DENDRO_103 * (DENDRO_100 * DENDRO_215 + DENDRO_102 * DENDRO_233 +
                             DENDRO_196 * DENDRO_99 + DENDRO_397) -
              DENDRO_121 * (DENDRO_413 - 1.0 * DENDRO_415) +
              DENDRO_121 * (DENDRO_420 - 1.0 * DENDRO_783) -
              DENDRO_121 * (DENDRO_416 * DENDRO_748 + DENDRO_418) +
              DENDRO_138 * DENDRO_436 + DENDRO_153 * DENDRO_435 +
              DENDRO_158 * (-DENDRO_410 + DENDRO_411) +
              DENDRO_158 * (DENDRO_148 * DENDRO_737 + DENDRO_153 * DENDRO_407) +
              DENDRO_158 * (DENDRO_407 * DENDRO_748 + DENDRO_784) +
              DENDRO_158 * (DENDRO_63 * DENDRO_786 + DENDRO_785) +
              DENDRO_164 * (DENDRO_423 + DENDRO_743) +
              DENDRO_164 * (DENDRO_425 + DENDRO_742) +
              DENDRO_164 * (DENDRO_427 - 1.0 * DENDRO_787) +
              DENDRO_238 * (DENDRO_426 - DENDRO_787) +
              DENDRO_238 * (-DENDRO_138 * grad_2_gt3[pp] + DENDRO_432) +
              DENDRO_238 * (-DENDRO_153 * grad_0_gt3[pp] + DENDRO_429) -
              DENDRO_262 * DENDRO_748 * (DENDRO_186 + DENDRO_354) +
              6.0 * DENDRO_267 * DENDRO_741 * grad_1_gt3[pp] +
              DENDRO_365 * (DENDRO_419 - DENDRO_783) +
              DENDRO_365 * (DENDRO_430 - DENDRO_737 * grad_0_gt3[pp]) +
              DENDRO_365 * (DENDRO_431 - DENDRO_748 * grad_2_gt3[pp]) -
              DENDRO_395 * DENDRO_758 + DENDRO_434 * DENDRO_737 + DENDRO_441 -
              DENDRO_755 * grad_0_gt3[pp] - DENDRO_756 * grad_1_gt3[pp] -
              DENDRO_757 * grad_2_gt3[pp]) -
         DENDRO_71 * (DENDRO_391 - DENDRO_392 + DENDRO_393 + DENDRO_396) +
         DENDRO_721 * gt3[pp] +
         DENDRO_782 *
             (DENDRO_196 - DENDRO_781 * (-DENDRO_11 + DENDRO_15 + DENDRO_92)) +
         DENDRO_83 * (DENDRO_233 - DENDRO_732 * DENDRO_781) -
         12 * grad2_1_1_alpha[pp]) +
    DENDRO_778 * grad_1_beta2[pp] -
    alpha[pp] * (-At3[pp] * K[pp] + DENDRO_46 * DENDRO_729 +
                 DENDRO_726 * DENDRO_779 + DENDRO_731 * DENDRO_780) +
    beta0[pp] * grad_0_At3[pp] + beta1[pp] * grad_1_At3[pp] +
    beta2[pp] * grad_2_At3[pp];
//--
At_rhs12[pp] =
    At1[pp] * grad_2_beta0[pp] + At2[pp] * grad_1_beta0[pp] +
    At3[pp] * grad_2_beta1[pp] - At4[pp] * DENDRO_34 +
    At5[pp] * grad_1_beta2[pp] +
    DENDRO_722 *
        (DENDRO_306 *
             (-DENDRO_103 * (DENDRO_177 * DENDRO_524 + DENDRO_206 * DENDRO_526 +
                             DENDRO_224 * DENDRO_529 + DENDRO_623) +
              DENDRO_121 * (DENDRO_187 * DENDRO_741 + DENDRO_645) +
              DENDRO_121 * (-DENDRO_194 * DENDRO_767 + DENDRO_791) +
              DENDRO_121 * (DENDRO_334 * DENDRO_748 + DENDRO_791) -
              DENDRO_121 * (DENDRO_104 * DENDRO_769 + DENDRO_499 * DENDRO_737 +
                            DENDRO_644) +
              DENDRO_121 * (-DENDRO_513 * DENDRO_737 + DENDRO_604 -
                            DENDRO_618 * DENDRO_737) +
              DENDRO_157 * (DENDRO_736 + DENDRO_738 + DENDRO_777) -
              DENDRO_158 * (DENDRO_143 * DENDRO_187 + DENDRO_776 - DENDRO_790) -
              DENDRO_158 * (-DENDRO_165 * DENDRO_737 - DENDRO_63 * DENDRO_769 +
                            DENDRO_775) +
              DENDRO_158 * (DENDRO_407 * DENDRO_767 + DENDRO_771 + DENDRO_790) +
              DENDRO_164 * (-DENDRO_142 * DENDRO_381 + DENDRO_792) +
              DENDRO_164 * (-DENDRO_357 * DENDRO_741 + DENDRO_647) +
              DENDRO_164 * (-DENDRO_618 * DENDRO_748 + DENDRO_792) +
              DENDRO_164 *
                  (-DENDRO_147 * DENDRO_285 + DENDRO_648 - 0.25 * DENDRO_765) +
              DENDRO_164 *
                  (-DENDRO_281 * DENDRO_737 + DENDRO_612 - DENDRO_785) +
              DENDRO_164 *
                  (-DENDRO_407 * DENDRO_741 + DENDRO_473 + DENDRO_614) -
              DENDRO_262 * (DENDRO_187 * DENDRO_767 - 0.5 * DENDRO_789) -
              DENDRO_262 * (-DENDRO_139 * DENDRO_737 + DENDRO_338 * DENDRO_737 -
                            DENDRO_407 * DENDRO_769) -
              DENDRO_268 * (-DENDRO_381 * DENDRO_748 + DENDRO_640) -
              DENDRO_268 * (-DENDRO_194 * DENDRO_741 - DENDRO_195 * DENDRO_741 +
                            DENDRO_420) -
              DENDRO_268 *
                  (-DENDRO_285 * DENDRO_737 + DENDRO_590 + DENDRO_641) +
              DENDRO_272 * (DENDRO_122 * DENDRO_138 + DENDRO_747) +
              DENDRO_272 * (DENDRO_155 + DENDRO_744 + DENDRO_774) +
              DENDRO_365 * (DENDRO_440 - DENDRO_741 * grad_1_gt5[pp]) -
              DENDRO_628 * DENDRO_758 + DENDRO_629 + DENDRO_653 -
              DENDRO_755 * grad_0_gt4[pp] - DENDRO_756 * grad_1_gt4[pp] -
              DENDRO_757 * grad_2_gt4[pp]) +
         DENDRO_628 * DENDRO_720 + DENDRO_632 * DENDRO_735 -
         DENDRO_733 * (-DENDRO_633 + DENDRO_634 + DENDRO_635 + DENDRO_636) -
         DENDRO_734 * DENDRO_9 * (-DENDRO_638 + DENDRO_737) -
         12 * grad2_1_2_alpha[pp]) +
    DENDRO_788 * grad_1_beta1[pp] + DENDRO_788 * grad_2_beta2[pp] -
    alpha[pp] * (-At4[pp] * K[pp] + DENDRO_46 * DENDRO_763 +
                 DENDRO_761 * DENDRO_779 + DENDRO_764 * DENDRO_780) +
    beta0[pp] * grad_0_At4[pp] + beta1[pp] * grad_1_At4[pp] +
    beta2[pp] * grad_2_At4[pp];
//--
At_rhs22[pp] =
    -At5[pp] * DENDRO_33 - At5[pp] * DENDRO_34 + At5[pp] * DENDRO_38 +
    DENDRO_42 * grad_2_beta0[pp] +
    DENDRO_722 *
        (DENDRO_306 *
             (-DENDRO_103 * (DENDRO_100 * DENDRO_210 + DENDRO_102 * DENDRO_230 +
                             DENDRO_189 * DENDRO_99 + DENDRO_321) -
              DENDRO_121 * (DENDRO_344 - 1.0 * DENDRO_346) -
              DENDRO_121 * (0.25 * DENDRO_789 + 1.0 * DENDRO_793) +
              DENDRO_142 * DENDRO_373 * grad_0_gt5[pp] +
              DENDRO_147 * DENDRO_371 +
              DENDRO_158 * (DENDRO_335 - 1.0 * DENDRO_337) +
              DENDRO_158 * (0.25 * DENDRO_766 + 1.0 * DENDRO_794) -
              DENDRO_158 * (-DENDRO_340 * DENDRO_769 + DENDRO_768) -
              DENDRO_164 * (DENDRO_122 * DENDRO_748 + DENDRO_770) -
              DENDRO_164 * (DENDRO_142 * DENDRO_184 + DENDRO_790) -
              DENDRO_164 * (DENDRO_147 * DENDRO_357 + DENDRO_154 * DENDRO_737) -
              DENDRO_164 * (DENDRO_66 * DENDRO_786 + 0.25 * DENDRO_739) +
              DENDRO_169 * (DENDRO_766 + DENDRO_794) +
              DENDRO_169 *
                  (DENDRO_147 * grad_0_gt5[pp] + DENDRO_769 * grad_2_gt0[pp]) +
              6.0 * DENDRO_261 * DENDRO_767 * grad_2_gt5[pp] -
              DENDRO_262 * DENDRO_769 * (DENDRO_124 + DENDRO_323) +
              DENDRO_274 * DENDRO_748 * grad_1_gt5[pp] -
              DENDRO_318 * DENDRO_758 - DENDRO_365 * (DENDRO_789 + DENDRO_793) -
              DENDRO_365 *
                  (DENDRO_141 * DENDRO_769 + DENDRO_737 * grad_0_gt5[pp]) +
              DENDRO_370 * DENDRO_737 + DENDRO_388 -
              DENDRO_755 * grad_0_gt5[pp] - DENDRO_756 * grad_1_gt5[pp] -
              DENDRO_757 * grad_2_gt5[pp]) +
         DENDRO_308 * DENDRO_72 + DENDRO_721 * gt5[pp] -
         DENDRO_782 * (-DENDRO_311 + DENDRO_769) -
         DENDRO_82 * (DENDRO_314 - DENDRO_315 + DENDRO_316 + DENDRO_319) -
         12 * grad2_2_2_alpha[pp]) +
    DENDRO_778 * grad_2_beta1[pp] -
    alpha[pp] * (At5[pp] * DENDRO_51 * DENDRO_764 - At5[pp] * K[pp] +
                 DENDRO_55 * DENDRO_763 + DENDRO_761 * DENDRO_780) +
    beta0[pp] * grad_0_At5[pp] + beta1[pp] * grad_1_At5[pp] +
    beta2[pp] * grad_2_At5[pp];
//--
K_rhs[pp] =
    -DENDRO_289 * chi[pp] *
        (DENDRO_799 * (DENDRO_444 + DENDRO_60 * DENDRO_803) +
         DENDRO_802 * (DENDRO_445 + DENDRO_60 * DENDRO_801) -
         grad2_0_0_alpha[pp] +
         grad_0_alpha[pp] *
             (DENDRO_442 * DENDRO_9 + DENDRO_443 * DENDRO_9 -
              DENDRO_59 * (-0.5 * DENDRO_795 * DENDRO_91 + DENDRO_90) +
              DENDRO_89)) -
    DENDRO_291 * chi[pp] *
        (DENDRO_797 * (DENDRO_399 + DENDRO_781 * DENDRO_795) +
         DENDRO_802 * (DENDRO_402 + DENDRO_781 * DENDRO_800) -
         grad2_1_1_alpha[pp] +
         grad_1_alpha[pp] *
             (DENDRO_392 + DENDRO_400 * DENDRO_9 + DENDRO_401 * DENDRO_9 -
              DENDRO_59 * (DENDRO_394 - DENDRO_395 * DENDRO_803))) -
    DENDRO_293 * chi[pp] *
        (DENDRO_797 * (DENDRO_327 + DENDRO_795 * DENDRO_796) +
         DENDRO_799 * (DENDRO_329 + DENDRO_796 * DENDRO_798) -
         grad2_2_2_alpha[pp] +
         grad_2_alpha[pp] *
             (DENDRO_331 * DENDRO_9 + DENDRO_332 * DENDRO_9 +
              DENDRO_333 * DENDRO_9 -
              DENDRO_59 * (DENDRO_317 - DENDRO_318 * DENDRO_801))) +
    DENDRO_806 * chi[pp] *
        (0.5 * DENDRO_797 * (DENDRO_624 + DENDRO_637 * DENDRO_795) +
         DENDRO_804 *
             (-DENDRO_59 * (-DENDRO_628 * DENDRO_798 + grad_2_chi[pp]) +
              DENDRO_625 * DENDRO_9 + DENDRO_631) +
         DENDRO_805 *
             (-DENDRO_59 * (-DENDRO_628 * DENDRO_800 + grad_1_chi[pp]) +
              DENDRO_626 * DENDRO_9 + DENDRO_627 * DENDRO_9 + DENDRO_633) -
         grad2_1_2_alpha[pp]) -
    DENDRO_808 * chi[pp] *
        (0.5 * DENDRO_799 * (DENDRO_525 + DENDRO_562 * DENDRO_798) +
         DENDRO_805 *
             (DENDRO_527 * DENDRO_9 + DENDRO_528 * DENDRO_9 + DENDRO_558 -
              DENDRO_59 * (-DENDRO_549 * DENDRO_800 + grad_0_chi[pp])) +
         DENDRO_807 *
             (DENDRO_522 * DENDRO_9 + DENDRO_523 * DENDRO_9 + DENDRO_553 -
              DENDRO_59 * (-DENDRO_549 * DENDRO_795 + grad_2_chi[pp])) -
         grad2_0_2_alpha[pp]) +
    DENDRO_809 * chi[pp] *
        (0.5 * DENDRO_802 * (DENDRO_687 + DENDRO_713 * DENDRO_800) +
         DENDRO_804 *
             (-DENDRO_59 * (-DENDRO_702 * DENDRO_798 + grad_0_chi[pp]) +
              DENDRO_686 * DENDRO_9 + DENDRO_712) +
         DENDRO_807 *
             (-DENDRO_59 * (-DENDRO_702 * DENDRO_795 + grad_1_chi[pp]) +
              DENDRO_684 * DENDRO_9 + DENDRO_685 * DENDRO_9 + DENDRO_704) -
         grad2_0_1_alpha[pp]) +
    DENDRO_842 *
        (At0[pp] * DENDRO_813 * DENDRO_814 + At1[pp] * DENDRO_821 * DENDRO_831 +
         At2[pp] * DENDRO_820 * DENDRO_821 + At3[pp] * DENDRO_814 * DENDRO_817 +
         At4[pp] * DENDRO_821 * DENDRO_841 + At5[pp] * DENDRO_814 * DENDRO_819 +
         pow(K[pp], 2)) +
    beta0[pp] * grad_0_K[pp] + beta1[pp] * grad_1_K[pp] +
    beta2[pp] * grad_2_K[pp];
//--
Gt_rhs0[pp] =
    DENDRO_172 * DENDRO_863 + DENDRO_177 * DENDRO_865 +
    DENDRO_179 * DENDRO_864 + DENDRO_189 * DENDRO_855 +
    DENDRO_196 * DENDRO_854 - 4.0 / 3.0 * DENDRO_289 * grad2_0_0_beta0[pp] -
    DENDRO_291 * grad2_1_1_beta0[pp] - DENDRO_293 * grad2_2_2_beta0[pp] +
    DENDRO_358 * DENDRO_870 - DENDRO_358 * grad_0_beta0[pp] -
    DENDRO_360 * grad_1_beta0[pp] - DENDRO_362 * grad_2_beta0[pp] +
    DENDRO_806 * grad2_1_2_beta0[pp] - DENDRO_813 * DENDRO_850 -
    DENDRO_820 * DENDRO_860 - DENDRO_831 * DENDRO_861 +
    DENDRO_844 * grad2_0_1_beta0[pp] - DENDRO_845 * grad2_0_2_beta0[pp] +
    DENDRO_846 * grad2_1_1_beta1[pp] + DENDRO_846 * grad2_1_2_beta2[pp] -
    DENDRO_847 * grad2_1_2_beta1[pp] - DENDRO_847 * grad2_2_2_beta2[pp] -
    DENDRO_848 * grad2_0_1_beta1[pp] - DENDRO_848 * grad2_0_2_beta2[pp] +
    DENDRO_853 * DENDRO_98 -
    DENDRO_859 * (DENDRO_10 * DENDRO_868 + DENDRO_831 * DENDRO_869) -
    DENDRO_859 * (-DENDRO_14 * DENDRO_866 + DENDRO_820 * DENDRO_867) -
    DENDRO_859 * (-DENDRO_19 * DENDRO_856 + DENDRO_813 * DENDRO_858) +
    beta0[pp] * grad_0_Gt0[pp] + beta1[pp] * grad_1_Gt0[pp] +
    beta2[pp] * grad_2_Gt0[pp];
//--
Gt_rhs1[pp] = DENDRO_109 * DENDRO_863 - DENDRO_116 * DENDRO_862 * DENDRO_876 +
              DENDRO_199 * grad_0_beta1[pp] -
              DENDRO_206 * DENDRO_862 * DENDRO_877 + DENDRO_210 * DENDRO_855 -
              DENDRO_218 * DENDRO_870 + DENDRO_218 * grad_1_beta1[pp] +
              DENDRO_236 * grad_2_beta1[pp] - DENDRO_289 * grad2_0_0_beta1[pp] -
              4.0 / 3.0 * DENDRO_291 * grad2_1_1_beta1[pp] -
              DENDRO_293 * grad2_2_2_beta1[pp] + DENDRO_69 * DENDRO_853 -
              DENDRO_741 * DENDRO_854 - DENDRO_808 * grad2_0_2_beta1[pp] -
              DENDRO_817 * DENDRO_861 + DENDRO_844 * grad2_0_1_beta1[pp] +
              DENDRO_846 * grad2_0_0_beta0[pp] +
              DENDRO_846 * grad2_0_2_beta2[pp] + DENDRO_850 * DENDRO_876 -
              DENDRO_859 * (DENDRO_10 * DENDRO_856 - DENDRO_858 * DENDRO_876) -
              DENDRO_859 * (DENDRO_16 * DENDRO_866 - DENDRO_867 * DENDRO_877) +
              DENDRO_859 * (DENDRO_22 * DENDRO_868 - DENDRO_817 * DENDRO_869) +
              DENDRO_860 * DENDRO_877 + DENDRO_871 * DENDRO_872 +
              DENDRO_873 * grad2_1_2_beta1[pp] +
              DENDRO_874 * grad2_2_2_beta2[pp] -
              DENDRO_875 * grad2_0_1_beta0[pp] -
              DENDRO_875 * grad2_1_2_beta2[pp] + beta0[pp] * grad_0_Gt1[pp] +
              beta1[pp] * grad_1_Gt1[pp] + beta2[pp] * grad_2_Gt1[pp];
//--
Gt_rhs2[pp] =
    DENDRO_101 * DENDRO_853 + DENDRO_129 * DENDRO_864 +
    DENDRO_134 * DENDRO_863 + DENDRO_224 * DENDRO_865 +
    DENDRO_230 * DENDRO_855 + DENDRO_233 * DENDRO_854 -
    DENDRO_289 * grad2_0_0_beta2[pp] - DENDRO_291 * grad2_1_1_beta2[pp] -
    DENDRO_293 * DENDRO_872 - 1.0 / 3.0 * DENDRO_293 * grad2_1_2_beta1[pp] -
    4.0 / 3.0 * DENDRO_293 * grad2_2_2_beta2[pp] -
    DENDRO_358 * grad_0_beta2[pp] - DENDRO_360 * grad_1_beta2[pp] +
    DENDRO_362 * DENDRO_870 - DENDRO_362 * grad_2_beta2[pp] +
    DENDRO_809 * grad2_0_1_beta2[pp] - DENDRO_819 * DENDRO_860 -
    DENDRO_820 * DENDRO_850 - DENDRO_841 * DENDRO_861 -
    DENDRO_845 * grad2_0_2_beta2[pp] - DENDRO_847 * grad2_0_0_beta0[pp] -
    DENDRO_847 * grad2_0_1_beta1[pp] -
    DENDRO_859 * (-DENDRO_14 * DENDRO_856 + DENDRO_820 * DENDRO_858) -
    DENDRO_859 * (DENDRO_16 * DENDRO_868 + DENDRO_841 * DENDRO_869) -
    DENDRO_859 * (-DENDRO_25 * DENDRO_866 + DENDRO_819 * DENDRO_867) +
    DENDRO_873 * grad2_1_2_beta2[pp] + DENDRO_874 * grad2_0_1_beta0[pp] +
    DENDRO_874 * grad2_1_1_beta1[pp] + beta0[pp] * grad_0_Gt2[pp] +
    beta1[pp] * grad_1_Gt2[pp] + beta2[pp] * grad_2_Gt2[pp];
// Dendro: reduced ops: 3913
// Dendro: }}}
