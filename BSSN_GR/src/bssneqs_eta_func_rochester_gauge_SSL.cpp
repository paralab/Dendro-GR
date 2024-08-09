// CODEGEN: SSL was enabled, adding term to gauge condition!
// Codgen: generating unstage version
// Codgen: using rochester gauge
// Codgen: using eta func damping
// Dendro: {{{
// Dendro: original ops: 607517
// Dendro: printing temp variables
const double DENDRO_0 = 2 * alpha[pp];
const double DENDRO_1 = sqrt(chi[pp]);
const double DENDRO_2 = (3.0 / 4.0) * BSSN_XI[2];
const double DENDRO_3 = pow(gt4[pp], 2);
const double DENDRO_4 = pow(gt1[pp], 2);
const double DENDRO_5 = pow(gt2[pp], 2);
const double DENDRO_6 = gt3[pp] * gt5[pp];
const double DENDRO_7 = gt1[pp] * gt4[pp];
const double DENDRO_8 = 2 * gt2[pp];
const double DENDRO_9 = DENDRO_3 * gt0[pp] + DENDRO_4 * gt5[pp] +
                        DENDRO_5 * gt3[pp] - DENDRO_6 * gt0[pp] -
                        DENDRO_7 * DENDRO_8;
const double DENDRO_10 = 1.0 / DENDRO_9;
const double DENDRO_11 = gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp];
const double DENDRO_12 = DENDRO_7 - gt2[pp] * gt3[pp];
const double DENDRO_13 = DENDRO_12 * grad_2_chi[pp];
const double DENDRO_14 = gt0[pp] * gt4[pp] - gt1[pp] * gt2[pp];
const double DENDRO_15 = pow(grad_0_chi[pp], 2);
const double DENDRO_16 = -DENDRO_3 + DENDRO_6;
const double DENDRO_17 = pow(grad_1_chi[pp], 2);
const double DENDRO_18 = -DENDRO_5 + gt0[pp] * gt5[pp];
const double DENDRO_19 = pow(grad_2_chi[pp], 2);
const double DENDRO_20 = -DENDRO_4 + gt0[pp] * gt3[pp];
const double DENDRO_21 =
    BSSN_ETA_R0 *
    sqrt(DENDRO_10 * (2 * DENDRO_11 * grad_0_chi[pp] * grad_1_chi[pp] -
                      2 * DENDRO_13 * grad_0_chi[pp] +
                      2 * DENDRO_14 * grad_1_chi[pp] * grad_2_chi[pp] -
                      DENDRO_15 * DENDRO_16 - DENDRO_17 * DENDRO_18 -
                      DENDRO_19 * DENDRO_20)) *
    pow(1 - pow(chi[pp], BSSN_ETA_POWER[0]), -BSSN_ETA_POWER[1]);
const double DENDRO_22 = (4.0 / 3.0) * grad_0_beta0[pp];
const double DENDRO_23 = 2 * gt1[pp];
const double DENDRO_24 = (2.0 / 3.0) * gt0[pp];
const double DENDRO_25 = (1.0 / 3.0) * gt1[pp];
const double DENDRO_26 = (2.0 / 3.0) * grad_2_beta2[pp];
const double DENDRO_27 = (1.0 / 3.0) * gt2[pp];
const double DENDRO_28 = (2.0 / 3.0) * grad_1_beta1[pp];
const double DENDRO_29 = (2.0 / 3.0) * grad_0_beta0[pp];
const double DENDRO_30 = (4.0 / 3.0) * grad_1_beta1[pp];
const double DENDRO_31 = 2 * gt4[pp];
const double DENDRO_32 = (1.0 / 3.0) * gt4[pp];
const double DENDRO_33 = (4.0 / 3.0) * grad_2_beta2[pp];
const double DENDRO_34 = (2.0 / 3.0) * chi[pp];
const double DENDRO_35 = grad_0_beta0[pp] + grad_1_beta1[pp] + grad_2_beta2[pp];
const double DENDRO_36 = 2 * At1[pp];
const double DENDRO_37 = 2 * At2[pp];
const double DENDRO_38 = At0[pp] * DENDRO_11;
const double DENDRO_39 = At2[pp] * DENDRO_14;
const double DENDRO_40 = DENDRO_10 * DENDRO_36;
const double DENDRO_41 = At2[pp] * DENDRO_12;
const double DENDRO_42 = At0[pp] * DENDRO_16;
const double DENDRO_43 = 2 * DENDRO_10;
const double DENDRO_44 = At0[pp] * DENDRO_43;
const double DENDRO_45 = At2[pp] * DENDRO_20;
const double DENDRO_46 = DENDRO_10 * DENDRO_37;
const double DENDRO_47 =
    DENDRO_11 * grad_0_chi[pp] + DENDRO_14 * grad_2_chi[pp];
const double DENDRO_48 = -DENDRO_18 * grad_1_chi[pp] + DENDRO_47;
const double DENDRO_49 = 0.5 * DENDRO_48;
const double DENDRO_50 = 1.0 / chi[pp];
const double DENDRO_51 = DENDRO_50 * gt0[pp];
const double DENDRO_52 = 1.0 * grad_0_gt1[pp];
const double DENDRO_53 = 0.5 * grad_1_gt0[pp];
const double DENDRO_54 = DENDRO_52 - DENDRO_53;
const double DENDRO_55 = 1.0 * grad_0_gt2[pp];
const double DENDRO_56 = 0.5 * grad_2_gt0[pp];
const double DENDRO_57 = DENDRO_55 - DENDRO_56;
const double DENDRO_58 = 0.5 * grad_0_gt0[pp];
const double DENDRO_59 = DENDRO_11 * DENDRO_58 + DENDRO_14 * DENDRO_57;
const double DENDRO_60 = -DENDRO_18 * DENDRO_54 + DENDRO_59;
const double DENDRO_61 = DENDRO_49 * DENDRO_51 + DENDRO_60;
const double DENDRO_62 = 12 * grad_1_alpha[pp];
const double DENDRO_63 = DENDRO_10 * DENDRO_62;
const double DENDRO_64 = DENDRO_12 * grad_0_chi[pp] -
                         DENDRO_14 * grad_1_chi[pp] +
                         DENDRO_20 * grad_2_chi[pp];
const double DENDRO_65 = -DENDRO_64;
const double DENDRO_66 =
    DENDRO_12 * DENDRO_58 - DENDRO_14 * DENDRO_54 + DENDRO_20 * DENDRO_57;
const double DENDRO_67 = -0.5 * DENDRO_50 * DENDRO_65 * gt0[pp] + DENDRO_66;
const double DENDRO_68 = 12 * grad_2_alpha[pp];
const double DENDRO_69 = DENDRO_10 * DENDRO_68;
const double DENDRO_70 = DENDRO_16 * DENDRO_58;
const double DENDRO_71 = DENDRO_12 * DENDRO_57;
const double DENDRO_72 = DENDRO_11 * DENDRO_54;
const double DENDRO_73 = DENDRO_10 * DENDRO_72;
const double DENDRO_74 = 1.0 * grad_0_chi[pp];
const double DENDRO_75 = DENDRO_10 * gt0[pp];
const double DENDRO_76 =
    -DENDRO_11 * grad_1_chi[pp] + DENDRO_13 + DENDRO_16 * grad_0_chi[pp];
const double DENDRO_77 = -DENDRO_76;
const double DENDRO_78 = 0.5 * DENDRO_77;
const double DENDRO_79 = DENDRO_10 * DENDRO_70 + DENDRO_10 * DENDRO_71 +
                         DENDRO_50 * (DENDRO_74 - DENDRO_75 * DENDRO_78) -
                         DENDRO_73;
const double DENDRO_80  = 12 * grad_0_alpha[pp];
const double DENDRO_81  = pow(chi[pp], -2);
const double DENDRO_82  = DENDRO_15 * DENDRO_81;
const double DENDRO_83  = 4.0 * DENDRO_10;
const double DENDRO_84  = DENDRO_11 * DENDRO_83;
const double DENDRO_85  = DENDRO_84 * grad2_0_1_gt0[pp];
const double DENDRO_86  = DENDRO_14 * DENDRO_83;
const double DENDRO_87  = DENDRO_86 * grad2_1_2_gt0[pp];
const double DENDRO_88  = pow(DENDRO_9, -2);
const double DENDRO_89  = DENDRO_11 * grad_0_gt3[pp];
const double DENDRO_90  = DENDRO_16 * grad_1_gt0[pp];
const double DENDRO_91  = grad_0_gt4[pp] + grad_1_gt2[pp] - grad_2_gt1[pp];
const double DENDRO_92  = DENDRO_12 * DENDRO_91;
const double DENDRO_93  = -DENDRO_89 + DENDRO_90 + DENDRO_92;
const double DENDRO_94  = DENDRO_12 * grad_0_gt5[pp];
const double DENDRO_95  = DENDRO_16 * grad_2_gt0[pp];
const double DENDRO_96  = grad_0_gt4[pp] - grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_97  = DENDRO_11 * DENDRO_96;
const double DENDRO_98  = DENDRO_94 + DENDRO_95 - DENDRO_97;
const double DENDRO_99  = -DENDRO_11 * DENDRO_54 + DENDRO_70 + DENDRO_71;
const double DENDRO_100 = DENDRO_18 * grad_0_gt3[pp];
const double DENDRO_101 = DENDRO_11 * grad_1_gt0[pp];
const double DENDRO_102 = DENDRO_14 * DENDRO_91;
const double DENDRO_103 = DENDRO_101 + DENDRO_102;
const double DENDRO_104 = -DENDRO_100 + DENDRO_103;
const double DENDRO_105 = 0.25 * grad_0_gt3[pp];
const double DENDRO_106 = 1.0 * grad_1_gt1[pp];
const double DENDRO_107 = -DENDRO_106;
const double DENDRO_108 = DENDRO_18 * DENDRO_88;
const double DENDRO_109 = 4 * DENDRO_108;
const double DENDRO_110 = DENDRO_104 * DENDRO_109 * (-DENDRO_105 - DENDRO_107);
const double DENDRO_111 = DENDRO_12 * grad_2_gt0[pp];
const double DENDRO_112 = DENDRO_20 * grad_0_gt5[pp];
const double DENDRO_113 = DENDRO_111 + DENDRO_112 - DENDRO_14 * DENDRO_96;
const double DENDRO_114 = 0.25 * grad_0_gt5[pp];
const double DENDRO_115 = 1.0 * grad_2_gt2[pp];
const double DENDRO_116 = -DENDRO_115;
const double DENDRO_117 = DENDRO_114 + DENDRO_116;
const double DENDRO_118 = DENDRO_20 * DENDRO_88;
const double DENDRO_119 = 4 * DENDRO_118;
const double DENDRO_120 =
    DENDRO_11 * grad_2_gt0[pp] + DENDRO_14 * grad_0_gt5[pp];
const double DENDRO_121 = DENDRO_120 - DENDRO_18 * DENDRO_96;
const double DENDRO_122 = 0.25 * grad_0_gt4[pp];
const double DENDRO_123 = -DENDRO_122;
const double DENDRO_124 = 0.25 * grad_1_gt2[pp];
const double DENDRO_125 = 0.75 * grad_2_gt1[pp];
const double DENDRO_126 =
    DENDRO_119 * DENDRO_121 * (DENDRO_123 + DENDRO_124 + DENDRO_125);
const double DENDRO_127 = 0.75 * grad_1_gt2[pp];
const double DENDRO_128 = 0.25 * grad_2_gt1[pp];
const double DENDRO_129 = DENDRO_123 + DENDRO_127 + DENDRO_128;
const double DENDRO_130 = DENDRO_12 * grad_1_gt0[pp];
const double DENDRO_131 = DENDRO_20 * DENDRO_91;
const double DENDRO_132 = DENDRO_130 + DENDRO_131 - DENDRO_14 * grad_0_gt3[pp];
const double DENDRO_133 = DENDRO_16 * DENDRO_60;
const double DENDRO_134 = 4 * DENDRO_88;
const double DENDRO_135 = DENDRO_133 * DENDRO_134 * (DENDRO_52 + DENDRO_53);
const double DENDRO_136 = DENDRO_55 + DENDRO_56;
const double DENDRO_137 = 0.25 * grad_1_gt0[pp];
const double DENDRO_138 = DENDRO_137 * DENDRO_98;
const double DENDRO_139 = DENDRO_14 * DENDRO_88;
const double DENDRO_140 = 4 * DENDRO_139;
const double DENDRO_141 = 0.25 * grad_2_gt0[pp];
const double DENDRO_142 = DENDRO_141 * DENDRO_93;
const double DENDRO_143 = DENDRO_93 * grad_0_gt0[pp];
const double DENDRO_144 = DENDRO_99 * grad_1_gt0[pp];
const double DENDRO_145 = DENDRO_11 * DENDRO_88;
const double DENDRO_146 = 2.0 * DENDRO_145;
const double DENDRO_147 = DENDRO_98 * grad_0_gt0[pp];
const double DENDRO_148 = DENDRO_99 * grad_2_gt0[pp];
const double DENDRO_149 = DENDRO_104 * grad_1_gt0[pp];
const double DENDRO_150 = DENDRO_60 * grad_0_gt3[pp];
const double DENDRO_151 = DENDRO_149 + DENDRO_150;
const double DENDRO_152 = DENDRO_113 * grad_2_gt0[pp];
const double DENDRO_153 = DENDRO_66 * grad_0_gt5[pp];
const double DENDRO_154 = DENDRO_105 * DENDRO_121;
const double DENDRO_155 = -grad_0_gt4[pp] + grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_156 = 0.5 * DENDRO_155;
const double DENDRO_157 = DENDRO_104 * DENDRO_156 + DENDRO_154;
const double DENDRO_158 = DENDRO_114 * DENDRO_132;
const double DENDRO_159 = 0.5 * DENDRO_113;
const double DENDRO_160 = 0.25 * DENDRO_143;
const double DENDRO_161 = 4 * DENDRO_145;
const double DENDRO_162 = 0.25 * DENDRO_147;
const double DENDRO_163 = DENDRO_121 * grad_1_gt0[pp];
const double DENDRO_164 = DENDRO_60 * DENDRO_96;
const double DENDRO_165 = DENDRO_12 * DENDRO_88;
const double DENDRO_166 = 2.0 * DENDRO_165;
const double DENDRO_167 = DENDRO_166 * (DENDRO_163 + DENDRO_164);
const double DENDRO_168 = DENDRO_132 * grad_2_gt0[pp];
const double DENDRO_169 = DENDRO_66 * DENDRO_91;
const double DENDRO_170 = 0.5 * grad_0_gt3[pp];
const double DENDRO_171 = DENDRO_107 + DENDRO_170;
const double DENDRO_172 = DENDRO_121 * DENDRO_171;
const double DENDRO_173 = DENDRO_104 * DENDRO_96;
const double DENDRO_174 = 0.25 * DENDRO_173;
const double DENDRO_175 = -DENDRO_174;
const double DENDRO_176 = 0.5 * grad_0_gt5[pp];
const double DENDRO_177 = DENDRO_116 + DENDRO_176;
const double DENDRO_178 = -DENDRO_132;
const double DENDRO_179 = DENDRO_177 * DENDRO_178;
const double DENDRO_180 = -DENDRO_113;
const double DENDRO_181 = DENDRO_180 * DENDRO_91;
const double DENDRO_182 = -0.25 * DENDRO_181;
const double DENDRO_183 = 0.5 * DENDRO_54;
const double DENDRO_184 = DENDRO_121 * DENDRO_183;
const double DENDRO_185 = 4 * DENDRO_165;
const double DENDRO_186 = DENDRO_185 * (DENDRO_155 * DENDRO_60 + DENDRO_184);
const double DENDRO_187 = 0.5 * DENDRO_57;
const double DENDRO_188 = DENDRO_104 * DENDRO_183;
const double DENDRO_189 = -2 * DENDRO_171 * DENDRO_60 + DENDRO_188;
const double DENDRO_190 = -grad2_0_0_chi[pp];
const double DENDRO_191 = -DENDRO_99;
const double DENDRO_192 = DENDRO_10 * grad_0_chi[pp];
const double DENDRO_193 = DENDRO_10 * grad_1_chi[pp];
const double DENDRO_194 = -DENDRO_66;
const double DENDRO_195 = DENDRO_10 * grad_2_chi[pp];
const double DENDRO_196 = 2 * DENDRO_50;
const double DENDRO_197 = DENDRO_12 * DENDRO_121;
const double DENDRO_198 = DENDRO_18 * grad_2_gt3[pp];
const double DENDRO_199 = DENDRO_14 * grad_1_gt5[pp];
const double DENDRO_200 = DENDRO_11 * DENDRO_155;
const double DENDRO_201 = DENDRO_199 + DENDRO_200;
const double DENDRO_202 = -DENDRO_198 + DENDRO_201;
const double DENDRO_203 = 0.5 * grad_2_gt5[pp];
const double DENDRO_204 = DENDRO_14 * DENDRO_203;
const double DENDRO_205 = 0.5 * grad_1_gt5[pp];
const double DENDRO_206 = 1.0 * grad_2_gt4[pp];
const double DENDRO_207 = -DENDRO_206;
const double DENDRO_208 = DENDRO_205 + DENDRO_207;
const double DENDRO_209 =
    -DENDRO_11 * DENDRO_177 + DENDRO_18 * DENDRO_208 + DENDRO_204;
const double DENDRO_210 = DENDRO_20 * DENDRO_209;
const double DENDRO_211 = 0.5 * grad_1_gt3[pp];
const double DENDRO_212 = DENDRO_18 * DENDRO_211;
const double DENDRO_213 = DENDRO_11 * DENDRO_171;
const double DENDRO_214 = 1.0 * grad_1_gt4[pp];
const double DENDRO_215 = 0.5 * grad_2_gt3[pp];
const double DENDRO_216 = DENDRO_214 - DENDRO_215;
const double DENDRO_217 = -DENDRO_14 * DENDRO_216 + DENDRO_212 + DENDRO_213;
const double DENDRO_218 = -DENDRO_217;
const double DENDRO_219 = DENDRO_18 * DENDRO_218;
const double DENDRO_220 = DENDRO_133 + DENDRO_210 + DENDRO_219;
const double DENDRO_221 =
    DENDRO_88 * (-1.0 * DENDRO_104 * DENDRO_11 - 1.0 * DENDRO_14 * DENDRO_202 +
                 DENDRO_197 + DENDRO_220);
const double DENDRO_222 = 2.0 * grad_1_gt0[pp];
const double DENDRO_223 = -DENDRO_98;
const double DENDRO_224 = DENDRO_12 * DENDRO_223;
const double DENDRO_225 = -DENDRO_11 * grad_2_gt3[pp] +
                          DENDRO_12 * grad_1_gt5[pp] + DENDRO_155 * DENDRO_16;
const double DENDRO_226 = -DENDRO_225;
const double DENDRO_227 = -DENDRO_93;
const double DENDRO_228 =
    DENDRO_11 * DENDRO_208 + DENDRO_12 * DENDRO_203 - DENDRO_16 * DENDRO_177;
const double DENDRO_229 = -DENDRO_228;
const double DENDRO_230 = DENDRO_20 * DENDRO_229;
const double DENDRO_231 = DENDRO_11 * DENDRO_211;
const double DENDRO_232 =
    -DENDRO_12 * DENDRO_216 + DENDRO_16 * DENDRO_171 + DENDRO_231;
const double DENDRO_233 = DENDRO_18 * DENDRO_232;
const double DENDRO_234 = DENDRO_16 * DENDRO_191;
const double DENDRO_235 = DENDRO_230 + DENDRO_233 + DENDRO_234;
const double DENDRO_236 =
    DENDRO_88 * (-1.0 * DENDRO_11 * DENDRO_227 - 1.0 * DENDRO_14 * DENDRO_226 +
                 DENDRO_224 + DENDRO_235);
const double DENDRO_237 = 2.0 * grad_0_gt0[pp];
const double DENDRO_238 = DENDRO_12 * DENDRO_180;
const double DENDRO_239 = DENDRO_20 * grad_1_gt5[pp];
const double DENDRO_240 = DENDRO_12 * DENDRO_155;
const double DENDRO_241 = -DENDRO_14 * grad_2_gt3[pp] + DENDRO_239 + DENDRO_240;
const double DENDRO_242 = -DENDRO_241;
const double DENDRO_243 = DENDRO_20 * DENDRO_203;
const double DENDRO_244 = DENDRO_14 * DENDRO_208;
const double DENDRO_245 = -DENDRO_12 * DENDRO_177 + DENDRO_243 + DENDRO_244;
const double DENDRO_246 = -DENDRO_245;
const double DENDRO_247 = DENDRO_20 * DENDRO_246;
const double DENDRO_248 = DENDRO_14 * DENDRO_211;
const double DENDRO_249 =
    DENDRO_12 * DENDRO_171 - DENDRO_20 * DENDRO_216 + DENDRO_248;
const double DENDRO_250 = DENDRO_18 * DENDRO_249;
const double DENDRO_251 = DENDRO_16 * DENDRO_194;
const double DENDRO_252 = DENDRO_247 + DENDRO_250 + DENDRO_251;
const double DENDRO_253 =
    DENDRO_88 * (-1.0 * DENDRO_11 * DENDRO_178 - 1.0 * DENDRO_14 * DENDRO_242 +
                 DENDRO_238 + DENDRO_252);
const double DENDRO_254 = 2.0 * grad_2_gt0[pp];
const double DENDRO_255 = 3 * DENDRO_50;
const double DENDRO_256 = grad_0_chi[pp] * grad_1_chi[pp];
const double DENDRO_257 = 2 * DENDRO_11;
const double DENDRO_258 = DENDRO_255 * grad_2_chi[pp];
const double DENDRO_259 = 2 * DENDRO_12;
const double DENDRO_260 = 2 * DENDRO_14;
const double DENDRO_261 = DENDRO_104 * DENDRO_11 + DENDRO_14 * DENDRO_202 -
                          1.0 * DENDRO_197 - DENDRO_220;
const double DENDRO_262 = DENDRO_11 * DENDRO_227 + DENDRO_14 * DENDRO_226 -
                          1.0 * DENDRO_224 - DENDRO_235;
const double DENDRO_263 = DENDRO_11 * DENDRO_178 + DENDRO_14 * DENDRO_242 -
                          1.0 * DENDRO_238 - DENDRO_252;
const double DENDRO_264 =
    DENDRO_16 * (-DENDRO_15 * DENDRO_255 + 2 * grad2_0_0_chi[pp]) +
    DENDRO_18 * (-DENDRO_17 * DENDRO_255 + 2 * grad2_1_1_chi[pp]) +
    2 * DENDRO_192 * DENDRO_262 + 2 * DENDRO_193 * DENDRO_261 +
    2 * DENDRO_195 * DENDRO_263 +
    DENDRO_20 * (-DENDRO_19 * DENDRO_255 + 2 * grad2_2_2_chi[pp]) -
    DENDRO_257 * (-DENDRO_255 * DENDRO_256 + 2 * grad2_0_1_chi[pp]) +
    DENDRO_259 * (-DENDRO_258 * grad_0_chi[pp] + 2 * grad2_0_2_chi[pp]) -
    DENDRO_260 * (-DENDRO_258 * grad_1_chi[pp] + 2 * grad2_1_2_chi[pp]);
const double DENDRO_265 = DENDRO_50 * DENDRO_75;
const double DENDRO_266 = 3 * alpha[pp];
const double DENDRO_267 = DENDRO_50 * gt5[pp];
const double DENDRO_268 = DENDRO_209 + DENDRO_267 * DENDRO_49;
const double DENDRO_269 = 4 * grad_1_alpha[pp];
const double DENDRO_270 = DENDRO_10 * DENDRO_269;
const double DENDRO_271 = DENDRO_228 - 0.5 * DENDRO_50 * DENDRO_77 * gt5[pp];
const double DENDRO_272 = 4 * grad_0_alpha[pp];
const double DENDRO_273 = DENDRO_10 * DENDRO_272;
const double DENDRO_274 = 1.0 * grad_2_chi[pp];
const double DENDRO_275 = DENDRO_10 * gt5[pp];
const double DENDRO_276 = 0.5 * DENDRO_65;
const double DENDRO_277 = -DENDRO_10 * DENDRO_12 * DENDRO_177 +
                          DENDRO_10 * DENDRO_243 + DENDRO_10 * DENDRO_244 +
                          DENDRO_50 * (DENDRO_274 - DENDRO_275 * DENDRO_276);
const double DENDRO_278 = 4 * grad_2_alpha[pp];
const double DENDRO_279 = 4 * gt2[pp];
const double DENDRO_280 = 4 * gt4[pp];
const double DENDRO_281 = DENDRO_19 * DENDRO_81;
const double DENDRO_282 = DENDRO_84 * grad2_0_1_gt5[pp];
const double DENDRO_283 = DENDRO_86 * grad2_1_2_gt5[pp];
const double DENDRO_284 = DENDRO_10 * DENDRO_12;
const double DENDRO_285 = 4 * DENDRO_284;
const double DENDRO_286 = DENDRO_10 * DENDRO_16;
const double DENDRO_287 = 2.0 * DENDRO_286;
const double DENDRO_288 = DENDRO_10 * DENDRO_18;
const double DENDRO_289 = 2.0 * DENDRO_288;
const double DENDRO_290 = DENDRO_10 * DENDRO_20;
const double DENDRO_291 = 2.0 * DENDRO_290;
const double DENDRO_292 = DENDRO_180 * grad_0_gt5[pp];
const double DENDRO_293 = DENDRO_16 * DENDRO_88;
const double DENDRO_294 = 3.0 * DENDRO_293;
const double DENDRO_295 = DENDRO_242 * grad_1_gt5[pp];
const double DENDRO_296 = 3.0 * DENDRO_108;
const double DENDRO_297 = 6.0 * DENDRO_88;
const double DENDRO_298 = 0.25 * grad_2_gt3[pp];
const double DENDRO_299 = DENDRO_109 * DENDRO_202 * (DENDRO_214 - DENDRO_298);
const double DENDRO_300 = -DENDRO_141 + DENDRO_55;
const double DENDRO_301 = 4 * DENDRO_293;
const double DENDRO_302 = 0.75 * grad_0_gt4[pp];
const double DENDRO_303 = -DENDRO_128;
const double DENDRO_304 =
    DENDRO_121 * DENDRO_301 * (DENDRO_124 + DENDRO_302 + DENDRO_303);
const double DENDRO_305 = DENDRO_122 + DENDRO_127 + DENDRO_303;
const double DENDRO_306 = DENDRO_115 + DENDRO_176;
const double DENDRO_307 = DENDRO_134 * DENDRO_210 * (DENDRO_205 + DENDRO_206);
const double DENDRO_308 = DENDRO_114 * DENDRO_242;
const double DENDRO_309 = 0.25 * grad_1_gt5[pp];
const double DENDRO_310 = DENDRO_180 * DENDRO_309;
const double DENDRO_311 = DENDRO_202 * grad_1_gt5[pp];
const double DENDRO_312 = DENDRO_209 * grad_2_gt3[pp];
const double DENDRO_313 = DENDRO_311 + DENDRO_312;
const double DENDRO_314 = 2.0 * DENDRO_139;
const double DENDRO_315 = DENDRO_229 * grad_2_gt0[pp];
const double DENDRO_316 = DENDRO_180 * grad_2_gt5[pp];
const double DENDRO_317 = DENDRO_246 * grad_0_gt5[pp];
const double DENDRO_318 = DENDRO_242 * grad_2_gt5[pp];
const double DENDRO_319 = DENDRO_246 * grad_1_gt5[pp];
const double DENDRO_320 = DENDRO_121 * DENDRO_298;
const double DENDRO_321 = 0.5 * DENDRO_91;
const double DENDRO_322 = DENDRO_202 * DENDRO_321 + DENDRO_320;
const double DENDRO_323 = 0.25 * DENDRO_316;
const double DENDRO_324 = 0.25 * DENDRO_318;
const double DENDRO_325 = DENDRO_141 * DENDRO_226;
const double DENDRO_326 = DENDRO_121 * grad_1_gt5[pp];
const double DENDRO_327 = DENDRO_209 * DENDRO_96;
const double DENDRO_328 = DENDRO_166 * (DENDRO_326 + DENDRO_327);
const double DENDRO_329 = DENDRO_121 * DENDRO_216;
const double DENDRO_330 = DENDRO_202 * DENDRO_96;
const double DENDRO_331 = 0.25 * DENDRO_330;
const double DENDRO_332 = DENDRO_329 + DENDRO_331;
const double DENDRO_333 = DENDRO_155 * DENDRO_229;
const double DENDRO_334 = DENDRO_226 * DENDRO_57;
const double DENDRO_335 = DENDRO_155 * DENDRO_223;
const double DENDRO_336 = 0.25 * DENDRO_335;
const double DENDRO_337 = 0.5 * DENDRO_208;
const double DENDRO_338 = DENDRO_121 * DENDRO_337;
const double DENDRO_339 = 0.5 * DENDRO_177;
const double DENDRO_340 = DENDRO_226 * DENDRO_339;
const double DENDRO_341 = -DENDRO_202 * DENDRO_337;
const double DENDRO_342 = 2 * DENDRO_209 * DENDRO_216 + DENDRO_341;
const double DENDRO_343 = -DENDRO_223 * DENDRO_339;
const double DENDRO_344 = 2 * DENDRO_57;
const double DENDRO_345 = -grad2_2_2_chi[pp];
const double DENDRO_346 = -DENDRO_12;
const double DENDRO_347 = -DENDRO_177;
const double DENDRO_348 = -DENDRO_16;
const double DENDRO_349 = -DENDRO_208;
const double DENDRO_350 =
    DENDRO_11 * DENDRO_349 + DENDRO_203 * DENDRO_346 + DENDRO_347 * DENDRO_348;
const double DENDRO_351 = -DENDRO_18;
const double DENDRO_352 =
    DENDRO_11 * DENDRO_347 + DENDRO_204 + DENDRO_349 * DENDRO_351;
const double DENDRO_353 = -DENDRO_20;
const double DENDRO_354 = DENDRO_203 * DENDRO_353;
const double DENDRO_355 = DENDRO_346 * DENDRO_347;
const double DENDRO_356 = DENDRO_14 * DENDRO_349;
const double DENDRO_357 = DENDRO_261 * DENDRO_88;
const double DENDRO_358 = 2.0 * DENDRO_357;
const double DENDRO_359 = DENDRO_262 * DENDRO_88;
const double DENDRO_360 = 2.0 * DENDRO_359;
const double DENDRO_361 = DENDRO_263 * DENDRO_88;
const double DENDRO_362 = 2.0 * DENDRO_361;
const double DENDRO_363 = -DENDRO_264;
const double DENDRO_364 = DENDRO_363 * DENDRO_50;
const double DENDRO_365 = DENDRO_50 * gt3[pp];
const double DENDRO_366 = DENDRO_10 * DENDRO_278;
const double DENDRO_367 = 1.0 * grad_1_chi[pp];
const double DENDRO_368 = DENDRO_10 * gt3[pp];
const double DENDRO_369 = -DENDRO_10 * DENDRO_14 * DENDRO_216 +
                          DENDRO_10 * DENDRO_212 + DENDRO_10 * DENDRO_213 +
                          DENDRO_50 * (DENDRO_367 - DENDRO_368 * DENDRO_49);
const double DENDRO_370 = -grad2_1_1_chi[pp];
const double DENDRO_371 = -DENDRO_171;
const double DENDRO_372 =
    DENDRO_216 * DENDRO_346 + DENDRO_231 + DENDRO_348 * DENDRO_371;
const double DENDRO_373 = DENDRO_211 * DENDRO_351;
const double DENDRO_374 = DENDRO_11 * DENDRO_371;
const double DENDRO_375 = DENDRO_14 * DENDRO_216;
const double DENDRO_376 =
    DENDRO_216 * DENDRO_353 + DENDRO_248 + DENDRO_346 * DENDRO_371;
const double DENDRO_377 = DENDRO_226 * DENDRO_54;
const double DENDRO_378 = DENDRO_155 * DENDRO_227;
const double DENDRO_379 = 0.25 * DENDRO_378;
const double DENDRO_380 = DENDRO_137 * DENDRO_226;
const double DENDRO_381 = 0.5 * DENDRO_96;
const double DENDRO_382 = DENDRO_178 * DENDRO_309;
const double DENDRO_383 = DENDRO_178 * DENDRO_208;
const double DENDRO_384 = DENDRO_242 * DENDRO_91;
const double DENDRO_385 = -0.25 * DENDRO_384;
const double DENDRO_386 = 0.5 * DENDRO_171;
const double DENDRO_387 = DENDRO_226 * DENDRO_386;
const double DENDRO_388 = 0.5 * DENDRO_216;
const double DENDRO_389 = DENDRO_242 * DENDRO_388;
const double DENDRO_390 = 2 * DENDRO_208 * DENDRO_249;
const double DENDRO_391 = DENDRO_202 * grad_1_gt3[pp];
const double DENDRO_392 = 0.25 * DENDRO_391;
const double DENDRO_393 = DENDRO_218 * grad_2_gt3[pp];
const double DENDRO_394 = DENDRO_178 * DENDRO_388;
const double DENDRO_395 = -DENDRO_227 * DENDRO_386;
const double DENDRO_396 = 2 * DENDRO_232 * DENDRO_54;
const double DENDRO_397 = DENDRO_104 * grad_1_gt3[pp];
const double DENDRO_398 = 0.25 * DENDRO_397;
const double DENDRO_399 = DENDRO_218 * grad_0_gt3[pp];
const double DENDRO_400 = DENDRO_232 * grad_1_gt0[pp];
const double DENDRO_401 = DENDRO_155 * DENDRO_232;
const double DENDRO_402 = DENDRO_249 * grad_1_gt5[pp];
const double DENDRO_403 = DENDRO_249 * DENDRO_91;
const double DENDRO_404 = DENDRO_207 + DENDRO_309;
const double DENDRO_405 = -DENDRO_124;
const double DENDRO_406 = DENDRO_119 * (DENDRO_122 + DENDRO_125 + DENDRO_405);
const double DENDRO_407 = DENDRO_301 * (-DENDRO_137 + DENDRO_52);
const double DENDRO_408 = DENDRO_301 * (DENDRO_128 + DENDRO_302 + DENDRO_405);
const double DENDRO_409 = 4 * gt1[pp];
const double DENDRO_410 = DENDRO_105 * DENDRO_202;
const double DENDRO_411 = DENDRO_104 * DENDRO_298;
const double DENDRO_412 = DENDRO_104 * grad_0_gt3[pp];
const double DENDRO_413 = DENDRO_202 * grad_2_gt3[pp];
const double DENDRO_414 = 3.0 * DENDRO_118;
const double DENDRO_415 =
    -DENDRO_134 * DENDRO_233 * (DENDRO_106 + DENDRO_170) -
    DENDRO_134 * DENDRO_250 * (DENDRO_214 + DENDRO_215) -
    DENDRO_17 * DENDRO_81 -
    DENDRO_185 * (DENDRO_104 * DENDRO_215 + DENDRO_410) -
    DENDRO_185 * (DENDRO_170 * DENDRO_202 + DENDRO_411) +
    DENDRO_280 * grad_1_Gt2[pp] + DENDRO_285 * grad2_0_2_gt3[pp] +
    DENDRO_287 * grad2_0_0_gt3[pp] + DENDRO_289 * grad2_1_1_gt3[pp] +
    DENDRO_291 * grad2_2_2_gt3[pp] - DENDRO_294 * DENDRO_412 +
    DENDRO_409 * grad_1_Gt0[pp] - DENDRO_413 * DENDRO_414 -
    DENDRO_84 * grad2_0_1_gt3[pp] - DENDRO_86 * grad2_1_2_gt3[pp] +
    4 * grad_1_Gt1[pp] * gt3[pp];
const double DENDRO_416 = DENDRO_227 * grad_1_gt0[pp];
const double DENDRO_417 = DENDRO_223 * grad_2_gt0[pp];
const double DENDRO_418 = DENDRO_137 * DENDRO_223;
const double DENDRO_419 = DENDRO_141 * DENDRO_227;
const double DENDRO_420 = DENDRO_227 * grad_0_gt0[pp];
const double DENDRO_421 = DENDRO_191 * grad_1_gt0[pp];
const double DENDRO_422 = DENDRO_223 * grad_0_gt0[pp];
const double DENDRO_423 = DENDRO_191 * grad_2_gt0[pp];
const double DENDRO_424 = DENDRO_194 * grad_0_gt5[pp];
const double DENDRO_425 = 0.25 * DENDRO_420;
const double DENDRO_426 = 0.25 * DENDRO_422;
const double DENDRO_427 = DENDRO_114 * DENDRO_178;
const double DENDRO_428 = DENDRO_194 * DENDRO_91;
const double DENDRO_429 = DENDRO_178 * DENDRO_187;
const double DENDRO_430 = DENDRO_180 * DENDRO_187;
const double DENDRO_431 = DENDRO_348 * DENDRO_58;
const double DENDRO_432 = DENDRO_346 * DENDRO_57;
const double DENDRO_433 = DENDRO_351 * DENDRO_54 + DENDRO_59;
const double DENDRO_434 =
    DENDRO_14 * DENDRO_54 + DENDRO_346 * DENDRO_58 + DENDRO_353 * DENDRO_57;
const double DENDRO_435 = DENDRO_202 * grad_1_gt0[pp];
const double DENDRO_436 = DENDRO_121 * grad_0_gt3[pp];
const double DENDRO_437 = DENDRO_104 * DENDRO_91;
const double DENDRO_438 = DENDRO_436 + DENDRO_437;
const double DENDRO_439 = DENDRO_242 * grad_2_gt0[pp];
const double DENDRO_440 = DENDRO_178 * grad_0_gt5[pp];
const double DENDRO_441 = DENDRO_181 + DENDRO_440;
const double DENDRO_442 = DENDRO_114 * DENDRO_223;
const double DENDRO_443 = -DENDRO_177 * DENDRO_246;
const double DENDRO_444 = DENDRO_156 * DENDRO_209 + 0.25 * DENDRO_326;
const double DENDRO_445 = DENDRO_227 * DENDRO_91;
const double DENDRO_446 = 0.25 * DENDRO_445;
const double DENDRO_447 = DENDRO_104 * DENDRO_388;
const double DENDRO_448 = -DENDRO_202 * DENDRO_386;
const double DENDRO_449 = DENDRO_447 + DENDRO_448;
const double DENDRO_450 = DENDRO_191 * DENDRO_57;
const double DENDRO_451 = 0.25 * DENDRO_163 + DENDRO_321 * DENDRO_60;
const double DENDRO_452 = DENDRO_141 * DENDRO_180;
const double DENDRO_453 = 0.25 * DENDRO_121;
const double DENDRO_454 = DENDRO_155 * DENDRO_453 + DENDRO_205 * DENDRO_60;
const double DENDRO_455 = DENDRO_194 * DENDRO_203;
const double DENDRO_456 = -DENDRO_180 * DENDRO_339;
const double DENDRO_457 = DENDRO_141 * DENDRO_223;
const double DENDRO_458 = DENDRO_229 * DENDRO_58;
const double DENDRO_459 = DENDRO_187 * DENDRO_223 + DENDRO_458;
const double DENDRO_460 = DENDRO_453 * DENDRO_96;
const double DENDRO_461 = DENDRO_209 * DENDRO_53 + DENDRO_453 * DENDRO_91;
const double DENDRO_462 = DENDRO_171 * DENDRO_209;
const double DENDRO_463 = DENDRO_229 * DENDRO_53;
const double DENDRO_464 = DENDRO_114 * DENDRO_227 + DENDRO_463;
const double DENDRO_465 = DENDRO_156 * DENDRO_246;
const double DENDRO_466 = DENDRO_308 + DENDRO_310;
const double DENDRO_467 = DENDRO_155 * DENDRO_202;
const double DENDRO_468 = 0.25 * DENDRO_467;
const double DENDRO_469 = DENDRO_170 * DENDRO_209;
const double DENDRO_470 = DENDRO_104 * DENDRO_309 + DENDRO_469;
const double DENDRO_471 = 0.25 * DENDRO_91;
const double DENDRO_472 = DENDRO_223 * DENDRO_471 + DENDRO_463;
const double DENDRO_473 = DENDRO_242 * DENDRO_339;
const double DENDRO_474 = -DENDRO_473;
const double DENDRO_475 = 0.25 * DENDRO_178;
const double DENDRO_476 = DENDRO_475 * grad_2_gt5[pp];
const double DENDRO_477 = DENDRO_246 * DENDRO_321 + DENDRO_476;
const double DENDRO_478 = DENDRO_183 * DENDRO_202;
const double DENDRO_479 = -0.5 * DENDRO_172 + DENDRO_216 * DENDRO_60;
const double DENDRO_480 = 0.25 * DENDRO_226;
const double DENDRO_481 = DENDRO_480 * grad_0_gt0[pp];
const double DENDRO_482 = DENDRO_187 * DENDRO_227;
const double DENDRO_483 = DENDRO_481 + DENDRO_482;
const double DENDRO_484 = DENDRO_191 * DENDRO_321 + DENDRO_418;
const double DENDRO_485 = DENDRO_194 * DENDRO_205;
const double DENDRO_486 = 0.25 * DENDRO_155;
const double DENDRO_487 = DENDRO_180 * DENDRO_486;
const double DENDRO_488 = DENDRO_155 * DENDRO_242;
const double DENDRO_489 = DENDRO_178 * grad_1_gt5[pp];
const double DENDRO_490 = DENDRO_384 + DENDRO_489;
const double DENDRO_491 = 1.0 * DENDRO_108;
const double DENDRO_492 = -grad2_0_2_chi[pp];
const double DENDRO_493 = DENDRO_346 * grad_0_gt5[pp];
const double DENDRO_494 = DENDRO_348 * grad_2_gt0[pp];
const double DENDRO_495 = 0.5 * DENDRO_192;
const double DENDRO_496 = DENDRO_120 + DENDRO_351 * DENDRO_96;
const double DENDRO_497 = 0.5 * DENDRO_193;
const double DENDRO_498 = DENDRO_353 * grad_0_gt5[pp];
const double DENDRO_499 = DENDRO_346 * grad_2_gt0[pp];
const double DENDRO_500 = DENDRO_14 * DENDRO_96;
const double DENDRO_501 = 0.5 * DENDRO_195;
const double DENDRO_502 = 2.0 * gt2[pp];
const double DENDRO_503 = 2.0 * gt4[pp];
const double DENDRO_504 = 2.0 * gt5[pp];
const double DENDRO_505 = 2.0 * grad_2_Gt0[pp];
const double DENDRO_506 = 2.0 * grad_2_Gt1[pp];
const double DENDRO_507 = DENDRO_81 * grad_2_chi[pp];
const double DENDRO_508 = DENDRO_507 * grad_0_chi[pp];
const double DENDRO_509 = DENDRO_84 * grad2_0_1_gt2[pp];
const double DENDRO_510 = DENDRO_86 * grad2_1_2_gt2[pp];
const double DENDRO_511 = DENDRO_10 * gt2[pp];
const double DENDRO_512 =
    -DENDRO_196 *
        (DENDRO_492 + DENDRO_495 * (DENDRO_493 + DENDRO_494 + DENDRO_97) +
         DENDRO_496 * DENDRO_497 +
         DENDRO_501 * (DENDRO_498 + DENDRO_499 + DENDRO_500)) +
    DENDRO_285 * grad2_0_2_gt2[pp] + DENDRO_287 * grad2_0_0_gt2[pp] +
    DENDRO_289 * grad2_1_1_gt2[pp] + DENDRO_291 * grad2_2_2_gt2[pp] +
    DENDRO_358 * grad_1_gt2[pp] + DENDRO_360 * grad_0_gt2[pp] +
    DENDRO_362 * grad_2_gt2[pp] + DENDRO_364 * DENDRO_511 +
    DENDRO_502 * grad_0_Gt0[pp] + DENDRO_502 * grad_2_Gt2[pp] +
    DENDRO_503 * grad_0_Gt1[pp] + DENDRO_504 * grad_0_Gt2[pp] +
    DENDRO_505 * gt0[pp] + DENDRO_506 * gt1[pp] - DENDRO_508 - DENDRO_509 -
    DENDRO_510;
const double DENDRO_513 =
    -DENDRO_10 * DENDRO_11 * DENDRO_96 + DENDRO_10 * DENDRO_94 +
    DENDRO_10 * DENDRO_95 +
    DENDRO_50 * (-DENDRO_511 * DENDRO_77 + grad_2_chi[pp]);
const double DENDRO_514 = 2.0 * grad_0_alpha[pp];
const double DENDRO_515 =
    DENDRO_10 * DENDRO_111 + DENDRO_10 * DENDRO_112 -
    DENDRO_10 * DENDRO_14 * DENDRO_96 +
    DENDRO_50 * (-DENDRO_511 * DENDRO_65 + grad_0_chi[pp]);
const double DENDRO_516 = 2.0 * grad_2_alpha[pp];
const double DENDRO_517 = 2.0 * grad_1_alpha[pp];
const double DENDRO_518 = DENDRO_50 * gt2[pp];
const double DENDRO_519 = DENDRO_10 * (DENDRO_121 + DENDRO_48 * DENDRO_518);
const double DENDRO_520 = -DENDRO_513 * DENDRO_514 - DENDRO_515 * DENDRO_516 +
                          DENDRO_517 * DENDRO_519 - 4 * grad2_0_2_alpha[pp];
const double DENDRO_521 = DENDRO_226 * grad_2_gt0[pp];
const double DENDRO_522 = DENDRO_227 * grad_0_gt5[pp] + DENDRO_521;
const double DENDRO_523 = DENDRO_121 * grad_2_gt3[pp];
const double DENDRO_524 = DENDRO_104 * grad_1_gt5[pp] + DENDRO_523;
const double DENDRO_525 = DENDRO_467 + DENDRO_524;
const double DENDRO_526 = DENDRO_119 * (-DENDRO_338 + DENDRO_444);
const double DENDRO_527 = 0.25 * DENDRO_488;
const double DENDRO_528 = DENDRO_109 * (DENDRO_411 + DENDRO_449);
const double DENDRO_529 = DENDRO_301 * (0.5 * DENDRO_164 + DENDRO_451);
const double DENDRO_530 = DENDRO_185 * (-DENDRO_208 * DENDRO_60 + DENDRO_461);
const double DENDRO_531 = DENDRO_114 * DENDRO_180 + DENDRO_455;
const double DENDRO_532 = DENDRO_185 * (DENDRO_454 + DENDRO_460);
const double DENDRO_533 = DENDRO_104 * DENDRO_337;
const double DENDRO_534 =
    -0.5 * DENDRO_121 * DENDRO_216 + DENDRO_462 + DENDRO_533;
const double DENDRO_535 = DENDRO_310 + DENDRO_476;
const double DENDRO_536 = -DENDRO_227 * DENDRO_339;
const double DENDRO_537 = DENDRO_174 + DENDRO_479;
const double DENDRO_538 = DENDRO_156 * DENDRO_191;
const double DENDRO_539 = DENDRO_141 * DENDRO_242;
const double DENDRO_540 = DENDRO_427 + DENDRO_485;
const double DENDRO_541 = 0.25 * DENDRO_437;
const double DENDRO_542 = DENDRO_215 * DENDRO_60;
const double DENDRO_543 = DENDRO_137 * DENDRO_202 + DENDRO_542;
const double DENDRO_544 = DENDRO_541 + DENDRO_543;
const double DENDRO_545 = DENDRO_226 * grad_1_gt0[pp];
const double DENDRO_546 = DENDRO_445 + DENDRO_545;
const double DENDRO_547 = DENDRO_223 * grad_0_gt3[pp];
const double DENDRO_548 = DENDRO_180 * grad_2_gt3[pp];
const double DENDRO_549 = 0.25 * DENDRO_311;
const double DENDRO_550 = -DENDRO_208 * DENDRO_246;
const double DENDRO_551 = DENDRO_114 * DENDRO_226 + DENDRO_229 * DENDRO_381;
const double DENDRO_552 = DENDRO_216 * DENDRO_218;
const double DENDRO_553 = DENDRO_232 * DENDRO_321;
const double DENDRO_554 = DENDRO_105 * DENDRO_226 + DENDRO_553;
const double DENDRO_555 = DENDRO_242 * DENDRO_298;
const double DENDRO_556 = DENDRO_183 * DENDRO_223;
const double DENDRO_557 = DENDRO_482 + DENDRO_556;
const double DENDRO_558 = DENDRO_229 * DENDRO_54 + 0.5 * DENDRO_334;
const double DENDRO_559 = DENDRO_223 * DENDRO_96;
const double DENDRO_560 = 0.25 * DENDRO_559;
const double DENDRO_561 = DENDRO_246 * DENDRO_381;
const double DENDRO_562 = DENDRO_202 * DENDRO_471 + DENDRO_469;
const double DENDRO_563 = DENDRO_180 * DENDRO_337;
const double DENDRO_564 = -DENDRO_563;
const double DENDRO_565 = DENDRO_203 * DENDRO_249;
const double DENDRO_566 = -DENDRO_242 * DENDRO_337 + DENDRO_565;
const double DENDRO_567 = DENDRO_176 * DENDRO_232 + DENDRO_480 * DENDRO_96;
const double DENDRO_568 = DENDRO_209 * DENDRO_211;
const double DENDRO_569 = DENDRO_202 * DENDRO_298 + DENDRO_568;
const double DENDRO_570 = DENDRO_202 * DENDRO_388;
const double DENDRO_571 = DENDRO_155 * DENDRO_480;
const double DENDRO_572 = DENDRO_170 * DENDRO_229 + DENDRO_226 * DENDRO_471;
const double DENDRO_573 = -DENDRO_223 * DENDRO_386;
const double DENDRO_574 = DENDRO_232 * DENDRO_57 + 0.5 * DENDRO_377;
const double DENDRO_575 = DENDRO_453 * grad_1_gt3[pp];
const double DENDRO_576 = DENDRO_410 + DENDRO_575;
const double DENDRO_577 = DENDRO_218 * DENDRO_321;
const double DENDRO_578 = 0.25 * DENDRO_96;
const double DENDRO_579 = DENDRO_176 * DENDRO_249;
const double DENDRO_580 = DENDRO_242 * DENDRO_578 + DENDRO_579;
const double DENDRO_581 = DENDRO_180 * DENDRO_96;
const double DENDRO_582 = 1.0 * DENDRO_293;
const double DENDRO_583 = -grad2_1_2_chi[pp];
const double DENDRO_584 = DENDRO_11 * grad_2_gt3[pp] + DENDRO_155 * DENDRO_348 +
                          DENDRO_346 * grad_1_gt5[pp];
const double DENDRO_585 = DENDRO_351 * grad_2_gt3[pp];
const double DENDRO_586 = DENDRO_353 * grad_1_gt5[pp];
const double DENDRO_587 = DENDRO_14 * grad_2_gt3[pp];
const double DENDRO_588 = DENDRO_155 * DENDRO_346;
const double DENDRO_589 = DENDRO_10 * gt4[pp];
const double DENDRO_590 =
    DENDRO_285 * grad2_0_2_gt4[pp] + DENDRO_287 * grad2_0_0_gt4[pp] +
    DENDRO_289 * grad2_1_1_gt4[pp] + DENDRO_291 * grad2_2_2_gt4[pp] +
    DENDRO_502 * grad_1_Gt0[pp] + DENDRO_503 * grad_1_Gt1[pp] +
    DENDRO_503 * grad_2_Gt2[pp] + DENDRO_504 * grad_1_Gt2[pp] +
    DENDRO_505 * gt1[pp] + DENDRO_506 * gt3[pp] - DENDRO_507 * grad_1_chi[pp] -
    DENDRO_84 * grad2_0_1_gt4[pp] - DENDRO_86 * grad2_1_2_gt4[pp];
const double DENDRO_591 =
    -DENDRO_196 *
        (DENDRO_495 * DENDRO_584 + DENDRO_497 * (DENDRO_201 + DENDRO_585) +
         DENDRO_501 * (DENDRO_586 + DENDRO_587 + DENDRO_588) + DENDRO_583) +
    DENDRO_358 * grad_1_gt4[pp] + DENDRO_360 * grad_0_gt4[pp] +
    DENDRO_362 * grad_2_gt4[pp] + DENDRO_364 * DENDRO_589 + DENDRO_590;
const double DENDRO_592 = DENDRO_10 * DENDRO_199 + DENDRO_10 * DENDRO_200;
const double DENDRO_593 =
    -DENDRO_10 * DENDRO_198 -
    DENDRO_50 * (-DENDRO_48 * DENDRO_589 + grad_2_chi[pp]) + DENDRO_592;
const double DENDRO_594 =
    -DENDRO_10 * DENDRO_14 * grad_2_gt3[pp] + DENDRO_10 * DENDRO_239 +
    DENDRO_10 * DENDRO_240 +
    DENDRO_50 * (-DENDRO_589 * DENDRO_65 + grad_1_chi[pp]);
const double DENDRO_595 = DENDRO_225 - DENDRO_50 * DENDRO_77 * gt4[pp];
const double DENDRO_596 = -DENDRO_10 * DENDRO_514 * DENDRO_595 -
                          DENDRO_516 * DENDRO_594 + DENDRO_517 * DENDRO_593 -
                          4 * grad2_1_2_alpha[pp];
const double DENDRO_597 = 1.0 * DENDRO_402;
const double DENDRO_598 = 0.5 * DENDRO_401;
const double DENDRO_599 = 0.25 * DENDRO_581;
const double DENDRO_600 = DENDRO_308 + DENDRO_476;
const double DENDRO_601 = DENDRO_177 * DENDRO_232;
const double DENDRO_602 = DENDRO_568 + DENDRO_570;
const double DENDRO_603 = DENDRO_242 * DENDRO_309;
const double DENDRO_604 = DENDRO_410 + DENDRO_411;
const double DENDRO_605 = DENDRO_232 * DENDRO_56;
const double DENDRO_606 = DENDRO_105 * DENDRO_223 + DENDRO_605;
const double DENDRO_607 = DENDRO_218 * DENDRO_381 + DENDRO_575;
const double DENDRO_608 = DENDRO_180 * DENDRO_298 + DENDRO_579;
const double DENDRO_609 = 1.0 * DENDRO_165;
const double DENDRO_610 =
    -DENDRO_119 * (DENDRO_209 * DENDRO_215 + DENDRO_341 + DENDRO_549) -
    DENDRO_185 * (-DENDRO_533 + DENDRO_562) -
    DENDRO_582 * (DENDRO_173 + DENDRO_438) -
    DENDRO_609 * (DENDRO_330 + DENDRO_524);
const double DENDRO_611 = DENDRO_473 + DENDRO_563;
const double DENDRO_612 = DENDRO_105 * DENDRO_227;
const double DENDRO_613 = -DENDRO_171 * DENDRO_218;
const double DENDRO_614 = DENDRO_156 * DENDRO_249 + DENDRO_178 * DENDRO_298;
const double DENDRO_615 = DENDRO_191 * DENDRO_54;
const double DENDRO_616 = 0.25 * DENDRO_149;
const double DENDRO_617 = DENDRO_141 * DENDRO_178 + DENDRO_194 * DENDRO_381;
const double DENDRO_618 = 0.5 * DENDRO_179 + DENDRO_194 * DENDRO_208;
const double DENDRO_619 = DENDRO_481 + DENDRO_556;
const double DENDRO_620 = DENDRO_191 * DENDRO_381 + DENDRO_419;
const double DENDRO_621 = DENDRO_104 * DENDRO_486 + DENDRO_542;
const double DENDRO_622 = DENDRO_177 * DENDRO_249 + 0.5 * DENDRO_383;
const double DENDRO_623 = DENDRO_156 * DENDRO_218;
const double DENDRO_624 = DENDRO_227 * DENDRO_578 + DENDRO_605;
const double DENDRO_625 = -DENDRO_104 * DENDRO_386;
const double DENDRO_626 = DENDRO_211 * DENDRO_60;
const double DENDRO_627 = DENDRO_155 * DENDRO_475 + DENDRO_194 * DENDRO_215;
const double DENDRO_628 = DENDRO_137 * DENDRO_227;
const double DENDRO_629 = DENDRO_232 * DENDRO_58;
const double DENDRO_630 = DENDRO_183 * DENDRO_227 + DENDRO_629;
const double DENDRO_631 = DENDRO_475 * DENDRO_91;
const double DENDRO_632 = DENDRO_249 * DENDRO_56 + DENDRO_475 * DENDRO_96;
const double DENDRO_633 = 1.0 * DENDRO_118;
const double DENDRO_634 = -grad2_0_1_chi[pp];
const double DENDRO_635 = DENDRO_348 * grad_1_gt0[pp];
const double DENDRO_636 = DENDRO_346 * DENDRO_91;
const double DENDRO_637 = DENDRO_351 * grad_0_gt3[pp];
const double DENDRO_638 = DENDRO_14 * grad_0_gt3[pp];
const double DENDRO_639 =
    DENDRO_346 * grad_1_gt0[pp] + DENDRO_353 * DENDRO_91 + DENDRO_638;
const double DENDRO_640 = 2.0 * gt1[pp];
const double DENDRO_641 = DENDRO_640 * grad_0_Gt0[pp];
const double DENDRO_642 = 2.0 * grad_0_Gt1[pp] * gt3[pp];
const double DENDRO_643 = DENDRO_503 * grad_0_Gt2[pp];
const double DENDRO_644 = 2.0 * grad_1_Gt0[pp] * gt0[pp];
const double DENDRO_645 = DENDRO_640 * grad_1_Gt1[pp];
const double DENDRO_646 = DENDRO_502 * grad_1_Gt2[pp];
const double DENDRO_647 = -DENDRO_256 * DENDRO_81;
const double DENDRO_648 = DENDRO_285 * grad2_0_2_gt1[pp];
const double DENDRO_649 = DENDRO_287 * grad2_0_0_gt1[pp];
const double DENDRO_650 = DENDRO_289 * grad2_1_1_gt1[pp];
const double DENDRO_651 = DENDRO_291 * grad2_2_2_gt1[pp];
const double DENDRO_652 = -DENDRO_84 * grad2_0_1_gt1[pp];
const double DENDRO_653 = -DENDRO_86 * grad2_1_2_gt1[pp];
const double DENDRO_654 = DENDRO_10 * gt1[pp];
const double DENDRO_655 =
    -DENDRO_196 * (DENDRO_495 * (DENDRO_635 + DENDRO_636 + DENDRO_89) +
                   DENDRO_497 * (DENDRO_103 + DENDRO_637) +
                   DENDRO_501 * DENDRO_639 + DENDRO_634) +
    DENDRO_358 * grad_1_gt1[pp] + DENDRO_360 * grad_0_gt1[pp] +
    DENDRO_362 * grad_2_gt1[pp] + DENDRO_364 * DENDRO_654 + DENDRO_641 +
    DENDRO_642 + DENDRO_643 + DENDRO_644 + DENDRO_645 + DENDRO_646 +
    DENDRO_647 + DENDRO_648 + DENDRO_649 + DENDRO_650 + DENDRO_651 +
    DENDRO_652 + DENDRO_653;
const double DENDRO_656 =
    -DENDRO_10 * DENDRO_11 * grad_0_gt3[pp] + DENDRO_10 * DENDRO_90 +
    DENDRO_10 * DENDRO_92 +
    DENDRO_50 * (-DENDRO_654 * DENDRO_77 + grad_1_chi[pp]);
const double DENDRO_657 =
    DENDRO_10 * DENDRO_100 - DENDRO_10 * DENDRO_11 * grad_1_gt0[pp] -
    DENDRO_10 * DENDRO_14 * DENDRO_91 +
    DENDRO_50 * (-DENDRO_48 * DENDRO_654 + grad_0_chi[pp]);
const double DENDRO_658 = DENDRO_50 * gt1[pp];
const double DENDRO_659 =
    DENDRO_10 * DENDRO_516 *
        (-DENDRO_130 - DENDRO_131 + DENDRO_638 + DENDRO_65 * DENDRO_658) -
    DENDRO_514 * DENDRO_656 - DENDRO_517 * DENDRO_657 - 4 * grad2_0_1_alpha[pp];
const double DENDRO_660 = -0.25 * DENDRO_180 * grad_1_gt5[pp] + DENDRO_611;
const double DENDRO_661 = 0.5 * DENDRO_397;
const double DENDRO_662 = DENDRO_232 * DENDRO_53;
const double DENDRO_663 = DENDRO_182 + DENDRO_618;
const double DENDRO_664 = -0.5 * DENDRO_180 * DENDRO_216 + DENDRO_622;
const double DENDRO_665 = DENDRO_411 + DENDRO_575;
const double DENDRO_666 = DENDRO_104 * DENDRO_105 + DENDRO_626;
const double DENDRO_667 = -DENDRO_119 * (DENDRO_121 * DENDRO_215 + DENDRO_468) +
                          DENDRO_140 * (DENDRO_448 + DENDRO_665) +
                          DENDRO_161 * (DENDRO_625 + DENDRO_666) -
                          DENDRO_185 * (DENDRO_154 + DENDRO_543) -
                          DENDRO_185 * (DENDRO_154 + DENDRO_621) -
                          DENDRO_301 * (1.0 * DENDRO_150 + DENDRO_616);
const double DENDRO_668 =
    -DENDRO_11 *
        (DENDRO_659 +
         alpha[pp] *
             (-DENDRO_109 * (DENDRO_394 + DENDRO_614) -
              DENDRO_109 * (DENDRO_613 + DENDRO_661) -
              DENDRO_109 * (DENDRO_395 + DENDRO_612 + DENDRO_662) +
              DENDRO_119 * DENDRO_660 +
              DENDRO_139 * (DENDRO_378 + DENDRO_545 + DENDRO_547) +
              DENDRO_139 * (DENDRO_488 + DENDRO_489 + DENDRO_548) -
              DENDRO_140 * DENDRO_664 + DENDRO_140 * (DENDRO_573 + DENDRO_624) +
              DENDRO_140 * (DENDRO_623 + DENDRO_665) +
              DENDRO_146 * (DENDRO_191 * grad_0_gt3[pp] + DENDRO_416) +
              DENDRO_161 * (DENDRO_627 + DENDRO_631) +
              DENDRO_161 * (-DENDRO_171 * DENDRO_191 + DENDRO_630) +
              DENDRO_161 * (DENDRO_194 * DENDRO_216 + DENDRO_632) +
              DENDRO_161 * (DENDRO_218 * DENDRO_53 + DENDRO_666) +
              DENDRO_185 * DENDRO_663 - DENDRO_185 * (DENDRO_418 + DENDRO_620) -
              DENDRO_185 * (DENDRO_538 + DENDRO_619) -
              DENDRO_185 * (DENDRO_485 + DENDRO_539 + DENDRO_599) -
              DENDRO_301 * (0.5 * DENDRO_428 + DENDRO_617) -
              DENDRO_301 * (DENDRO_191 * DENDRO_53 + DENDRO_425 + DENDRO_615) -
              DENDRO_633 * (DENDRO_335 + DENDRO_521 + DENDRO_559) + DENDRO_655 +
              DENDRO_667)) -
    DENDRO_11 *
        (DENDRO_659 +
         alpha[pp] *
             (-DENDRO_109 * (1.0 * DENDRO_400 + DENDRO_612) -
              DENDRO_109 * (0.5 * DENDRO_403 + DENDRO_614) -
              DENDRO_109 * (DENDRO_170 * DENDRO_218 + DENDRO_398 + DENDRO_613) -
              DENDRO_119 * (DENDRO_226 * DENDRO_56 + DENDRO_560) -
              DENDRO_119 * (0.25 * DENDRO_242 * grad_0_gt5[pp] - DENDRO_611) +
              DENDRO_140 * (DENDRO_380 + DENDRO_606) +
              DENDRO_140 * (DENDRO_380 + DENDRO_624) +
              DENDRO_140 * (-DENDRO_385 - DENDRO_622) +
              DENDRO_140 * (DENDRO_448 + DENDRO_607) +
              DENDRO_140 * (DENDRO_527 + DENDRO_608) +
              DENDRO_140 * (DENDRO_604 + DENDRO_623) +
              DENDRO_146 * (DENDRO_218 * grad_1_gt0[pp] + DENDRO_412) +
              DENDRO_161 * (DENDRO_628 + DENDRO_630) +
              DENDRO_161 * (DENDRO_631 + DENDRO_632) +
              DENDRO_161 * (DENDRO_249 * DENDRO_57 + DENDRO_627) +
              DENDRO_161 * (DENDRO_170 * DENDRO_191 + DENDRO_628 + DENDRO_629) +
              DENDRO_161 * (DENDRO_218 * DENDRO_54 + DENDRO_625 + DENDRO_626) -
              DENDRO_185 * (DENDRO_419 + DENDRO_619) -
              DENDRO_185 * (DENDRO_478 + DENDRO_621) -
              DENDRO_185 * (DENDRO_481 + DENDRO_620) -
              DENDRO_185 * (0.5 * DENDRO_242 * DENDRO_57 - DENDRO_618) -
              DENDRO_301 * (0.5 * DENDRO_420 + DENDRO_615) -
              DENDRO_301 * (DENDRO_429 + DENDRO_617) -
              DENDRO_301 * (DENDRO_170 * DENDRO_60 + DENDRO_188 + DENDRO_616) -
              DENDRO_609 * (DENDRO_173 + DENDRO_435 + DENDRO_436) -
              DENDRO_609 * (DENDRO_439 + DENDRO_440 + DENDRO_581) -
              DENDRO_633 * (DENDRO_330 + DENDRO_467 + DENDRO_523) +
              DENDRO_655)) +
    DENDRO_12 *
        (DENDRO_520 +
         alpha[pp] *
             (-DENDRO_109 * (DENDRO_410 + DENDRO_449) -
              DENDRO_109 * (DENDRO_226 * DENDRO_53 + DENDRO_446) -
              DENDRO_119 * (1.0 * DENDRO_315 + DENDRO_442) -
              DENDRO_119 * (0.5 * DENDRO_327 + DENDRO_444) -
              DENDRO_119 * (DENDRO_176 * DENDRO_246 + DENDRO_323 + DENDRO_443) +
              DENDRO_140 * (DENDRO_325 + DENDRO_464) +
              DENDRO_140 * (DENDRO_325 + DENDRO_472) +
              DENDRO_140 * (DENDRO_465 + DENDRO_466) +
              DENDRO_140 * (DENDRO_468 + DENDRO_470) +
              DENDRO_140 * (DENDRO_474 + DENDRO_477) +
              DENDRO_140 * (0.5 * DENDRO_329 + DENDRO_331 - DENDRO_462) +
              DENDRO_145 * (DENDRO_435 + DENDRO_438) +
              DENDRO_145 * (DENDRO_439 + DENDRO_441) +
              DENDRO_161 * (DENDRO_418 + DENDRO_483) +
              DENDRO_161 * (DENDRO_478 + DENDRO_479) +
              DENDRO_161 * (DENDRO_481 + DENDRO_484) +
              DENDRO_161 * (DENDRO_187 * DENDRO_242 + DENDRO_485 + DENDRO_487) -
              DENDRO_166 * (DENDRO_246 * grad_2_gt0[pp] + DENDRO_292) -
              DENDRO_185 * (DENDRO_457 + DENDRO_459) -
              DENDRO_185 * (DENDRO_460 + DENDRO_461) -
              DENDRO_185 * (DENDRO_209 * DENDRO_54 + DENDRO_454) -
              DENDRO_185 * (DENDRO_176 * DENDRO_191 + DENDRO_457 + DENDRO_458) -
              DENDRO_185 * (DENDRO_246 * DENDRO_57 + DENDRO_455 + DENDRO_456) -
              DENDRO_301 * (DENDRO_184 + DENDRO_451) -
              DENDRO_301 * (0.5 * DENDRO_422 + DENDRO_450) -
              DENDRO_301 * (DENDRO_176 * DENDRO_194 + DENDRO_430 + DENDRO_452) -
              DENDRO_491 * (DENDRO_488 + DENDRO_490) + DENDRO_512)) +
    DENDRO_12 *
        (DENDRO_520 +
         alpha[pp] *
             (-DENDRO_109 * (DENDRO_178 * DENDRO_205 + DENDRO_527) -
              DENDRO_119 * (0.5 * DENDRO_316 + DENDRO_443) -
              DENDRO_119 * (DENDRO_229 * DENDRO_56 + DENDRO_343 + DENDRO_442) +
              DENDRO_139 * DENDRO_525 + DENDRO_139 * (DENDRO_335 + DENDRO_522) -
              DENDRO_140 * DENDRO_534 + DENDRO_140 * (DENDRO_465 + DENDRO_535) +
              DENDRO_140 * (DENDRO_472 + DENDRO_536) +
              DENDRO_140 * (DENDRO_474 + DENDRO_535) + DENDRO_161 * DENDRO_537 +
              DENDRO_161 * DENDRO_544 + DENDRO_161 * (DENDRO_419 + DENDRO_484) +
              DENDRO_161 * (DENDRO_483 + DENDRO_538) +
              DENDRO_161 * (DENDRO_487 + DENDRO_540) +
              DENDRO_161 * (DENDRO_539 + DENDRO_540) -
              DENDRO_166 * (DENDRO_191 * grad_0_gt5[pp] + DENDRO_417) -
              DENDRO_185 * (DENDRO_456 + DENDRO_531) -
              DENDRO_185 * (-DENDRO_177 * DENDRO_191 + DENDRO_459) -
              DENDRO_185 * (DENDRO_246 * DENDRO_56 + DENDRO_531) -
              DENDRO_301 * (1.0 * DENDRO_424 + DENDRO_452) -
              DENDRO_301 * (DENDRO_191 * DENDRO_56 + DENDRO_426 + DENDRO_450) -
              DENDRO_491 * (DENDRO_378 + DENDRO_546) + DENDRO_512 - DENDRO_526 -
              DENDRO_528 - DENDRO_529 - DENDRO_530 - DENDRO_532)) -
    DENDRO_14 *
        (DENDRO_596 +
         alpha[pp] *
             (-DENDRO_109 * (DENDRO_554 + DENDRO_598) -
              DENDRO_109 * (DENDRO_555 + DENDRO_597) -
              DENDRO_109 * (DENDRO_215 * DENDRO_218 + DENDRO_392 + DENDRO_552) -
              DENDRO_119 * (0.5 * DENDRO_318 + DENDRO_550) -
              DENDRO_119 * (-DENDRO_340 + DENDRO_551) +
              DENDRO_140 * (DENDRO_566 + DENDRO_603) +
              DENDRO_140 * (DENDRO_567 + DENDRO_571) +
              DENDRO_140 * (DENDRO_572 - DENDRO_601) +
              DENDRO_140 * (-DENDRO_208 * DENDRO_218 + DENDRO_602) +
              DENDRO_140 * (DENDRO_215 * DENDRO_246 + DENDRO_565 + DENDRO_603) +
              DENDRO_161 * (DENDRO_379 + DENDRO_574) +
              DENDRO_161 * (DENDRO_382 + DENDRO_580) +
              DENDRO_161 * (DENDRO_382 + DENDRO_608) +
              DENDRO_161 * (DENDRO_446 + DENDRO_606) +
              DENDRO_161 * (DENDRO_447 + DENDRO_607) +
              DENDRO_161 * (DENDRO_577 + DENDRO_604) -
              DENDRO_185 * (DENDRO_536 + DENDRO_558) -
              DENDRO_185 * (DENDRO_561 + DENDRO_600) -
              DENDRO_185 * (DENDRO_564 + DENDRO_600) -
              DENDRO_301 * (DENDRO_419 + DENDRO_557) -
              DENDRO_301 * (DENDRO_176 * DENDRO_178 + DENDRO_599) +
              DENDRO_314 * (DENDRO_218 * grad_1_gt5[pp] + DENDRO_413) +
              DENDRO_591 - DENDRO_609 * (DENDRO_522 + DENDRO_559) +
              DENDRO_610)) -
    DENDRO_14 *
        (DENDRO_596 +
         alpha[pp] *
             (-DENDRO_109 * (-DENDRO_387 + DENDRO_554) -
              DENDRO_109 * (0.5 * DENDRO_391 + DENDRO_552) -
              DENDRO_109 * (DENDRO_205 * DENDRO_249 + DENDRO_389 + DENDRO_555) -
              DENDRO_119 * (1.0 * DENDRO_312 + DENDRO_549) -
              DENDRO_119 * (0.5 * DENDRO_333 + DENDRO_551) -
              DENDRO_119 * (DENDRO_205 * DENDRO_246 + DENDRO_324 + DENDRO_550) +
              DENDRO_140 * (DENDRO_569 + DENDRO_570) +
              DENDRO_140 * (DENDRO_571 + DENDRO_572) +
              DENDRO_140 * (-DENDRO_171 * DENDRO_229 + DENDRO_567) +
              DENDRO_140 * (DENDRO_205 * DENDRO_218 + DENDRO_569) +
              DENDRO_140 * (DENDRO_216 * DENDRO_246 + DENDRO_566) +
              DENDRO_145 * (DENDRO_490 + DENDRO_548) +
              DENDRO_145 * (DENDRO_546 + DENDRO_547) +
              DENDRO_161 * (DENDRO_447 + DENDRO_576) +
              DENDRO_161 * (DENDRO_573 + DENDRO_574) +
              DENDRO_161 * (DENDRO_576 + DENDRO_577) +
              DENDRO_161 * (DENDRO_180 * DENDRO_388 + DENDRO_580) -
              DENDRO_185 * (DENDRO_320 + DENDRO_470) -
              DENDRO_185 * (DENDRO_320 + DENDRO_562) -
              DENDRO_185 * (DENDRO_336 + DENDRO_558) -
              DENDRO_185 * (DENDRO_464 + DENDRO_560) -
              DENDRO_185 * (DENDRO_466 + DENDRO_561) -
              DENDRO_185 * (DENDRO_477 + DENDRO_564) -
              DENDRO_301 * (DENDRO_418 + DENDRO_557) -
              DENDRO_301 * (DENDRO_121 * DENDRO_170 + DENDRO_541) +
              DENDRO_314 * (DENDRO_246 * grad_2_gt3[pp] + DENDRO_295) -
              DENDRO_582 * (DENDRO_441 + DENDRO_581) + DENDRO_591)) +
    DENDRO_16 *
        (DENDRO_270 * DENDRO_61 - DENDRO_272 * DENDRO_79 -
         DENDRO_366 * DENDRO_67 +
         alpha[pp] *
             (-DENDRO_109 * DENDRO_129 * DENDRO_178 - DENDRO_110 +
              DENDRO_117 * DENDRO_119 * DENDRO_180 - DENDRO_126 -
              DENDRO_134 * DENDRO_136 * DENDRO_251 - DENDRO_135 +
              DENDRO_140 * DENDRO_157 +
              DENDRO_140 * (-1.0 * DENDRO_172 - DENDRO_175) +
              DENDRO_140 * (-1.0 * DENDRO_179 - DENDRO_182) +
              DENDRO_140 * (DENDRO_156 * DENDRO_180 + DENDRO_427) +
              DENDRO_140 * (DENDRO_223 * DENDRO_53 + DENDRO_419) +
              DENDRO_140 * (DENDRO_227 * DENDRO_56 + DENDRO_418) +
              DENDRO_146 * DENDRO_151 + DENDRO_146 * (DENDRO_420 + DENDRO_421) +
              DENDRO_146 * (DENDRO_178 * grad_2_gt0[pp] + DENDRO_428) +
              DENDRO_161 * DENDRO_189 +
              DENDRO_161 * (1.0 * DENDRO_421 + DENDRO_425) +
              DENDRO_161 * (DENDRO_155 * DENDRO_194 + DENDRO_429) -
              DENDRO_166 * (DENDRO_422 + DENDRO_423) -
              DENDRO_166 * (DENDRO_180 * grad_2_gt0[pp] + DENDRO_424) -
              DENDRO_167 - DENDRO_185 * (1.0 * DENDRO_423 + DENDRO_426) -
              DENDRO_185 * (-2 * DENDRO_177 * DENDRO_194 + DENDRO_430) -
              DENDRO_186 -
              DENDRO_196 * (DENDRO_190 +
                            DENDRO_192 * (DENDRO_431 + DENDRO_432 + DENDRO_72) +
                            DENDRO_193 * DENDRO_433 + DENDRO_195 * DENDRO_434) +
              DENDRO_222 * DENDRO_357 -
              6.0 * DENDRO_234 * DENDRO_88 * grad_0_gt0[pp] +
              DENDRO_237 * DENDRO_359 + DENDRO_254 * DENDRO_361 +
              DENDRO_265 * DENDRO_363 + DENDRO_279 * grad_0_Gt2[pp] +
              DENDRO_285 * grad2_0_2_gt0[pp] + DENDRO_287 * grad2_0_0_gt0[pp] +
              DENDRO_289 * grad2_1_1_gt0[pp] + DENDRO_291 * grad2_2_2_gt0[pp] -
              DENDRO_296 * DENDRO_416 + DENDRO_409 * grad_0_Gt1[pp] -
              DENDRO_414 * DENDRO_417 - DENDRO_82 - DENDRO_85 - DENDRO_87 +
              4 * grad_0_Gt0[pp] * gt0[pp]) -
         4 * grad2_0_0_alpha[pp]) +
    DENDRO_18 *
        (-DENDRO_269 * DENDRO_369 +
         DENDRO_273 * (DENDRO_232 + DENDRO_365 * DENDRO_78) +
         DENDRO_366 * (DENDRO_249 + DENDRO_276 * DENDRO_365) +
         alpha[pp] *
             (DENDRO_119 * DENDRO_242 * DENDRO_404 +
              DENDRO_140 * (DENDRO_389 - DENDRO_390) +
              DENDRO_140 * (DENDRO_392 + 1.0 * DENDRO_393) +
              DENDRO_140 * (DENDRO_232 * DENDRO_96 - DENDRO_387) +
              DENDRO_146 * (DENDRO_397 + DENDRO_399) +
              DENDRO_146 * (DENDRO_178 * grad_2_gt3[pp] + DENDRO_403) +
              DENDRO_146 * (DENDRO_227 * grad_0_gt3[pp] + DENDRO_400) +
              DENDRO_161 * (DENDRO_395 + DENDRO_396) +
              DENDRO_161 * (DENDRO_398 + 1.0 * DENDRO_399) +
              DENDRO_161 * (DENDRO_249 * DENDRO_96 + DENDRO_394) -
              DENDRO_178 * DENDRO_408 - DENDRO_185 * (DENDRO_377 + DENDRO_379) -
              DENDRO_185 * (-1.0 * DENDRO_383 - DENDRO_385) -
              DENDRO_185 * (DENDRO_227 * DENDRO_381 + DENDRO_380) -
              DENDRO_185 * (DENDRO_242 * DENDRO_381 + DENDRO_382) -
              DENDRO_196 *
                  (DENDRO_192 * DENDRO_372 +
                   DENDRO_193 * (DENDRO_373 + DENDRO_374 + DENDRO_375) +
                   DENDRO_195 * DENDRO_376 + DENDRO_370) -
              DENDRO_219 * DENDRO_297 * grad_1_gt3[pp] -
              DENDRO_226 * DENDRO_406 - DENDRO_227 * DENDRO_407 +
              DENDRO_314 * (DENDRO_391 + DENDRO_393) +
              DENDRO_314 * (DENDRO_226 * grad_0_gt3[pp] + DENDRO_401) +
              DENDRO_314 * (DENDRO_242 * grad_2_gt3[pp] + DENDRO_402) +
              DENDRO_358 * grad_1_gt3[pp] + DENDRO_360 * grad_0_gt3[pp] +
              DENDRO_362 * grad_2_gt3[pp] + DENDRO_364 * DENDRO_368 +
              DENDRO_415) -
         4 * grad2_1_1_alpha[pp]) +
    DENDRO_20 *
        (DENDRO_268 * DENDRO_270 - DENDRO_271 * DENDRO_273 -
         DENDRO_277 * DENDRO_278 +
         alpha[pp] *
             (-DENDRO_109 * DENDRO_226 * DENDRO_305 -
              DENDRO_134 * DENDRO_230 * DENDRO_306 + DENDRO_140 * DENDRO_342 +
              DENDRO_140 * (1.0 * DENDRO_319 + DENDRO_324) +
              DENDRO_140 * (DENDRO_229 * DENDRO_91 - DENDRO_340) +
              DENDRO_161 * DENDRO_322 + DENDRO_161 * DENDRO_332 +
              DENDRO_161 * (DENDRO_334 + DENDRO_336) +
              DENDRO_161 * (DENDRO_176 * DENDRO_242 + DENDRO_310) +
              DENDRO_161 * (DENDRO_180 * DENDRO_205 + DENDRO_308) +
              DENDRO_161 * (DENDRO_223 * DENDRO_321 + DENDRO_325) -
              DENDRO_166 * (DENDRO_316 + DENDRO_317) -
              DENDRO_166 * (DENDRO_223 * grad_0_gt5[pp] + DENDRO_315) -
              DENDRO_185 * (1.0 * DENDRO_317 + DENDRO_323) -
              DENDRO_185 * (DENDRO_209 * DENDRO_91 - DENDRO_338) -
              DENDRO_185 * (DENDRO_229 * DENDRO_344 + DENDRO_343) -
              DENDRO_196 *
                  (DENDRO_192 * DENDRO_350 + DENDRO_193 * DENDRO_352 +
                   DENDRO_195 * (DENDRO_354 + DENDRO_355 + DENDRO_356) +
                   DENDRO_345) -
              DENDRO_223 * DENDRO_300 * DENDRO_301 -
              DENDRO_247 * DENDRO_297 * grad_2_gt5[pp] +
              DENDRO_275 * DENDRO_364 + DENDRO_279 * grad_2_Gt0[pp] +
              DENDRO_280 * grad_2_Gt1[pp] - DENDRO_281 - DENDRO_282 -
              DENDRO_283 + DENDRO_285 * grad2_0_2_gt5[pp] +
              DENDRO_287 * grad2_0_0_gt5[pp] + DENDRO_289 * grad2_1_1_gt5[pp] +
              DENDRO_291 * grad2_2_2_gt5[pp] - DENDRO_292 * DENDRO_294 -
              DENDRO_295 * DENDRO_296 - DENDRO_299 - DENDRO_304 - DENDRO_307 +
              DENDRO_313 * DENDRO_314 + DENDRO_314 * (DENDRO_318 + DENDRO_319) +
              DENDRO_314 * (DENDRO_226 * grad_0_gt5[pp] + DENDRO_333) -
              DENDRO_328 + DENDRO_358 * grad_1_gt5[pp] +
              DENDRO_360 * grad_0_gt5[pp] + DENDRO_362 * grad_2_gt5[pp] +
              4 * grad_2_Gt2[pp] * gt5[pp]) -
         4 * grad2_2_2_alpha[pp]);
const double DENDRO_669 = DENDRO_10 * DENDRO_668;
const double DENDRO_670 = (1.0 / 12.0) * chi[pp];
const double DENDRO_671 = (1.0 / 3.0) * At1[pp];
const double DENDRO_672 = At1[pp] * DENDRO_11;
const double DENDRO_673 = At4[pp] * DENDRO_14;
const double DENDRO_674 = At3[pp] * DENDRO_18;
const double DENDRO_675 = DENDRO_672 + DENDRO_673 - DENDRO_674;
const double DENDRO_676 = At4[pp] * DENDRO_12;
const double DENDRO_677 =
    -At1[pp] * DENDRO_16 + At3[pp] * DENDRO_11 - DENDRO_676;
const double DENDRO_678 = At4[pp] * DENDRO_20;
const double DENDRO_679 =
    -At1[pp] * DENDRO_12 + At3[pp] * DENDRO_14 - DENDRO_678;
const double DENDRO_680 = 6.0 * grad_2_alpha[pp];
const double DENDRO_681 = 6.0 * grad_0_alpha[pp];
const double DENDRO_682 = 6.0 * grad_1_alpha[pp];
const double DENDRO_683 = DENDRO_96 * DENDRO_98;
const double DENDRO_684 = DENDRO_225 * grad_2_gt0[pp];
const double DENDRO_685 = DENDRO_155 * DENDRO_98;
const double DENDRO_686 = DENDRO_684 + DENDRO_685;
const double DENDRO_687 = DENDRO_386 * DENDRO_93;
const double DENDRO_688 = DENDRO_132 * DENDRO_388;
const double DENDRO_689 = DENDRO_183 * DENDRO_98;
const double DENDRO_690 =
    DENDRO_156 * DENDRO_99 + 0.25 * DENDRO_225 * grad_0_gt0[pp];
const double DENDRO_691 = DENDRO_138 + DENDRO_142;
const double DENDRO_692 = DENDRO_113 * DENDRO_578;
const double DENDRO_693 = DENDRO_205 * DENDRO_66;
const double DENDRO_694 = DENDRO_141 * DENDRO_241 + DENDRO_693;
const double DENDRO_695 = DENDRO_155 * DENDRO_93;
const double DENDRO_696 = DENDRO_225 * grad_1_gt0[pp] + DENDRO_695;
const double DENDRO_697 = 1.0 * DENDRO_139;
const double DENDRO_698 = DENDRO_155 * DENDRO_241;
const double DENDRO_699 = 2.0 * DENDRO_236;
const double DENDRO_700 = 2.0 * DENDRO_221;
const double DENDRO_701 = 2.0 * DENDRO_253;
const double DENDRO_702 = DENDRO_264 * DENDRO_50;
const double DENDRO_703 = (1.0 / 3.0) * At2[pp];
const double DENDRO_704 = At5[pp] * DENDRO_14;
const double DENDRO_705 =
    At2[pp] * DENDRO_11 - At4[pp] * DENDRO_18 + DENDRO_704;
const double DENDRO_706 = At5[pp] * DENDRO_12;
const double DENDRO_707 =
    -At2[pp] * DENDRO_16 + At4[pp] * DENDRO_11 - DENDRO_706;
const double DENDRO_708 = At4[pp] * DENDRO_14 - At5[pp] * DENDRO_20 - DENDRO_41;
const double DENDRO_709 = DENDRO_113 * grad_2_gt5[pp];
const double DENDRO_710 = DENDRO_93 * grad_0_gt5[pp];
const double DENDRO_711 = DENDRO_113 * DENDRO_309;
const double DENDRO_712 = 0.25 * DENDRO_132 * grad_2_gt5[pp];
const double DENDRO_713 = DENDRO_711 + DENDRO_712;
const double DENDRO_714 = DENDRO_91 * DENDRO_93;
const double DENDRO_715 = DENDRO_113 * DENDRO_114 + DENDRO_203 * DENDRO_66;
const double DENDRO_716 = -0.5 * DENDRO_177 * DENDRO_98;
const double DENDRO_717 = -0.5 * DENDRO_177 * DENDRO_93;
const double DENDRO_718 = DENDRO_187 * DENDRO_93;
const double DENDRO_719 = 2 * At4[pp];
const double DENDRO_720 = At3[pp] * DENDRO_43;
const double DENDRO_721 = DENDRO_10 * DENDRO_719;
const double DENDRO_722 = 0.5 * DENDRO_365;
const double DENDRO_723 = DENDRO_10 * DENDRO_80;
const double DENDRO_724 = DENDRO_217 * grad_2_gt3[pp];
const double DENDRO_725 = DENDRO_132 * DENDRO_309;
const double DENDRO_726 = 1.0 * DENDRO_225;
const double DENDRO_727 = 0.25 * DENDRO_695;
const double DENDRO_728 = DENDRO_217 * grad_0_gt3[pp];
const double DENDRO_729 = (1.0 / 3.0) * At4[pp];
const double DENDRO_730 = DENDRO_241 * grad_2_gt5[pp];
const double DENDRO_731 = DENDRO_114 * DENDRO_241;
const double DENDRO_732 = DENDRO_712 + DENDRO_731;
const double DENDRO_733 = DENDRO_241 * DENDRO_309;
const double DENDRO_734 = -0.5 * DENDRO_249 * grad_0_gt5[pp] + DENDRO_725;
const double DENDRO_735 = DENDRO_245 * grad_0_gt5[pp];
const double DENDRO_736 = DENDRO_245 * grad_1_gt5[pp];
const double DENDRO_737 = DENDRO_11 * grad_1_chi[pp] +
                          DENDRO_346 * grad_2_chi[pp] +
                          DENDRO_348 * grad_0_chi[pp];
const double DENDRO_738 = 0.5 * DENDRO_267;
const double DENDRO_739 = DENDRO_10 * grad_0_alpha[pp];
const double DENDRO_740 = DENDRO_351 * grad_1_chi[pp] + DENDRO_47;
const double DENDRO_741 = DENDRO_10 * grad_1_alpha[pp];
const double DENDRO_742 = DENDRO_14 * grad_1_chi[pp] +
                          DENDRO_346 * grad_0_chi[pp] +
                          DENDRO_353 * grad_2_chi[pp];
const double DENDRO_743 = 0.5 * DENDRO_742;
const double DENDRO_744 = DENDRO_10 * grad_2_alpha[pp];
const double DENDRO_745 = 0.5 * DENDRO_740;
const double DENDRO_746 = 0.5 * grad_1_alpha[pp];
const double DENDRO_747 = 0.5 * grad_2_alpha[pp];
const double DENDRO_748 = 0.5 * grad_0_alpha[pp];
const double DENDRO_749 = DENDRO_10 * DENDRO_259;
const double DENDRO_750 = (DENDRO_11 * DENDRO_11);
const double DENDRO_751 = (DENDRO_12 * DENDRO_12);
const double DENDRO_752 = 2 * DENDRO_16;
const double DENDRO_753 = At0[pp] * (DENDRO_16 * DENDRO_16) +
                          At3[pp] * DENDRO_750 + At5[pp] * DENDRO_751 -
                          DENDRO_257 * DENDRO_676 + DENDRO_41 * DENDRO_752 -
                          DENDRO_672 * DENDRO_752;
const double DENDRO_754 = 3 * DENDRO_88;
const double DENDRO_755 = (DENDRO_14 * DENDRO_14);
const double DENDRO_756 = 2 * DENDRO_18;
const double DENDRO_757 = At0[pp] * DENDRO_750 +
                          At3[pp] * (DENDRO_18 * DENDRO_18) +
                          At5[pp] * DENDRO_755 + DENDRO_257 * DENDRO_39 -
                          DENDRO_672 * DENDRO_756 - DENDRO_673 * DENDRO_756;
const double DENDRO_758 = At1[pp] * DENDRO_14;
const double DENDRO_759 = 2 * DENDRO_20;
const double DENDRO_760 = At0[pp] * DENDRO_751 + At3[pp] * DENDRO_755 +
                          At5[pp] * (DENDRO_20 * DENDRO_20) -
                          DENDRO_259 * DENDRO_758 + DENDRO_41 * DENDRO_759 -
                          DENDRO_673 * DENDRO_759;
const double DENDRO_761 = At3[pp] * DENDRO_11;
const double DENDRO_762 =
    At2[pp] * DENDRO_751 - DENDRO_11 * DENDRO_678 + DENDRO_12 * DENDRO_42 -
    DENDRO_12 * DENDRO_672 - DENDRO_14 * DENDRO_676 + DENDRO_14 * DENDRO_761 +
    DENDRO_16 * DENDRO_45 - DENDRO_16 * DENDRO_758 + DENDRO_20 * DENDRO_706;
const double DENDRO_763 = 6 * DENDRO_88;
const double DENDRO_764 =
    -At1[pp] * DENDRO_16 * DENDRO_18 - At1[pp] * DENDRO_750 -
    At4[pp] * DENDRO_11 * DENDRO_14 - At4[pp] * DENDRO_12 * DENDRO_18 +
    DENDRO_11 * DENDRO_41 + DENDRO_11 * DENDRO_42 + DENDRO_14 * DENDRO_706 +
    DENDRO_16 * DENDRO_39 + DENDRO_18 * DENDRO_761;
const double DENDRO_765 = -DENDRO_764;
const double DENDRO_766 =
    -At1[pp] * DENDRO_11 * DENDRO_14 - At1[pp] * DENDRO_12 * DENDRO_18 -
    At4[pp] * DENDRO_18 * DENDRO_20 - At4[pp] * DENDRO_755 +
    DENDRO_11 * DENDRO_45 + DENDRO_12 * DENDRO_38 + DENDRO_14 * DENDRO_41 +
    DENDRO_14 * DENDRO_674 + DENDRO_20 * DENDRO_704;
const double DENDRO_767 = -DENDRO_766;
const double DENDRO_768 = (1.0 / 3.0) * alpha[pp];
const double DENDRO_769 = (7.0 / 3.0) * DENDRO_284;
const double DENDRO_770 = (1.0 / 3.0) * DENDRO_284;
const double DENDRO_771 = (1.0 / 3.0) * DENDRO_286;
const double DENDRO_772 = 2 * DENDRO_88;
const double DENDRO_773 = DENDRO_772 * grad_0_alpha[pp];
const double DENDRO_774 = pow(DENDRO_9, -3);
const double DENDRO_775 = 4 * grad_0_K[pp];
const double DENDRO_776 = DENDRO_10 * DENDRO_768;
const double DENDRO_777 = DENDRO_772 * grad_2_alpha[pp];
const double DENDRO_778 = DENDRO_772 * grad_1_alpha[pp];
const double DENDRO_779 = 4 * grad_2_K[pp];
const double DENDRO_780 = 4 * grad_1_K[pp];
const double DENDRO_781 = 9 * DENDRO_50;
const double DENDRO_782 = DENDRO_193 * DENDRO_781;
const double DENDRO_783 = DENDRO_10 * DENDRO_11;
const double DENDRO_784 = (1.0 / 3.0) * DENDRO_783;
const double DENDRO_785 = DENDRO_10 * DENDRO_14;
const double DENDRO_786 = (1.0 / 3.0) * grad2_0_2_beta0[pp];
const double DENDRO_787 = (1.0 / 3.0) * DENDRO_288;
const double DENDRO_788 = DENDRO_0 * DENDRO_774;
const double DENDRO_789 = 2.0 * DENDRO_774 * alpha[pp];

// Dendro: printing variables
//--
a_rhs[pp] =
    -DENDRO_0 * K[pp] -
    0.59999999999999998 * DENDRO_1 * (-DENDRO_1 + alpha[pp]) *
        exp(-1.0 / 800.0 * (t * t)) +
    lambda[0] * (beta0[pp] * grad_0_alpha[pp] + beta1[pp] * grad_1_alpha[pp] +
                 beta2[pp] * grad_2_alpha[pp]);
//--
b_rhs0[pp] =
    BSSN_XI[1] * (beta0[pp] * grad_0_beta0[pp] + beta1[pp] * grad_1_beta0[pp] +
                  beta2[pp] * grad_2_beta0[pp]) +
    DENDRO_2 * Gt0[pp] - DENDRO_21 * beta0[pp];
//--
b_rhs1[pp] =
    BSSN_XI[1] * (beta0[pp] * grad_0_beta1[pp] + beta1[pp] * grad_1_beta1[pp] +
                  beta2[pp] * grad_2_beta1[pp]) +
    DENDRO_2 * Gt1[pp] - DENDRO_21 * beta1[pp];
//--
b_rhs2[pp] =
    BSSN_XI[1] * (beta0[pp] * grad_0_beta2[pp] + beta1[pp] * grad_1_beta2[pp] +
                  beta2[pp] * grad_2_beta2[pp]) +
    DENDRO_2 * Gt2[pp] - DENDRO_21 * beta2[pp];
//--
gt_rhs00[pp] = -At0[pp] * DENDRO_0 + DENDRO_22 * gt0[pp] +
               DENDRO_23 * grad_0_beta1[pp] - DENDRO_24 * grad_1_beta1[pp] -
               DENDRO_24 * grad_2_beta2[pp] + DENDRO_8 * grad_0_beta2[pp] +
               beta0[pp] * grad_0_gt0[pp] + beta1[pp] * grad_1_gt0[pp] +
               beta2[pp] * grad_2_gt0[pp];
//--
gt_rhs01[pp] = -At1[pp] * DENDRO_0 + DENDRO_25 * grad_0_beta0[pp] +
               DENDRO_25 * grad_1_beta1[pp] - DENDRO_26 * gt1[pp] +
               beta0[pp] * grad_0_gt1[pp] + beta1[pp] * grad_1_gt1[pp] +
               beta2[pp] * grad_2_gt1[pp] + grad_0_beta1[pp] * gt3[pp] +
               grad_0_beta2[pp] * gt4[pp] + grad_1_beta0[pp] * gt0[pp] +
               grad_1_beta2[pp] * gt2[pp];
//--
gt_rhs02[pp] = -At2[pp] * DENDRO_0 + DENDRO_27 * grad_0_beta0[pp] +
               DENDRO_27 * grad_2_beta2[pp] - DENDRO_28 * gt2[pp] +
               beta0[pp] * grad_0_gt2[pp] + beta1[pp] * grad_1_gt2[pp] +
               beta2[pp] * grad_2_gt2[pp] + grad_0_beta1[pp] * gt4[pp] +
               grad_0_beta2[pp] * gt5[pp] + grad_2_beta0[pp] * gt0[pp] +
               grad_2_beta1[pp] * gt1[pp];
//--
gt_rhs11[pp] = -At3[pp] * DENDRO_0 + DENDRO_23 * grad_1_beta0[pp] -
               DENDRO_26 * gt3[pp] - DENDRO_29 * gt3[pp] + DENDRO_30 * gt3[pp] +
               DENDRO_31 * grad_1_beta2[pp] + beta0[pp] * grad_0_gt3[pp] +
               beta1[pp] * grad_1_gt3[pp] + beta2[pp] * grad_2_gt3[pp];
//--
gt_rhs12[pp] = -At4[pp] * DENDRO_0 - DENDRO_29 * gt4[pp] +
               DENDRO_32 * grad_1_beta1[pp] + DENDRO_32 * grad_2_beta2[pp] +
               beta0[pp] * grad_0_gt4[pp] + beta1[pp] * grad_1_gt4[pp] +
               beta2[pp] * grad_2_gt4[pp] + grad_1_beta0[pp] * gt2[pp] +
               grad_1_beta2[pp] * gt5[pp] + grad_2_beta0[pp] * gt1[pp] +
               grad_2_beta1[pp] * gt3[pp];
//--
gt_rhs22[pp] = -At5[pp] * DENDRO_0 - DENDRO_28 * gt5[pp] - DENDRO_29 * gt5[pp] +
               DENDRO_31 * grad_2_beta1[pp] + DENDRO_33 * gt5[pp] +
               DENDRO_8 * grad_2_beta0[pp] + beta0[pp] * grad_0_gt5[pp] +
               beta1[pp] * grad_1_gt5[pp] + beta2[pp] * grad_2_gt5[pp];
//--
chi_rhs[pp] = -DENDRO_34 * DENDRO_35 + DENDRO_34 * K[pp] * alpha[pp] +
              beta0[pp] * grad_0_chi[pp] + beta1[pp] * grad_1_chi[pp] +
              beta2[pp] * grad_2_chi[pp];
//--
At_rhs00[pp] =
    At0[pp] * DENDRO_22 - At0[pp] * DENDRO_26 - At0[pp] * DENDRO_28 +
    DENDRO_36 * grad_0_beta1[pp] + DENDRO_37 * grad_0_beta2[pp] +
    DENDRO_670 *
        (DENDRO_266 *
             (4 * DENDRO_10 * DENDRO_12 * grad2_0_2_gt0[pp] +
              2.0 * DENDRO_10 * DENDRO_16 * grad2_0_0_gt0[pp] +
              2.0 * DENDRO_10 * DENDRO_18 * grad2_1_1_gt0[pp] +
              2.0 * DENDRO_10 * DENDRO_20 * grad2_2_2_gt0[pp] +
              2.0 * DENDRO_11 * DENDRO_151 * DENDRO_88 +
              4 * DENDRO_11 * DENDRO_189 * DENDRO_88 - DENDRO_110 -
              DENDRO_113 * DENDRO_117 * DENDRO_119 +
              2.0 * DENDRO_12 * DENDRO_88 * (DENDRO_147 + DENDRO_148) +
              4 * DENDRO_12 * DENDRO_88 * (1.0 * DENDRO_148 + DENDRO_162) +
              2.0 * DENDRO_12 * DENDRO_88 * (DENDRO_152 + DENDRO_153) -
              DENDRO_126 + 4 * DENDRO_129 * DENDRO_132 * DENDRO_18 * DENDRO_88 -
              DENDRO_135 + 4 * DENDRO_136 * DENDRO_16 * DENDRO_66 * DENDRO_88 +
              4 * DENDRO_14 * DENDRO_157 * DENDRO_88 -
              DENDRO_140 * (DENDRO_138 + DENDRO_56 * DENDRO_93) -
              DENDRO_140 * (DENDRO_142 + DENDRO_53 * DENDRO_98) -
              DENDRO_140 * (DENDRO_172 + DENDRO_175) -
              DENDRO_140 * (DENDRO_179 + DENDRO_182) -
              DENDRO_140 * (DENDRO_155 * DENDRO_159 + DENDRO_158) -
              DENDRO_146 * (DENDRO_143 + DENDRO_144) -
              DENDRO_146 * (DENDRO_168 + DENDRO_169) +
              6.0 * DENDRO_16 * DENDRO_88 * DENDRO_99 * grad_0_gt0[pp] -
              DENDRO_161 * (1.0 * DENDRO_144 + DENDRO_160) -
              DENDRO_161 *
                  (DENDRO_132 * DENDRO_187 + 1.0 * DENDRO_155 * DENDRO_66) -
              DENDRO_167 +
              3.0 * DENDRO_18 * DENDRO_88 * DENDRO_93 * grad_1_gt0[pp] -
              DENDRO_185 *
                  (-DENDRO_159 * DENDRO_57 + 2 * DENDRO_177 * DENDRO_66) -
              DENDRO_186 -
              DENDRO_196 * (DENDRO_190 + DENDRO_191 * DENDRO_192 +
                            DENDRO_193 * DENDRO_60 + DENDRO_194 * DENDRO_195) +
              3.0 * DENDRO_20 * DENDRO_88 * DENDRO_98 * grad_2_gt0[pp] -
              DENDRO_221 * DENDRO_222 - DENDRO_236 * DENDRO_237 -
              DENDRO_253 * DENDRO_254 - DENDRO_264 * DENDRO_265 - DENDRO_82 -
              DENDRO_85 - DENDRO_87 + 4 * grad_0_Gt0[pp] * gt0[pp] +
              4 * grad_0_Gt1[pp] * gt1[pp] + 4 * grad_0_Gt2[pp] * gt2[pp]) +
         DENDRO_61 * DENDRO_63 + DENDRO_669 * gt0[pp] - DENDRO_67 * DENDRO_69 -
         DENDRO_79 * DENDRO_80 - 12 * grad2_0_0_alpha[pp]) -
    alpha[pp] *
        (-At0[pp] * K[pp] +
         DENDRO_40 * (-At1[pp] * DENDRO_18 + DENDRO_38 + DENDRO_39) +
         DENDRO_44 * (At1[pp] * DENDRO_11 - DENDRO_41 - DENDRO_42) +
         DENDRO_46 * (-At0[pp] * DENDRO_12 + At1[pp] * DENDRO_14 - DENDRO_45)) +
    beta0[pp] * grad_0_At0[pp] + beta1[pp] * grad_1_At0[pp] +
    beta2[pp] * grad_2_At0[pp];
//--
At_rhs01[pp] =
    At0[pp] * grad_1_beta0[pp] - At1[pp] * DENDRO_26 +
    At2[pp] * grad_1_beta2[pp] + At3[pp] * grad_0_beta1[pp] +
    At4[pp] * grad_0_beta2[pp] +
    DENDRO_670 *
        (DENDRO_10 * DENDRO_680 * (-DENDRO_132 - DENDRO_64 * DENDRO_658) +
         DENDRO_266 *
             (-DENDRO_109 * (DENDRO_171 * DENDRO_217 + DENDRO_661) -
              DENDRO_109 * (-DENDRO_105 * DENDRO_93 + DENDRO_662 + DENDRO_687) -
              DENDRO_109 * (-DENDRO_132 * DENDRO_298 +
                            0.5 * DENDRO_155 * DENDRO_249 - DENDRO_688) +
              DENDRO_118 * (DENDRO_683 + DENDRO_686) + DENDRO_119 * DENDRO_660 -
              DENDRO_140 * DENDRO_664 +
              DENDRO_140 * (-DENDRO_156 * DENDRO_217 + DENDRO_665) +
              DENDRO_140 * (DENDRO_386 * DENDRO_98 - DENDRO_578 * DENDRO_93 +
                            DENDRO_605) -
              DENDRO_146 *
                  (DENDRO_93 * grad_1_gt0[pp] + DENDRO_99 * grad_0_gt3[pp]) +
              DENDRO_161 * (-DENDRO_217 * DENDRO_53 + DENDRO_666) -
              DENDRO_161 * (DENDRO_132 * DENDRO_471 + DENDRO_132 * DENDRO_486 +
                            DENDRO_215 * DENDRO_66) +
              DENDRO_161 * (-DENDRO_132 * DENDRO_578 - DENDRO_216 * DENDRO_66 +
                            0.5 * DENDRO_249 * grad_2_gt0[pp]) +
              DENDRO_161 * (DENDRO_171 * DENDRO_99 - DENDRO_183 * DENDRO_93 +
                            DENDRO_629) +
              DENDRO_185 * DENDRO_663 + DENDRO_185 * (DENDRO_689 + DENDRO_690) +
              DENDRO_185 * (DENDRO_692 + DENDRO_694) +
              DENDRO_185 * (DENDRO_381 * DENDRO_99 + DENDRO_691) -
              DENDRO_196 * (DENDRO_104 * DENDRO_497 + DENDRO_178 * DENDRO_501 +
                            DENDRO_227 * DENDRO_495 + DENDRO_634) +
              DENDRO_301 *
                  (DENDRO_160 + DENDRO_53 * DENDRO_99 + DENDRO_54 * DENDRO_99) +
              DENDRO_301 * (0.25 * DENDRO_168 + 0.5 * DENDRO_169 +
                            DENDRO_381 * DENDRO_66) +
              DENDRO_641 + DENDRO_642 + DENDRO_643 + DENDRO_644 + DENDRO_645 +
              DENDRO_646 + DENDRO_647 + DENDRO_648 + DENDRO_649 + DENDRO_650 +
              DENDRO_651 + DENDRO_652 + DENDRO_653 - DENDRO_654 * DENDRO_702 +
              DENDRO_667 -
              DENDRO_697 * (DENDRO_696 + DENDRO_98 * grad_0_gt3[pp]) -
              DENDRO_697 * (DENDRO_113 * grad_2_gt3[pp] +
                            DENDRO_132 * grad_1_gt5[pp] + DENDRO_698) -
              DENDRO_699 * grad_0_gt1[pp] - DENDRO_700 * grad_1_gt1[pp] -
              DENDRO_701 * grad_2_gt1[pp]) +
         DENDRO_654 * DENDRO_668 - DENDRO_656 * DENDRO_681 -
         DENDRO_657 * DENDRO_682 - 12 * grad2_0_1_alpha[pp]) +
    DENDRO_671 * grad_0_beta0[pp] + DENDRO_671 * grad_1_beta1[pp] -
    alpha[pp] * (-At1[pp] * K[pp] + DENDRO_40 * DENDRO_675 +
                 DENDRO_44 * DENDRO_677 + DENDRO_46 * DENDRO_679) +
    beta0[pp] * grad_0_At1[pp] + beta1[pp] * grad_1_At1[pp] +
    beta2[pp] * grad_2_At1[pp];
//--
At_rhs02[pp] =
    At0[pp] * grad_2_beta0[pp] + At1[pp] * grad_2_beta1[pp] -
    At2[pp] * DENDRO_28 + At4[pp] * grad_0_beta1[pp] +
    At5[pp] * grad_0_beta2[pp] +
    DENDRO_670 *
        (DENDRO_266 *
             (4 * DENDRO_10 * DENDRO_12 * grad2_0_2_gt2[pp] +
              2.0 * DENDRO_10 * DENDRO_16 * grad2_0_0_gt2[pp] +
              2.0 * DENDRO_10 * DENDRO_18 * grad2_1_1_gt2[pp] +
              2.0 * DENDRO_10 * DENDRO_20 * grad2_2_2_gt2[pp] +
              4 * DENDRO_11 * DENDRO_537 * DENDRO_88 +
              4 * DENDRO_11 * DENDRO_544 * DENDRO_88 -
              DENDRO_119 * (DENDRO_177 * DENDRO_245 - 0.5 * DENDRO_709) -
              DENDRO_119 * (-DENDRO_114 * DENDRO_98 - DENDRO_228 * DENDRO_56 -
                            DENDRO_716) +
              4 * DENDRO_12 * DENDRO_88 *
                  (DENDRO_245 * DENDRO_56 + DENDRO_715) +
              2.0 * DENDRO_12 * DENDRO_88 *
                  (DENDRO_98 * grad_2_gt0[pp] + DENDRO_99 * grad_0_gt5[pp]) +
              DENDRO_14 * DENDRO_525 * DENDRO_88 +
              4 * DENDRO_14 * DENDRO_88 *
                  (0.5 * DENDRO_177 * DENDRO_241 - DENDRO_713) +
              4 * DENDRO_14 * DENDRO_88 *
                  (-DENDRO_228 * DENDRO_53 - DENDRO_471 * DENDRO_98 -
                   DENDRO_717) -
              DENDRO_140 * DENDRO_534 -
              DENDRO_140 * (DENDRO_156 * DENDRO_245 + DENDRO_713) +
              4 * DENDRO_16 * DENDRO_88 *
                  (0.25 * DENDRO_152 + 1.0 * DENDRO_153) +
              4 * DENDRO_16 * DENDRO_88 *
                  (DENDRO_162 + DENDRO_56 * DENDRO_99 + DENDRO_57 * DENDRO_99) -
              DENDRO_161 * (DENDRO_158 + DENDRO_694) -
              DENDRO_161 * (DENDRO_690 + DENDRO_718) -
              DENDRO_161 * (DENDRO_321 * DENDRO_99 + DENDRO_691) -
              DENDRO_161 * (DENDRO_113 * DENDRO_486 + DENDRO_158 + DENDRO_693) +
              DENDRO_18 * DENDRO_88 * (DENDRO_696 + DENDRO_714) +
              4 * DENDRO_18 * DENDRO_88 *
                  (DENDRO_132 * DENDRO_205 + 0.25 * DENDRO_698) -
              DENDRO_185 * (0.5 * DENDRO_113 * DENDRO_177 - DENDRO_715) -
              DENDRO_185 * (DENDRO_177 * DENDRO_99 - DENDRO_187 * DENDRO_98 -
                            DENDRO_228 * DENDRO_58) -
              DENDRO_196 * (DENDRO_121 * DENDRO_497 + DENDRO_180 * DENDRO_501 +
                            DENDRO_223 * DENDRO_495 + DENDRO_492) -
              DENDRO_508 - DENDRO_509 - DENDRO_510 - DENDRO_511 * DENDRO_702 -
              DENDRO_526 - DENDRO_528 - DENDRO_529 - DENDRO_530 - DENDRO_532 -
              DENDRO_697 * (DENDRO_686 + DENDRO_710) -
              DENDRO_699 * grad_0_gt2[pp] - DENDRO_700 * grad_1_gt2[pp] -
              DENDRO_701 * grad_2_gt2[pp] + 2.0 * grad_0_Gt0[pp] * gt2[pp] +
              2.0 * grad_0_Gt1[pp] * gt4[pp] + 2.0 * grad_0_Gt2[pp] * gt5[pp] +
              2.0 * grad_2_Gt0[pp] * gt0[pp] + 2.0 * grad_2_Gt1[pp] * gt1[pp] +
              2.0 * grad_2_Gt2[pp] * gt2[pp]) +
         DENDRO_511 * DENDRO_668 - DENDRO_513 * DENDRO_681 -
         DENDRO_515 * DENDRO_680 + DENDRO_519 * DENDRO_682 -
         12 * grad2_0_2_alpha[pp]) +
    DENDRO_703 * grad_0_beta0[pp] + DENDRO_703 * grad_2_beta2[pp] -
    alpha[pp] * (-At2[pp] * K[pp] + DENDRO_40 * DENDRO_705 +
                 DENDRO_44 * DENDRO_707 + DENDRO_46 * DENDRO_708) +
    beta0[pp] * grad_0_At2[pp] + beta1[pp] * grad_1_At2[pp] +
    beta2[pp] * grad_2_At2[pp];
//--
At_rhs11[pp] =
    -At3[pp] * DENDRO_26 - At3[pp] * DENDRO_29 + At3[pp] * DENDRO_30 +
    DENDRO_36 * grad_1_beta0[pp] +
    DENDRO_670 *
        (DENDRO_266 *
             (6.0 * DENDRO_108 * DENDRO_217 * grad_1_gt3[pp] -
              DENDRO_119 * DENDRO_241 * DENDRO_404 + DENDRO_132 * DENDRO_408 +
              DENDRO_140 * (DENDRO_392 - 1.0 * DENDRO_724) -
              DENDRO_140 * (-1.0 * DENDRO_232 * DENDRO_96 + DENDRO_387) -
              DENDRO_140 * (DENDRO_241 * DENDRO_388 + DENDRO_390) +
              DENDRO_146 * (DENDRO_104 * grad_1_gt3[pp] - DENDRO_728) +
              DENDRO_146 *
                  (-DENDRO_132 * grad_2_gt3[pp] + DENDRO_249 * DENDRO_91) +
              DENDRO_146 *
                  (DENDRO_232 * grad_1_gt0[pp] - DENDRO_93 * grad_0_gt3[pp]) +
              DENDRO_161 * (DENDRO_396 + DENDRO_687) +
              DENDRO_161 *
                  (0.25 * DENDRO_104 * grad_1_gt3[pp] - 1.0 * DENDRO_728) +
              DENDRO_161 * (DENDRO_249 * DENDRO_96 - DENDRO_688) +
              DENDRO_185 * (DENDRO_383 + DENDRO_385) +
              DENDRO_185 * (DENDRO_137 * DENDRO_225 + DENDRO_381 * DENDRO_93) +
              DENDRO_185 * (DENDRO_241 * DENDRO_381 + DENDRO_725) +
              DENDRO_185 * (DENDRO_54 * DENDRO_726 + DENDRO_727) -
              DENDRO_196 * (DENDRO_192 * DENDRO_232 + DENDRO_193 * DENDRO_218 +
                            DENDRO_195 * DENDRO_249 + DENDRO_370) +
              DENDRO_225 * DENDRO_406 + DENDRO_314 * (DENDRO_391 - DENDRO_724) +
              DENDRO_314 *
                  (DENDRO_155 * DENDRO_232 - DENDRO_225 * grad_0_gt3[pp]) +
              DENDRO_314 * (-DENDRO_241 * grad_2_gt3[pp] + DENDRO_402) -
              DENDRO_368 * DENDRO_702 + DENDRO_407 * DENDRO_93 + DENDRO_415 -
              DENDRO_699 * grad_0_gt3[pp] - DENDRO_700 * grad_1_gt3[pp] -
              DENDRO_701 * grad_2_gt3[pp]) -
         DENDRO_369 * DENDRO_62 + DENDRO_669 * gt3[pp] +
         DENDRO_69 * (DENDRO_249 - DENDRO_64 * DENDRO_722) +
         DENDRO_723 * (DENDRO_232 - DENDRO_722 * DENDRO_76) -
         12 * grad2_1_1_alpha[pp]) +
    DENDRO_719 * grad_1_beta2[pp] -
    alpha[pp] * (-At3[pp] * K[pp] + DENDRO_40 * DENDRO_677 +
                 DENDRO_675 * DENDRO_720 + DENDRO_679 * DENDRO_721) +
    beta0[pp] * grad_0_At3[pp] + beta1[pp] * grad_1_At3[pp] +
    beta2[pp] * grad_2_At3[pp];
//--
At_rhs12[pp] =
    At1[pp] * grad_2_beta0[pp] + At2[pp] * grad_1_beta0[pp] +
    At3[pp] * grad_2_beta1[pp] - At4[pp] * DENDRO_29 +
    At5[pp] * grad_1_beta2[pp] +
    DENDRO_670 *
        (-DENDRO_10 * DENDRO_595 * DENDRO_681 +
         DENDRO_266 *
             (-DENDRO_109 * (-DENDRO_241 * DENDRO_298 + DENDRO_597) -
              DENDRO_109 *
                  (-DENDRO_105 * DENDRO_225 + DENDRO_553 + DENDRO_598) -
              DENDRO_109 * (0.25 * DENDRO_202 * grad_1_gt3[pp] -
                            DENDRO_215 * DENDRO_217 - DENDRO_216 * DENDRO_217) -
              DENDRO_119 * (DENDRO_208 * DENDRO_245 - 0.5 * DENDRO_730) -
              DENDRO_119 *
                  (-DENDRO_114 * DENDRO_225 + 0.5 * DENDRO_177 * DENDRO_225 -
                   DENDRO_228 * DENDRO_381) +
              DENDRO_140 * (DENDRO_208 * DENDRO_217 + DENDRO_602) -
              DENDRO_140 * (DENDRO_170 * DENDRO_228 + DENDRO_225 * DENDRO_471 +
                            DENDRO_601) +
              DENDRO_140 * (-DENDRO_215 * DENDRO_245 +
                            0.5 * DENDRO_249 * grad_2_gt5[pp] - DENDRO_733) +
              DENDRO_140 * (-DENDRO_225 * DENDRO_486 - DENDRO_225 * DENDRO_578 +
                            0.5 * DENDRO_232 * grad_0_gt5[pp]) +
              DENDRO_140 * (DENDRO_241 * DENDRO_337 + DENDRO_565 - DENDRO_733) +
              DENDRO_161 * (-DENDRO_113 * DENDRO_298 - DENDRO_734) +
              DENDRO_161 * (-DENDRO_217 * DENDRO_321 + DENDRO_604) +
              DENDRO_161 * (-DENDRO_241 * DENDRO_578 - DENDRO_734) +
              DENDRO_161 *
                  (-DENDRO_105 * DENDRO_98 + 0.5 * DENDRO_232 * grad_2_gt0[pp] -
                   0.25 * DENDRO_714) +
              DENDRO_161 * (-DENDRO_183 * DENDRO_225 + DENDRO_232 * DENDRO_57 -
                            DENDRO_727) +
              DENDRO_161 *
                  (-DENDRO_217 * DENDRO_381 + DENDRO_447 + DENDRO_575) +
              DENDRO_165 * (DENDRO_683 + DENDRO_684 + DENDRO_710) -
              DENDRO_185 * (0.5 * DENDRO_113 * DENDRO_208 - DENDRO_732) +
              DENDRO_185 * (DENDRO_245 * DENDRO_381 + DENDRO_732) -
              DENDRO_185 * (-DENDRO_187 * DENDRO_225 - DENDRO_228 * DENDRO_54 -
                            DENDRO_717) -
              DENDRO_196 * (DENDRO_202 * DENDRO_497 + DENDRO_226 * DENDRO_495 +
                            DENDRO_242 * DENDRO_501 + DENDRO_583) +
              DENDRO_301 * (DENDRO_132 * DENDRO_176 + DENDRO_692) +
              DENDRO_301 * (DENDRO_142 + DENDRO_689 + DENDRO_718) +
              DENDRO_314 *
                  (DENDRO_202 * grad_2_gt3[pp] - DENDRO_217 * grad_1_gt5[pp]) -
              DENDRO_589 * DENDRO_702 + DENDRO_590 + DENDRO_610 -
              DENDRO_699 * grad_0_gt4[pp] - DENDRO_700 * grad_1_gt4[pp] -
              DENDRO_701 * grad_2_gt4[pp]) +
         DENDRO_589 * DENDRO_668 + DENDRO_593 * DENDRO_682 -
         DENDRO_594 * DENDRO_680 - 12 * grad2_1_2_alpha[pp]) +
    DENDRO_729 * grad_1_beta1[pp] + DENDRO_729 * grad_2_beta2[pp] -
    alpha[pp] * (-At4[pp] * K[pp] + DENDRO_40 * DENDRO_707 +
                 DENDRO_705 * DENDRO_720 + DENDRO_708 * DENDRO_721) +
    beta0[pp] * grad_0_At4[pp] + beta1[pp] * grad_1_At4[pp] +
    beta2[pp] * grad_2_At4[pp];
//--
At_rhs22[pp] =
    -At5[pp] * DENDRO_28 - At5[pp] * DENDRO_29 + At5[pp] * DENDRO_33 +
    DENDRO_37 * grad_2_beta0[pp] +
    DENDRO_670 *
        (DENDRO_266 *
             (4 * DENDRO_10 * DENDRO_12 * grad2_0_2_gt5[pp] +
              2.0 * DENDRO_10 * DENDRO_16 * grad2_0_0_gt5[pp] +
              2.0 * DENDRO_10 * DENDRO_18 * grad2_1_1_gt5[pp] +
              2.0 * DENDRO_10 * DENDRO_20 * grad2_2_2_gt5[pp] +
              4 * DENDRO_11 * DENDRO_322 * DENDRO_88 +
              4 * DENDRO_11 * DENDRO_332 * DENDRO_88 +
              3.0 * DENDRO_113 * DENDRO_16 * DENDRO_88 * grad_0_gt5[pp] +
              DENDRO_119 * DENDRO_228 * DENDRO_306 +
              4 * DENDRO_12 * DENDRO_88 *
                  (0.25 * DENDRO_709 + 1.0 * DENDRO_735) +
              2.0 * DENDRO_12 * DENDRO_88 * (DENDRO_709 + DENDRO_735) +
              4 * DENDRO_12 * DENDRO_88 *
                  (-1.0 * DENDRO_209 * DENDRO_91 + DENDRO_338) +
              2.0 * DENDRO_12 * DENDRO_88 *
                  (DENDRO_228 * grad_2_gt0[pp] + DENDRO_98 * grad_0_gt5[pp]) +
              2.0 * DENDRO_14 * DENDRO_313 * DENDRO_88 +
              4 * DENDRO_14 * DENDRO_342 * DENDRO_88 -
              DENDRO_140 * (0.25 * DENDRO_730 + 1.0 * DENDRO_736) -
              DENDRO_140 * (-1.0 * DENDRO_229 * DENDRO_91 + DENDRO_340) +
              4 * DENDRO_16 * DENDRO_300 * DENDRO_88 * DENDRO_98 -
              DENDRO_161 * (DENDRO_113 * DENDRO_205 + DENDRO_731) -
              DENDRO_161 * (DENDRO_141 * DENDRO_225 + DENDRO_321 * DENDRO_98) -
              DENDRO_161 * (DENDRO_176 * DENDRO_241 + DENDRO_711) -
              DENDRO_161 * (DENDRO_57 * DENDRO_726 + 0.25 * DENDRO_685) +
              4 * DENDRO_18 * DENDRO_225 * DENDRO_305 * DENDRO_88 +
              3.0 * DENDRO_18 * DENDRO_241 * DENDRO_88 * grad_1_gt5[pp] -
              DENDRO_185 * (-DENDRO_228 * DENDRO_344 - DENDRO_716) -
              DENDRO_196 * (DENDRO_192 * DENDRO_229 + DENDRO_193 * DENDRO_209 +
                            DENDRO_195 * DENDRO_246 + DENDRO_345) +
              6.0 * DENDRO_20 * DENDRO_245 * DENDRO_88 * grad_2_gt5[pp] -
              DENDRO_275 * DENDRO_702 - DENDRO_281 - DENDRO_282 - DENDRO_283 -
              DENDRO_299 - DENDRO_304 - DENDRO_307 -
              DENDRO_314 * (DENDRO_730 + DENDRO_736) -
              DENDRO_314 *
                  (DENDRO_155 * DENDRO_228 + DENDRO_225 * grad_0_gt5[pp]) -
              DENDRO_328 - DENDRO_699 * grad_0_gt5[pp] -
              DENDRO_700 * grad_1_gt5[pp] - DENDRO_701 * grad_2_gt5[pp] +
              4 * grad_2_Gt0[pp] * gt2[pp] + 4 * grad_2_Gt1[pp] * gt4[pp] +
              4 * grad_2_Gt2[pp] * gt5[pp]) +
         DENDRO_268 * DENDRO_63 - DENDRO_271 * DENDRO_723 -
         DENDRO_277 * DENDRO_68 + DENDRO_669 * gt5[pp] -
         12 * grad2_2_2_alpha[pp]) +
    DENDRO_719 * grad_2_beta1[pp] -
    alpha[pp] * (At5[pp] * DENDRO_43 * DENDRO_708 - At5[pp] * K[pp] +
                 DENDRO_46 * DENDRO_707 + DENDRO_705 * DENDRO_721) +
    beta0[pp] * grad_0_At5[pp] + beta1[pp] * grad_1_At5[pp] +
    beta2[pp] * grad_2_At5[pp];
//--
K_rhs[pp] =
    DENDRO_10 * DENDRO_257 * chi[pp] *
        (0.5 * DENDRO_744 * (DENDRO_639 + DENDRO_658 * DENDRO_742) +
         DENDRO_746 *
             (DENDRO_10 * DENDRO_101 + DENDRO_10 * DENDRO_102 +
              DENDRO_10 * DENDRO_637 -
              DENDRO_50 * (-DENDRO_654 * DENDRO_740 + grad_0_chi[pp])) +
         DENDRO_748 *
             (DENDRO_10 * DENDRO_635 + DENDRO_10 * DENDRO_636 +
              DENDRO_10 * DENDRO_89 -
              DENDRO_50 * (-DENDRO_654 * DENDRO_737 + grad_1_chi[pp])) -
         grad2_0_1_alpha[pp]) +
    DENDRO_10 * DENDRO_260 * chi[pp] *
        (0.5 * DENDRO_739 * (DENDRO_50 * DENDRO_737 * gt4[pp] + DENDRO_584) +
         DENDRO_746 * (DENDRO_10 * DENDRO_585 -
                       DENDRO_50 * (-DENDRO_589 * DENDRO_740 + grad_2_chi[pp]) +
                       DENDRO_592) +
         DENDRO_747 *
             (DENDRO_10 * DENDRO_586 + DENDRO_10 * DENDRO_587 +
              DENDRO_10 * DENDRO_588 -
              DENDRO_50 * (-DENDRO_589 * DENDRO_742 + grad_1_chi[pp])) -
         grad2_1_2_alpha[pp]) -
    DENDRO_286 * chi[pp] *
        (DENDRO_741 * (DENDRO_433 + DENDRO_51 * DENDRO_745) +
         DENDRO_744 * (DENDRO_434 + DENDRO_51 * DENDRO_743) -
         grad2_0_0_alpha[pp] +
         grad_0_alpha[pp] *
             (DENDRO_10 * DENDRO_431 + DENDRO_10 * DENDRO_432 -
              DENDRO_50 * (-0.5 * DENDRO_737 * DENDRO_75 + DENDRO_74) +
              DENDRO_73)) -
    DENDRO_288 * chi[pp] *
        (DENDRO_739 * (DENDRO_372 + DENDRO_722 * DENDRO_737) +
         DENDRO_744 * (DENDRO_376 + DENDRO_722 * DENDRO_742) -
         grad2_1_1_alpha[pp] +
         grad_1_alpha[pp] *
             (DENDRO_10 * DENDRO_373 + DENDRO_10 * DENDRO_374 +
              DENDRO_10 * DENDRO_375 -
              DENDRO_50 * (DENDRO_367 - DENDRO_368 * DENDRO_745))) -
    DENDRO_290 * chi[pp] *
        (DENDRO_739 * (DENDRO_350 + DENDRO_737 * DENDRO_738) +
         DENDRO_741 * (DENDRO_352 + DENDRO_738 * DENDRO_740) -
         grad2_2_2_alpha[pp] +
         grad_2_alpha[pp] *
             (DENDRO_10 * DENDRO_354 + DENDRO_10 * DENDRO_355 +
              DENDRO_10 * DENDRO_356 -
              DENDRO_50 * (DENDRO_274 - DENDRO_275 * DENDRO_743))) -
    DENDRO_749 * chi[pp] *
        (0.5 * DENDRO_741 * (DENDRO_496 + DENDRO_518 * DENDRO_740) +
         DENDRO_747 *
             (DENDRO_10 * DENDRO_498 + DENDRO_10 * DENDRO_499 +
              DENDRO_10 * DENDRO_500 -
              DENDRO_50 * (-DENDRO_511 * DENDRO_742 + grad_0_chi[pp])) +
         DENDRO_748 *
             (DENDRO_10 * DENDRO_493 + DENDRO_10 * DENDRO_494 +
              DENDRO_10 * DENDRO_97 -
              DENDRO_50 * (-DENDRO_511 * DENDRO_737 + grad_2_chi[pp])) -
         grad2_0_2_alpha[pp]) +
    DENDRO_768 *
        (At0[pp] * DENDRO_753 * DENDRO_754 + At1[pp] * DENDRO_763 * DENDRO_765 +
         At2[pp] * DENDRO_762 * DENDRO_763 + At3[pp] * DENDRO_754 * DENDRO_757 +
         At4[pp] * DENDRO_763 * DENDRO_767 + At5[pp] * DENDRO_754 * DENDRO_760 +
         pow(K[pp], 2)) +
    beta0[pp] * grad_0_K[pp] + beta1[pp] * grad_1_K[pp] +
    beta2[pp] * grad_2_K[pp];
//--
Gt_rhs0[pp] =
    (7.0 / 3.0) * DENDRO_10 * DENDRO_11 * grad2_0_1_beta0[pp] +
    (1.0 / 3.0) * DENDRO_10 * DENDRO_11 * grad2_1_1_beta1[pp] +
    (1.0 / 3.0) * DENDRO_10 * DENDRO_11 * grad2_1_2_beta2[pp] +
    2 * DENDRO_10 * DENDRO_14 * grad2_1_2_beta0[pp] +
    2 * DENDRO_191 * DENDRO_753 * DENDRO_774 * alpha[pp] +
    2.0 * DENDRO_223 * DENDRO_762 * DENDRO_774 * alpha[pp] +
    2.0 * DENDRO_226 * DENDRO_767 * DENDRO_774 * alpha[pp] +
    2.0 * DENDRO_227 * DENDRO_765 * DENDRO_774 * alpha[pp] +
    2 * DENDRO_229 * DENDRO_760 * DENDRO_774 * alpha[pp] +
    2 * DENDRO_232 * DENDRO_757 * DENDRO_774 * alpha[pp] +
    (2.0 / 3.0) * DENDRO_262 * DENDRO_35 * DENDRO_88 -
    4.0 / 3.0 * DENDRO_286 * grad2_0_0_beta0[pp] -
    DENDRO_288 * grad2_1_1_beta0[pp] - DENDRO_290 * grad2_2_2_beta0[pp] -
    DENDRO_357 * grad_1_beta0[pp] - DENDRO_359 * grad_0_beta0[pp] -
    DENDRO_361 * grad_2_beta0[pp] - DENDRO_753 * DENDRO_773 -
    DENDRO_762 * DENDRO_777 - DENDRO_765 * DENDRO_778 -
    DENDRO_769 * grad2_0_2_beta0[pp] - DENDRO_770 * grad2_1_2_beta1[pp] -
    DENDRO_770 * grad2_2_2_beta2[pp] - DENDRO_771 * grad2_0_1_beta1[pp] -
    DENDRO_771 * grad2_0_2_beta2[pp] -
    DENDRO_776 * (DENDRO_11 * DENDRO_780 + DENDRO_765 * DENDRO_782) -
    DENDRO_776 * (9 * DENDRO_10 * DENDRO_50 * DENDRO_753 * grad_0_chi[pp] -
                  DENDRO_16 * DENDRO_775) -
    DENDRO_776 * (9 * DENDRO_10 * DENDRO_50 * DENDRO_762 * grad_2_chi[pp] -
                  DENDRO_12 * DENDRO_779) +
    beta0[pp] * grad_0_Gt0[pp] + beta1[pp] * grad_1_Gt0[pp] +
    beta2[pp] * grad_2_Gt0[pp];
//--
Gt_rhs1[pp] =
    -DENDRO_104 * DENDRO_764 * DENDRO_789 +
    DENDRO_121 * DENDRO_762 * DENDRO_789 -
    DENDRO_202 * DENDRO_766 * DENDRO_789 +
    DENDRO_209 * DENDRO_760 * DENDRO_788 -
    DENDRO_217 * DENDRO_757 * DENDRO_788 - 2.0 / 3.0 * DENDRO_221 * DENDRO_35 +
    DENDRO_221 * grad_1_beta1[pp] + DENDRO_236 * grad_0_beta1[pp] +
    DENDRO_253 * grad_2_beta1[pp] - DENDRO_286 * grad2_0_0_beta1[pp] -
    4.0 / 3.0 * DENDRO_288 * grad2_1_1_beta1[pp] -
    DENDRO_290 * grad2_2_2_beta1[pp] + DENDRO_60 * DENDRO_753 * DENDRO_788 -
    DENDRO_749 * grad2_0_2_beta1[pp] - DENDRO_757 * DENDRO_778 +
    DENDRO_764 * DENDRO_773 + DENDRO_766 * DENDRO_777 -
    DENDRO_776 *
        (DENDRO_11 * DENDRO_775 - DENDRO_192 * DENDRO_764 * DENDRO_781) -
    DENDRO_776 *
        (DENDRO_14 * DENDRO_779 - DENDRO_195 * DENDRO_766 * DENDRO_781) +
    DENDRO_776 * (DENDRO_18 * DENDRO_780 - DENDRO_757 * DENDRO_782) +
    (7.0 / 3.0) * DENDRO_783 * grad2_0_1_beta1[pp] +
    DENDRO_784 * grad2_0_0_beta0[pp] + DENDRO_784 * grad2_0_2_beta2[pp] +
    DENDRO_785 * DENDRO_786 + (7.0 / 3.0) * DENDRO_785 * grad2_1_2_beta1[pp] +
    (1.0 / 3.0) * DENDRO_785 * grad2_2_2_beta2[pp] -
    DENDRO_787 * grad2_0_1_beta0[pp] - DENDRO_787 * grad2_1_2_beta2[pp] +
    beta0[pp] * grad_0_Gt1[pp] + beta1[pp] * grad_1_Gt1[pp] +
    beta2[pp] * grad_2_Gt1[pp];
//--
Gt_rhs2[pp] =
    2 * DENDRO_10 * DENDRO_11 * grad2_0_1_beta2[pp] +
    (1.0 / 3.0) * DENDRO_10 * DENDRO_14 * grad2_0_1_beta0[pp] +
    (1.0 / 3.0) * DENDRO_10 * DENDRO_14 * grad2_1_1_beta1[pp] +
    (7.0 / 3.0) * DENDRO_10 * DENDRO_14 * grad2_1_2_beta2[pp] +
    2.0 * DENDRO_178 * DENDRO_765 * DENDRO_774 * alpha[pp] +
    2.0 * DENDRO_180 * DENDRO_762 * DENDRO_774 * alpha[pp] +
    2 * DENDRO_194 * DENDRO_753 * DENDRO_774 * alpha[pp] +
    2.0 * DENDRO_242 * DENDRO_767 * DENDRO_774 * alpha[pp] +
    2 * DENDRO_246 * DENDRO_760 * DENDRO_774 * alpha[pp] +
    2 * DENDRO_249 * DENDRO_757 * DENDRO_774 * alpha[pp] +
    (2.0 / 3.0) * DENDRO_263 * DENDRO_35 * DENDRO_88 -
    DENDRO_286 * grad2_0_0_beta2[pp] - DENDRO_288 * grad2_1_1_beta2[pp] -
    DENDRO_290 * DENDRO_786 - 1.0 / 3.0 * DENDRO_290 * grad2_1_2_beta1[pp] -
    4.0 / 3.0 * DENDRO_290 * grad2_2_2_beta2[pp] -
    DENDRO_357 * grad_1_beta2[pp] - DENDRO_359 * grad_0_beta2[pp] -
    DENDRO_361 * grad_2_beta2[pp] - DENDRO_760 * DENDRO_777 -
    DENDRO_762 * DENDRO_773 - DENDRO_767 * DENDRO_778 -
    DENDRO_769 * grad2_0_2_beta2[pp] - DENDRO_770 * grad2_0_0_beta0[pp] -
    DENDRO_770 * grad2_0_1_beta1[pp] -
    DENDRO_776 * (DENDRO_14 * DENDRO_780 + DENDRO_767 * DENDRO_782) -
    DENDRO_776 * (9 * DENDRO_10 * DENDRO_50 * DENDRO_760 * grad_2_chi[pp] -
                  DENDRO_20 * DENDRO_779) -
    DENDRO_776 * (9 * DENDRO_10 * DENDRO_50 * DENDRO_762 * grad_0_chi[pp] -
                  DENDRO_12 * DENDRO_775) +
    beta0[pp] * grad_0_Gt2[pp] + beta1[pp] * grad_1_Gt2[pp] +
    beta2[pp] * grad_2_Gt2[pp];
// Dendro: reduced ops: 4194
// Dendro: }}}
