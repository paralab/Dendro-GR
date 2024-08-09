// CODEGEN: SSL was enabled, adding term to gauge condition!
// Codgen: generating unstage version
// Codgen: using standard gauge
// Codgen: using eta func damping
// Dendro: {{{
// Dendro: original ops: 623741
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
const double DENDRO_71 = DENDRO_33 * grad_2_chi[pp];
const double DENDRO_72 =
    -DENDRO_21 * grad_1_chi[pp] + DENDRO_35 * grad_0_chi[pp] + DENDRO_71;
const double DENDRO_73 = -DENDRO_72;
const double DENDRO_74 = 0.5 * DENDRO_73;
const double DENDRO_75 = DENDRO_31 * DENDRO_65 + DENDRO_31 * DENDRO_66 +
                         DENDRO_45 * (DENDRO_69 - DENDRO_70 * DENDRO_74) -
                         DENDRO_68;
const double DENDRO_76  = 12 * grad_0_alpha[pp];
const double DENDRO_77  = pow(chi[pp], -2);
const double DENDRO_78  = pow(grad_0_chi[pp], 2);
const double DENDRO_79  = DENDRO_77 * DENDRO_78;
const double DENDRO_80  = 4.0 * DENDRO_31;
const double DENDRO_81  = DENDRO_21 * DENDRO_80;
const double DENDRO_82  = DENDRO_81 * grad2_0_1_gt0[pp];
const double DENDRO_83  = DENDRO_23 * DENDRO_80;
const double DENDRO_84  = DENDRO_83 * grad2_1_2_gt0[pp];
const double DENDRO_85  = pow(DENDRO_30, -2);
const double DENDRO_86  = DENDRO_21 * grad_0_gt3[pp];
const double DENDRO_87  = DENDRO_35 * grad_1_gt0[pp];
const double DENDRO_88  = grad_0_gt4[pp] + grad_1_gt2[pp] - grad_2_gt1[pp];
const double DENDRO_89  = DENDRO_33 * DENDRO_88;
const double DENDRO_90  = -DENDRO_86 + DENDRO_87 + DENDRO_89;
const double DENDRO_91  = DENDRO_33 * grad_0_gt5[pp];
const double DENDRO_92  = DENDRO_35 * grad_2_gt0[pp];
const double DENDRO_93  = grad_0_gt4[pp] - grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_94  = DENDRO_21 * DENDRO_93;
const double DENDRO_95  = DENDRO_91 + DENDRO_92 - DENDRO_94;
const double DENDRO_96  = -DENDRO_21 * DENDRO_49 + DENDRO_65 + DENDRO_66;
const double DENDRO_97  = DENDRO_26 * grad_0_gt3[pp];
const double DENDRO_98  = DENDRO_21 * grad_1_gt0[pp];
const double DENDRO_99  = DENDRO_23 * DENDRO_88;
const double DENDRO_100 = DENDRO_98 + DENDRO_99;
const double DENDRO_101 = DENDRO_100 - DENDRO_97;
const double DENDRO_102 = 0.25 * grad_0_gt3[pp];
const double DENDRO_103 = 1.0 * grad_1_gt1[pp];
const double DENDRO_104 = -DENDRO_103;
const double DENDRO_105 = DENDRO_26 * DENDRO_85;
const double DENDRO_106 = 4 * DENDRO_105;
const double DENDRO_107 = DENDRO_101 * DENDRO_106 * (-DENDRO_102 - DENDRO_104);
const double DENDRO_108 = DENDRO_33 * grad_2_gt0[pp];
const double DENDRO_109 = DENDRO_39 * grad_0_gt5[pp];
const double DENDRO_110 = DENDRO_108 + DENDRO_109 - DENDRO_23 * DENDRO_93;
const double DENDRO_111 = 0.25 * grad_0_gt5[pp];
const double DENDRO_112 = 1.0 * grad_2_gt2[pp];
const double DENDRO_113 = -DENDRO_112;
const double DENDRO_114 = DENDRO_111 + DENDRO_113;
const double DENDRO_115 = DENDRO_39 * DENDRO_85;
const double DENDRO_116 = 4 * DENDRO_115;
const double DENDRO_117 =
    DENDRO_21 * grad_2_gt0[pp] + DENDRO_23 * grad_0_gt5[pp];
const double DENDRO_118 = DENDRO_117 - DENDRO_26 * DENDRO_93;
const double DENDRO_119 = 0.25 * grad_0_gt4[pp];
const double DENDRO_120 = -DENDRO_119;
const double DENDRO_121 = 0.25 * grad_1_gt2[pp];
const double DENDRO_122 = 0.75 * grad_2_gt1[pp];
const double DENDRO_123 =
    DENDRO_116 * DENDRO_118 * (DENDRO_120 + DENDRO_121 + DENDRO_122);
const double DENDRO_124 = 0.75 * grad_1_gt2[pp];
const double DENDRO_125 = 0.25 * grad_2_gt1[pp];
const double DENDRO_126 = DENDRO_120 + DENDRO_124 + DENDRO_125;
const double DENDRO_127 = DENDRO_33 * grad_1_gt0[pp];
const double DENDRO_128 = DENDRO_39 * DENDRO_88;
const double DENDRO_129 = DENDRO_127 + DENDRO_128 - DENDRO_23 * grad_0_gt3[pp];
const double DENDRO_130 = DENDRO_35 * DENDRO_55;
const double DENDRO_131 = 4 * DENDRO_85;
const double DENDRO_132 = DENDRO_130 * DENDRO_131 * (DENDRO_47 + DENDRO_48);
const double DENDRO_133 = DENDRO_50 + DENDRO_51;
const double DENDRO_134 = 0.25 * grad_1_gt0[pp];
const double DENDRO_135 = DENDRO_134 * DENDRO_95;
const double DENDRO_136 = DENDRO_23 * DENDRO_85;
const double DENDRO_137 = 4 * DENDRO_136;
const double DENDRO_138 = 0.25 * grad_2_gt0[pp];
const double DENDRO_139 = DENDRO_138 * DENDRO_90;
const double DENDRO_140 = DENDRO_90 * grad_0_gt0[pp];
const double DENDRO_141 = DENDRO_96 * grad_1_gt0[pp];
const double DENDRO_142 = DENDRO_21 * DENDRO_85;
const double DENDRO_143 = 2.0 * DENDRO_142;
const double DENDRO_144 = DENDRO_95 * grad_0_gt0[pp];
const double DENDRO_145 = DENDRO_96 * grad_2_gt0[pp];
const double DENDRO_146 = DENDRO_101 * grad_1_gt0[pp];
const double DENDRO_147 = DENDRO_55 * grad_0_gt3[pp];
const double DENDRO_148 = DENDRO_146 + DENDRO_147;
const double DENDRO_149 = DENDRO_110 * grad_2_gt0[pp];
const double DENDRO_150 = DENDRO_61 * grad_0_gt5[pp];
const double DENDRO_151 = DENDRO_102 * DENDRO_118;
const double DENDRO_152 = -grad_0_gt4[pp] + grad_1_gt2[pp] + grad_2_gt1[pp];
const double DENDRO_153 = 0.5 * DENDRO_152;
const double DENDRO_154 = DENDRO_101 * DENDRO_153 + DENDRO_151;
const double DENDRO_155 = DENDRO_111 * DENDRO_129;
const double DENDRO_156 = 0.5 * DENDRO_110;
const double DENDRO_157 = 0.25 * DENDRO_140;
const double DENDRO_158 = 4 * DENDRO_142;
const double DENDRO_159 = 0.25 * DENDRO_144;
const double DENDRO_160 = DENDRO_118 * grad_1_gt0[pp];
const double DENDRO_161 = DENDRO_55 * DENDRO_93;
const double DENDRO_162 = DENDRO_33 * DENDRO_85;
const double DENDRO_163 = 2.0 * DENDRO_162;
const double DENDRO_164 = DENDRO_163 * (DENDRO_160 + DENDRO_161);
const double DENDRO_165 = DENDRO_129 * grad_2_gt0[pp];
const double DENDRO_166 = DENDRO_61 * DENDRO_88;
const double DENDRO_167 = 0.5 * grad_0_gt3[pp];
const double DENDRO_168 = DENDRO_104 + DENDRO_167;
const double DENDRO_169 = DENDRO_118 * DENDRO_168;
const double DENDRO_170 = DENDRO_101 * DENDRO_93;
const double DENDRO_171 = 0.25 * DENDRO_170;
const double DENDRO_172 = -DENDRO_171;
const double DENDRO_173 = 0.5 * grad_0_gt5[pp];
const double DENDRO_174 = DENDRO_113 + DENDRO_173;
const double DENDRO_175 = -DENDRO_129;
const double DENDRO_176 = DENDRO_174 * DENDRO_175;
const double DENDRO_177 = -DENDRO_110;
const double DENDRO_178 = DENDRO_177 * DENDRO_88;
const double DENDRO_179 = -0.25 * DENDRO_178;
const double DENDRO_180 = 0.5 * DENDRO_49;
const double DENDRO_181 = DENDRO_118 * DENDRO_180;
const double DENDRO_182 = 4 * DENDRO_162;
const double DENDRO_183 = DENDRO_182 * (DENDRO_152 * DENDRO_55 + DENDRO_181);
const double DENDRO_184 = 0.5 * DENDRO_52;
const double DENDRO_185 = DENDRO_101 * DENDRO_180;
const double DENDRO_186 = -2 * DENDRO_168 * DENDRO_55 + DENDRO_185;
const double DENDRO_187 = -grad2_0_0_chi[pp];
const double DENDRO_188 = -DENDRO_96;
const double DENDRO_189 = DENDRO_31 * grad_0_chi[pp];
const double DENDRO_190 = DENDRO_31 * grad_1_chi[pp];
const double DENDRO_191 = -DENDRO_61;
const double DENDRO_192 = DENDRO_31 * grad_2_chi[pp];
const double DENDRO_193 = 2 * DENDRO_45;
const double DENDRO_194 = DENDRO_118 * DENDRO_33;
const double DENDRO_195 = DENDRO_26 * grad_2_gt3[pp];
const double DENDRO_196 = DENDRO_23 * grad_1_gt5[pp];
const double DENDRO_197 = DENDRO_152 * DENDRO_21;
const double DENDRO_198 = DENDRO_196 + DENDRO_197;
const double DENDRO_199 = -DENDRO_195 + DENDRO_198;
const double DENDRO_200 = 0.5 * grad_2_gt5[pp];
const double DENDRO_201 = DENDRO_200 * DENDRO_23;
const double DENDRO_202 = 0.5 * grad_1_gt5[pp];
const double DENDRO_203 = 1.0 * grad_2_gt4[pp];
const double DENDRO_204 = -DENDRO_203;
const double DENDRO_205 = DENDRO_202 + DENDRO_204;
const double DENDRO_206 =
    -DENDRO_174 * DENDRO_21 + DENDRO_201 + DENDRO_205 * DENDRO_26;
const double DENDRO_207 = DENDRO_206 * DENDRO_39;
const double DENDRO_208 = 0.5 * grad_1_gt3[pp];
const double DENDRO_209 = DENDRO_208 * DENDRO_26;
const double DENDRO_210 = DENDRO_168 * DENDRO_21;
const double DENDRO_211 = 1.0 * grad_1_gt4[pp];
const double DENDRO_212 = 0.5 * grad_2_gt3[pp];
const double DENDRO_213 = DENDRO_211 - DENDRO_212;
const double DENDRO_214 = DENDRO_209 + DENDRO_210 - DENDRO_213 * DENDRO_23;
const double DENDRO_215 = -DENDRO_214;
const double DENDRO_216 = DENDRO_215 * DENDRO_26;
const double DENDRO_217 = DENDRO_130 + DENDRO_207 + DENDRO_216;
const double DENDRO_218 =
    DENDRO_85 * (-1.0 * DENDRO_101 * DENDRO_21 + DENDRO_194 -
                 1.0 * DENDRO_199 * DENDRO_23 + DENDRO_217);
const double DENDRO_219 = 2.0 * grad_1_gt0[pp];
const double DENDRO_220 = -DENDRO_95;
const double DENDRO_221 = DENDRO_220 * DENDRO_33;
const double DENDRO_222 = DENDRO_152 * DENDRO_35 - DENDRO_21 * grad_2_gt3[pp] +
                          DENDRO_33 * grad_1_gt5[pp];
const double DENDRO_223 = -DENDRO_222;
const double DENDRO_224 = -DENDRO_90;
const double DENDRO_225 =
    -DENDRO_174 * DENDRO_35 + DENDRO_200 * DENDRO_33 + DENDRO_205 * DENDRO_21;
const double DENDRO_226 = -DENDRO_225;
const double DENDRO_227 = DENDRO_226 * DENDRO_39;
const double DENDRO_228 = DENDRO_208 * DENDRO_21;
const double DENDRO_229 =
    DENDRO_168 * DENDRO_35 - DENDRO_213 * DENDRO_33 + DENDRO_228;
const double DENDRO_230 = DENDRO_229 * DENDRO_26;
const double DENDRO_231 = DENDRO_188 * DENDRO_35;
const double DENDRO_232 = DENDRO_227 + DENDRO_230 + DENDRO_231;
const double DENDRO_233 =
    DENDRO_85 * (-1.0 * DENDRO_21 * DENDRO_224 + DENDRO_221 -
                 1.0 * DENDRO_223 * DENDRO_23 + DENDRO_232);
const double DENDRO_234 = 2.0 * grad_0_gt0[pp];
const double DENDRO_235 = DENDRO_177 * DENDRO_33;
const double DENDRO_236 = DENDRO_39 * grad_1_gt5[pp];
const double DENDRO_237 = DENDRO_152 * DENDRO_33;
const double DENDRO_238 = -DENDRO_23 * grad_2_gt3[pp] + DENDRO_236 + DENDRO_237;
const double DENDRO_239 = -DENDRO_238;
const double DENDRO_240 = DENDRO_200 * DENDRO_39;
const double DENDRO_241 = DENDRO_205 * DENDRO_23;
const double DENDRO_242 = -DENDRO_174 * DENDRO_33 + DENDRO_240 + DENDRO_241;
const double DENDRO_243 = -DENDRO_242;
const double DENDRO_244 = DENDRO_243 * DENDRO_39;
const double DENDRO_245 = DENDRO_208 * DENDRO_23;
const double DENDRO_246 =
    DENDRO_168 * DENDRO_33 - DENDRO_213 * DENDRO_39 + DENDRO_245;
const double DENDRO_247 = DENDRO_246 * DENDRO_26;
const double DENDRO_248 = DENDRO_191 * DENDRO_35;
const double DENDRO_249 = DENDRO_244 + DENDRO_247 + DENDRO_248;
const double DENDRO_250 =
    DENDRO_85 * (-1.0 * DENDRO_175 * DENDRO_21 - 1.0 * DENDRO_23 * DENDRO_239 +
                 DENDRO_235 + DENDRO_249);
const double DENDRO_251 = 2.0 * grad_2_gt0[pp];
const double DENDRO_252 = 3 * DENDRO_45;
const double DENDRO_253 = grad_0_chi[pp] * grad_1_chi[pp];
const double DENDRO_254 = 2 * DENDRO_21;
const double DENDRO_255 = DENDRO_252 * grad_2_chi[pp];
const double DENDRO_256 = 2 * DENDRO_33;
const double DENDRO_257 = 2 * DENDRO_23;
const double DENDRO_258 = pow(grad_1_chi[pp], 2);
const double DENDRO_259 = pow(grad_2_chi[pp], 2);
const double DENDRO_260 = DENDRO_101 * DENDRO_21 - 1.0 * DENDRO_194 +
                          DENDRO_199 * DENDRO_23 - DENDRO_217;
const double DENDRO_261 = DENDRO_21 * DENDRO_224 - 1.0 * DENDRO_221 +
                          DENDRO_223 * DENDRO_23 - DENDRO_232;
const double DENDRO_262 = DENDRO_175 * DENDRO_21 + DENDRO_23 * DENDRO_239 -
                          1.0 * DENDRO_235 - DENDRO_249;
const double DENDRO_263 =
    2 * DENDRO_189 * DENDRO_261 + 2 * DENDRO_190 * DENDRO_260 +
    2 * DENDRO_192 * DENDRO_262 -
    DENDRO_254 * (-DENDRO_252 * DENDRO_253 + 2 * grad2_0_1_chi[pp]) +
    DENDRO_256 * (-DENDRO_255 * grad_0_chi[pp] + 2 * grad2_0_2_chi[pp]) -
    DENDRO_257 * (-DENDRO_255 * grad_1_chi[pp] + 2 * grad2_1_2_chi[pp]) +
    DENDRO_26 * (-DENDRO_252 * DENDRO_258 + 2 * grad2_1_1_chi[pp]) +
    DENDRO_35 * (-DENDRO_252 * DENDRO_78 + 2 * grad2_0_0_chi[pp]) +
    DENDRO_39 * (-DENDRO_252 * DENDRO_259 + 2 * grad2_2_2_chi[pp]);
const double DENDRO_264 = DENDRO_45 * DENDRO_70;
const double DENDRO_265 = 3 * alpha[pp];
const double DENDRO_266 = DENDRO_45 * gt5[pp];
const double DENDRO_267 = DENDRO_206 + DENDRO_266 * DENDRO_44;
const double DENDRO_268 = 4 * grad_1_alpha[pp];
const double DENDRO_269 = DENDRO_268 * DENDRO_31;
const double DENDRO_270 = DENDRO_225 - 0.5 * DENDRO_45 * DENDRO_73 * gt5[pp];
const double DENDRO_271 = 4 * grad_0_alpha[pp];
const double DENDRO_272 = DENDRO_271 * DENDRO_31;
const double DENDRO_273 = 1.0 * grad_2_chi[pp];
const double DENDRO_274 = DENDRO_31 * gt5[pp];
const double DENDRO_275 = 0.5 * DENDRO_60;
const double DENDRO_276 = -DENDRO_174 * DENDRO_31 * DENDRO_33 +
                          DENDRO_240 * DENDRO_31 + DENDRO_241 * DENDRO_31 +
                          DENDRO_45 * (DENDRO_273 - DENDRO_274 * DENDRO_275);
const double DENDRO_277 = 4 * grad_2_alpha[pp];
const double DENDRO_278 = 4 * gt2[pp];
const double DENDRO_279 = 4 * gt4[pp];
const double DENDRO_280 = DENDRO_259 * DENDRO_77;
const double DENDRO_281 = DENDRO_81 * grad2_0_1_gt5[pp];
const double DENDRO_282 = DENDRO_83 * grad2_1_2_gt5[pp];
const double DENDRO_283 = DENDRO_31 * DENDRO_33;
const double DENDRO_284 = 4 * DENDRO_283;
const double DENDRO_285 = DENDRO_31 * DENDRO_35;
const double DENDRO_286 = 2.0 * DENDRO_285;
const double DENDRO_287 = DENDRO_26 * DENDRO_31;
const double DENDRO_288 = 2.0 * DENDRO_287;
const double DENDRO_289 = DENDRO_31 * DENDRO_39;
const double DENDRO_290 = 2.0 * DENDRO_289;
const double DENDRO_291 = DENDRO_177 * grad_0_gt5[pp];
const double DENDRO_292 = DENDRO_35 * DENDRO_85;
const double DENDRO_293 = 3.0 * DENDRO_292;
const double DENDRO_294 = DENDRO_239 * grad_1_gt5[pp];
const double DENDRO_295 = 3.0 * DENDRO_105;
const double DENDRO_296 = 6.0 * DENDRO_85;
const double DENDRO_297 = 0.25 * grad_2_gt3[pp];
const double DENDRO_298 = DENDRO_106 * DENDRO_199 * (DENDRO_211 - DENDRO_297);
const double DENDRO_299 = -DENDRO_138 + DENDRO_50;
const double DENDRO_300 = 4 * DENDRO_292;
const double DENDRO_301 = 0.75 * grad_0_gt4[pp];
const double DENDRO_302 = -DENDRO_125;
const double DENDRO_303 =
    DENDRO_118 * DENDRO_300 * (DENDRO_121 + DENDRO_301 + DENDRO_302);
const double DENDRO_304 = DENDRO_119 + DENDRO_124 + DENDRO_302;
const double DENDRO_305 = DENDRO_112 + DENDRO_173;
const double DENDRO_306 = DENDRO_131 * DENDRO_207 * (DENDRO_202 + DENDRO_203);
const double DENDRO_307 = DENDRO_111 * DENDRO_239;
const double DENDRO_308 = 0.25 * grad_1_gt5[pp];
const double DENDRO_309 = DENDRO_177 * DENDRO_308;
const double DENDRO_310 = DENDRO_199 * grad_1_gt5[pp];
const double DENDRO_311 = DENDRO_206 * grad_2_gt3[pp];
const double DENDRO_312 = DENDRO_310 + DENDRO_311;
const double DENDRO_313 = 2.0 * DENDRO_136;
const double DENDRO_314 = DENDRO_226 * grad_2_gt0[pp];
const double DENDRO_315 = DENDRO_177 * grad_2_gt5[pp];
const double DENDRO_316 = DENDRO_243 * grad_0_gt5[pp];
const double DENDRO_317 = DENDRO_239 * grad_2_gt5[pp];
const double DENDRO_318 = DENDRO_243 * grad_1_gt5[pp];
const double DENDRO_319 = DENDRO_118 * DENDRO_297;
const double DENDRO_320 = 0.5 * DENDRO_88;
const double DENDRO_321 = DENDRO_199 * DENDRO_320 + DENDRO_319;
const double DENDRO_322 = 0.25 * DENDRO_315;
const double DENDRO_323 = 0.25 * DENDRO_317;
const double DENDRO_324 = DENDRO_138 * DENDRO_223;
const double DENDRO_325 = DENDRO_118 * grad_1_gt5[pp];
const double DENDRO_326 = DENDRO_206 * DENDRO_93;
const double DENDRO_327 = DENDRO_163 * (DENDRO_325 + DENDRO_326);
const double DENDRO_328 = DENDRO_118 * DENDRO_213;
const double DENDRO_329 = DENDRO_199 * DENDRO_93;
const double DENDRO_330 = 0.25 * DENDRO_329;
const double DENDRO_331 = DENDRO_328 + DENDRO_330;
const double DENDRO_332 = DENDRO_152 * DENDRO_226;
const double DENDRO_333 = DENDRO_223 * DENDRO_52;
const double DENDRO_334 = DENDRO_152 * DENDRO_220;
const double DENDRO_335 = 0.25 * DENDRO_334;
const double DENDRO_336 = 0.5 * DENDRO_205;
const double DENDRO_337 = DENDRO_118 * DENDRO_336;
const double DENDRO_338 = 0.5 * DENDRO_174;
const double DENDRO_339 = DENDRO_223 * DENDRO_338;
const double DENDRO_340 = -DENDRO_199 * DENDRO_336;
const double DENDRO_341 = 2 * DENDRO_206 * DENDRO_213 + DENDRO_340;
const double DENDRO_342 = -DENDRO_220 * DENDRO_338;
const double DENDRO_343 = 2 * DENDRO_52;
const double DENDRO_344 = -grad2_2_2_chi[pp];
const double DENDRO_345 = -DENDRO_33;
const double DENDRO_346 = -DENDRO_174;
const double DENDRO_347 = -DENDRO_35;
const double DENDRO_348 = -DENDRO_205;
const double DENDRO_349 =
    DENDRO_200 * DENDRO_345 + DENDRO_21 * DENDRO_348 + DENDRO_346 * DENDRO_347;
const double DENDRO_350 = -DENDRO_26;
const double DENDRO_351 =
    DENDRO_201 + DENDRO_21 * DENDRO_346 + DENDRO_348 * DENDRO_350;
const double DENDRO_352 = -DENDRO_39;
const double DENDRO_353 = DENDRO_200 * DENDRO_352;
const double DENDRO_354 = DENDRO_345 * DENDRO_346;
const double DENDRO_355 = DENDRO_23 * DENDRO_348;
const double DENDRO_356 = DENDRO_260 * DENDRO_85;
const double DENDRO_357 = 2.0 * DENDRO_356;
const double DENDRO_358 = DENDRO_261 * DENDRO_85;
const double DENDRO_359 = 2.0 * DENDRO_358;
const double DENDRO_360 = DENDRO_262 * DENDRO_85;
const double DENDRO_361 = 2.0 * DENDRO_360;
const double DENDRO_362 = -DENDRO_263;
const double DENDRO_363 = DENDRO_362 * DENDRO_45;
const double DENDRO_364 = DENDRO_45 * gt3[pp];
const double DENDRO_365 = DENDRO_277 * DENDRO_31;
const double DENDRO_366 = 1.0 * grad_1_chi[pp];
const double DENDRO_367 = DENDRO_31 * gt3[pp];
const double DENDRO_368 = DENDRO_209 * DENDRO_31 + DENDRO_210 * DENDRO_31 -
                          DENDRO_213 * DENDRO_23 * DENDRO_31 +
                          DENDRO_45 * (DENDRO_366 - DENDRO_367 * DENDRO_44);
const double DENDRO_369 = -grad2_1_1_chi[pp];
const double DENDRO_370 = -DENDRO_168;
const double DENDRO_371 =
    DENDRO_213 * DENDRO_345 + DENDRO_228 + DENDRO_347 * DENDRO_370;
const double DENDRO_372 = DENDRO_208 * DENDRO_350;
const double DENDRO_373 = DENDRO_21 * DENDRO_370;
const double DENDRO_374 = DENDRO_213 * DENDRO_23;
const double DENDRO_375 =
    DENDRO_213 * DENDRO_352 + DENDRO_245 + DENDRO_345 * DENDRO_370;
const double DENDRO_376 = DENDRO_223 * DENDRO_49;
const double DENDRO_377 = DENDRO_152 * DENDRO_224;
const double DENDRO_378 = 0.25 * DENDRO_377;
const double DENDRO_379 = DENDRO_134 * DENDRO_223;
const double DENDRO_380 = 0.5 * DENDRO_93;
const double DENDRO_381 = DENDRO_175 * DENDRO_308;
const double DENDRO_382 = DENDRO_175 * DENDRO_205;
const double DENDRO_383 = DENDRO_239 * DENDRO_88;
const double DENDRO_384 = -0.25 * DENDRO_383;
const double DENDRO_385 = 0.5 * DENDRO_168;
const double DENDRO_386 = DENDRO_223 * DENDRO_385;
const double DENDRO_387 = 0.5 * DENDRO_213;
const double DENDRO_388 = DENDRO_239 * DENDRO_387;
const double DENDRO_389 = 2 * DENDRO_205 * DENDRO_246;
const double DENDRO_390 = DENDRO_199 * grad_1_gt3[pp];
const double DENDRO_391 = 0.25 * DENDRO_390;
const double DENDRO_392 = DENDRO_215 * grad_2_gt3[pp];
const double DENDRO_393 = DENDRO_175 * DENDRO_387;
const double DENDRO_394 = -DENDRO_224 * DENDRO_385;
const double DENDRO_395 = 2 * DENDRO_229 * DENDRO_49;
const double DENDRO_396 = DENDRO_101 * grad_1_gt3[pp];
const double DENDRO_397 = 0.25 * DENDRO_396;
const double DENDRO_398 = DENDRO_215 * grad_0_gt3[pp];
const double DENDRO_399 = DENDRO_229 * grad_1_gt0[pp];
const double DENDRO_400 = DENDRO_152 * DENDRO_229;
const double DENDRO_401 = DENDRO_246 * grad_1_gt5[pp];
const double DENDRO_402 = DENDRO_246 * DENDRO_88;
const double DENDRO_403 = DENDRO_204 + DENDRO_308;
const double DENDRO_404 = -DENDRO_121;
const double DENDRO_405 = DENDRO_116 * (DENDRO_119 + DENDRO_122 + DENDRO_404);
const double DENDRO_406 = DENDRO_300 * (-DENDRO_134 + DENDRO_47);
const double DENDRO_407 = DENDRO_300 * (DENDRO_125 + DENDRO_301 + DENDRO_404);
const double DENDRO_408 = 4 * gt1[pp];
const double DENDRO_409 = DENDRO_102 * DENDRO_199;
const double DENDRO_410 = DENDRO_101 * DENDRO_297;
const double DENDRO_411 = DENDRO_101 * grad_0_gt3[pp];
const double DENDRO_412 = DENDRO_199 * grad_2_gt3[pp];
const double DENDRO_413 = 3.0 * DENDRO_115;
const double DENDRO_414 =
    -DENDRO_131 * DENDRO_230 * (DENDRO_103 + DENDRO_167) -
    DENDRO_131 * DENDRO_247 * (DENDRO_211 + DENDRO_212) -
    DENDRO_182 * (DENDRO_101 * DENDRO_212 + DENDRO_409) -
    DENDRO_182 * (DENDRO_167 * DENDRO_199 + DENDRO_410) -
    DENDRO_258 * DENDRO_77 + DENDRO_279 * grad_1_Gt2[pp] +
    DENDRO_284 * grad2_0_2_gt3[pp] + DENDRO_286 * grad2_0_0_gt3[pp] +
    DENDRO_288 * grad2_1_1_gt3[pp] + DENDRO_290 * grad2_2_2_gt3[pp] -
    DENDRO_293 * DENDRO_411 + DENDRO_408 * grad_1_Gt0[pp] -
    DENDRO_412 * DENDRO_413 - DENDRO_81 * grad2_0_1_gt3[pp] -
    DENDRO_83 * grad2_1_2_gt3[pp] + 4 * grad_1_Gt1[pp] * gt3[pp];
const double DENDRO_415 = DENDRO_224 * grad_1_gt0[pp];
const double DENDRO_416 = DENDRO_220 * grad_2_gt0[pp];
const double DENDRO_417 = DENDRO_134 * DENDRO_220;
const double DENDRO_418 = DENDRO_138 * DENDRO_224;
const double DENDRO_419 = DENDRO_224 * grad_0_gt0[pp];
const double DENDRO_420 = DENDRO_188 * grad_1_gt0[pp];
const double DENDRO_421 = DENDRO_220 * grad_0_gt0[pp];
const double DENDRO_422 = DENDRO_188 * grad_2_gt0[pp];
const double DENDRO_423 = DENDRO_191 * grad_0_gt5[pp];
const double DENDRO_424 = 0.25 * DENDRO_419;
const double DENDRO_425 = 0.25 * DENDRO_421;
const double DENDRO_426 = DENDRO_111 * DENDRO_175;
const double DENDRO_427 = DENDRO_191 * DENDRO_88;
const double DENDRO_428 = DENDRO_175 * DENDRO_184;
const double DENDRO_429 = DENDRO_177 * DENDRO_184;
const double DENDRO_430 = DENDRO_347 * DENDRO_53;
const double DENDRO_431 = DENDRO_345 * DENDRO_52;
const double DENDRO_432 = DENDRO_350 * DENDRO_49 + DENDRO_54;
const double DENDRO_433 =
    DENDRO_23 * DENDRO_49 + DENDRO_345 * DENDRO_53 + DENDRO_352 * DENDRO_52;
const double DENDRO_434 = DENDRO_199 * grad_1_gt0[pp];
const double DENDRO_435 = DENDRO_118 * grad_0_gt3[pp];
const double DENDRO_436 = DENDRO_101 * DENDRO_88;
const double DENDRO_437 = DENDRO_435 + DENDRO_436;
const double DENDRO_438 = DENDRO_239 * grad_2_gt0[pp];
const double DENDRO_439 = DENDRO_175 * grad_0_gt5[pp];
const double DENDRO_440 = DENDRO_178 + DENDRO_439;
const double DENDRO_441 = DENDRO_111 * DENDRO_220;
const double DENDRO_442 = -DENDRO_174 * DENDRO_243;
const double DENDRO_443 = DENDRO_153 * DENDRO_206 + 0.25 * DENDRO_325;
const double DENDRO_444 = DENDRO_224 * DENDRO_88;
const double DENDRO_445 = 0.25 * DENDRO_444;
const double DENDRO_446 = DENDRO_101 * DENDRO_387;
const double DENDRO_447 = -DENDRO_199 * DENDRO_385;
const double DENDRO_448 = DENDRO_446 + DENDRO_447;
const double DENDRO_449 = DENDRO_188 * DENDRO_52;
const double DENDRO_450 = 0.25 * DENDRO_160 + DENDRO_320 * DENDRO_55;
const double DENDRO_451 = DENDRO_138 * DENDRO_177;
const double DENDRO_452 = 0.25 * DENDRO_118;
const double DENDRO_453 = DENDRO_152 * DENDRO_452 + DENDRO_202 * DENDRO_55;
const double DENDRO_454 = DENDRO_191 * DENDRO_200;
const double DENDRO_455 = -DENDRO_177 * DENDRO_338;
const double DENDRO_456 = DENDRO_138 * DENDRO_220;
const double DENDRO_457 = DENDRO_226 * DENDRO_53;
const double DENDRO_458 = DENDRO_184 * DENDRO_220 + DENDRO_457;
const double DENDRO_459 = DENDRO_452 * DENDRO_93;
const double DENDRO_460 = DENDRO_206 * DENDRO_48 + DENDRO_452 * DENDRO_88;
const double DENDRO_461 = DENDRO_168 * DENDRO_206;
const double DENDRO_462 = DENDRO_226 * DENDRO_48;
const double DENDRO_463 = DENDRO_111 * DENDRO_224 + DENDRO_462;
const double DENDRO_464 = DENDRO_153 * DENDRO_243;
const double DENDRO_465 = DENDRO_307 + DENDRO_309;
const double DENDRO_466 = DENDRO_152 * DENDRO_199;
const double DENDRO_467 = 0.25 * DENDRO_466;
const double DENDRO_468 = DENDRO_167 * DENDRO_206;
const double DENDRO_469 = DENDRO_101 * DENDRO_308 + DENDRO_468;
const double DENDRO_470 = 0.25 * DENDRO_88;
const double DENDRO_471 = DENDRO_220 * DENDRO_470 + DENDRO_462;
const double DENDRO_472 = DENDRO_239 * DENDRO_338;
const double DENDRO_473 = -DENDRO_472;
const double DENDRO_474 = 0.25 * DENDRO_175;
const double DENDRO_475 = DENDRO_474 * grad_2_gt5[pp];
const double DENDRO_476 = DENDRO_243 * DENDRO_320 + DENDRO_475;
const double DENDRO_477 = DENDRO_180 * DENDRO_199;
const double DENDRO_478 = -0.5 * DENDRO_169 + DENDRO_213 * DENDRO_55;
const double DENDRO_479 = 0.25 * DENDRO_223;
const double DENDRO_480 = DENDRO_479 * grad_0_gt0[pp];
const double DENDRO_481 = DENDRO_184 * DENDRO_224;
const double DENDRO_482 = DENDRO_480 + DENDRO_481;
const double DENDRO_483 = DENDRO_188 * DENDRO_320 + DENDRO_417;
const double DENDRO_484 = DENDRO_191 * DENDRO_202;
const double DENDRO_485 = 0.25 * DENDRO_152;
const double DENDRO_486 = DENDRO_177 * DENDRO_485;
const double DENDRO_487 = DENDRO_152 * DENDRO_239;
const double DENDRO_488 = DENDRO_175 * grad_1_gt5[pp];
const double DENDRO_489 = DENDRO_383 + DENDRO_488;
const double DENDRO_490 = 1.0 * DENDRO_105;
const double DENDRO_491 = -grad2_0_2_chi[pp];
const double DENDRO_492 = DENDRO_345 * grad_0_gt5[pp];
const double DENDRO_493 = DENDRO_347 * grad_2_gt0[pp];
const double DENDRO_494 = 0.5 * DENDRO_189;
const double DENDRO_495 = DENDRO_117 + DENDRO_350 * DENDRO_93;
const double DENDRO_496 = 0.5 * DENDRO_190;
const double DENDRO_497 = DENDRO_352 * grad_0_gt5[pp];
const double DENDRO_498 = DENDRO_345 * grad_2_gt0[pp];
const double DENDRO_499 = DENDRO_23 * DENDRO_93;
const double DENDRO_500 = 0.5 * DENDRO_192;
const double DENDRO_501 = 2.0 * gt2[pp];
const double DENDRO_502 = 2.0 * gt4[pp];
const double DENDRO_503 = 2.0 * gt5[pp];
const double DENDRO_504 = 2.0 * grad_2_Gt0[pp];
const double DENDRO_505 = 2.0 * grad_2_Gt1[pp];
const double DENDRO_506 = DENDRO_77 * grad_2_chi[pp];
const double DENDRO_507 = DENDRO_506 * grad_0_chi[pp];
const double DENDRO_508 = DENDRO_81 * grad2_0_1_gt2[pp];
const double DENDRO_509 = DENDRO_83 * grad2_1_2_gt2[pp];
const double DENDRO_510 = DENDRO_31 * gt2[pp];
const double DENDRO_511 =
    -DENDRO_193 *
        (DENDRO_491 + DENDRO_494 * (DENDRO_492 + DENDRO_493 + DENDRO_94) +
         DENDRO_495 * DENDRO_496 +
         DENDRO_500 * (DENDRO_497 + DENDRO_498 + DENDRO_499)) +
    DENDRO_284 * grad2_0_2_gt2[pp] + DENDRO_286 * grad2_0_0_gt2[pp] +
    DENDRO_288 * grad2_1_1_gt2[pp] + DENDRO_290 * grad2_2_2_gt2[pp] +
    DENDRO_357 * grad_1_gt2[pp] + DENDRO_359 * grad_0_gt2[pp] +
    DENDRO_361 * grad_2_gt2[pp] + DENDRO_363 * DENDRO_510 +
    DENDRO_501 * grad_0_Gt0[pp] + DENDRO_501 * grad_2_Gt2[pp] +
    DENDRO_502 * grad_0_Gt1[pp] + DENDRO_503 * grad_0_Gt2[pp] +
    DENDRO_504 * gt0[pp] + DENDRO_505 * gt1[pp] - DENDRO_507 - DENDRO_508 -
    DENDRO_509;
const double DENDRO_512 =
    -DENDRO_21 * DENDRO_31 * DENDRO_93 + DENDRO_31 * DENDRO_91 +
    DENDRO_31 * DENDRO_92 +
    DENDRO_45 * (-DENDRO_510 * DENDRO_73 + grad_2_chi[pp]);
const double DENDRO_513 = 2.0 * grad_0_alpha[pp];
const double DENDRO_514 =
    DENDRO_108 * DENDRO_31 + DENDRO_109 * DENDRO_31 -
    DENDRO_23 * DENDRO_31 * DENDRO_93 +
    DENDRO_45 * (-DENDRO_510 * DENDRO_60 + grad_0_chi[pp]);
const double DENDRO_515 = 2.0 * grad_2_alpha[pp];
const double DENDRO_516 = 2.0 * grad_1_alpha[pp];
const double DENDRO_517 = DENDRO_45 * gt2[pp];
const double DENDRO_518 = DENDRO_31 * (DENDRO_118 + DENDRO_43 * DENDRO_517);
const double DENDRO_519 = -DENDRO_512 * DENDRO_513 - DENDRO_514 * DENDRO_515 +
                          DENDRO_516 * DENDRO_518 - 4 * grad2_0_2_alpha[pp];
const double DENDRO_520 = DENDRO_223 * grad_2_gt0[pp];
const double DENDRO_521 = DENDRO_224 * grad_0_gt5[pp] + DENDRO_520;
const double DENDRO_522 = DENDRO_118 * grad_2_gt3[pp];
const double DENDRO_523 = DENDRO_101 * grad_1_gt5[pp] + DENDRO_522;
const double DENDRO_524 = DENDRO_466 + DENDRO_523;
const double DENDRO_525 = DENDRO_116 * (-DENDRO_337 + DENDRO_443);
const double DENDRO_526 = 0.25 * DENDRO_487;
const double DENDRO_527 = DENDRO_106 * (DENDRO_410 + DENDRO_448);
const double DENDRO_528 = DENDRO_300 * (0.5 * DENDRO_161 + DENDRO_450);
const double DENDRO_529 = DENDRO_182 * (-DENDRO_205 * DENDRO_55 + DENDRO_460);
const double DENDRO_530 = DENDRO_111 * DENDRO_177 + DENDRO_454;
const double DENDRO_531 = DENDRO_182 * (DENDRO_453 + DENDRO_459);
const double DENDRO_532 = DENDRO_101 * DENDRO_336;
const double DENDRO_533 =
    -0.5 * DENDRO_118 * DENDRO_213 + DENDRO_461 + DENDRO_532;
const double DENDRO_534 = DENDRO_309 + DENDRO_475;
const double DENDRO_535 = -DENDRO_224 * DENDRO_338;
const double DENDRO_536 = DENDRO_171 + DENDRO_478;
const double DENDRO_537 = DENDRO_153 * DENDRO_188;
const double DENDRO_538 = DENDRO_138 * DENDRO_239;
const double DENDRO_539 = DENDRO_426 + DENDRO_484;
const double DENDRO_540 = 0.25 * DENDRO_436;
const double DENDRO_541 = DENDRO_212 * DENDRO_55;
const double DENDRO_542 = DENDRO_134 * DENDRO_199 + DENDRO_541;
const double DENDRO_543 = DENDRO_540 + DENDRO_542;
const double DENDRO_544 = DENDRO_223 * grad_1_gt0[pp];
const double DENDRO_545 = DENDRO_444 + DENDRO_544;
const double DENDRO_546 = DENDRO_220 * grad_0_gt3[pp];
const double DENDRO_547 = DENDRO_177 * grad_2_gt3[pp];
const double DENDRO_548 = 0.25 * DENDRO_310;
const double DENDRO_549 = -DENDRO_205 * DENDRO_243;
const double DENDRO_550 = DENDRO_111 * DENDRO_223 + DENDRO_226 * DENDRO_380;
const double DENDRO_551 = DENDRO_213 * DENDRO_215;
const double DENDRO_552 = DENDRO_229 * DENDRO_320;
const double DENDRO_553 = DENDRO_102 * DENDRO_223 + DENDRO_552;
const double DENDRO_554 = DENDRO_239 * DENDRO_297;
const double DENDRO_555 = DENDRO_180 * DENDRO_220;
const double DENDRO_556 = DENDRO_481 + DENDRO_555;
const double DENDRO_557 = DENDRO_226 * DENDRO_49 + 0.5 * DENDRO_333;
const double DENDRO_558 = DENDRO_220 * DENDRO_93;
const double DENDRO_559 = 0.25 * DENDRO_558;
const double DENDRO_560 = DENDRO_243 * DENDRO_380;
const double DENDRO_561 = DENDRO_199 * DENDRO_470 + DENDRO_468;
const double DENDRO_562 = DENDRO_177 * DENDRO_336;
const double DENDRO_563 = -DENDRO_562;
const double DENDRO_564 = DENDRO_200 * DENDRO_246;
const double DENDRO_565 = -DENDRO_239 * DENDRO_336 + DENDRO_564;
const double DENDRO_566 = DENDRO_173 * DENDRO_229 + DENDRO_479 * DENDRO_93;
const double DENDRO_567 = DENDRO_206 * DENDRO_208;
const double DENDRO_568 = DENDRO_199 * DENDRO_297 + DENDRO_567;
const double DENDRO_569 = DENDRO_199 * DENDRO_387;
const double DENDRO_570 = DENDRO_152 * DENDRO_479;
const double DENDRO_571 = DENDRO_167 * DENDRO_226 + DENDRO_223 * DENDRO_470;
const double DENDRO_572 = -DENDRO_220 * DENDRO_385;
const double DENDRO_573 = DENDRO_229 * DENDRO_52 + 0.5 * DENDRO_376;
const double DENDRO_574 = DENDRO_452 * grad_1_gt3[pp];
const double DENDRO_575 = DENDRO_409 + DENDRO_574;
const double DENDRO_576 = DENDRO_215 * DENDRO_320;
const double DENDRO_577 = 0.25 * DENDRO_93;
const double DENDRO_578 = DENDRO_173 * DENDRO_246;
const double DENDRO_579 = DENDRO_239 * DENDRO_577 + DENDRO_578;
const double DENDRO_580 = DENDRO_177 * DENDRO_93;
const double DENDRO_581 = 1.0 * DENDRO_292;
const double DENDRO_582 = -grad2_1_2_chi[pp];
const double DENDRO_583 = DENDRO_152 * DENDRO_347 + DENDRO_21 * grad_2_gt3[pp] +
                          DENDRO_345 * grad_1_gt5[pp];
const double DENDRO_584 = DENDRO_350 * grad_2_gt3[pp];
const double DENDRO_585 = DENDRO_352 * grad_1_gt5[pp];
const double DENDRO_586 = DENDRO_23 * grad_2_gt3[pp];
const double DENDRO_587 = DENDRO_152 * DENDRO_345;
const double DENDRO_588 = DENDRO_31 * gt4[pp];
const double DENDRO_589 =
    DENDRO_284 * grad2_0_2_gt4[pp] + DENDRO_286 * grad2_0_0_gt4[pp] +
    DENDRO_288 * grad2_1_1_gt4[pp] + DENDRO_290 * grad2_2_2_gt4[pp] +
    DENDRO_501 * grad_1_Gt0[pp] + DENDRO_502 * grad_1_Gt1[pp] +
    DENDRO_502 * grad_2_Gt2[pp] + DENDRO_503 * grad_1_Gt2[pp] +
    DENDRO_504 * gt1[pp] + DENDRO_505 * gt3[pp] - DENDRO_506 * grad_1_chi[pp] -
    DENDRO_81 * grad2_0_1_gt4[pp] - DENDRO_83 * grad2_1_2_gt4[pp];
const double DENDRO_590 =
    -DENDRO_193 *
        (DENDRO_494 * DENDRO_583 + DENDRO_496 * (DENDRO_198 + DENDRO_584) +
         DENDRO_500 * (DENDRO_585 + DENDRO_586 + DENDRO_587) + DENDRO_582) +
    DENDRO_357 * grad_1_gt4[pp] + DENDRO_359 * grad_0_gt4[pp] +
    DENDRO_361 * grad_2_gt4[pp] + DENDRO_363 * DENDRO_588 + DENDRO_589;
const double DENDRO_591 = DENDRO_196 * DENDRO_31 + DENDRO_197 * DENDRO_31;
const double DENDRO_592 =
    -DENDRO_195 * DENDRO_31 -
    DENDRO_45 * (-DENDRO_43 * DENDRO_588 + grad_2_chi[pp]) + DENDRO_591;
const double DENDRO_593 =
    -DENDRO_23 * DENDRO_31 * grad_2_gt3[pp] + DENDRO_236 * DENDRO_31 +
    DENDRO_237 * DENDRO_31 +
    DENDRO_45 * (-DENDRO_588 * DENDRO_60 + grad_1_chi[pp]);
const double DENDRO_594 = DENDRO_222 - DENDRO_45 * DENDRO_73 * gt4[pp];
const double DENDRO_595 = -DENDRO_31 * DENDRO_513 * DENDRO_594 -
                          DENDRO_515 * DENDRO_593 + DENDRO_516 * DENDRO_592 -
                          4 * grad2_1_2_alpha[pp];
const double DENDRO_596 = 1.0 * DENDRO_401;
const double DENDRO_597 = 0.5 * DENDRO_400;
const double DENDRO_598 = 0.25 * DENDRO_580;
const double DENDRO_599 = DENDRO_307 + DENDRO_475;
const double DENDRO_600 = DENDRO_174 * DENDRO_229;
const double DENDRO_601 = DENDRO_567 + DENDRO_569;
const double DENDRO_602 = DENDRO_239 * DENDRO_308;
const double DENDRO_603 = DENDRO_409 + DENDRO_410;
const double DENDRO_604 = DENDRO_229 * DENDRO_51;
const double DENDRO_605 = DENDRO_102 * DENDRO_220 + DENDRO_604;
const double DENDRO_606 = DENDRO_215 * DENDRO_380 + DENDRO_574;
const double DENDRO_607 = DENDRO_177 * DENDRO_297 + DENDRO_578;
const double DENDRO_608 = 1.0 * DENDRO_162;
const double DENDRO_609 =
    -DENDRO_116 * (DENDRO_206 * DENDRO_212 + DENDRO_340 + DENDRO_548) -
    DENDRO_182 * (-DENDRO_532 + DENDRO_561) -
    DENDRO_581 * (DENDRO_170 + DENDRO_437) -
    DENDRO_608 * (DENDRO_329 + DENDRO_523);
const double DENDRO_610 = DENDRO_472 + DENDRO_562;
const double DENDRO_611 = DENDRO_102 * DENDRO_224;
const double DENDRO_612 = -DENDRO_168 * DENDRO_215;
const double DENDRO_613 = DENDRO_153 * DENDRO_246 + DENDRO_175 * DENDRO_297;
const double DENDRO_614 = DENDRO_188 * DENDRO_49;
const double DENDRO_615 = 0.25 * DENDRO_146;
const double DENDRO_616 = DENDRO_138 * DENDRO_175 + DENDRO_191 * DENDRO_380;
const double DENDRO_617 = 0.5 * DENDRO_176 + DENDRO_191 * DENDRO_205;
const double DENDRO_618 = DENDRO_480 + DENDRO_555;
const double DENDRO_619 = DENDRO_188 * DENDRO_380 + DENDRO_418;
const double DENDRO_620 = DENDRO_101 * DENDRO_485 + DENDRO_541;
const double DENDRO_621 = DENDRO_174 * DENDRO_246 + 0.5 * DENDRO_382;
const double DENDRO_622 = DENDRO_153 * DENDRO_215;
const double DENDRO_623 = DENDRO_224 * DENDRO_577 + DENDRO_604;
const double DENDRO_624 = -DENDRO_101 * DENDRO_385;
const double DENDRO_625 = DENDRO_208 * DENDRO_55;
const double DENDRO_626 = DENDRO_152 * DENDRO_474 + DENDRO_191 * DENDRO_212;
const double DENDRO_627 = DENDRO_134 * DENDRO_224;
const double DENDRO_628 = DENDRO_229 * DENDRO_53;
const double DENDRO_629 = DENDRO_180 * DENDRO_224 + DENDRO_628;
const double DENDRO_630 = DENDRO_474 * DENDRO_88;
const double DENDRO_631 = DENDRO_246 * DENDRO_51 + DENDRO_474 * DENDRO_93;
const double DENDRO_632 = 1.0 * DENDRO_115;
const double DENDRO_633 = -grad2_0_1_chi[pp];
const double DENDRO_634 = DENDRO_347 * grad_1_gt0[pp];
const double DENDRO_635 = DENDRO_345 * DENDRO_88;
const double DENDRO_636 = DENDRO_350 * grad_0_gt3[pp];
const double DENDRO_637 = DENDRO_23 * grad_0_gt3[pp];
const double DENDRO_638 =
    DENDRO_345 * grad_1_gt0[pp] + DENDRO_352 * DENDRO_88 + DENDRO_637;
const double DENDRO_639 = 2.0 * gt1[pp];
const double DENDRO_640 = DENDRO_639 * grad_0_Gt0[pp];
const double DENDRO_641 = 2.0 * grad_0_Gt1[pp] * gt3[pp];
const double DENDRO_642 = DENDRO_502 * grad_0_Gt2[pp];
const double DENDRO_643 = 2.0 * grad_1_Gt0[pp] * gt0[pp];
const double DENDRO_644 = DENDRO_639 * grad_1_Gt1[pp];
const double DENDRO_645 = DENDRO_501 * grad_1_Gt2[pp];
const double DENDRO_646 = -DENDRO_253 * DENDRO_77;
const double DENDRO_647 = DENDRO_284 * grad2_0_2_gt1[pp];
const double DENDRO_648 = DENDRO_286 * grad2_0_0_gt1[pp];
const double DENDRO_649 = DENDRO_288 * grad2_1_1_gt1[pp];
const double DENDRO_650 = DENDRO_290 * grad2_2_2_gt1[pp];
const double DENDRO_651 = -DENDRO_81 * grad2_0_1_gt1[pp];
const double DENDRO_652 = -DENDRO_83 * grad2_1_2_gt1[pp];
const double DENDRO_653 = DENDRO_31 * gt1[pp];
const double DENDRO_654 =
    -DENDRO_193 * (DENDRO_494 * (DENDRO_634 + DENDRO_635 + DENDRO_86) +
                   DENDRO_496 * (DENDRO_100 + DENDRO_636) +
                   DENDRO_500 * DENDRO_638 + DENDRO_633) +
    DENDRO_357 * grad_1_gt1[pp] + DENDRO_359 * grad_0_gt1[pp] +
    DENDRO_361 * grad_2_gt1[pp] + DENDRO_363 * DENDRO_653 + DENDRO_640 +
    DENDRO_641 + DENDRO_642 + DENDRO_643 + DENDRO_644 + DENDRO_645 +
    DENDRO_646 + DENDRO_647 + DENDRO_648 + DENDRO_649 + DENDRO_650 +
    DENDRO_651 + DENDRO_652;
const double DENDRO_655 =
    -DENDRO_21 * DENDRO_31 * grad_0_gt3[pp] + DENDRO_31 * DENDRO_87 +
    DENDRO_31 * DENDRO_89 +
    DENDRO_45 * (-DENDRO_653 * DENDRO_73 + grad_1_chi[pp]);
const double DENDRO_656 =
    -DENDRO_21 * DENDRO_31 * grad_1_gt0[pp] -
    DENDRO_23 * DENDRO_31 * DENDRO_88 + DENDRO_31 * DENDRO_97 +
    DENDRO_45 * (-DENDRO_43 * DENDRO_653 + grad_0_chi[pp]);
const double DENDRO_657 = DENDRO_45 * gt1[pp];
const double DENDRO_658 =
    DENDRO_31 * DENDRO_515 *
        (-DENDRO_127 - DENDRO_128 + DENDRO_60 * DENDRO_657 + DENDRO_637) -
    DENDRO_513 * DENDRO_655 - DENDRO_516 * DENDRO_656 - 4 * grad2_0_1_alpha[pp];
const double DENDRO_659 = -0.25 * DENDRO_177 * grad_1_gt5[pp] + DENDRO_610;
const double DENDRO_660 = 0.5 * DENDRO_396;
const double DENDRO_661 = DENDRO_229 * DENDRO_48;
const double DENDRO_662 = DENDRO_179 + DENDRO_617;
const double DENDRO_663 = -0.5 * DENDRO_177 * DENDRO_213 + DENDRO_621;
const double DENDRO_664 = DENDRO_410 + DENDRO_574;
const double DENDRO_665 = DENDRO_101 * DENDRO_102 + DENDRO_625;
const double DENDRO_666 = -DENDRO_116 * (DENDRO_118 * DENDRO_212 + DENDRO_467) +
                          DENDRO_137 * (DENDRO_447 + DENDRO_664) +
                          DENDRO_158 * (DENDRO_624 + DENDRO_665) -
                          DENDRO_182 * (DENDRO_151 + DENDRO_542) -
                          DENDRO_182 * (DENDRO_151 + DENDRO_620) -
                          DENDRO_300 * (1.0 * DENDRO_147 + DENDRO_615);
const double DENDRO_667 =
    -DENDRO_21 *
        (DENDRO_658 +
         alpha[pp] *
             (-DENDRO_106 * (DENDRO_393 + DENDRO_613) -
              DENDRO_106 * (DENDRO_612 + DENDRO_660) -
              DENDRO_106 * (DENDRO_394 + DENDRO_611 + DENDRO_661) +
              DENDRO_116 * DENDRO_659 +
              DENDRO_136 * (DENDRO_377 + DENDRO_544 + DENDRO_546) +
              DENDRO_136 * (DENDRO_487 + DENDRO_488 + DENDRO_547) -
              DENDRO_137 * DENDRO_663 + DENDRO_137 * (DENDRO_572 + DENDRO_623) +
              DENDRO_137 * (DENDRO_622 + DENDRO_664) +
              DENDRO_143 * (DENDRO_188 * grad_0_gt3[pp] + DENDRO_415) +
              DENDRO_158 * (DENDRO_626 + DENDRO_630) +
              DENDRO_158 * (-DENDRO_168 * DENDRO_188 + DENDRO_629) +
              DENDRO_158 * (DENDRO_191 * DENDRO_213 + DENDRO_631) +
              DENDRO_158 * (DENDRO_215 * DENDRO_48 + DENDRO_665) +
              DENDRO_182 * DENDRO_662 - DENDRO_182 * (DENDRO_417 + DENDRO_619) -
              DENDRO_182 * (DENDRO_537 + DENDRO_618) -
              DENDRO_182 * (DENDRO_484 + DENDRO_538 + DENDRO_598) -
              DENDRO_300 * (0.5 * DENDRO_427 + DENDRO_616) -
              DENDRO_300 * (DENDRO_188 * DENDRO_48 + DENDRO_424 + DENDRO_614) -
              DENDRO_632 * (DENDRO_334 + DENDRO_520 + DENDRO_558) + DENDRO_654 +
              DENDRO_666)) -
    DENDRO_21 *
        (DENDRO_658 +
         alpha[pp] *
             (-DENDRO_106 * (1.0 * DENDRO_399 + DENDRO_611) -
              DENDRO_106 * (0.5 * DENDRO_402 + DENDRO_613) -
              DENDRO_106 * (DENDRO_167 * DENDRO_215 + DENDRO_397 + DENDRO_612) -
              DENDRO_116 * (DENDRO_223 * DENDRO_51 + DENDRO_559) -
              DENDRO_116 * (0.25 * DENDRO_239 * grad_0_gt5[pp] - DENDRO_610) +
              DENDRO_137 * (DENDRO_379 + DENDRO_605) +
              DENDRO_137 * (DENDRO_379 + DENDRO_623) +
              DENDRO_137 * (-DENDRO_384 - DENDRO_621) +
              DENDRO_137 * (DENDRO_447 + DENDRO_606) +
              DENDRO_137 * (DENDRO_526 + DENDRO_607) +
              DENDRO_137 * (DENDRO_603 + DENDRO_622) +
              DENDRO_143 * (DENDRO_215 * grad_1_gt0[pp] + DENDRO_411) +
              DENDRO_158 * (DENDRO_627 + DENDRO_629) +
              DENDRO_158 * (DENDRO_630 + DENDRO_631) +
              DENDRO_158 * (DENDRO_246 * DENDRO_52 + DENDRO_626) +
              DENDRO_158 * (DENDRO_167 * DENDRO_188 + DENDRO_627 + DENDRO_628) +
              DENDRO_158 * (DENDRO_215 * DENDRO_49 + DENDRO_624 + DENDRO_625) -
              DENDRO_182 * (DENDRO_418 + DENDRO_618) -
              DENDRO_182 * (DENDRO_477 + DENDRO_620) -
              DENDRO_182 * (DENDRO_480 + DENDRO_619) -
              DENDRO_182 * (0.5 * DENDRO_239 * DENDRO_52 - DENDRO_617) -
              DENDRO_300 * (0.5 * DENDRO_419 + DENDRO_614) -
              DENDRO_300 * (DENDRO_428 + DENDRO_616) -
              DENDRO_300 * (DENDRO_167 * DENDRO_55 + DENDRO_185 + DENDRO_615) -
              DENDRO_608 * (DENDRO_170 + DENDRO_434 + DENDRO_435) -
              DENDRO_608 * (DENDRO_438 + DENDRO_439 + DENDRO_580) -
              DENDRO_632 * (DENDRO_329 + DENDRO_466 + DENDRO_522) +
              DENDRO_654)) -
    DENDRO_23 *
        (DENDRO_595 +
         alpha[pp] *
             (-DENDRO_106 * (DENDRO_553 + DENDRO_597) -
              DENDRO_106 * (DENDRO_554 + DENDRO_596) -
              DENDRO_106 * (DENDRO_212 * DENDRO_215 + DENDRO_391 + DENDRO_551) -
              DENDRO_116 * (0.5 * DENDRO_317 + DENDRO_549) -
              DENDRO_116 * (-DENDRO_339 + DENDRO_550) +
              DENDRO_137 * (DENDRO_565 + DENDRO_602) +
              DENDRO_137 * (DENDRO_566 + DENDRO_570) +
              DENDRO_137 * (DENDRO_571 - DENDRO_600) +
              DENDRO_137 * (-DENDRO_205 * DENDRO_215 + DENDRO_601) +
              DENDRO_137 * (DENDRO_212 * DENDRO_243 + DENDRO_564 + DENDRO_602) +
              DENDRO_158 * (DENDRO_378 + DENDRO_573) +
              DENDRO_158 * (DENDRO_381 + DENDRO_579) +
              DENDRO_158 * (DENDRO_381 + DENDRO_607) +
              DENDRO_158 * (DENDRO_445 + DENDRO_605) +
              DENDRO_158 * (DENDRO_446 + DENDRO_606) +
              DENDRO_158 * (DENDRO_576 + DENDRO_603) -
              DENDRO_182 * (DENDRO_535 + DENDRO_557) -
              DENDRO_182 * (DENDRO_560 + DENDRO_599) -
              DENDRO_182 * (DENDRO_563 + DENDRO_599) -
              DENDRO_300 * (DENDRO_418 + DENDRO_556) -
              DENDRO_300 * (DENDRO_173 * DENDRO_175 + DENDRO_598) +
              DENDRO_313 * (DENDRO_215 * grad_1_gt5[pp] + DENDRO_412) +
              DENDRO_590 - DENDRO_608 * (DENDRO_521 + DENDRO_558) +
              DENDRO_609)) -
    DENDRO_23 *
        (DENDRO_595 +
         alpha[pp] *
             (-DENDRO_106 * (-DENDRO_386 + DENDRO_553) -
              DENDRO_106 * (0.5 * DENDRO_390 + DENDRO_551) -
              DENDRO_106 * (DENDRO_202 * DENDRO_246 + DENDRO_388 + DENDRO_554) -
              DENDRO_116 * (1.0 * DENDRO_311 + DENDRO_548) -
              DENDRO_116 * (0.5 * DENDRO_332 + DENDRO_550) -
              DENDRO_116 * (DENDRO_202 * DENDRO_243 + DENDRO_323 + DENDRO_549) +
              DENDRO_137 * (DENDRO_568 + DENDRO_569) +
              DENDRO_137 * (DENDRO_570 + DENDRO_571) +
              DENDRO_137 * (-DENDRO_168 * DENDRO_226 + DENDRO_566) +
              DENDRO_137 * (DENDRO_202 * DENDRO_215 + DENDRO_568) +
              DENDRO_137 * (DENDRO_213 * DENDRO_243 + DENDRO_565) +
              DENDRO_142 * (DENDRO_489 + DENDRO_547) +
              DENDRO_142 * (DENDRO_545 + DENDRO_546) +
              DENDRO_158 * (DENDRO_446 + DENDRO_575) +
              DENDRO_158 * (DENDRO_572 + DENDRO_573) +
              DENDRO_158 * (DENDRO_575 + DENDRO_576) +
              DENDRO_158 * (DENDRO_177 * DENDRO_387 + DENDRO_579) -
              DENDRO_182 * (DENDRO_319 + DENDRO_469) -
              DENDRO_182 * (DENDRO_319 + DENDRO_561) -
              DENDRO_182 * (DENDRO_335 + DENDRO_557) -
              DENDRO_182 * (DENDRO_463 + DENDRO_559) -
              DENDRO_182 * (DENDRO_465 + DENDRO_560) -
              DENDRO_182 * (DENDRO_476 + DENDRO_563) -
              DENDRO_300 * (DENDRO_417 + DENDRO_556) -
              DENDRO_300 * (DENDRO_118 * DENDRO_167 + DENDRO_540) +
              DENDRO_313 * (DENDRO_243 * grad_2_gt3[pp] + DENDRO_294) -
              DENDRO_581 * (DENDRO_440 + DENDRO_580) + DENDRO_590)) +
    DENDRO_26 *
        (-DENDRO_268 * DENDRO_368 +
         DENDRO_272 * (DENDRO_229 + DENDRO_364 * DENDRO_74) +
         DENDRO_365 * (DENDRO_246 + DENDRO_275 * DENDRO_364) +
         alpha[pp] *
             (DENDRO_116 * DENDRO_239 * DENDRO_403 +
              DENDRO_137 * (DENDRO_388 - DENDRO_389) +
              DENDRO_137 * (DENDRO_391 + 1.0 * DENDRO_392) +
              DENDRO_137 * (DENDRO_229 * DENDRO_93 - DENDRO_386) +
              DENDRO_143 * (DENDRO_396 + DENDRO_398) +
              DENDRO_143 * (DENDRO_175 * grad_2_gt3[pp] + DENDRO_402) +
              DENDRO_143 * (DENDRO_224 * grad_0_gt3[pp] + DENDRO_399) +
              DENDRO_158 * (DENDRO_394 + DENDRO_395) +
              DENDRO_158 * (DENDRO_397 + 1.0 * DENDRO_398) +
              DENDRO_158 * (DENDRO_246 * DENDRO_93 + DENDRO_393) -
              DENDRO_175 * DENDRO_407 - DENDRO_182 * (DENDRO_376 + DENDRO_378) -
              DENDRO_182 * (-1.0 * DENDRO_382 - DENDRO_384) -
              DENDRO_182 * (DENDRO_224 * DENDRO_380 + DENDRO_379) -
              DENDRO_182 * (DENDRO_239 * DENDRO_380 + DENDRO_381) -
              DENDRO_193 *
                  (DENDRO_189 * DENDRO_371 +
                   DENDRO_190 * (DENDRO_372 + DENDRO_373 + DENDRO_374) +
                   DENDRO_192 * DENDRO_375 + DENDRO_369) -
              DENDRO_216 * DENDRO_296 * grad_1_gt3[pp] -
              DENDRO_223 * DENDRO_405 - DENDRO_224 * DENDRO_406 +
              DENDRO_313 * (DENDRO_390 + DENDRO_392) +
              DENDRO_313 * (DENDRO_223 * grad_0_gt3[pp] + DENDRO_400) +
              DENDRO_313 * (DENDRO_239 * grad_2_gt3[pp] + DENDRO_401) +
              DENDRO_357 * grad_1_gt3[pp] + DENDRO_359 * grad_0_gt3[pp] +
              DENDRO_361 * grad_2_gt3[pp] + DENDRO_363 * DENDRO_367 +
              DENDRO_414) -
         4 * grad2_1_1_alpha[pp]) +
    DENDRO_33 *
        (DENDRO_519 +
         alpha[pp] *
             (-DENDRO_106 * (DENDRO_409 + DENDRO_448) -
              DENDRO_106 * (DENDRO_223 * DENDRO_48 + DENDRO_445) -
              DENDRO_116 * (1.0 * DENDRO_314 + DENDRO_441) -
              DENDRO_116 * (0.5 * DENDRO_326 + DENDRO_443) -
              DENDRO_116 * (DENDRO_173 * DENDRO_243 + DENDRO_322 + DENDRO_442) +
              DENDRO_137 * (DENDRO_324 + DENDRO_463) +
              DENDRO_137 * (DENDRO_324 + DENDRO_471) +
              DENDRO_137 * (DENDRO_464 + DENDRO_465) +
              DENDRO_137 * (DENDRO_467 + DENDRO_469) +
              DENDRO_137 * (DENDRO_473 + DENDRO_476) +
              DENDRO_137 * (0.5 * DENDRO_328 + DENDRO_330 - DENDRO_461) +
              DENDRO_142 * (DENDRO_434 + DENDRO_437) +
              DENDRO_142 * (DENDRO_438 + DENDRO_440) +
              DENDRO_158 * (DENDRO_417 + DENDRO_482) +
              DENDRO_158 * (DENDRO_477 + DENDRO_478) +
              DENDRO_158 * (DENDRO_480 + DENDRO_483) +
              DENDRO_158 * (DENDRO_184 * DENDRO_239 + DENDRO_484 + DENDRO_486) -
              DENDRO_163 * (DENDRO_243 * grad_2_gt0[pp] + DENDRO_291) -
              DENDRO_182 * (DENDRO_456 + DENDRO_458) -
              DENDRO_182 * (DENDRO_459 + DENDRO_460) -
              DENDRO_182 * (DENDRO_206 * DENDRO_49 + DENDRO_453) -
              DENDRO_182 * (DENDRO_173 * DENDRO_188 + DENDRO_456 + DENDRO_457) -
              DENDRO_182 * (DENDRO_243 * DENDRO_52 + DENDRO_454 + DENDRO_455) -
              DENDRO_300 * (DENDRO_181 + DENDRO_450) -
              DENDRO_300 * (0.5 * DENDRO_421 + DENDRO_449) -
              DENDRO_300 * (DENDRO_173 * DENDRO_191 + DENDRO_429 + DENDRO_451) -
              DENDRO_490 * (DENDRO_487 + DENDRO_489) + DENDRO_511)) +
    DENDRO_33 *
        (DENDRO_519 +
         alpha[pp] *
             (-DENDRO_106 * (DENDRO_175 * DENDRO_202 + DENDRO_526) -
              DENDRO_116 * (0.5 * DENDRO_315 + DENDRO_442) -
              DENDRO_116 * (DENDRO_226 * DENDRO_51 + DENDRO_342 + DENDRO_441) +
              DENDRO_136 * DENDRO_524 + DENDRO_136 * (DENDRO_334 + DENDRO_521) -
              DENDRO_137 * DENDRO_533 + DENDRO_137 * (DENDRO_464 + DENDRO_534) +
              DENDRO_137 * (DENDRO_471 + DENDRO_535) +
              DENDRO_137 * (DENDRO_473 + DENDRO_534) + DENDRO_158 * DENDRO_536 +
              DENDRO_158 * DENDRO_543 + DENDRO_158 * (DENDRO_418 + DENDRO_483) +
              DENDRO_158 * (DENDRO_482 + DENDRO_537) +
              DENDRO_158 * (DENDRO_486 + DENDRO_539) +
              DENDRO_158 * (DENDRO_538 + DENDRO_539) -
              DENDRO_163 * (DENDRO_188 * grad_0_gt5[pp] + DENDRO_416) -
              DENDRO_182 * (DENDRO_455 + DENDRO_530) -
              DENDRO_182 * (-DENDRO_174 * DENDRO_188 + DENDRO_458) -
              DENDRO_182 * (DENDRO_243 * DENDRO_51 + DENDRO_530) -
              DENDRO_300 * (1.0 * DENDRO_423 + DENDRO_451) -
              DENDRO_300 * (DENDRO_188 * DENDRO_51 + DENDRO_425 + DENDRO_449) -
              DENDRO_490 * (DENDRO_377 + DENDRO_545) + DENDRO_511 - DENDRO_525 -
              DENDRO_527 - DENDRO_528 - DENDRO_529 - DENDRO_531)) +
    DENDRO_35 *
        (DENDRO_269 * DENDRO_56 - DENDRO_271 * DENDRO_75 -
         DENDRO_365 * DENDRO_62 +
         alpha[pp] *
             (-DENDRO_106 * DENDRO_126 * DENDRO_175 - DENDRO_107 +
              DENDRO_114 * DENDRO_116 * DENDRO_177 - DENDRO_123 -
              DENDRO_131 * DENDRO_133 * DENDRO_248 - DENDRO_132 +
              DENDRO_137 * DENDRO_154 +
              DENDRO_137 * (-1.0 * DENDRO_169 - DENDRO_172) +
              DENDRO_137 * (-1.0 * DENDRO_176 - DENDRO_179) +
              DENDRO_137 * (DENDRO_153 * DENDRO_177 + DENDRO_426) +
              DENDRO_137 * (DENDRO_220 * DENDRO_48 + DENDRO_418) +
              DENDRO_137 * (DENDRO_224 * DENDRO_51 + DENDRO_417) +
              DENDRO_143 * DENDRO_148 + DENDRO_143 * (DENDRO_419 + DENDRO_420) +
              DENDRO_143 * (DENDRO_175 * grad_2_gt0[pp] + DENDRO_427) +
              DENDRO_158 * DENDRO_186 +
              DENDRO_158 * (1.0 * DENDRO_420 + DENDRO_424) +
              DENDRO_158 * (DENDRO_152 * DENDRO_191 + DENDRO_428) -
              DENDRO_163 * (DENDRO_421 + DENDRO_422) -
              DENDRO_163 * (DENDRO_177 * grad_2_gt0[pp] + DENDRO_423) -
              DENDRO_164 - DENDRO_182 * (1.0 * DENDRO_422 + DENDRO_425) -
              DENDRO_182 * (-2 * DENDRO_174 * DENDRO_191 + DENDRO_429) -
              DENDRO_183 -
              DENDRO_193 * (DENDRO_187 +
                            DENDRO_189 * (DENDRO_430 + DENDRO_431 + DENDRO_67) +
                            DENDRO_190 * DENDRO_432 + DENDRO_192 * DENDRO_433) +
              DENDRO_219 * DENDRO_356 -
              6.0 * DENDRO_231 * DENDRO_85 * grad_0_gt0[pp] +
              DENDRO_234 * DENDRO_358 + DENDRO_251 * DENDRO_360 +
              DENDRO_264 * DENDRO_362 + DENDRO_278 * grad_0_Gt2[pp] +
              DENDRO_284 * grad2_0_2_gt0[pp] + DENDRO_286 * grad2_0_0_gt0[pp] +
              DENDRO_288 * grad2_1_1_gt0[pp] + DENDRO_290 * grad2_2_2_gt0[pp] -
              DENDRO_295 * DENDRO_415 + DENDRO_408 * grad_0_Gt1[pp] -
              DENDRO_413 * DENDRO_416 - DENDRO_79 - DENDRO_82 - DENDRO_84 +
              4 * grad_0_Gt0[pp] * gt0[pp]) -
         4 * grad2_0_0_alpha[pp]) +
    DENDRO_39 *
        (DENDRO_267 * DENDRO_269 - DENDRO_270 * DENDRO_272 -
         DENDRO_276 * DENDRO_277 +
         alpha[pp] *
             (-DENDRO_106 * DENDRO_223 * DENDRO_304 -
              DENDRO_131 * DENDRO_227 * DENDRO_305 + DENDRO_137 * DENDRO_341 +
              DENDRO_137 * (1.0 * DENDRO_318 + DENDRO_323) +
              DENDRO_137 * (DENDRO_226 * DENDRO_88 - DENDRO_339) +
              DENDRO_158 * DENDRO_321 + DENDRO_158 * DENDRO_331 +
              DENDRO_158 * (DENDRO_333 + DENDRO_335) +
              DENDRO_158 * (DENDRO_173 * DENDRO_239 + DENDRO_309) +
              DENDRO_158 * (DENDRO_177 * DENDRO_202 + DENDRO_307) +
              DENDRO_158 * (DENDRO_220 * DENDRO_320 + DENDRO_324) -
              DENDRO_163 * (DENDRO_315 + DENDRO_316) -
              DENDRO_163 * (DENDRO_220 * grad_0_gt5[pp] + DENDRO_314) -
              DENDRO_182 * (1.0 * DENDRO_316 + DENDRO_322) -
              DENDRO_182 * (DENDRO_206 * DENDRO_88 - DENDRO_337) -
              DENDRO_182 * (DENDRO_226 * DENDRO_343 + DENDRO_342) -
              DENDRO_193 *
                  (DENDRO_189 * DENDRO_349 + DENDRO_190 * DENDRO_351 +
                   DENDRO_192 * (DENDRO_353 + DENDRO_354 + DENDRO_355) +
                   DENDRO_344) -
              DENDRO_220 * DENDRO_299 * DENDRO_300 -
              DENDRO_244 * DENDRO_296 * grad_2_gt5[pp] +
              DENDRO_274 * DENDRO_363 + DENDRO_278 * grad_2_Gt0[pp] +
              DENDRO_279 * grad_2_Gt1[pp] - DENDRO_280 - DENDRO_281 -
              DENDRO_282 + DENDRO_284 * grad2_0_2_gt5[pp] +
              DENDRO_286 * grad2_0_0_gt5[pp] + DENDRO_288 * grad2_1_1_gt5[pp] +
              DENDRO_290 * grad2_2_2_gt5[pp] - DENDRO_291 * DENDRO_293 -
              DENDRO_294 * DENDRO_295 - DENDRO_298 - DENDRO_303 - DENDRO_306 +
              DENDRO_312 * DENDRO_313 + DENDRO_313 * (DENDRO_317 + DENDRO_318) +
              DENDRO_313 * (DENDRO_223 * grad_0_gt5[pp] + DENDRO_332) -
              DENDRO_327 + DENDRO_357 * grad_1_gt5[pp] +
              DENDRO_359 * grad_0_gt5[pp] + DENDRO_361 * grad_2_gt5[pp] +
              4 * grad_2_Gt2[pp] * gt5[pp]) -
         4 * grad2_2_2_alpha[pp]);
const double DENDRO_668 = DENDRO_31 * DENDRO_667;
const double DENDRO_669 = (1.0 / 12.0) * chi[pp];
const double DENDRO_670 = (1.0 / 3.0) * At1[pp];
const double DENDRO_671 = At1[pp] * DENDRO_21;
const double DENDRO_672 = At4[pp] * DENDRO_23;
const double DENDRO_673 = At3[pp] * DENDRO_26;
const double DENDRO_674 = DENDRO_671 + DENDRO_672 - DENDRO_673;
const double DENDRO_675 = At4[pp] * DENDRO_33;
const double DENDRO_676 =
    -At1[pp] * DENDRO_35 + At3[pp] * DENDRO_21 - DENDRO_675;
const double DENDRO_677 = At4[pp] * DENDRO_39;
const double DENDRO_678 =
    -At1[pp] * DENDRO_33 + At3[pp] * DENDRO_23 - DENDRO_677;
const double DENDRO_679 = 6.0 * grad_2_alpha[pp];
const double DENDRO_680 = 6.0 * grad_0_alpha[pp];
const double DENDRO_681 = 6.0 * grad_1_alpha[pp];
const double DENDRO_682 = DENDRO_93 * DENDRO_95;
const double DENDRO_683 = DENDRO_222 * grad_2_gt0[pp];
const double DENDRO_684 = DENDRO_152 * DENDRO_95;
const double DENDRO_685 = DENDRO_683 + DENDRO_684;
const double DENDRO_686 = DENDRO_385 * DENDRO_90;
const double DENDRO_687 = DENDRO_129 * DENDRO_387;
const double DENDRO_688 = DENDRO_180 * DENDRO_95;
const double DENDRO_689 =
    DENDRO_153 * DENDRO_96 + 0.25 * DENDRO_222 * grad_0_gt0[pp];
const double DENDRO_690 = DENDRO_135 + DENDRO_139;
const double DENDRO_691 = DENDRO_110 * DENDRO_577;
const double DENDRO_692 = DENDRO_202 * DENDRO_61;
const double DENDRO_693 = DENDRO_138 * DENDRO_238 + DENDRO_692;
const double DENDRO_694 = DENDRO_152 * DENDRO_90;
const double DENDRO_695 = DENDRO_222 * grad_1_gt0[pp] + DENDRO_694;
const double DENDRO_696 = 1.0 * DENDRO_136;
const double DENDRO_697 = DENDRO_152 * DENDRO_238;
const double DENDRO_698 = 2.0 * DENDRO_233;
const double DENDRO_699 = 2.0 * DENDRO_218;
const double DENDRO_700 = 2.0 * DENDRO_250;
const double DENDRO_701 = DENDRO_263 * DENDRO_45;
const double DENDRO_702 = (1.0 / 3.0) * At2[pp];
const double DENDRO_703 = At5[pp] * DENDRO_23;
const double DENDRO_704 =
    At2[pp] * DENDRO_21 - At4[pp] * DENDRO_26 + DENDRO_703;
const double DENDRO_705 = At5[pp] * DENDRO_33;
const double DENDRO_706 =
    -At2[pp] * DENDRO_35 + At4[pp] * DENDRO_21 - DENDRO_705;
const double DENDRO_707 = At4[pp] * DENDRO_23 - At5[pp] * DENDRO_39 - DENDRO_34;
const double DENDRO_708 = DENDRO_110 * grad_2_gt5[pp];
const double DENDRO_709 = DENDRO_90 * grad_0_gt5[pp];
const double DENDRO_710 = DENDRO_110 * DENDRO_308;
const double DENDRO_711 = 0.25 * DENDRO_129 * grad_2_gt5[pp];
const double DENDRO_712 = DENDRO_710 + DENDRO_711;
const double DENDRO_713 = DENDRO_88 * DENDRO_90;
const double DENDRO_714 = DENDRO_110 * DENDRO_111 + DENDRO_200 * DENDRO_61;
const double DENDRO_715 = -0.5 * DENDRO_174 * DENDRO_95;
const double DENDRO_716 = -0.5 * DENDRO_174 * DENDRO_90;
const double DENDRO_717 = DENDRO_184 * DENDRO_90;
const double DENDRO_718 = 2 * At4[pp];
const double DENDRO_719 = At3[pp] * DENDRO_37;
const double DENDRO_720 = DENDRO_31 * DENDRO_718;
const double DENDRO_721 = 0.5 * DENDRO_364;
const double DENDRO_722 = DENDRO_31 * DENDRO_76;
const double DENDRO_723 = DENDRO_214 * grad_2_gt3[pp];
const double DENDRO_724 = DENDRO_129 * DENDRO_308;
const double DENDRO_725 = 1.0 * DENDRO_222;
const double DENDRO_726 = 0.25 * DENDRO_694;
const double DENDRO_727 = DENDRO_214 * grad_0_gt3[pp];
const double DENDRO_728 = (1.0 / 3.0) * At4[pp];
const double DENDRO_729 = DENDRO_238 * grad_2_gt5[pp];
const double DENDRO_730 = DENDRO_111 * DENDRO_238;
const double DENDRO_731 = DENDRO_711 + DENDRO_730;
const double DENDRO_732 = DENDRO_238 * DENDRO_308;
const double DENDRO_733 = -0.5 * DENDRO_246 * grad_0_gt5[pp] + DENDRO_724;
const double DENDRO_734 = DENDRO_242 * grad_0_gt5[pp];
const double DENDRO_735 = DENDRO_242 * grad_1_gt5[pp];
const double DENDRO_736 = DENDRO_21 * grad_1_chi[pp] +
                          DENDRO_345 * grad_2_chi[pp] +
                          DENDRO_347 * grad_0_chi[pp];
const double DENDRO_737 = 0.5 * DENDRO_266;
const double DENDRO_738 = DENDRO_31 * grad_0_alpha[pp];
const double DENDRO_739 = DENDRO_350 * grad_1_chi[pp] + DENDRO_42;
const double DENDRO_740 = DENDRO_31 * grad_1_alpha[pp];
const double DENDRO_741 = DENDRO_23 * grad_1_chi[pp] +
                          DENDRO_345 * grad_0_chi[pp] +
                          DENDRO_352 * grad_2_chi[pp];
const double DENDRO_742 = 0.5 * DENDRO_741;
const double DENDRO_743 = DENDRO_31 * grad_2_alpha[pp];
const double DENDRO_744 = 0.5 * DENDRO_739;
const double DENDRO_745 = 0.5 * grad_1_alpha[pp];
const double DENDRO_746 = 0.5 * grad_2_alpha[pp];
const double DENDRO_747 = 0.5 * grad_0_alpha[pp];
const double DENDRO_748 = DENDRO_256 * DENDRO_31;
const double DENDRO_749 = (DENDRO_21 * DENDRO_21);
const double DENDRO_750 = (DENDRO_33 * DENDRO_33);
const double DENDRO_751 = 2 * DENDRO_35;
const double DENDRO_752 = At0[pp] * (DENDRO_35 * DENDRO_35) +
                          At3[pp] * DENDRO_749 + At5[pp] * DENDRO_750 -
                          DENDRO_254 * DENDRO_675 + DENDRO_34 * DENDRO_751 -
                          DENDRO_671 * DENDRO_751;
const double DENDRO_753 = 3 * DENDRO_85;
const double DENDRO_754 = (DENDRO_23 * DENDRO_23);
const double DENDRO_755 = 2 * DENDRO_26;
const double DENDRO_756 = At0[pp] * DENDRO_749 +
                          At3[pp] * (DENDRO_26 * DENDRO_26) +
                          At5[pp] * DENDRO_754 + DENDRO_24 * DENDRO_254 -
                          DENDRO_671 * DENDRO_755 - DENDRO_672 * DENDRO_755;
const double DENDRO_757 = At1[pp] * DENDRO_23;
const double DENDRO_758 = 2 * DENDRO_39;
const double DENDRO_759 = At0[pp] * DENDRO_750 + At3[pp] * DENDRO_754 +
                          At5[pp] * (DENDRO_39 * DENDRO_39) -
                          DENDRO_256 * DENDRO_757 + DENDRO_34 * DENDRO_758 -
                          DENDRO_672 * DENDRO_758;
const double DENDRO_760 = At3[pp] * DENDRO_21;
const double DENDRO_761 =
    At2[pp] * DENDRO_750 - DENDRO_21 * DENDRO_677 - DENDRO_23 * DENDRO_675 +
    DENDRO_23 * DENDRO_760 + DENDRO_33 * DENDRO_36 - DENDRO_33 * DENDRO_671 +
    DENDRO_35 * DENDRO_40 - DENDRO_35 * DENDRO_757 + DENDRO_39 * DENDRO_705;
const double DENDRO_762 = 6 * DENDRO_85;
const double DENDRO_763 =
    -At1[pp] * DENDRO_26 * DENDRO_35 - At1[pp] * DENDRO_749 -
    At4[pp] * DENDRO_21 * DENDRO_23 - At4[pp] * DENDRO_26 * DENDRO_33 +
    DENDRO_21 * DENDRO_34 + DENDRO_21 * DENDRO_36 + DENDRO_23 * DENDRO_705 +
    DENDRO_24 * DENDRO_35 + DENDRO_26 * DENDRO_760;
const double DENDRO_764 = -DENDRO_763;
const double DENDRO_765 =
    -At1[pp] * DENDRO_21 * DENDRO_23 - At1[pp] * DENDRO_26 * DENDRO_33 -
    At4[pp] * DENDRO_26 * DENDRO_39 - At4[pp] * DENDRO_754 +
    DENDRO_21 * DENDRO_40 + DENDRO_22 * DENDRO_33 + DENDRO_23 * DENDRO_34 +
    DENDRO_23 * DENDRO_673 + DENDRO_39 * DENDRO_703;
const double DENDRO_766 = -DENDRO_765;
const double DENDRO_767 = (1.0 / 3.0) * alpha[pp];
const double DENDRO_768 = (7.0 / 3.0) * DENDRO_283;
const double DENDRO_769 = (1.0 / 3.0) * DENDRO_283;
const double DENDRO_770 = (1.0 / 3.0) * DENDRO_285;
const double DENDRO_771 = 2 * DENDRO_85;
const double DENDRO_772 = DENDRO_771 * grad_0_alpha[pp];
const double DENDRO_773 = pow(DENDRO_30, -3);
const double DENDRO_774 = 4 * grad_0_K[pp];
const double DENDRO_775 = DENDRO_31 * DENDRO_767;
const double DENDRO_776 = DENDRO_771 * grad_2_alpha[pp];
const double DENDRO_777 = DENDRO_771 * grad_1_alpha[pp];
const double DENDRO_778 = 4 * grad_2_K[pp];
const double DENDRO_779 = 4 * grad_1_K[pp];
const double DENDRO_780 = 9 * DENDRO_45;
const double DENDRO_781 = DENDRO_190 * DENDRO_780;
const double DENDRO_782 =
    -2.0 / 3.0 * DENDRO_17 * DENDRO_261 * DENDRO_85 -
    2 * DENDRO_188 * DENDRO_752 * DENDRO_773 * alpha[pp] -
    7.0 / 3.0 * DENDRO_21 * DENDRO_31 * grad2_0_1_beta0[pp] -
    1.0 / 3.0 * DENDRO_21 * DENDRO_31 * grad2_1_1_beta1[pp] -
    1.0 / 3.0 * DENDRO_21 * DENDRO_31 * grad2_1_2_beta2[pp] -
    2.0 * DENDRO_220 * DENDRO_761 * DENDRO_773 * alpha[pp] -
    2.0 * DENDRO_223 * DENDRO_766 * DENDRO_773 * alpha[pp] -
    2.0 * DENDRO_224 * DENDRO_764 * DENDRO_773 * alpha[pp] -
    2 * DENDRO_226 * DENDRO_759 * DENDRO_773 * alpha[pp] -
    2 * DENDRO_229 * DENDRO_756 * DENDRO_773 * alpha[pp] -
    2 * DENDRO_23 * DENDRO_31 * grad2_1_2_beta0[pp] +
    (4.0 / 3.0) * DENDRO_285 * grad2_0_0_beta0[pp] +
    DENDRO_287 * grad2_1_1_beta0[pp] + DENDRO_289 * grad2_2_2_beta0[pp] +
    DENDRO_356 * grad_1_beta0[pp] + DENDRO_358 * grad_0_beta0[pp] +
    DENDRO_360 * grad_2_beta0[pp] + DENDRO_752 * DENDRO_772 +
    DENDRO_761 * DENDRO_776 + DENDRO_764 * DENDRO_777 +
    DENDRO_768 * grad2_0_2_beta0[pp] + DENDRO_769 * grad2_1_2_beta1[pp] +
    DENDRO_769 * grad2_2_2_beta2[pp] + DENDRO_770 * grad2_0_1_beta1[pp] +
    DENDRO_770 * grad2_0_2_beta2[pp] +
    DENDRO_775 * (DENDRO_21 * DENDRO_779 + DENDRO_764 * DENDRO_781) +
    DENDRO_775 * (9 * DENDRO_31 * DENDRO_45 * DENDRO_752 * grad_0_chi[pp] -
                  DENDRO_35 * DENDRO_774) +
    DENDRO_775 * (9 * DENDRO_31 * DENDRO_45 * DENDRO_761 * grad_2_chi[pp] -
                  DENDRO_33 * DENDRO_778) -
    beta0[pp] * grad_0_Gt0[pp] - beta1[pp] * grad_1_Gt0[pp] -
    beta2[pp] * grad_2_Gt0[pp];
const double DENDRO_783 = DENDRO_285 * grad2_0_0_beta1[pp];
const double DENDRO_784 = DENDRO_289 * grad2_2_2_beta1[pp];
const double DENDRO_785 = DENDRO_748 * grad2_0_2_beta1[pp];
const double DENDRO_786 = DENDRO_756 * DENDRO_777;
const double DENDRO_787 = (4.0 / 3.0) * DENDRO_287 * grad2_1_1_beta1[pp];
const double DENDRO_788 = DENDRO_21 * DENDRO_774;
const double DENDRO_789 = DENDRO_189 * DENDRO_780;
const double DENDRO_790 = DENDRO_23 * DENDRO_778;
const double DENDRO_791 = DENDRO_192 * DENDRO_780;
const double DENDRO_792 = (1.0 / 3.0) * DENDRO_287;
const double DENDRO_793 = DENDRO_792 * grad2_0_1_beta0[pp];
const double DENDRO_794 = DENDRO_792 * grad2_1_2_beta2[pp];
const double DENDRO_795 = DENDRO_26 * DENDRO_779 - DENDRO_756 * DENDRO_781;
const double DENDRO_796 = DENDRO_21 * DENDRO_31;
const double DENDRO_797 = (1.0 / 3.0) * DENDRO_796;
const double DENDRO_798 = DENDRO_23 * DENDRO_31;
const double DENDRO_799 = (1.0 / 3.0) * grad2_0_2_beta0[pp];
const double DENDRO_800 = DENDRO_0 * DENDRO_773;
const double DENDRO_801 = 2.0 * DENDRO_773 * alpha[pp];
const double DENDRO_802 = beta0[pp] * grad_0_Gt1[pp] +
                          beta1[pp] * grad_1_Gt1[pp] +
                          beta2[pp] * grad_2_Gt1[pp];
const double DENDRO_803 =
    -2.0 / 3.0 * DENDRO_17 * DENDRO_262 * DENDRO_85 -
    2.0 * DENDRO_175 * DENDRO_764 * DENDRO_773 * alpha[pp] -
    2.0 * DENDRO_177 * DENDRO_761 * DENDRO_773 * alpha[pp] -
    2 * DENDRO_191 * DENDRO_752 * DENDRO_773 * alpha[pp] -
    2 * DENDRO_21 * DENDRO_31 * grad2_0_1_beta2[pp] -
    1.0 / 3.0 * DENDRO_23 * DENDRO_31 * grad2_0_1_beta0[pp] -
    1.0 / 3.0 * DENDRO_23 * DENDRO_31 * grad2_1_1_beta1[pp] -
    7.0 / 3.0 * DENDRO_23 * DENDRO_31 * grad2_1_2_beta2[pp] -
    2.0 * DENDRO_239 * DENDRO_766 * DENDRO_773 * alpha[pp] -
    2 * DENDRO_243 * DENDRO_759 * DENDRO_773 * alpha[pp] -
    2 * DENDRO_246 * DENDRO_756 * DENDRO_773 * alpha[pp] +
    DENDRO_285 * grad2_0_0_beta2[pp] + DENDRO_287 * grad2_1_1_beta2[pp] +
    DENDRO_289 * DENDRO_799 + (1.0 / 3.0) * DENDRO_289 * grad2_1_2_beta1[pp] +
    (4.0 / 3.0) * DENDRO_289 * grad2_2_2_beta2[pp] +
    DENDRO_356 * grad_1_beta2[pp] + DENDRO_358 * grad_0_beta2[pp] +
    DENDRO_360 * grad_2_beta2[pp] + DENDRO_759 * DENDRO_776 +
    DENDRO_761 * DENDRO_772 + DENDRO_766 * DENDRO_777 +
    DENDRO_768 * grad2_0_2_beta2[pp] + DENDRO_769 * grad2_0_0_beta0[pp] +
    DENDRO_769 * grad2_0_1_beta1[pp] +
    DENDRO_775 * (DENDRO_23 * DENDRO_779 + DENDRO_766 * DENDRO_781) +
    DENDRO_775 * (9 * DENDRO_31 * DENDRO_45 * DENDRO_759 * grad_2_chi[pp] -
                  DENDRO_39 * DENDRO_778) +
    DENDRO_775 * (9 * DENDRO_31 * DENDRO_45 * DENDRO_761 * grad_0_chi[pp] -
                  DENDRO_33 * DENDRO_774) -
    beta0[pp] * grad_0_Gt2[pp] - beta1[pp] * grad_1_Gt2[pp] -
    beta2[pp] * grad_2_Gt2[pp];
const double DENDRO_804 =
    BSSN_ETA_R0 *
    sqrt(DENDRO_31 * (2 * DENDRO_21 * grad_0_chi[pp] * grad_1_chi[pp] +
                      2 * DENDRO_23 * grad_1_chi[pp] * grad_2_chi[pp] -
                      DENDRO_258 * DENDRO_26 - DENDRO_259 * DENDRO_39 -
                      DENDRO_35 * DENDRO_78 - 2 * DENDRO_71 * grad_0_chi[pp])) *
    pow(1 - pow(chi[pp], BSSN_ETA_POWER[0]), -BSSN_ETA_POWER[1]);

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
chi_rhs[pp] = -DENDRO_16 * DENDRO_17 + DENDRO_16 * K[pp] * alpha[pp] +
              beta0[pp] * grad_0_chi[pp] + beta1[pp] * grad_1_chi[pp] +
              beta2[pp] * grad_2_chi[pp];
//--
At_rhs00[pp] =
    -At0[pp] * DENDRO_10 + At0[pp] * DENDRO_3 - At0[pp] * DENDRO_8 +
    DENDRO_18 * grad_0_beta1[pp] + DENDRO_19 * grad_0_beta2[pp] +
    DENDRO_669 *
        (DENDRO_265 *
             (-DENDRO_107 - DENDRO_110 * DENDRO_114 * DENDRO_116 - DENDRO_123 +
              4 * DENDRO_126 * DENDRO_129 * DENDRO_26 * DENDRO_85 - DENDRO_132 +
              4 * DENDRO_133 * DENDRO_35 * DENDRO_61 * DENDRO_85 -
              DENDRO_137 * (DENDRO_135 + DENDRO_51 * DENDRO_90) -
              DENDRO_137 * (DENDRO_139 + DENDRO_48 * DENDRO_95) -
              DENDRO_137 * (DENDRO_169 + DENDRO_172) -
              DENDRO_137 * (DENDRO_176 + DENDRO_179) -
              DENDRO_137 * (DENDRO_152 * DENDRO_156 + DENDRO_155) -
              DENDRO_143 * (DENDRO_140 + DENDRO_141) -
              DENDRO_143 * (DENDRO_165 + DENDRO_166) +
              2.0 * DENDRO_148 * DENDRO_21 * DENDRO_85 +
              4 * DENDRO_154 * DENDRO_23 * DENDRO_85 -
              DENDRO_158 * (1.0 * DENDRO_141 + DENDRO_157) -
              DENDRO_158 *
                  (DENDRO_129 * DENDRO_184 + 1.0 * DENDRO_152 * DENDRO_61) -
              DENDRO_164 -
              DENDRO_182 *
                  (-DENDRO_156 * DENDRO_52 + 2 * DENDRO_174 * DENDRO_61) -
              DENDRO_183 + 4 * DENDRO_186 * DENDRO_21 * DENDRO_85 -
              DENDRO_193 * (DENDRO_187 + DENDRO_188 * DENDRO_189 +
                            DENDRO_190 * DENDRO_55 + DENDRO_191 * DENDRO_192) -
              DENDRO_218 * DENDRO_219 - DENDRO_233 * DENDRO_234 -
              DENDRO_250 * DENDRO_251 +
              2.0 * DENDRO_26 * DENDRO_31 * grad2_1_1_gt0[pp] +
              3.0 * DENDRO_26 * DENDRO_85 * DENDRO_90 * grad_1_gt0[pp] -
              DENDRO_263 * DENDRO_264 +
              4 * DENDRO_31 * DENDRO_33 * grad2_0_2_gt0[pp] +
              2.0 * DENDRO_31 * DENDRO_35 * grad2_0_0_gt0[pp] +
              2.0 * DENDRO_31 * DENDRO_39 * grad2_2_2_gt0[pp] +
              2.0 * DENDRO_33 * DENDRO_85 * (DENDRO_144 + DENDRO_145) +
              4 * DENDRO_33 * DENDRO_85 * (1.0 * DENDRO_145 + DENDRO_159) +
              2.0 * DENDRO_33 * DENDRO_85 * (DENDRO_149 + DENDRO_150) +
              6.0 * DENDRO_35 * DENDRO_85 * DENDRO_96 * grad_0_gt0[pp] +
              3.0 * DENDRO_39 * DENDRO_85 * DENDRO_95 * grad_2_gt0[pp] -
              DENDRO_79 - DENDRO_82 - DENDRO_84 + 4 * grad_0_Gt0[pp] * gt0[pp] +
              4 * grad_0_Gt1[pp] * gt1[pp] + 4 * grad_0_Gt2[pp] * gt2[pp]) +
         DENDRO_56 * DENDRO_58 - DENDRO_62 * DENDRO_64 + DENDRO_668 * gt0[pp] -
         DENDRO_75 * DENDRO_76 - 12 * grad2_0_0_alpha[pp]) -
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
    DENDRO_669 *
        (DENDRO_265 *
             (-DENDRO_106 * (DENDRO_168 * DENDRO_214 + DENDRO_660) -
              DENDRO_106 * (-DENDRO_102 * DENDRO_90 + DENDRO_661 + DENDRO_686) -
              DENDRO_106 * (-DENDRO_129 * DENDRO_297 +
                            0.5 * DENDRO_152 * DENDRO_246 - DENDRO_687) +
              DENDRO_115 * (DENDRO_682 + DENDRO_685) + DENDRO_116 * DENDRO_659 -
              DENDRO_137 * DENDRO_663 +
              DENDRO_137 * (-DENDRO_153 * DENDRO_214 + DENDRO_664) +
              DENDRO_137 * (DENDRO_385 * DENDRO_95 - DENDRO_577 * DENDRO_90 +
                            DENDRO_604) -
              DENDRO_143 *
                  (DENDRO_90 * grad_1_gt0[pp] + DENDRO_96 * grad_0_gt3[pp]) +
              DENDRO_158 * (-DENDRO_214 * DENDRO_48 + DENDRO_665) -
              DENDRO_158 * (DENDRO_129 * DENDRO_470 + DENDRO_129 * DENDRO_485 +
                            DENDRO_212 * DENDRO_61) +
              DENDRO_158 * (-DENDRO_129 * DENDRO_577 - DENDRO_213 * DENDRO_61 +
                            0.5 * DENDRO_246 * grad_2_gt0[pp]) +
              DENDRO_158 * (DENDRO_168 * DENDRO_96 - DENDRO_180 * DENDRO_90 +
                            DENDRO_628) +
              DENDRO_182 * DENDRO_662 + DENDRO_182 * (DENDRO_688 + DENDRO_689) +
              DENDRO_182 * (DENDRO_691 + DENDRO_693) +
              DENDRO_182 * (DENDRO_380 * DENDRO_96 + DENDRO_690) -
              DENDRO_193 * (DENDRO_101 * DENDRO_496 + DENDRO_175 * DENDRO_500 +
                            DENDRO_224 * DENDRO_494 + DENDRO_633) +
              DENDRO_300 *
                  (DENDRO_157 + DENDRO_48 * DENDRO_96 + DENDRO_49 * DENDRO_96) +
              DENDRO_300 * (0.25 * DENDRO_165 + 0.5 * DENDRO_166 +
                            DENDRO_380 * DENDRO_61) +
              DENDRO_640 + DENDRO_641 + DENDRO_642 + DENDRO_643 + DENDRO_644 +
              DENDRO_645 + DENDRO_646 + DENDRO_647 + DENDRO_648 + DENDRO_649 +
              DENDRO_650 + DENDRO_651 + DENDRO_652 - DENDRO_653 * DENDRO_701 +
              DENDRO_666 -
              DENDRO_696 * (DENDRO_695 + DENDRO_95 * grad_0_gt3[pp]) -
              DENDRO_696 * (DENDRO_110 * grad_2_gt3[pp] +
                            DENDRO_129 * grad_1_gt5[pp] + DENDRO_697) -
              DENDRO_698 * grad_0_gt1[pp] - DENDRO_699 * grad_1_gt1[pp] -
              DENDRO_700 * grad_2_gt1[pp]) +
         DENDRO_31 * DENDRO_679 * (-DENDRO_129 - DENDRO_59 * DENDRO_657) +
         DENDRO_653 * DENDRO_667 - DENDRO_655 * DENDRO_680 -
         DENDRO_656 * DENDRO_681 - 12 * grad2_0_1_alpha[pp]) +
    DENDRO_670 * grad_0_beta0[pp] + DENDRO_670 * grad_1_beta1[pp] -
    alpha[pp] * (-At1[pp] * K[pp] + DENDRO_32 * DENDRO_674 +
                 DENDRO_38 * DENDRO_676 + DENDRO_41 * DENDRO_678) +
    beta0[pp] * grad_0_At1[pp] + beta1[pp] * grad_1_At1[pp] +
    beta2[pp] * grad_2_At1[pp];
//--
At_rhs02[pp] =
    At0[pp] * grad_2_beta0[pp] + At1[pp] * grad_2_beta1[pp] -
    At2[pp] * DENDRO_10 + At4[pp] * grad_0_beta1[pp] +
    At5[pp] * grad_0_beta2[pp] +
    DENDRO_669 *
        (DENDRO_265 *
             (-DENDRO_116 * (DENDRO_174 * DENDRO_242 - 0.5 * DENDRO_708) -
              DENDRO_116 * (-DENDRO_111 * DENDRO_95 - DENDRO_225 * DENDRO_51 -
                            DENDRO_715) -
              DENDRO_137 * DENDRO_533 -
              DENDRO_137 * (DENDRO_153 * DENDRO_242 + DENDRO_712) -
              DENDRO_158 * (DENDRO_155 + DENDRO_693) -
              DENDRO_158 * (DENDRO_689 + DENDRO_717) -
              DENDRO_158 * (DENDRO_320 * DENDRO_96 + DENDRO_690) -
              DENDRO_158 * (DENDRO_110 * DENDRO_485 + DENDRO_155 + DENDRO_692) -
              DENDRO_182 * (0.5 * DENDRO_110 * DENDRO_174 - DENDRO_714) -
              DENDRO_182 * (DENDRO_174 * DENDRO_96 - DENDRO_184 * DENDRO_95 -
                            DENDRO_225 * DENDRO_53) -
              DENDRO_193 * (DENDRO_118 * DENDRO_496 + DENDRO_177 * DENDRO_500 +
                            DENDRO_220 * DENDRO_494 + DENDRO_491) +
              4 * DENDRO_21 * DENDRO_536 * DENDRO_85 +
              4 * DENDRO_21 * DENDRO_543 * DENDRO_85 +
              DENDRO_23 * DENDRO_524 * DENDRO_85 +
              4 * DENDRO_23 * DENDRO_85 *
                  (0.5 * DENDRO_174 * DENDRO_238 - DENDRO_712) +
              4 * DENDRO_23 * DENDRO_85 *
                  (-DENDRO_225 * DENDRO_48 - DENDRO_470 * DENDRO_95 -
                   DENDRO_716) +
              2.0 * DENDRO_26 * DENDRO_31 * grad2_1_1_gt2[pp] +
              DENDRO_26 * DENDRO_85 * (DENDRO_695 + DENDRO_713) +
              4 * DENDRO_26 * DENDRO_85 *
                  (DENDRO_129 * DENDRO_202 + 0.25 * DENDRO_697) +
              4 * DENDRO_31 * DENDRO_33 * grad2_0_2_gt2[pp] +
              2.0 * DENDRO_31 * DENDRO_35 * grad2_0_0_gt2[pp] +
              2.0 * DENDRO_31 * DENDRO_39 * grad2_2_2_gt2[pp] +
              4 * DENDRO_33 * DENDRO_85 *
                  (DENDRO_242 * DENDRO_51 + DENDRO_714) +
              2.0 * DENDRO_33 * DENDRO_85 *
                  (DENDRO_95 * grad_2_gt0[pp] + DENDRO_96 * grad_0_gt5[pp]) +
              4 * DENDRO_35 * DENDRO_85 *
                  (0.25 * DENDRO_149 + 1.0 * DENDRO_150) +
              4 * DENDRO_35 * DENDRO_85 *
                  (DENDRO_159 + DENDRO_51 * DENDRO_96 + DENDRO_52 * DENDRO_96) -
              DENDRO_507 - DENDRO_508 - DENDRO_509 - DENDRO_510 * DENDRO_701 -
              DENDRO_525 - DENDRO_527 - DENDRO_528 - DENDRO_529 - DENDRO_531 -
              DENDRO_696 * (DENDRO_685 + DENDRO_709) -
              DENDRO_698 * grad_0_gt2[pp] - DENDRO_699 * grad_1_gt2[pp] -
              DENDRO_700 * grad_2_gt2[pp] + 2.0 * grad_0_Gt0[pp] * gt2[pp] +
              2.0 * grad_0_Gt1[pp] * gt4[pp] + 2.0 * grad_0_Gt2[pp] * gt5[pp] +
              2.0 * grad_2_Gt0[pp] * gt0[pp] + 2.0 * grad_2_Gt1[pp] * gt1[pp] +
              2.0 * grad_2_Gt2[pp] * gt2[pp]) +
         DENDRO_510 * DENDRO_667 - DENDRO_512 * DENDRO_680 -
         DENDRO_514 * DENDRO_679 + DENDRO_518 * DENDRO_681 -
         12 * grad2_0_2_alpha[pp]) +
    DENDRO_702 * grad_0_beta0[pp] + DENDRO_702 * grad_2_beta2[pp] -
    alpha[pp] * (-At2[pp] * K[pp] + DENDRO_32 * DENDRO_704 +
                 DENDRO_38 * DENDRO_706 + DENDRO_41 * DENDRO_707) +
    beta0[pp] * grad_0_At2[pp] + beta1[pp] * grad_1_At2[pp] +
    beta2[pp] * grad_2_At2[pp];
//--
At_rhs11[pp] =
    -At3[pp] * DENDRO_11 + At3[pp] * DENDRO_12 - At3[pp] * DENDRO_8 +
    DENDRO_18 * grad_1_beta0[pp] +
    DENDRO_669 *
        (DENDRO_265 *
             (6.0 * DENDRO_105 * DENDRO_214 * grad_1_gt3[pp] -
              DENDRO_116 * DENDRO_238 * DENDRO_403 + DENDRO_129 * DENDRO_407 +
              DENDRO_137 * (DENDRO_391 - 1.0 * DENDRO_723) -
              DENDRO_137 * (-1.0 * DENDRO_229 * DENDRO_93 + DENDRO_386) -
              DENDRO_137 * (DENDRO_238 * DENDRO_387 + DENDRO_389) +
              DENDRO_143 * (DENDRO_101 * grad_1_gt3[pp] - DENDRO_727) +
              DENDRO_143 *
                  (-DENDRO_129 * grad_2_gt3[pp] + DENDRO_246 * DENDRO_88) +
              DENDRO_143 *
                  (DENDRO_229 * grad_1_gt0[pp] - DENDRO_90 * grad_0_gt3[pp]) +
              DENDRO_158 * (DENDRO_395 + DENDRO_686) +
              DENDRO_158 *
                  (0.25 * DENDRO_101 * grad_1_gt3[pp] - 1.0 * DENDRO_727) +
              DENDRO_158 * (DENDRO_246 * DENDRO_93 - DENDRO_687) +
              DENDRO_182 * (DENDRO_382 + DENDRO_384) +
              DENDRO_182 * (DENDRO_134 * DENDRO_222 + DENDRO_380 * DENDRO_90) +
              DENDRO_182 * (DENDRO_238 * DENDRO_380 + DENDRO_724) +
              DENDRO_182 * (DENDRO_49 * DENDRO_725 + DENDRO_726) -
              DENDRO_193 * (DENDRO_189 * DENDRO_229 + DENDRO_190 * DENDRO_215 +
                            DENDRO_192 * DENDRO_246 + DENDRO_369) +
              DENDRO_222 * DENDRO_405 + DENDRO_313 * (DENDRO_390 - DENDRO_723) +
              DENDRO_313 *
                  (DENDRO_152 * DENDRO_229 - DENDRO_222 * grad_0_gt3[pp]) +
              DENDRO_313 * (-DENDRO_238 * grad_2_gt3[pp] + DENDRO_401) -
              DENDRO_367 * DENDRO_701 + DENDRO_406 * DENDRO_90 + DENDRO_414 -
              DENDRO_698 * grad_0_gt3[pp] - DENDRO_699 * grad_1_gt3[pp] -
              DENDRO_700 * grad_2_gt3[pp]) -
         DENDRO_368 * DENDRO_57 +
         DENDRO_64 * (DENDRO_246 - DENDRO_59 * DENDRO_721) +
         DENDRO_668 * gt3[pp] +
         DENDRO_722 * (DENDRO_229 - DENDRO_72 * DENDRO_721) -
         12 * grad2_1_1_alpha[pp]) +
    DENDRO_718 * grad_1_beta2[pp] -
    alpha[pp] * (-At3[pp] * K[pp] + DENDRO_32 * DENDRO_676 +
                 DENDRO_674 * DENDRO_719 + DENDRO_678 * DENDRO_720) +
    beta0[pp] * grad_0_At3[pp] + beta1[pp] * grad_1_At3[pp] +
    beta2[pp] * grad_2_At3[pp];
//--
At_rhs12[pp] =
    At1[pp] * grad_2_beta0[pp] + At2[pp] * grad_1_beta0[pp] +
    At3[pp] * grad_2_beta1[pp] - At4[pp] * DENDRO_11 +
    At5[pp] * grad_1_beta2[pp] +
    DENDRO_669 *
        (DENDRO_265 *
             (-DENDRO_106 * (-DENDRO_238 * DENDRO_297 + DENDRO_596) -
              DENDRO_106 *
                  (-DENDRO_102 * DENDRO_222 + DENDRO_552 + DENDRO_597) -
              DENDRO_106 * (0.25 * DENDRO_199 * grad_1_gt3[pp] -
                            DENDRO_212 * DENDRO_214 - DENDRO_213 * DENDRO_214) -
              DENDRO_116 * (DENDRO_205 * DENDRO_242 - 0.5 * DENDRO_729) -
              DENDRO_116 *
                  (-DENDRO_111 * DENDRO_222 + 0.5 * DENDRO_174 * DENDRO_222 -
                   DENDRO_225 * DENDRO_380) +
              DENDRO_137 * (DENDRO_205 * DENDRO_214 + DENDRO_601) -
              DENDRO_137 * (DENDRO_167 * DENDRO_225 + DENDRO_222 * DENDRO_470 +
                            DENDRO_600) +
              DENDRO_137 * (-DENDRO_212 * DENDRO_242 +
                            0.5 * DENDRO_246 * grad_2_gt5[pp] - DENDRO_732) +
              DENDRO_137 * (-DENDRO_222 * DENDRO_485 - DENDRO_222 * DENDRO_577 +
                            0.5 * DENDRO_229 * grad_0_gt5[pp]) +
              DENDRO_137 * (DENDRO_238 * DENDRO_336 + DENDRO_564 - DENDRO_732) +
              DENDRO_158 * (-DENDRO_110 * DENDRO_297 - DENDRO_733) +
              DENDRO_158 * (-DENDRO_214 * DENDRO_320 + DENDRO_603) +
              DENDRO_158 * (-DENDRO_238 * DENDRO_577 - DENDRO_733) +
              DENDRO_158 *
                  (-DENDRO_102 * DENDRO_95 + 0.5 * DENDRO_229 * grad_2_gt0[pp] -
                   0.25 * DENDRO_713) +
              DENDRO_158 * (-DENDRO_180 * DENDRO_222 + DENDRO_229 * DENDRO_52 -
                            DENDRO_726) +
              DENDRO_158 *
                  (-DENDRO_214 * DENDRO_380 + DENDRO_446 + DENDRO_574) +
              DENDRO_162 * (DENDRO_682 + DENDRO_683 + DENDRO_709) -
              DENDRO_182 * (0.5 * DENDRO_110 * DENDRO_205 - DENDRO_731) +
              DENDRO_182 * (DENDRO_242 * DENDRO_380 + DENDRO_731) -
              DENDRO_182 * (-DENDRO_184 * DENDRO_222 - DENDRO_225 * DENDRO_49 -
                            DENDRO_716) -
              DENDRO_193 * (DENDRO_199 * DENDRO_496 + DENDRO_223 * DENDRO_494 +
                            DENDRO_239 * DENDRO_500 + DENDRO_582) +
              DENDRO_300 * (DENDRO_129 * DENDRO_173 + DENDRO_691) +
              DENDRO_300 * (DENDRO_139 + DENDRO_688 + DENDRO_717) +
              DENDRO_313 *
                  (DENDRO_199 * grad_2_gt3[pp] - DENDRO_214 * grad_1_gt5[pp]) -
              DENDRO_588 * DENDRO_701 + DENDRO_589 + DENDRO_609 -
              DENDRO_698 * grad_0_gt4[pp] - DENDRO_699 * grad_1_gt4[pp] -
              DENDRO_700 * grad_2_gt4[pp]) -
         DENDRO_31 * DENDRO_594 * DENDRO_680 + DENDRO_588 * DENDRO_667 +
         DENDRO_592 * DENDRO_681 - DENDRO_593 * DENDRO_679 -
         12 * grad2_1_2_alpha[pp]) +
    DENDRO_728 * grad_1_beta1[pp] + DENDRO_728 * grad_2_beta2[pp] -
    alpha[pp] * (-At4[pp] * K[pp] + DENDRO_32 * DENDRO_706 +
                 DENDRO_704 * DENDRO_719 + DENDRO_707 * DENDRO_720) +
    beta0[pp] * grad_0_At4[pp] + beta1[pp] * grad_1_At4[pp] +
    beta2[pp] * grad_2_At4[pp];
//--
At_rhs22[pp] =
    -At5[pp] * DENDRO_10 - At5[pp] * DENDRO_11 + At5[pp] * DENDRO_15 +
    DENDRO_19 * grad_2_beta0[pp] +
    DENDRO_669 *
        (DENDRO_265 *
             (3.0 * DENDRO_110 * DENDRO_35 * DENDRO_85 * grad_0_gt5[pp] +
              DENDRO_116 * DENDRO_225 * DENDRO_305 -
              DENDRO_137 * (0.25 * DENDRO_729 + 1.0 * DENDRO_735) -
              DENDRO_137 * (-1.0 * DENDRO_226 * DENDRO_88 + DENDRO_339) -
              DENDRO_158 * (DENDRO_110 * DENDRO_202 + DENDRO_730) -
              DENDRO_158 * (DENDRO_138 * DENDRO_222 + DENDRO_320 * DENDRO_95) -
              DENDRO_158 * (DENDRO_173 * DENDRO_238 + DENDRO_710) -
              DENDRO_158 * (DENDRO_52 * DENDRO_725 + 0.25 * DENDRO_684) -
              DENDRO_182 * (-DENDRO_225 * DENDRO_343 - DENDRO_715) -
              DENDRO_193 * (DENDRO_189 * DENDRO_226 + DENDRO_190 * DENDRO_206 +
                            DENDRO_192 * DENDRO_243 + DENDRO_344) +
              4 * DENDRO_21 * DENDRO_321 * DENDRO_85 +
              4 * DENDRO_21 * DENDRO_331 * DENDRO_85 +
              4 * DENDRO_222 * DENDRO_26 * DENDRO_304 * DENDRO_85 +
              2.0 * DENDRO_23 * DENDRO_312 * DENDRO_85 +
              4 * DENDRO_23 * DENDRO_341 * DENDRO_85 +
              3.0 * DENDRO_238 * DENDRO_26 * DENDRO_85 * grad_1_gt5[pp] +
              6.0 * DENDRO_242 * DENDRO_39 * DENDRO_85 * grad_2_gt5[pp] +
              2.0 * DENDRO_26 * DENDRO_31 * grad2_1_1_gt5[pp] -
              DENDRO_274 * DENDRO_701 - DENDRO_280 - DENDRO_281 - DENDRO_282 -
              DENDRO_298 + 4 * DENDRO_299 * DENDRO_35 * DENDRO_85 * DENDRO_95 -
              DENDRO_303 - DENDRO_306 +
              4 * DENDRO_31 * DENDRO_33 * grad2_0_2_gt5[pp] +
              2.0 * DENDRO_31 * DENDRO_35 * grad2_0_0_gt5[pp] +
              2.0 * DENDRO_31 * DENDRO_39 * grad2_2_2_gt5[pp] -
              DENDRO_313 * (DENDRO_729 + DENDRO_735) -
              DENDRO_313 *
                  (DENDRO_152 * DENDRO_225 + DENDRO_222 * grad_0_gt5[pp]) -
              DENDRO_327 +
              4 * DENDRO_33 * DENDRO_85 *
                  (0.25 * DENDRO_708 + 1.0 * DENDRO_734) +
              2.0 * DENDRO_33 * DENDRO_85 * (DENDRO_708 + DENDRO_734) +
              4 * DENDRO_33 * DENDRO_85 *
                  (-1.0 * DENDRO_206 * DENDRO_88 + DENDRO_337) +
              2.0 * DENDRO_33 * DENDRO_85 *
                  (DENDRO_225 * grad_2_gt0[pp] + DENDRO_95 * grad_0_gt5[pp]) -
              DENDRO_698 * grad_0_gt5[pp] - DENDRO_699 * grad_1_gt5[pp] -
              DENDRO_700 * grad_2_gt5[pp] + 4 * grad_2_Gt0[pp] * gt2[pp] +
              4 * grad_2_Gt1[pp] * gt4[pp] + 4 * grad_2_Gt2[pp] * gt5[pp]) +
         DENDRO_267 * DENDRO_58 - DENDRO_270 * DENDRO_722 -
         DENDRO_276 * DENDRO_63 + DENDRO_668 * gt5[pp] -
         12 * grad2_2_2_alpha[pp]) +
    DENDRO_718 * grad_2_beta1[pp] -
    alpha[pp] * (At5[pp] * DENDRO_37 * DENDRO_707 - At5[pp] * K[pp] +
                 DENDRO_41 * DENDRO_706 + DENDRO_704 * DENDRO_720) +
    beta0[pp] * grad_0_At5[pp] + beta1[pp] * grad_1_At5[pp] +
    beta2[pp] * grad_2_At5[pp];
//--
K_rhs[pp] =
    DENDRO_254 * DENDRO_31 * chi[pp] *
        (0.5 * DENDRO_743 * (DENDRO_638 + DENDRO_657 * DENDRO_741) +
         DENDRO_745 *
             (DENDRO_31 * DENDRO_636 + DENDRO_31 * DENDRO_98 +
              DENDRO_31 * DENDRO_99 -
              DENDRO_45 * (-DENDRO_653 * DENDRO_739 + grad_0_chi[pp])) +
         DENDRO_747 *
             (DENDRO_31 * DENDRO_634 + DENDRO_31 * DENDRO_635 +
              DENDRO_31 * DENDRO_86 -
              DENDRO_45 * (-DENDRO_653 * DENDRO_736 + grad_1_chi[pp])) -
         grad2_0_1_alpha[pp]) +
    DENDRO_257 * DENDRO_31 * chi[pp] *
        (0.5 * DENDRO_738 * (DENDRO_45 * DENDRO_736 * gt4[pp] + DENDRO_583) +
         DENDRO_745 * (DENDRO_31 * DENDRO_584 -
                       DENDRO_45 * (-DENDRO_588 * DENDRO_739 + grad_2_chi[pp]) +
                       DENDRO_591) +
         DENDRO_746 *
             (DENDRO_31 * DENDRO_585 + DENDRO_31 * DENDRO_586 +
              DENDRO_31 * DENDRO_587 -
              DENDRO_45 * (-DENDRO_588 * DENDRO_741 + grad_1_chi[pp])) -
         grad2_1_2_alpha[pp]) -
    DENDRO_285 * chi[pp] *
        (DENDRO_740 * (DENDRO_432 + DENDRO_46 * DENDRO_744) +
         DENDRO_743 * (DENDRO_433 + DENDRO_46 * DENDRO_742) -
         grad2_0_0_alpha[pp] +
         grad_0_alpha[pp] *
             (DENDRO_31 * DENDRO_430 + DENDRO_31 * DENDRO_431 -
              DENDRO_45 * (DENDRO_69 - 0.5 * DENDRO_70 * DENDRO_736) +
              DENDRO_68)) -
    DENDRO_287 * chi[pp] *
        (DENDRO_738 * (DENDRO_371 + DENDRO_721 * DENDRO_736) +
         DENDRO_743 * (DENDRO_375 + DENDRO_721 * DENDRO_741) -
         grad2_1_1_alpha[pp] +
         grad_1_alpha[pp] *
             (DENDRO_31 * DENDRO_372 + DENDRO_31 * DENDRO_373 +
              DENDRO_31 * DENDRO_374 -
              DENDRO_45 * (DENDRO_366 - DENDRO_367 * DENDRO_744))) -
    DENDRO_289 * chi[pp] *
        (DENDRO_738 * (DENDRO_349 + DENDRO_736 * DENDRO_737) +
         DENDRO_740 * (DENDRO_351 + DENDRO_737 * DENDRO_739) -
         grad2_2_2_alpha[pp] +
         grad_2_alpha[pp] *
             (DENDRO_31 * DENDRO_353 + DENDRO_31 * DENDRO_354 +
              DENDRO_31 * DENDRO_355 -
              DENDRO_45 * (DENDRO_273 - DENDRO_274 * DENDRO_742))) -
    DENDRO_748 * chi[pp] *
        (0.5 * DENDRO_740 * (DENDRO_495 + DENDRO_517 * DENDRO_739) +
         DENDRO_746 *
             (DENDRO_31 * DENDRO_497 + DENDRO_31 * DENDRO_498 +
              DENDRO_31 * DENDRO_499 -
              DENDRO_45 * (-DENDRO_510 * DENDRO_741 + grad_0_chi[pp])) +
         DENDRO_747 *
             (DENDRO_31 * DENDRO_492 + DENDRO_31 * DENDRO_493 +
              DENDRO_31 * DENDRO_94 -
              DENDRO_45 * (-DENDRO_510 * DENDRO_736 + grad_2_chi[pp])) -
         grad2_0_2_alpha[pp]) +
    DENDRO_767 *
        (At0[pp] * DENDRO_752 * DENDRO_753 + At1[pp] * DENDRO_762 * DENDRO_764 +
         At2[pp] * DENDRO_761 * DENDRO_762 + At3[pp] * DENDRO_753 * DENDRO_756 +
         At4[pp] * DENDRO_762 * DENDRO_766 + At5[pp] * DENDRO_753 * DENDRO_759 +
         pow(K[pp], 2)) +
    beta0[pp] * grad_0_K[pp] + beta1[pp] * grad_1_K[pp] +
    beta2[pp] * grad_2_K[pp];
//--
Gt_rhs0[pp] = -DENDRO_782;
//--
Gt_rhs1[pp] =
    -DENDRO_101 * DENDRO_763 * DENDRO_801 +
    DENDRO_118 * DENDRO_761 * DENDRO_801 - 2.0 / 3.0 * DENDRO_17 * DENDRO_218 -
    DENDRO_199 * DENDRO_765 * DENDRO_801 +
    DENDRO_206 * DENDRO_759 * DENDRO_800 -
    DENDRO_214 * DENDRO_756 * DENDRO_800 + DENDRO_218 * grad_1_beta1[pp] +
    DENDRO_233 * grad_0_beta1[pp] + DENDRO_250 * grad_2_beta1[pp] +
    DENDRO_55 * DENDRO_752 * DENDRO_800 + DENDRO_763 * DENDRO_772 +
    DENDRO_765 * DENDRO_776 + DENDRO_775 * DENDRO_795 -
    DENDRO_775 * (-DENDRO_763 * DENDRO_789 + DENDRO_788) -
    DENDRO_775 * (-DENDRO_765 * DENDRO_791 + DENDRO_790) - DENDRO_783 -
    DENDRO_784 - DENDRO_785 - DENDRO_786 - DENDRO_787 - DENDRO_793 -
    DENDRO_794 + (7.0 / 3.0) * DENDRO_796 * grad2_0_1_beta1[pp] +
    DENDRO_797 * grad2_0_0_beta0[pp] + DENDRO_797 * grad2_0_2_beta2[pp] +
    DENDRO_798 * DENDRO_799 + (7.0 / 3.0) * DENDRO_798 * grad2_1_2_beta1[pp] +
    (1.0 / 3.0) * DENDRO_798 * grad2_2_2_beta2[pp] + DENDRO_802;
//--
Gt_rhs2[pp] = -DENDRO_803;
//--
B_rhs0[pp] =
    -B0[pp] * DENDRO_804 - DENDRO_782 +
    lambda[2] * (beta0[pp] * grad_0_B0[pp] + beta1[pp] * grad_1_B0[pp] +
                 beta2[pp] * grad_2_B0[pp]) -
    lambda[3] * (beta0[pp] * grad_0_Gt0[pp] + beta1[pp] * grad_1_Gt0[pp] +
                 beta2[pp] * grad_2_Gt0[pp]);
//--
B_rhs1[pp] =
    -B1[pp] * DENDRO_804 +
    2.0 * DENDRO_101 * DENDRO_764 * DENDRO_773 * alpha[pp] +
    2.0 * DENDRO_118 * DENDRO_761 * DENDRO_773 * alpha[pp] +
    (2.0 / 3.0) * DENDRO_17 * DENDRO_260 * DENDRO_85 +
    2.0 * DENDRO_199 * DENDRO_766 * DENDRO_773 * alpha[pp] +
    2 * DENDRO_206 * DENDRO_759 * DENDRO_773 * alpha[pp] +
    (1.0 / 3.0) * DENDRO_21 * DENDRO_31 * grad2_0_0_beta0[pp] +
    (7.0 / 3.0) * DENDRO_21 * DENDRO_31 * grad2_0_1_beta1[pp] +
    (1.0 / 3.0) * DENDRO_21 * DENDRO_31 * grad2_0_2_beta2[pp] +
    2 * DENDRO_215 * DENDRO_756 * DENDRO_773 * alpha[pp] +
    (1.0 / 3.0) * DENDRO_23 * DENDRO_31 * grad2_0_2_beta0[pp] +
    (7.0 / 3.0) * DENDRO_23 * DENDRO_31 * grad2_1_2_beta1[pp] +
    (1.0 / 3.0) * DENDRO_23 * DENDRO_31 * grad2_2_2_beta2[pp] -
    DENDRO_356 * grad_1_beta1[pp] - DENDRO_358 * grad_0_beta1[pp] -
    DENDRO_360 * grad_2_beta1[pp] +
    2 * DENDRO_55 * DENDRO_752 * DENDRO_773 * alpha[pp] -
    DENDRO_764 * DENDRO_772 - DENDRO_766 * DENDRO_776 +
    DENDRO_775 * DENDRO_795 -
    DENDRO_775 * (DENDRO_764 * DENDRO_789 + DENDRO_788) -
    DENDRO_775 * (DENDRO_766 * DENDRO_791 + DENDRO_790) - DENDRO_783 -
    DENDRO_784 - DENDRO_785 - DENDRO_786 - DENDRO_787 - DENDRO_793 -
    DENDRO_794 - DENDRO_802 * lambda[3] + beta0[pp] * grad_0_Gt1[pp] +
    beta1[pp] * grad_1_Gt1[pp] + beta2[pp] * grad_2_Gt1[pp] +
    lambda[2] * (beta0[pp] * grad_0_B1[pp] + beta1[pp] * grad_1_B1[pp] +
                 beta2[pp] * grad_2_B1[pp]);
//--
B_rhs2[pp] =
    -B2[pp] * DENDRO_804 - DENDRO_803 +
    lambda[2] * (beta0[pp] * grad_0_B2[pp] + beta1[pp] * grad_1_B2[pp] +
                 beta2[pp] * grad_2_B2[pp]) -
    lambda[3] * (beta0[pp] * grad_0_Gt2[pp] + beta1[pp] * grad_1_Gt2[pp] +
                 beta2[pp] * grad_2_Gt2[pp]);
// Dendro: reduced ops: 4330
// Dendro: }}}
