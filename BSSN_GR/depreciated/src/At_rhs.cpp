bssn::timer::t_rhs.start();
for (unsigned int k = PW; k < nz - PW; k++) {
    z = pmin[2] + k * hz;
    for (unsigned int j = PW; j < ny - PW; j++) {
        y = pmin[1] + j * hy;
        for (unsigned int i = PW; i < nx - PW; i++) {
            x       = pmin[0] + i * hx;
            pp      = i + nx * (j + ny * k);
            r_coord = sqrt(x * x + y * y + z * z);
            eta     = ETA_CONST;
            if (r_coord >= ETA_R0) {
                eta *= pow((ETA_R0 / r_coord), ETA_DAMPING_EXP);
            }
            // Dendro: {{{
            // Dendro: original ops:  630012
            // Dendro: printing temp variables
            double DENDRO_0  = grad_0_beta0[pp];
            double DENDRO_1  = (2.0L / 3.0L) * At0[pp];
            double DENDRO_2  = grad_1_beta1[pp];
            double DENDRO_3  = grad_2_beta2[pp];
            double DENDRO_4  = 2 * At1[pp];
            double DENDRO_5  = grad_0_beta1[pp];
            double DENDRO_6  = 2 * At2[pp];
            double DENDRO_7  = grad_0_beta2[pp];
            double DENDRO_8  = pow(gt4[pp], 2);
            double DENDRO_9  = pow(gt1[pp], 2);
            double DENDRO_10 = pow(gt2[pp], 2);
            double DENDRO_11 = gt0[pp] * gt3[pp];
            double DENDRO_12 = gt1[pp] * gt2[pp];
            double DENDRO_13 = DENDRO_10 * gt3[pp] - DENDRO_11 * gt5[pp] -
                               2 * DENDRO_12 * gt4[pp] + DENDRO_8 * gt0[pp] +
                               DENDRO_9 * gt5[pp];
            double DENDRO_14 = 1.0 / DENDRO_13;
            double DENDRO_15 = 2 * At1[pp] * DENDRO_14;
            double DENDRO_16 = gt1[pp] * gt5[pp] - gt2[pp] * gt4[pp];
            double DENDRO_17 = -DENDRO_12 + gt0[pp] * gt4[pp];
            double DENDRO_18 = gt0[pp] * gt5[pp];
            double DENDRO_19 = -DENDRO_10 + DENDRO_18;
            double DENDRO_20 = 2 * At0[pp] * DENDRO_14;
            double DENDRO_21 = At1[pp] * DENDRO_16;
            double DENDRO_22 = gt1[pp] * gt4[pp];
            double DENDRO_23 = gt2[pp] * gt3[pp];
            double DENDRO_24 = DENDRO_22 - DENDRO_23;
            double DENDRO_25 = -At2[pp] * DENDRO_24;
            double DENDRO_26 = gt3[pp] * gt5[pp];
            double DENDRO_27 = DENDRO_26 - DENDRO_8;
            double DENDRO_28 = 2 * At2[pp] * DENDRO_14;
            double DENDRO_29 = DENDRO_11 - DENDRO_9;
            double DENDRO_30 = (1.0L / 12.0L) * chi[pp];
            double DENDRO_31 = grad2_0_0_alpha[pp];
            double DENDRO_32 = grad_1_alpha[pp];
            double DENDRO_33 = 1.0 / chi[pp];
            double DENDRO_34 = grad_1_chi[pp];
            double DENDRO_35 = grad_2_chi[pp];
            double DENDRO_36 = grad_0_chi[pp];
            double DENDRO_37 = DENDRO_16 * DENDRO_36 + DENDRO_17 * DENDRO_35;
            double DENDRO_38 = -DENDRO_19 * DENDRO_34 + DENDRO_37;
            double DENDRO_39 = 0.5 * DENDRO_33 * DENDRO_38;
            double DENDRO_40 = grad_0_gt1[pp];
            double DENDRO_41 = 1.0 * DENDRO_40;
            double DENDRO_42 = grad_1_gt0[pp];
            double DENDRO_43 = 0.5 * DENDRO_42;
            double DENDRO_44 = DENDRO_41 - DENDRO_43;
            double DENDRO_45 = grad_0_gt0[pp];
            double DENDRO_46 = 0.5 * DENDRO_45;
            double DENDRO_47 = grad_0_gt2[pp];
            double DENDRO_48 = 1.0 * DENDRO_47;
            double DENDRO_49 = grad_2_gt0[pp];
            double DENDRO_50 = 0.5 * DENDRO_49;
            double DENDRO_51 = DENDRO_48 - DENDRO_50;
            double DENDRO_52 = DENDRO_16 * DENDRO_46 + DENDRO_17 * DENDRO_51;
            double DENDRO_53 = -DENDRO_19 * DENDRO_44 + DENDRO_52;
            double DENDRO_54 =
                DENDRO_14 * DENDRO_32 * (DENDRO_39 * gt0[pp] + DENDRO_53);
            double DENDRO_55 = grad_2_alpha[pp];
            double DENDRO_56 = 12 * DENDRO_14 * DENDRO_55;
            double DENDRO_57 = DENDRO_17 * DENDRO_34;
            double DENDRO_58 = DENDRO_24 * DENDRO_36;
            double DENDRO_59 = DENDRO_29 * DENDRO_35;
            double DENDRO_60 = DENDRO_57 - DENDRO_58 - DENDRO_59;
            double DENDRO_61 = 0.5 * DENDRO_33 * DENDRO_60;
            double DENDRO_62 = DENDRO_61 * gt0[pp];
            double DENDRO_63 = DENDRO_24 * DENDRO_46;
            double DENDRO_64 = DENDRO_29 * DENDRO_51;
            double DENDRO_65 = DENDRO_17 * DENDRO_44;
            double DENDRO_66 = DENDRO_63 + DENDRO_64 - DENDRO_65;
            double DENDRO_67 = grad_0_alpha[pp];
            double DENDRO_68 = DENDRO_27 * DENDRO_46;
            double DENDRO_69 = DENDRO_14 * DENDRO_68;
            double DENDRO_70 = DENDRO_24 * DENDRO_51;
            double DENDRO_71 = DENDRO_14 * DENDRO_70;
            double DENDRO_72 = DENDRO_16 * DENDRO_44;
            double DENDRO_73 = DENDRO_14 * DENDRO_72;
            double DENDRO_74 = -DENDRO_22 + DENDRO_23;
            double DENDRO_75 = DENDRO_16 * DENDRO_34;
            double DENDRO_76 = -DENDRO_26 + DENDRO_8;
            double DENDRO_77 =
                DENDRO_35 * DENDRO_74 + DENDRO_36 * DENDRO_76 + DENDRO_75;
            double DENDRO_78 =
                DENDRO_33 *
                (0.5 * DENDRO_14 * DENDRO_77 * gt0[pp] - 1.0 * DENDRO_36);
            double DENDRO_79  = 3 * alpha[pp];
            double DENDRO_80  = grad_0_Gt0[pp];
            double DENDRO_81  = 4 * gt1[pp];
            double DENDRO_82  = grad_0_Gt1[pp];
            double DENDRO_83  = 4 * gt2[pp];
            double DENDRO_84  = grad_0_Gt2[pp];
            double DENDRO_85  = pow(chi[pp], -2);
            double DENDRO_86  = pow(DENDRO_36, 2);
            double DENDRO_87  = 4.0 * DENDRO_14 * DENDRO_17;
            double DENDRO_88  = 4 * DENDRO_14 * DENDRO_24;
            double DENDRO_89  = 4.0 * DENDRO_14 * DENDRO_16;
            double DENDRO_90  = 2.0 * DENDRO_14 * DENDRO_29;
            double DENDRO_91  = 2.0 * DENDRO_14 * DENDRO_19;
            double DENDRO_92  = 2.0 * DENDRO_14 * DENDRO_27;
            double DENDRO_93  = pow(DENDRO_13, -2);
            double DENDRO_94  = 4 * DENDRO_29 * DENDRO_93;
            double DENDRO_95  = grad_0_gt4[pp];
            double DENDRO_96  = 0.25 * DENDRO_95;
            double DENDRO_97  = -DENDRO_96;
            double DENDRO_98  = grad_1_gt2[pp];
            double DENDRO_99  = 0.25 * DENDRO_98;
            double DENDRO_100 = grad_2_gt1[pp];
            double DENDRO_101 = 0.75 * DENDRO_100;
            double DENDRO_102 = DENDRO_100 + DENDRO_95 - DENDRO_98;
            double DENDRO_103 = grad_0_gt5[pp];
            double DENDRO_104 = DENDRO_103 * DENDRO_17 + DENDRO_16 * DENDRO_49;
            double DENDRO_105 = -DENDRO_102 * DENDRO_19 + DENDRO_104;
            double DENDRO_106 = 4 * DENDRO_93;
            double DENDRO_107 = 2.0 * DENDRO_16 * DENDRO_93;
            double DENDRO_108 = grad_0_gt3[pp];
            double DENDRO_109 = DENDRO_108 * DENDRO_19;
            double DENDRO_110 = DENDRO_16 * DENDRO_42;
            double DENDRO_111 = -DENDRO_100 + DENDRO_95 + DENDRO_98;
            double DENDRO_112 = DENDRO_111 * DENDRO_17;
            double DENDRO_113 = DENDRO_110 + DENDRO_112;
            double DENDRO_114 = -DENDRO_109 + DENDRO_113;
            double DENDRO_115 = DENDRO_108 * DENDRO_53;
            double DENDRO_116 = 4 * DENDRO_17 * DENDRO_93;
            double DENDRO_117 = 0.25 * DENDRO_108;
            double DENDRO_118 = DENDRO_105 * DENDRO_117;
            double DENDRO_119 = DENDRO_100 - DENDRO_95 + DENDRO_98;
            double DENDRO_120 = DENDRO_114 * DENDRO_119;
            double DENDRO_121 = 2.0 * DENDRO_24 * DENDRO_93;
            double DENDRO_122 = 4 * DENDRO_24 * DENDRO_93;
            double DENDRO_123 = 0.25 * DENDRO_42;
            double DENDRO_124 = -DENDRO_123;
            double DENDRO_125 = DENDRO_124 + 0.5 * DENDRO_40;
            double DENDRO_126 = DENDRO_105 * DENDRO_125;
            double DENDRO_127 = 4 * DENDRO_16 * DENDRO_93;
            double DENDRO_128 = DENDRO_114 * DENDRO_125;
            double DENDRO_129 = 0.5 * DENDRO_108;
            double DENDRO_130 = grad_1_gt1[pp];
            double DENDRO_131 = 1.0 * DENDRO_130;
            double DENDRO_132 = -DENDRO_131;
            double DENDRO_133 = DENDRO_129 + DENDRO_132;
            double DENDRO_134 =
                -DENDRO_105 * DENDRO_94 * (DENDRO_101 + DENDRO_97 + DENDRO_99) -
                DENDRO_106 * DENDRO_27 * DENDRO_53 * (DENDRO_41 + DENDRO_43) +
                DENDRO_107 * (DENDRO_114 * DENDRO_42 + DENDRO_115) +
                DENDRO_116 * (DENDRO_118 + 0.5 * DENDRO_120) -
                DENDRO_121 * (DENDRO_102 * DENDRO_53 + DENDRO_105 * DENDRO_42) -
                DENDRO_122 * (DENDRO_119 * DENDRO_53 + DENDRO_126) +
                DENDRO_127 * (DENDRO_128 - 2 * DENDRO_133 * DENDRO_53) +
                4 * DENDRO_80 * gt0[pp] + DENDRO_81 * DENDRO_82 +
                DENDRO_83 * DENDRO_84 - DENDRO_85 * DENDRO_86 -
                DENDRO_87 * grad2_1_2_gt0[pp] + DENDRO_88 * grad2_0_2_gt0[pp] -
                DENDRO_89 * grad2_0_1_gt0[pp] + DENDRO_90 * grad2_2_2_gt0[pp] +
                DENDRO_91 * grad2_1_1_gt0[pp] + DENDRO_92 * grad2_0_0_gt0[pp];
            double DENDRO_135 = 3.0 * DENDRO_29 * DENDRO_93;
            double DENDRO_136 = DENDRO_103 * DENDRO_24;
            double DENDRO_137 = DENDRO_27 * DENDRO_49;
            double DENDRO_138 = DENDRO_102 * DENDRO_16;
            double DENDRO_139 = DENDRO_136 + DENDRO_137 - DENDRO_138;
            double DENDRO_140 = DENDRO_139 * DENDRO_49;
            double DENDRO_141 = 3.0 * DENDRO_19 * DENDRO_93;
            double DENDRO_142 = DENDRO_108 * DENDRO_16;
            double DENDRO_143 = DENDRO_27 * DENDRO_42;
            double DENDRO_144 = DENDRO_111 * DENDRO_24;
            double DENDRO_145 = -DENDRO_142 + DENDRO_143 + DENDRO_144;
            double DENDRO_146 = DENDRO_145 * DENDRO_42;
            double DENDRO_147 = 4 * DENDRO_114 * DENDRO_19 * DENDRO_93;
            double DENDRO_148 = 6.0 * DENDRO_45 * DENDRO_93;
            double DENDRO_149 = DENDRO_68 + DENDRO_70 - DENDRO_72;
            double DENDRO_150 = 0.25 * DENDRO_103;
            double DENDRO_151 = grad_2_gt2[pp];
            double DENDRO_152 = 1.0 * DENDRO_151;
            double DENDRO_153 = -DENDRO_152;
            double DENDRO_154 = DENDRO_24 * DENDRO_49;
            double DENDRO_155 = DENDRO_103 * DENDRO_29;
            double DENDRO_156 = DENDRO_102 * DENDRO_17;
            double DENDRO_157 = -DENDRO_154 - DENDRO_155 + DENDRO_156;
            double DENDRO_158 = 4 * DENDRO_157 * DENDRO_29 * DENDRO_93;
            double DENDRO_159 = DENDRO_24 * DENDRO_42;
            double DENDRO_160 = DENDRO_108 * DENDRO_17;
            double DENDRO_161 = DENDRO_111 * DENDRO_29;
            double DENDRO_162 = DENDRO_159 - DENDRO_160 + DENDRO_161;
            double DENDRO_163 = 0.75 * DENDRO_98;
            double DENDRO_164 = 0.25 * DENDRO_100;
            double DENDRO_165 = 4 * DENDRO_19 * DENDRO_93 *
                                (DENDRO_163 + DENDRO_164 + DENDRO_97);
            double DENDRO_166 = 4 * DENDRO_27 * DENDRO_93;
            double DENDRO_167 = DENDRO_48 + DENDRO_50;
            double DENDRO_168 = 0.25 * DENDRO_49;
            double DENDRO_169 = DENDRO_145 * DENDRO_168;
            double DENDRO_170 = DENDRO_123 * DENDRO_139;
            double DENDRO_171 = DENDRO_154 + DENDRO_155 - DENDRO_156;
            double DENDRO_172 = DENDRO_103 * DENDRO_66;
            double DENDRO_173 = DENDRO_139 * DENDRO_45;
            double DENDRO_174 = DENDRO_149 * DENDRO_49;
            double DENDRO_175 = DENDRO_145 * DENDRO_45;
            double DENDRO_176 = DENDRO_149 * DENDRO_42;
            double DENDRO_177 = 0.25 * DENDRO_173;
            double DENDRO_178 = 0.25 * DENDRO_175;
            double DENDRO_179 = DENDRO_150 * DENDRO_162;
            double DENDRO_180 = DENDRO_119 * DENDRO_171;
            double DENDRO_181 = DENDRO_111 * DENDRO_66;
            double DENDRO_182 = DENDRO_105 * DENDRO_133;
            double DENDRO_183 = DENDRO_102 * DENDRO_114;
            double DENDRO_184 = 0.25 * DENDRO_183;
            double DENDRO_185 = 0.5 * DENDRO_103;
            double DENDRO_186 = DENDRO_153 + DENDRO_185;
            double DENDRO_187 = -DENDRO_159 + DENDRO_160 - DENDRO_161;
            double DENDRO_188 = DENDRO_186 * DENDRO_187;
            double DENDRO_189 = DENDRO_111 * DENDRO_157;
            double DENDRO_190 = 0.25 * DENDRO_189;
            double DENDRO_191 = -DENDRO_190;
            double DENDRO_192 = -DENDRO_168;
            double DENDRO_193 = DENDRO_192 + 0.5 * DENDRO_47;
            double DENDRO_194 = DENDRO_103 - 2.0 * DENDRO_151;
            double DENDRO_195 = 2 * DENDRO_33;
            double DENDRO_196 = grad2_0_0_chi[pp];
            double DENDRO_197 = -DENDRO_196;
            double DENDRO_198 = DENDRO_14 * DENDRO_35;
            double DENDRO_199 = -DENDRO_63 - DENDRO_64 + DENDRO_65;
            double DENDRO_200 = DENDRO_14 * DENDRO_34;
            double DENDRO_201 = DENDRO_14 * DENDRO_36;
            double DENDRO_202 = -DENDRO_68 - DENDRO_70 + DENDRO_72;
            double DENDRO_203 = DENDRO_105 * DENDRO_24;
            double DENDRO_204 = grad_2_gt3[pp];
            double DENDRO_205 = DENDRO_19 * DENDRO_204;
            double DENDRO_206 = grad_1_gt5[pp];
            double DENDRO_207 = DENDRO_17 * DENDRO_206;
            double DENDRO_208 = DENDRO_119 * DENDRO_16;
            double DENDRO_209 = DENDRO_207 + DENDRO_208;
            double DENDRO_210 = -DENDRO_205 + DENDRO_209;
            double DENDRO_211 = DENDRO_17 * DENDRO_210;
            double DENDRO_212 = DENDRO_114 * DENDRO_16;
            double DENDRO_213 = grad_2_gt5[pp];
            double DENDRO_214 = 0.5 * DENDRO_213;
            double DENDRO_215 = DENDRO_17 * DENDRO_214;
            double DENDRO_216 = 0.5 * DENDRO_206;
            double DENDRO_217 = grad_2_gt4[pp];
            double DENDRO_218 = 1.0 * DENDRO_217;
            double DENDRO_219 = -DENDRO_218;
            double DENDRO_220 = DENDRO_216 + DENDRO_219;
            double DENDRO_221 =
                -DENDRO_16 * DENDRO_186 + DENDRO_19 * DENDRO_220 + DENDRO_215;
            double DENDRO_222 = DENDRO_221 * DENDRO_29;
            double DENDRO_223 = grad_1_gt3[pp];
            double DENDRO_224 = 0.5 * DENDRO_223;
            double DENDRO_225 = DENDRO_19 * DENDRO_224;
            double DENDRO_226 = grad_1_gt4[pp];
            double DENDRO_227 = 1.0 * DENDRO_226;
            double DENDRO_228 = 0.5 * DENDRO_204;
            double DENDRO_229 = DENDRO_227 - DENDRO_228;
            double DENDRO_230 = DENDRO_17 * DENDRO_229;
            double DENDRO_231 = DENDRO_133 * DENDRO_16;
            double DENDRO_232 = -DENDRO_225 + DENDRO_230 - DENDRO_231;
            double DENDRO_233 = DENDRO_19 * DENDRO_232;
            double DENDRO_234 = DENDRO_27 * DENDRO_53;
            double DENDRO_235 =
                2.0 * DENDRO_93 *
                (DENDRO_203 - 1.0 * DENDRO_211 - 1.0 * DENDRO_212 + DENDRO_222 +
                 DENDRO_233 + DENDRO_234);
            double DENDRO_236 = DENDRO_157 * DENDRO_24;
            double DENDRO_237 = DENDRO_17 * DENDRO_204;
            double DENDRO_238 = DENDRO_206 * DENDRO_29;
            double DENDRO_239 = DENDRO_119 * DENDRO_24;
            double DENDRO_240 = DENDRO_237 - DENDRO_238 - DENDRO_239;
            double DENDRO_241 = DENDRO_17 * DENDRO_240;
            double DENDRO_242 = DENDRO_16 * DENDRO_187;
            double DENDRO_243 = DENDRO_214 * DENDRO_29;
            double DENDRO_244 = DENDRO_186 * DENDRO_24;
            double DENDRO_245 = DENDRO_17 * DENDRO_220;
            double DENDRO_246 = -DENDRO_243 + DENDRO_244 - DENDRO_245;
            double DENDRO_247 = DENDRO_246 * DENDRO_29;
            double DENDRO_248 = DENDRO_17 * DENDRO_224;
            double DENDRO_249 =
                DENDRO_133 * DENDRO_24 - DENDRO_229 * DENDRO_29 + DENDRO_248;
            double DENDRO_250 = DENDRO_19 * DENDRO_249;
            double DENDRO_251 = DENDRO_199 * DENDRO_27;
            double DENDRO_252 =
                2.0 * DENDRO_93 *
                (DENDRO_236 - 1.0 * DENDRO_241 - 1.0 * DENDRO_242 + DENDRO_247 +
                 DENDRO_250 + DENDRO_251);
            double DENDRO_253 = -DENDRO_136 - DENDRO_137 + DENDRO_138;
            double DENDRO_254 = DENDRO_24 * DENDRO_253;
            double DENDRO_255 = DENDRO_16 * DENDRO_204;
            double DENDRO_256 = DENDRO_206 * DENDRO_24;
            double DENDRO_257 = DENDRO_119 * DENDRO_27;
            double DENDRO_258 = DENDRO_255 - DENDRO_256 - DENDRO_257;
            double DENDRO_259 = DENDRO_17 * DENDRO_258;
            double DENDRO_260 = DENDRO_142 - DENDRO_143 - DENDRO_144;
            double DENDRO_261 = DENDRO_16 * DENDRO_260;
            double DENDRO_262 = DENDRO_214 * DENDRO_24;
            double DENDRO_263 = DENDRO_16 * DENDRO_220;
            double DENDRO_264 = DENDRO_186 * DENDRO_27;
            double DENDRO_265 = -DENDRO_262 - DENDRO_263 + DENDRO_264;
            double DENDRO_266 = DENDRO_265 * DENDRO_29;
            double DENDRO_267 = DENDRO_16 * DENDRO_224;
            double DENDRO_268 =
                DENDRO_133 * DENDRO_27 - DENDRO_229 * DENDRO_24 + DENDRO_267;
            double DENDRO_269 = DENDRO_19 * DENDRO_268;
            double DENDRO_270 = DENDRO_202 * DENDRO_27;
            double DENDRO_271 =
                2.0 * DENDRO_93 *
                (DENDRO_254 - 1.0 * DENDRO_259 - 1.0 * DENDRO_261 + DENDRO_266 +
                 DENDRO_269 + DENDRO_270);
            double DENDRO_272 = grad2_2_2_chi[pp];
            double DENDRO_273 = 3 * DENDRO_33;
            double DENDRO_274 = pow(DENDRO_35, 2);
            double DENDRO_275 =
                DENDRO_29 * (2 * DENDRO_272 - DENDRO_273 * DENDRO_274);
            double DENDRO_276 = grad2_1_1_chi[pp];
            double DENDRO_277 = pow(DENDRO_34, 2);
            double DENDRO_278 =
                DENDRO_19 * (-DENDRO_273 * DENDRO_277 + 2 * DENDRO_276);
            double DENDRO_279 =
                DENDRO_27 * (2 * DENDRO_196 - DENDRO_273 * DENDRO_86);
            double DENDRO_280 = grad2_1_2_chi[pp];
            double DENDRO_281 = DENDRO_34 * DENDRO_35;
            double DENDRO_282 =
                2 * DENDRO_17 * (-DENDRO_273 * DENDRO_281 + 2 * DENDRO_280);
            double DENDRO_283 = grad2_0_2_chi[pp];
            double DENDRO_284 = 3 * DENDRO_33 * DENDRO_36;
            double DENDRO_285 =
                2 * DENDRO_24 * (2 * DENDRO_283 - DENDRO_284 * DENDRO_35);
            double DENDRO_286 = grad2_0_1_chi[pp];
            double DENDRO_287 =
                2 * DENDRO_16 * (-DENDRO_284 * DENDRO_34 + 2 * DENDRO_286);
            double DENDRO_288 = 2 * DENDRO_14;
            double DENDRO_289 = -1.0 * DENDRO_203 + DENDRO_211 + DENDRO_212 -
                                DENDRO_222 - DENDRO_233 - DENDRO_234;
            double DENDRO_290 = DENDRO_288 * DENDRO_289 * DENDRO_34;
            double DENDRO_291 = -1.0 * DENDRO_236 + DENDRO_241 + DENDRO_242 -
                                DENDRO_247 - DENDRO_250 - DENDRO_251;
            double DENDRO_292 = DENDRO_288 * DENDRO_291 * DENDRO_35;
            double DENDRO_293 = -1.0 * DENDRO_254 + DENDRO_259 + DENDRO_261 -
                                DENDRO_266 - DENDRO_269 - DENDRO_270;
            double DENDRO_294 = DENDRO_288 * DENDRO_293 * DENDRO_36;
            double DENDRO_295 = DENDRO_275 + DENDRO_278 + DENDRO_279 -
                                DENDRO_282 + DENDRO_285 - DENDRO_287 +
                                DENDRO_290 + DENDRO_292 + DENDRO_294;
            double DENDRO_296 = DENDRO_14 * DENDRO_295 * DENDRO_33;
            double DENDRO_297 = grad2_2_2_alpha[pp];
            double DENDRO_298 =
                DENDRO_14 * DENDRO_32 * (DENDRO_221 + DENDRO_39 * gt5[pp]);
            double DENDRO_299 = 4 * DENDRO_14 * DENDRO_67;
            double DENDRO_300 = 0.5 * gt5[pp];
            double DENDRO_301 = DENDRO_24 * DENDRO_35;
            double DENDRO_302 = DENDRO_27 * DENDRO_36;
            double DENDRO_303 = -DENDRO_301 - DENDRO_302 + DENDRO_75;
            double DENDRO_304 = DENDRO_300 * DENDRO_303 * DENDRO_33;
            double DENDRO_305 = DENDRO_14 * DENDRO_243;
            double DENDRO_306 = DENDRO_14 * DENDRO_244;
            double DENDRO_307 = DENDRO_14 * DENDRO_245;
            double DENDRO_308 = -DENDRO_11 + DENDRO_9;
            double DENDRO_309 =
                DENDRO_308 * DENDRO_35 + DENDRO_36 * DENDRO_74 + DENDRO_57;
            double DENDRO_310 =
                DENDRO_33 *
                (DENDRO_14 * DENDRO_300 * DENDRO_309 - 1.0 * DENDRO_35);
            double DENDRO_311 = grad_2_Gt0[pp];
            double DENDRO_312 = 4 * gt4[pp];
            double DENDRO_313 = grad_2_Gt1[pp];
            double DENDRO_314 = grad_2_Gt2[pp];
            double DENDRO_315 = 4 * DENDRO_19 * DENDRO_93;
            double DENDRO_316 = 0.25 * DENDRO_204;
            double DENDRO_317 = -DENDRO_316;
            double DENDRO_318 = 0.75 * DENDRO_95;
            double DENDRO_319 = -DENDRO_164;
            double DENDRO_320 = 2.0 * DENDRO_17 * DENDRO_93;
            double DENDRO_321 = DENDRO_204 * DENDRO_221;
            double DENDRO_322 = DENDRO_105 * DENDRO_316;
            double DENDRO_323 = DENDRO_111 * DENDRO_210;
            double DENDRO_324 = DENDRO_105 * DENDRO_229;
            double DENDRO_325 = DENDRO_102 * DENDRO_210;
            double DENDRO_326 = 0.25 * DENDRO_325;
            double DENDRO_327 = 0.25 * DENDRO_206;
            double DENDRO_328 = -0.5 * DENDRO_217 + DENDRO_327;
            double DENDRO_329 = -DENDRO_210 * DENDRO_328;
            double DENDRO_330 =
                -DENDRO_105 * DENDRO_166 *
                    (DENDRO_318 + DENDRO_319 + DENDRO_99) +
                DENDRO_116 * (2 * DENDRO_221 * DENDRO_229 + DENDRO_329) -
                DENDRO_121 *
                    (DENDRO_102 * DENDRO_221 + DENDRO_105 * DENDRO_206) +
                DENDRO_127 * (DENDRO_322 + 0.5 * DENDRO_323) +
                DENDRO_127 * (DENDRO_324 + DENDRO_326) -
                DENDRO_210 * DENDRO_315 * (DENDRO_227 + DENDRO_317) -
                DENDRO_274 * DENDRO_85 + DENDRO_311 * DENDRO_83 +
                DENDRO_312 * DENDRO_313 + 4 * DENDRO_314 * gt5[pp] +
                DENDRO_320 * (DENDRO_206 * DENDRO_210 + DENDRO_321) -
                DENDRO_87 * grad2_1_2_gt5[pp] + DENDRO_88 * grad2_0_2_gt5[pp] -
                DENDRO_89 * grad2_0_1_gt5[pp] + DENDRO_90 * grad2_2_2_gt5[pp] +
                DENDRO_91 * grad2_1_1_gt5[pp] + DENDRO_92 * grad2_0_0_gt5[pp];
            double DENDRO_331 = DENDRO_206 * DENDRO_240;
            double DENDRO_332 = 3.0 * DENDRO_27 * DENDRO_93;
            double DENDRO_333 = DENDRO_103 * DENDRO_157;
            double DENDRO_334 = 6.0 * DENDRO_213 * DENDRO_93;
            double DENDRO_335 =
                4 * DENDRO_27 * DENDRO_93 * (DENDRO_192 + DENDRO_48);
            double DENDRO_336 = 4 * DENDRO_19 * DENDRO_93 *
                                (DENDRO_163 + DENDRO_319 + DENDRO_96);
            double DENDRO_337 = 4 * DENDRO_265 * DENDRO_29 * DENDRO_93;
            double DENDRO_338 = 4 * DENDRO_221 * DENDRO_29 * DENDRO_93;
            double DENDRO_339 = DENDRO_157 * DENDRO_327;
            double DENDRO_340 = DENDRO_150 * DENDRO_240;
            double DENDRO_341 = DENDRO_206 * DENDRO_246;
            double DENDRO_342 = DENDRO_103 * DENDRO_246;
            double DENDRO_343 = DENDRO_265 * DENDRO_49;
            double DENDRO_344 = 0.25 * DENDRO_213;
            double DENDRO_345 = DENDRO_240 * DENDRO_344;
            double DENDRO_346 = DENDRO_157 * DENDRO_344;
            double DENDRO_347 = DENDRO_168 * DENDRO_258;
            double DENDRO_348 = DENDRO_111 * DENDRO_253;
            double DENDRO_349 = DENDRO_119 * DENDRO_265;
            double DENDRO_350 = DENDRO_119 * DENDRO_253;
            double DENDRO_351 = 0.25 * DENDRO_350;
            double DENDRO_352 = DENDRO_105 * DENDRO_328;
            double DENDRO_353 = -DENDRO_352;
            double DENDRO_354 = DENDRO_111 * DENDRO_221;
            double DENDRO_355 = DENDRO_150 - 0.5 * DENDRO_151;
            double DENDRO_356 = DENDRO_258 * DENDRO_355;
            double DENDRO_357 = -DENDRO_356;
            double DENDRO_358 = DENDRO_111 * DENDRO_265;
            double DENDRO_359 = -DENDRO_253 * DENDRO_355;
            double DENDRO_360 = 2.0 * DENDRO_47 - 1.0 * DENDRO_49;
            double DENDRO_361 = -DENDRO_272;
            double DENDRO_362 = -DENDRO_216;
            double DENDRO_363 = DENDRO_218 + DENDRO_362;
            double DENDRO_364 = -DENDRO_185;
            double DENDRO_365 = DENDRO_152 + DENDRO_364;
            double DENDRO_366 = DENDRO_10 - DENDRO_18;
            double DENDRO_367 = 2.0 * DENDRO_289 * DENDRO_93;
            double DENDRO_368 = 2.0 * DENDRO_291 * DENDRO_93;
            double DENDRO_369 = 2.0 * DENDRO_293 * DENDRO_93;
            double DENDRO_370 = -DENDRO_275 - DENDRO_278 - DENDRO_279 +
                                DENDRO_282 - DENDRO_285 + DENDRO_287 -
                                DENDRO_290 - DENDRO_292 - DENDRO_294;
            double DENDRO_371 = DENDRO_14 * DENDRO_33 * DENDRO_370;
            double DENDRO_372 = grad2_1_1_alpha[pp];
            double DENDRO_373 = 4 * DENDRO_14 * DENDRO_55;
            double DENDRO_374 = 0.5 * DENDRO_33 * gt3[pp];
            double DENDRO_375 = -1.0 * DENDRO_34;
            double DENDRO_376 = 0.5 * DENDRO_14 * gt3[pp];
            double DENDRO_377 = DENDRO_34 * DENDRO_366 + DENDRO_37;
            double DENDRO_378 = -DENDRO_14 * DENDRO_225 +
                                DENDRO_14 * DENDRO_230 - DENDRO_14 * DENDRO_231;
            double DENDRO_379 = grad_1_Gt0[pp];
            double DENDRO_380 = grad_1_Gt1[pp];
            double DENDRO_381 = grad_1_Gt2[pp];
            double DENDRO_382 = DENDRO_204 * DENDRO_210;
            double DENDRO_383 = DENDRO_108 * DENDRO_114;
            double DENDRO_384 = DENDRO_114 * DENDRO_316;
            double DENDRO_385 = DENDRO_117 * DENDRO_210;
            double DENDRO_386 =
                -DENDRO_106 * DENDRO_19 * DENDRO_249 *
                    (DENDRO_227 + DENDRO_228) -
                DENDRO_122 * (DENDRO_114 * DENDRO_228 + DENDRO_385) -
                DENDRO_122 * (DENDRO_129 * DENDRO_210 + DENDRO_384) -
                DENDRO_135 * DENDRO_382 - DENDRO_277 * DENDRO_85 +
                DENDRO_312 * DENDRO_381 - DENDRO_332 * DENDRO_383 +
                DENDRO_379 * DENDRO_81 + 4 * DENDRO_380 * gt3[pp] -
                DENDRO_87 * grad2_1_2_gt3[pp] + DENDRO_88 * grad2_0_2_gt3[pp] -
                DENDRO_89 * grad2_0_1_gt3[pp] + DENDRO_90 * grad2_2_2_gt3[pp] +
                DENDRO_91 * grad2_1_1_gt3[pp] + DENDRO_92 * grad2_0_0_gt3[pp];
            double DENDRO_387 = 6.0 * DENDRO_223;
            double DENDRO_388 = 4 * DENDRO_240 * DENDRO_29 * DENDRO_93;
            double DENDRO_389 =
                4 * DENDRO_27 * DENDRO_93 * (DENDRO_124 + DENDRO_41);
            double DENDRO_390 = -DENDRO_99;
            double DENDRO_391 = 4 * DENDRO_29 * DENDRO_93 *
                                (DENDRO_101 + DENDRO_390 + DENDRO_96);
            double DENDRO_392 = 4 * DENDRO_27 * DENDRO_93 *
                                (DENDRO_164 + DENDRO_318 + DENDRO_390);
            double DENDRO_393 = 4 * DENDRO_19 * DENDRO_268 * DENDRO_93;
            double DENDRO_394 = DENDRO_210 * DENDRO_223;
            double DENDRO_395 = DENDRO_204 * DENDRO_232;
            double DENDRO_396 = DENDRO_114 * DENDRO_223;
            double DENDRO_397 = DENDRO_108 * DENDRO_232;
            double DENDRO_398 = DENDRO_206 * DENDRO_249;
            double DENDRO_399 = DENDRO_268 * DENDRO_42;
            double DENDRO_400 = 0.25 * DENDRO_394;
            double DENDRO_401 = 0.25 * DENDRO_396;
            double DENDRO_402 = DENDRO_187 * DENDRO_327;
            double DENDRO_403 = 0.5 * DENDRO_95;
            double DENDRO_404 = 0.5 * DENDRO_100;
            double DENDRO_405 = DENDRO_403 + DENDRO_404 - 0.5 * DENDRO_98;
            double DENDRO_406 = DENDRO_123 * DENDRO_258;
            double DENDRO_407 = DENDRO_119 * DENDRO_268;
            double DENDRO_408 = DENDRO_111 * DENDRO_249;
            double DENDRO_409 = DENDRO_119 * DENDRO_260;
            double DENDRO_410 = 0.25 * DENDRO_409;
            double DENDRO_411 = DENDRO_187 * DENDRO_220;
            double DENDRO_412 = DENDRO_111 * DENDRO_240;
            double DENDRO_413 = 0.25 * DENDRO_412;
            double DENDRO_414 = DENDRO_117 - 0.5 * DENDRO_130;
            double DENDRO_415 = DENDRO_258 * DENDRO_414;
            double DENDRO_416 = -DENDRO_415;
            double DENDRO_417 = DENDRO_102 * DENDRO_268;
            double DENDRO_418 = 0.5 * DENDRO_226 + DENDRO_317;
            double DENDRO_419 = DENDRO_240 * DENDRO_418;
            double DENDRO_420 = 2 * DENDRO_220 * DENDRO_249;
            double DENDRO_421 = -DENDRO_260 * DENDRO_414;
            double DENDRO_422 = 2 * DENDRO_268 * DENDRO_44;
            double DENDRO_423 = DENDRO_187 * DENDRO_418;
            double DENDRO_424 = DENDRO_102 * DENDRO_249;
            double DENDRO_425 = -DENDRO_276;
            double DENDRO_426 = -DENDRO_129;
            double DENDRO_427 = DENDRO_131 + DENDRO_426;
            double DENDRO_428 = DENDRO_253 * DENDRO_49;
            double DENDRO_429 = DENDRO_260 * DENDRO_42;
            double DENDRO_430 = DENDRO_168 * DENDRO_260;
            double DENDRO_431 = DENDRO_123 * DENDRO_253;
            double DENDRO_432 = DENDRO_103 * DENDRO_199;
            double DENDRO_433 = DENDRO_253 * DENDRO_45;
            double DENDRO_434 = DENDRO_202 * DENDRO_49;
            double DENDRO_435 = DENDRO_260 * DENDRO_45;
            double DENDRO_436 = DENDRO_202 * DENDRO_42;
            double DENDRO_437 = 0.25 * DENDRO_433;
            double DENDRO_438 = 0.25 * DENDRO_435;
            double DENDRO_439 = DENDRO_150 * DENDRO_187;
            double DENDRO_440 = DENDRO_119 * DENDRO_157;
            double DENDRO_441 = DENDRO_111 * DENDRO_199;
            double DENDRO_442 = DENDRO_157 * DENDRO_193;
            double DENDRO_443 = DENDRO_187 * DENDRO_193;
            double DENDRO_444 = grad2_0_2_alpha[pp];
            double DENDRO_445 = DENDRO_14 * DENDRO_32 *
                                (DENDRO_105 + DENDRO_33 * DENDRO_38 * gt2[pp]);
            double DENDRO_446 = 2.0 * DENDRO_55;
            double DENDRO_447 = DENDRO_14 * DENDRO_154;
            double DENDRO_448 = DENDRO_14 * DENDRO_155;
            double DENDRO_449 = DENDRO_14 * DENDRO_156;
            double DENDRO_450 = -DENDRO_36;
            double DENDRO_451 = DENDRO_14 * gt2[pp];
            double DENDRO_452 =
                DENDRO_33 * (DENDRO_309 * DENDRO_451 + DENDRO_450);
            double DENDRO_453 = 2.0 * DENDRO_67;
            double DENDRO_454 = DENDRO_136 * DENDRO_14;
            double DENDRO_455 = DENDRO_137 * DENDRO_14;
            double DENDRO_456 = DENDRO_138 * DENDRO_14;
            double DENDRO_457 = -DENDRO_35;
            double DENDRO_458 =
                DENDRO_33 * (DENDRO_451 * DENDRO_77 + DENDRO_457);
            double DENDRO_459 =
                -4 * DENDRO_444 + 2.0 * DENDRO_445 +
                DENDRO_446 *
                    (-DENDRO_447 - DENDRO_448 + DENDRO_449 + DENDRO_452) +
                DENDRO_453 *
                    (-DENDRO_454 - DENDRO_455 + DENDRO_456 + DENDRO_458);
            double DENDRO_460 = -DENDRO_283;
            double DENDRO_461 = 0.5 * DENDRO_14 * DENDRO_35;
            double DENDRO_462 = 0.5 * DENDRO_14 * DENDRO_34;
            double DENDRO_463 = 0.5 * DENDRO_14 * DENDRO_36;
            double DENDRO_464 = DENDRO_33 * DENDRO_370;
            double DENDRO_465 =
                DENDRO_151 * DENDRO_368 -
                DENDRO_195 *
                    (DENDRO_460 +
                     DENDRO_461 * (DENDRO_103 * DENDRO_308 + DENDRO_156 +
                                   DENDRO_49 * DENDRO_74) +
                     DENDRO_462 * (DENDRO_102 * DENDRO_366 + DENDRO_104) +
                     DENDRO_463 * (DENDRO_103 * DENDRO_74 + DENDRO_138 +
                                   DENDRO_49 * DENDRO_76)) +
                DENDRO_367 * DENDRO_98 + DENDRO_369 * DENDRO_47 +
                DENDRO_451 * DENDRO_464;
            double DENDRO_466 = 2.0 * gt0[pp];
            double DENDRO_467 = DENDRO_311 * DENDRO_466;
            double DENDRO_468 = 2.0 * gt1[pp];
            double DENDRO_469 = DENDRO_313 * DENDRO_468;
            double DENDRO_470 = 2.0 * gt2[pp];
            double DENDRO_471 = DENDRO_470 * DENDRO_80;
            double DENDRO_472 = DENDRO_314 * DENDRO_470;
            double DENDRO_473 = 2.0 * gt4[pp];
            double DENDRO_474 = DENDRO_473 * DENDRO_82;
            double DENDRO_475 = 2.0 * gt5[pp];
            double DENDRO_476 = DENDRO_475 * DENDRO_84;
            double DENDRO_477 = DENDRO_36 * DENDRO_85;
            double DENDRO_478 = -DENDRO_35 * DENDRO_477;
            double DENDRO_479 = -DENDRO_87 * grad2_1_2_gt2[pp];
            double DENDRO_480 = DENDRO_88 * grad2_0_2_gt2[pp];
            double DENDRO_481 = -DENDRO_89 * grad2_0_1_gt2[pp];
            double DENDRO_482 = DENDRO_90 * grad2_2_2_gt2[pp];
            double DENDRO_483 = DENDRO_91 * grad2_1_1_gt2[pp];
            double DENDRO_484 = DENDRO_92 * grad2_0_0_gt2[pp];
            double DENDRO_485 = DENDRO_150 * DENDRO_253;
            double DENDRO_486 = DENDRO_111 * DENDRO_260;
            double DENDRO_487 = 0.25 * DENDRO_486;
            double DENDRO_488 = DENDRO_202 * DENDRO_51;
            double DENDRO_489 = DENDRO_16 * DENDRO_93;
            double DENDRO_490 = DENDRO_111 * DENDRO_114;
            double DENDRO_491 = DENDRO_105 * DENDRO_108;
            double DENDRO_492 = DENDRO_210 * DENDRO_42 + DENDRO_491;
            double DENDRO_493 = DENDRO_103 * DENDRO_187;
            double DENDRO_494 = DENDRO_240 * DENDRO_49 + DENDRO_493;
            double DENDRO_495 = DENDRO_265 * DENDRO_43;
            double DENDRO_496 = DENDRO_150 * DENDRO_260 + DENDRO_495;
            double DENDRO_497 = 0.25 * DENDRO_45;
            double DENDRO_498 = DENDRO_258 * DENDRO_497;
            double DENDRO_499 = DENDRO_193 * DENDRO_260;
            double DENDRO_500 = DENDRO_431 + DENDRO_499;
            double DENDRO_501 = DENDRO_119 * DENDRO_210;
            double DENDRO_502 = 0.25 * DENDRO_501;
            double DENDRO_503 = DENDRO_129 * DENDRO_221;
            double DENDRO_504 = DENDRO_114 * DENDRO_327 + DENDRO_503;
            double DENDRO_505 = -DENDRO_210 * DENDRO_414;
            double DENDRO_506 = DENDRO_114 * DENDRO_418;
            double DENDRO_507 = DENDRO_385 + DENDRO_506;
            double DENDRO_508 = DENDRO_265 * DENDRO_46;
            double DENDRO_509 = DENDRO_168 * DENDRO_253 + DENDRO_508;
            double DENDRO_510 = DENDRO_193 * DENDRO_253;
            double DENDRO_511 = 0.5 * DENDRO_98;
            double DENDRO_512 = DENDRO_404 + DENDRO_511 - 0.5 * DENDRO_95;
            double DENDRO_513 = DENDRO_246 * DENDRO_512 + DENDRO_339;
            double DENDRO_514 = 0.25 * DENDRO_348 + DENDRO_495;
            double DENDRO_515 = 1.0 * DENDRO_19 * DENDRO_93;
            double DENDRO_516 = DENDRO_187 * DENDRO_206;
            double DENDRO_517 = DENDRO_119 * DENDRO_240;
            double DENDRO_518 = -0.5 * DENDRO_100 + DENDRO_403 + DENDRO_511;
            double DENDRO_519 = DENDRO_202 * DENDRO_518 + DENDRO_431;
            double DENDRO_520 = DENDRO_157 * DENDRO_168;
            double DENDRO_521 = DENDRO_186 * DENDRO_246;
            double DENDRO_522 = -DENDRO_521;
            double DENDRO_523 = -0.25 * DENDRO_102 * DENDRO_19 +
                                0.25 * DENDRO_103 * DENDRO_17 +
                                0.25 * DENDRO_16 * DENDRO_49;
            double DENDRO_524 = DENDRO_102 * DENDRO_523;
            double DENDRO_525 =
                DENDRO_111 * DENDRO_523 + DENDRO_221 * DENDRO_43;
            double DENDRO_526 =
                DENDRO_105 * DENDRO_123 + DENDRO_518 * DENDRO_53;
            double DENDRO_527 = DENDRO_240 * DENDRO_355;
            double DENDRO_528 = -DENDRO_527;
            double DENDRO_529 = DENDRO_187 * DENDRO_344;
            double DENDRO_530 = DENDRO_246 * DENDRO_518 + DENDRO_529;
            double DENDRO_531 = DENDRO_199 * DENDRO_216;
            double DENDRO_532 = DENDRO_193 * DENDRO_240;
            double DENDRO_533 = 0.25 * DENDRO_440;
            double DENDRO_534 = DENDRO_125 * DENDRO_210;
            double DENDRO_535 = -0.5 * DENDRO_182 + DENDRO_229 * DENDRO_53;
            double DENDRO_536 =
                DENDRO_119 * DENDRO_523 + DENDRO_216 * DENDRO_53;
            double DENDRO_537 = DENDRO_199 * DENDRO_214;
            double DENDRO_538 = -DENDRO_157 * DENDRO_355;
            double DENDRO_539 = 0.5 * DENDRO_324;
            double DENDRO_540 = DENDRO_133 * DENDRO_221;
            double DENDRO_541 = DENDRO_539 - DENDRO_540;
            double DENDRO_542 =
                DENDRO_105 * DENDRO_327 + DENDRO_221 * DENDRO_512;
            double DENDRO_543 = 0.25 * DENDRO_517;
            double DENDRO_544 = DENDRO_157 * DENDRO_214;
            double DENDRO_545 = DENDRO_17 * DENDRO_93;
            double DENDRO_546 = DENDRO_105 * DENDRO_204;
            double DENDRO_547 = DENDRO_114 * DENDRO_206 + DENDRO_546;
            double DENDRO_548 = DENDRO_545 * (DENDRO_501 + DENDRO_547);
            double DENDRO_549 = DENDRO_258 * DENDRO_49;
            double DENDRO_550 = DENDRO_103 * DENDRO_260 + DENDRO_549;
            double DENDRO_551 = DENDRO_168 * DENDRO_240;
            double DENDRO_552 = DENDRO_439 + DENDRO_531;
            double DENDRO_553 = DENDRO_339 + DENDRO_528;
            double DENDRO_554 = DENDRO_384 + DENDRO_505;
            double DENDRO_555 = -DENDRO_315 * (DENDRO_506 + DENDRO_554);
            double DENDRO_556 = DENDRO_150 * DENDRO_157 + DENDRO_537;
            double DENDRO_557 = DENDRO_123 * DENDRO_210;
            double DENDRO_558 = DENDRO_228 * DENDRO_53;
            double DENDRO_559 = 0.25 * DENDRO_490;
            double DENDRO_560 =
                DENDRO_127 * (DENDRO_557 + DENDRO_558 + DENDRO_559);
            double DENDRO_561 = DENDRO_258 * DENDRO_42;
            double DENDRO_562 = -DENDRO_94 * (DENDRO_353 + DENDRO_542);
            double DENDRO_563 = -DENDRO_122 * (DENDRO_524 + DENDRO_536);
            double DENDRO_564 = -DENDRO_260 * DENDRO_355;
            double DENDRO_565 = DENDRO_202 * DENDRO_512 + DENDRO_498;
            double DENDRO_566 = DENDRO_114 * DENDRO_328;
            double DENDRO_567 = -DENDRO_566;
            double DENDRO_568 = DENDRO_127 * (DENDRO_184 + DENDRO_535);
            double DENDRO_569 =
                -DENDRO_122 * (-DENDRO_220 * DENDRO_53 + DENDRO_525);
            double DENDRO_570 =
                -DENDRO_166 * (DENDRO_405 * DENDRO_53 + DENDRO_526);
            double DENDRO_571 = grad2_1_2_alpha[pp];
            double DENDRO_572 = DENDRO_14 * DENDRO_67;
            double DENDRO_573 = 2.0 * DENDRO_32;
            double DENDRO_574 = DENDRO_14 * gt4[pp];
            double DENDRO_575 = -DENDRO_14 * DENDRO_205 +
                                DENDRO_14 * DENDRO_207 + DENDRO_14 * DENDRO_208;
            double DENDRO_576 = DENDRO_14 * DENDRO_237;
            double DENDRO_577 = DENDRO_14 * DENDRO_238;
            double DENDRO_578 = DENDRO_14 * DENDRO_239;
            double DENDRO_579 = -DENDRO_34;
            double DENDRO_580 =
                DENDRO_33 * (DENDRO_309 * DENDRO_574 + DENDRO_579);
            double DENDRO_581 =
                DENDRO_446 *
                    (DENDRO_576 - DENDRO_577 - DENDRO_578 + DENDRO_580) -
                4 * DENDRO_571 +
                DENDRO_572 * (-2.0 * DENDRO_119 * DENDRO_27 +
                              2.0 * DENDRO_16 * DENDRO_204 -
                              2.0 * DENDRO_206 * DENDRO_24 +
                              2.0 * DENDRO_303 * DENDRO_33 * gt4[pp]) +
                DENDRO_573 *
                    (DENDRO_33 * (DENDRO_377 * DENDRO_574 + DENDRO_457) +
                     DENDRO_575);
            double DENDRO_582 = -DENDRO_280;
            double DENDRO_583 =
                -DENDRO_195 *
                    (DENDRO_461 * (DENDRO_119 * DENDRO_74 +
                                   DENDRO_206 * DENDRO_308 + DENDRO_237) +
                     DENDRO_462 * (DENDRO_204 * DENDRO_366 + DENDRO_209) +
                     DENDRO_463 * (DENDRO_119 * DENDRO_76 +
                                   DENDRO_206 * DENDRO_74 + DENDRO_255) +
                     DENDRO_582) +
                DENDRO_217 * DENDRO_368 + DENDRO_226 * DENDRO_367 +
                DENDRO_369 * DENDRO_95 + DENDRO_464 * DENDRO_574;
            double DENDRO_584 = DENDRO_311 * DENDRO_468;
            double DENDRO_585 = DENDRO_379 * DENDRO_470;
            double DENDRO_586 = 2.0 * gt3[pp];
            double DENDRO_587 = DENDRO_313 * DENDRO_586;
            double DENDRO_588 = DENDRO_380 * DENDRO_473;
            double DENDRO_589 = DENDRO_314 * DENDRO_473;
            double DENDRO_590 = DENDRO_381 * DENDRO_475;
            double DENDRO_591 = -DENDRO_281 * DENDRO_85;
            double DENDRO_592 = -DENDRO_87 * grad2_1_2_gt4[pp];
            double DENDRO_593 = DENDRO_88 * grad2_0_2_gt4[pp];
            double DENDRO_594 = -DENDRO_89 * grad2_0_1_gt4[pp];
            double DENDRO_595 = DENDRO_90 * grad2_2_2_gt4[pp];
            double DENDRO_596 = DENDRO_91 * grad2_1_1_gt4[pp];
            double DENDRO_597 = DENDRO_92 * grad2_0_0_gt4[pp];
            double DENDRO_598 = DENDRO_210 * DENDRO_327;
            double DENDRO_599 = DENDRO_229 * DENDRO_232;
            double DENDRO_600 = DENDRO_223 * DENDRO_523;
            double DENDRO_601 = DENDRO_157 * DENDRO_204 + DENDRO_516;
            double DENDRO_602 = DENDRO_108 * DENDRO_253 + DENDRO_561;
            double DENDRO_603 = DENDRO_210 * DENDRO_418;
            double DENDRO_604 = DENDRO_221 * DENDRO_224;
            double DENDRO_605 = DENDRO_210 * DENDRO_316 + DENDRO_604;
            double DENDRO_606 = 0.25 * DENDRO_323 + DENDRO_503;
            double DENDRO_607 = DENDRO_232 * DENDRO_518;
            double DENDRO_608 = DENDRO_240 * DENDRO_316;
            double DENDRO_609 = DENDRO_246 * DENDRO_405 + DENDRO_340;
            double DENDRO_610 = DENDRO_164 + DENDRO_390 + DENDRO_96;
            double DENDRO_611 = DENDRO_253 * DENDRO_610;
            double DENDRO_612 = 1.0 * DENDRO_27 * DENDRO_93;
            double DENDRO_613 = DENDRO_102 * DENDRO_157;
            double DENDRO_614 = DENDRO_125 * DENDRO_253;
            double DENDRO_615 = DENDRO_220 * DENDRO_246;
            double DENDRO_616 = -DENDRO_615;
            double DENDRO_617 = DENDRO_157 * DENDRO_328;
            double DENDRO_618 = -DENDRO_617;
            double DENDRO_619 = DENDRO_185 * DENDRO_249;
            double DENDRO_620 = DENDRO_157 * DENDRO_418;
            double DENDRO_621 = DENDRO_240 * DENDRO_610;
            double DENDRO_622 = -0.25 * DENDRO_119 * DENDRO_27 +
                                0.25 * DENDRO_16 * DENDRO_204 -
                                0.25 * DENDRO_206 * DENDRO_24;
            double DENDRO_623 = DENDRO_119 * DENDRO_622;
            double DENDRO_624 =
                DENDRO_111 * DENDRO_622 + DENDRO_129 * DENDRO_265;
            double DENDRO_625 = DENDRO_268 * DENDRO_518;
            double DENDRO_626 = DENDRO_117 * DENDRO_258 + DENDRO_625;
            double DENDRO_627 = DENDRO_214 * DENDRO_249;
            double DENDRO_628 = -DENDRO_240 * DENDRO_328;
            double DENDRO_629 = -DENDRO_253 * DENDRO_414;
            double DENDRO_630 = DENDRO_268 * DENDRO_51;
            double DENDRO_631 = DENDRO_125 * DENDRO_258 + DENDRO_630;
            double DENDRO_632 = DENDRO_185 * DENDRO_268;
            double DENDRO_633 = DENDRO_258 * DENDRO_610 + DENDRO_632;
            double DENDRO_634 =
                DENDRO_193 * DENDRO_258 + DENDRO_265 * DENDRO_44;
            double DENDRO_635 =
                DENDRO_150 * DENDRO_258 + DENDRO_265 * DENDRO_405;
            double DENDRO_636 = 1.0 * DENDRO_398;
            double DENDRO_637 = 0.25 * DENDRO_613;
            double DENDRO_638 = DENDRO_214 * DENDRO_240;
            double DENDRO_639 = 1.0 * DENDRO_24 * DENDRO_93;
            double DENDRO_640 = -DENDRO_639 * (DENDRO_325 + DENDRO_547);
            double DENDRO_641 = DENDRO_102 * DENDRO_253;
            double DENDRO_642 = DENDRO_157 * DENDRO_316;
            double DENDRO_643 = DENDRO_402 + DENDRO_619;
            double DENDRO_644 = DENDRO_340 + DENDRO_618;
            double DENDRO_645 =
                -DENDRO_612 * (DENDRO_183 + DENDRO_490 + DENDRO_491);
            double DENDRO_646 = -DENDRO_94 * (DENDRO_221 * DENDRO_228 +
                                              DENDRO_329 + DENDRO_598);
            double DENDRO_647 = DENDRO_240 * DENDRO_327 + DENDRO_627;
            double DENDRO_648 = DENDRO_384 + DENDRO_385;
            double DENDRO_649 = DENDRO_117 * DENDRO_253;
            double DENDRO_650 = DENDRO_268 * DENDRO_50;
            double DENDRO_651 = DENDRO_430 + DENDRO_614;
            double DENDRO_652 = -DENDRO_122 * (DENDRO_567 + DENDRO_606);
            double DENDRO_653 = DENDRO_232 * DENDRO_405 + DENDRO_600;
            double DENDRO_654 = DENDRO_603 + DENDRO_604;
            double DENDRO_655 = DENDRO_186 * DENDRO_268;
            double DENDRO_656 = 0.5 * DENDRO_407;
            double DENDRO_657 = grad2_0_1_alpha[pp];
            double DENDRO_658 = DENDRO_14 * DENDRO_55;
            double DENDRO_659 = DENDRO_14 * gt1[pp];
            double DENDRO_660 = -DENDRO_109 * DENDRO_14 +
                                DENDRO_110 * DENDRO_14 + DENDRO_112 * DENDRO_14;
            double DENDRO_661 = DENDRO_14 * DENDRO_142;
            double DENDRO_662 = DENDRO_14 * DENDRO_143;
            double DENDRO_663 = DENDRO_14 * DENDRO_144;
            double DENDRO_664 =
                DENDRO_33 * (DENDRO_579 + DENDRO_659 * DENDRO_77);
            double DENDRO_665 =
                DENDRO_453 *
                    (DENDRO_661 - DENDRO_662 - DENDRO_663 + DENDRO_664) +
                DENDRO_573 *
                    (DENDRO_33 * (DENDRO_377 * DENDRO_659 + DENDRO_450) +
                     DENDRO_660) -
                4 * DENDRO_657 +
                DENDRO_658 * (2.0 * DENDRO_108 * DENDRO_17 -
                              2.0 * DENDRO_111 * DENDRO_29 -
                              2.0 * DENDRO_24 * DENDRO_42 +
                              2.0 * DENDRO_33 * DENDRO_60 * gt1[pp]);
            double DENDRO_666 = -DENDRO_286;
            double DENDRO_667 =
                DENDRO_100 * DENDRO_368 + DENDRO_130 * DENDRO_367 -
                DENDRO_195 *
                    (DENDRO_461 * (DENDRO_111 * DENDRO_308 + DENDRO_160 +
                                   DENDRO_42 * DENDRO_74) +
                     DENDRO_462 * (DENDRO_108 * DENDRO_366 + DENDRO_113) +
                     DENDRO_463 * (DENDRO_111 * DENDRO_74 + DENDRO_142 +
                                   DENDRO_42 * DENDRO_76) +
                     DENDRO_666) +
                DENDRO_369 * DENDRO_40 + DENDRO_464 * DENDRO_659;
            double DENDRO_668 = DENDRO_379 * DENDRO_466;
            double DENDRO_669 = DENDRO_468 * DENDRO_80;
            double DENDRO_670 = DENDRO_380 * DENDRO_468;
            double DENDRO_671 = DENDRO_381 * DENDRO_470;
            double DENDRO_672 = DENDRO_586 * DENDRO_82;
            double DENDRO_673 = DENDRO_473 * DENDRO_84;
            double DENDRO_674 = -DENDRO_34 * DENDRO_477;
            double DENDRO_675 = -DENDRO_87 * grad2_1_2_gt1[pp];
            double DENDRO_676 = DENDRO_88 * grad2_0_2_gt1[pp];
            double DENDRO_677 = -DENDRO_89 * grad2_0_1_gt1[pp];
            double DENDRO_678 = DENDRO_90 * grad2_2_2_gt1[pp];
            double DENDRO_679 = DENDRO_91 * grad2_1_1_gt1[pp];
            double DENDRO_680 = DENDRO_92 * grad2_0_0_gt1[pp];
            double DENDRO_681 = DENDRO_114 * DENDRO_123;
            double DENDRO_682 = -DENDRO_166 * (1.0 * DENDRO_115 + DENDRO_681);
            double DENDRO_683 =
                -DENDRO_94 * (DENDRO_105 * DENDRO_228 + DENDRO_502);
            double DENDRO_684 = DENDRO_114 * DENDRO_224;
            double DENDRO_685 = DENDRO_133 * DENDRO_232;
            double DENDRO_686 = -DENDRO_685;
            double DENDRO_687 = DENDRO_116 * (DENDRO_554 + DENDRO_600);
            double DENDRO_688 = DENDRO_118 + DENDRO_558;
            double DENDRO_689 = -DENDRO_122 * (DENDRO_557 + DENDRO_688);
            double DENDRO_690 = -DENDRO_114 * DENDRO_414;
            double DENDRO_691 = DENDRO_224 * DENDRO_53;
            double DENDRO_692 = DENDRO_114 * DENDRO_117 + DENDRO_691;
            double DENDRO_693 = DENDRO_127 * (DENDRO_690 + DENDRO_692);
            double DENDRO_694 = DENDRO_232 * DENDRO_512;
            double DENDRO_695 = DENDRO_384 + DENDRO_600;
            double DENDRO_696 = 0.25 * DENDRO_120;
            double DENDRO_697 = -DENDRO_122 * (DENDRO_688 + DENDRO_696);
            double DENDRO_698 = 1.0 * DENDRO_29 * DENDRO_93;
            double DENDRO_699 = DENDRO_117 * DENDRO_260;
            double DENDRO_700 = DENDRO_268 * DENDRO_43;
            double DENDRO_701 = DENDRO_202 * DENDRO_405 + DENDRO_430;
            double DENDRO_702 = DENDRO_202 * DENDRO_44;
            double DENDRO_703 = DENDRO_260 * DENDRO_610;
            double DENDRO_704 = DENDRO_249 * DENDRO_512;
            double DENDRO_705 = DENDRO_187 * DENDRO_316 + DENDRO_704;
            double DENDRO_706 = 0.25 * DENDRO_108 * DENDRO_17 -
                                0.25 * DENDRO_111 * DENDRO_29 -
                                0.25 * DENDRO_24 * DENDRO_42;
            double DENDRO_707 = DENDRO_111 * DENDRO_706;
            double DENDRO_708 =
                DENDRO_119 * DENDRO_706 + DENDRO_199 * DENDRO_228;
            double DENDRO_709 = DENDRO_268 * DENDRO_46;
            double DENDRO_710 = DENDRO_125 * DENDRO_260;
            double DENDRO_711 = DENDRO_249 * DENDRO_50;
            double DENDRO_712 = DENDRO_187 * DENDRO_610 + DENDRO_711;
            double DENDRO_713 = 0.5 * DENDRO_411;
            double DENDRO_714 = DENDRO_186 * DENDRO_249;
            double DENDRO_715 = -DENDRO_713 - DENDRO_714;
            double DENDRO_716 = 0.5 * DENDRO_188;
            double DENDRO_717 = DENDRO_199 * DENDRO_220;
            double DENDRO_718 = -DENDRO_716 - DENDRO_717;
            double DENDRO_719 =
                DENDRO_168 * DENDRO_187 + DENDRO_199 * DENDRO_405;
            double DENDRO_720 = DENDRO_406 + DENDRO_650;
            double DENDRO_721 = DENDRO_123 * DENDRO_260 + DENDRO_709;
            double DENDRO_722 =
                -DENDRO_16 *
                    (DENDRO_665 +
                     alpha[pp] *
                         (DENDRO_107 * (DENDRO_108 * DENDRO_202 + DENDRO_429) +
                          DENDRO_116 * (DENDRO_620 + DENDRO_715) +
                          DENDRO_116 * (DENDRO_694 + DENDRO_695) +
                          DENDRO_116 * (DENDRO_629 + DENDRO_650 + DENDRO_703) -
                          DENDRO_122 * (DENDRO_190 + DENDRO_718) -
                          DENDRO_122 * (DENDRO_431 + DENDRO_701) -
                          DENDRO_122 * (DENDRO_565 + DENDRO_614) -
                          DENDRO_122 * (DENDRO_531 + DENDRO_551 + DENDRO_637) +
                          DENDRO_127 * (DENDRO_707 + DENDRO_708) +
                          DENDRO_127 * (DENDRO_199 * DENDRO_229 + DENDRO_712) +
                          DENDRO_127 * (DENDRO_232 * DENDRO_43 + DENDRO_692) +
                          DENDRO_127 * (-DENDRO_133 * DENDRO_202 + DENDRO_709 +
                                        DENDRO_710) -
                          DENDRO_166 * (0.5 * DENDRO_441 + DENDRO_719) -
                          DENDRO_166 * (DENDRO_202 * DENDRO_43 + DENDRO_438 +
                                        DENDRO_702) -
                          DENDRO_315 * (DENDRO_423 + DENDRO_705) -
                          DENDRO_315 * (DENDRO_684 + DENDRO_686) -
                          DENDRO_315 * (DENDRO_421 + DENDRO_699 + DENDRO_700) +
                          DENDRO_545 * (DENDRO_409 + DENDRO_602) +
                          DENDRO_545 * (DENDRO_517 + DENDRO_601) + DENDRO_667 +
                          DENDRO_668 + DENDRO_669 + DENDRO_670 + DENDRO_671 +
                          DENDRO_672 + DENDRO_673 + DENDRO_674 + DENDRO_675 +
                          DENDRO_676 + DENDRO_677 + DENDRO_678 + DENDRO_679 +
                          DENDRO_680 + DENDRO_682 + DENDRO_683 + DENDRO_687 +
                          DENDRO_689 + DENDRO_693 + DENDRO_697 -
                          DENDRO_698 * (DENDRO_350 + DENDRO_549 + DENDRO_641) -
                          DENDRO_94 * (DENDRO_553 + DENDRO_618))) -
                DENDRO_16 *
                    (DENDRO_665 +
                     alpha[pp] *
                         (DENDRO_107 * (DENDRO_232 * DENDRO_42 + DENDRO_383) +
                          DENDRO_116 * (DENDRO_413 + DENDRO_715) +
                          DENDRO_116 * (DENDRO_505 + DENDRO_653) +
                          DENDRO_116 * (DENDRO_648 + DENDRO_694) +
                          DENDRO_116 * (DENDRO_649 + DENDRO_720) +
                          DENDRO_116 * (DENDRO_703 + DENDRO_720) +
                          DENDRO_116 * (DENDRO_543 + DENDRO_619 + DENDRO_642) -
                          DENDRO_122 * (DENDRO_498 + DENDRO_651) -
                          DENDRO_122 * (DENDRO_498 + DENDRO_701) -
                          DENDRO_122 * (DENDRO_532 + DENDRO_718) -
                          DENDRO_122 * (DENDRO_534 + DENDRO_558 + DENDRO_696) +
                          DENDRO_127 * (DENDRO_707 + DENDRO_712) +
                          DENDRO_127 * (DENDRO_710 + DENDRO_721) +
                          DENDRO_127 * (DENDRO_129 * DENDRO_202 + DENDRO_721) +
                          DENDRO_127 * (DENDRO_249 * DENDRO_51 + DENDRO_708) +
                          DENDRO_127 * (DENDRO_232 * DENDRO_44 + DENDRO_690 +
                                        DENDRO_691) -
                          DENDRO_166 * (DENDRO_443 + DENDRO_719) -
                          DENDRO_166 * (DENDRO_260 * DENDRO_46 + DENDRO_702) -
                          DENDRO_166 * (DENDRO_128 + DENDRO_129 * DENDRO_53 +
                                        DENDRO_681) -
                          DENDRO_315 * (1.0 * DENDRO_399 + DENDRO_699) -
                          DENDRO_315 * (0.5 * DENDRO_408 + DENDRO_705) -
                          DENDRO_315 * (DENDRO_129 * DENDRO_232 + DENDRO_401 +
                                        DENDRO_686) -
                          DENDRO_639 * (DENDRO_183 + DENDRO_492) -
                          DENDRO_639 * (DENDRO_494 + DENDRO_613) + DENDRO_667 +
                          DENDRO_668 + DENDRO_669 + DENDRO_670 + DENDRO_671 +
                          DENDRO_672 + DENDRO_673 + DENDRO_674 + DENDRO_675 +
                          DENDRO_676 + DENDRO_677 + DENDRO_678 + DENDRO_679 +
                          DENDRO_680 -
                          DENDRO_698 * (DENDRO_325 + DENDRO_501 + DENDRO_546) -
                          DENDRO_94 * (DENDRO_528 + DENDRO_644) -
                          DENDRO_94 * (DENDRO_258 * DENDRO_50 + DENDRO_611))) -
                DENDRO_17 *
                    (DENDRO_581 +
                     alpha[pp] *
                         (DENDRO_116 * (DENDRO_603 + DENDRO_605) +
                          DENDRO_116 * (DENDRO_623 + DENDRO_624) +
                          DENDRO_116 * (-DENDRO_133 * DENDRO_265 + DENDRO_633) +
                          DENDRO_116 * (DENDRO_216 * DENDRO_232 + DENDRO_605) +
                          DENDRO_116 * (DENDRO_229 * DENDRO_246 + DENDRO_627 +
                                        DENDRO_628) -
                          DENDRO_122 * (DENDRO_322 + DENDRO_504) -
                          DENDRO_122 * (DENDRO_322 + DENDRO_606) -
                          DENDRO_122 * (DENDRO_339 + DENDRO_609) -
                          DENDRO_122 * (DENDRO_351 + DENDRO_634) -
                          DENDRO_122 * (DENDRO_496 + DENDRO_611) -
                          DENDRO_122 * (DENDRO_530 + DENDRO_618) +
                          DENDRO_127 * (DENDRO_507 + DENDRO_600) +
                          DENDRO_127 * (DENDRO_629 + DENDRO_631) +
                          DENDRO_127 * (DENDRO_385 + DENDRO_600 + DENDRO_607) +
                          DENDRO_127 * (DENDRO_619 + DENDRO_620 + DENDRO_621) -
                          DENDRO_166 * (DENDRO_500 + DENDRO_614) -
                          DENDRO_166 * (DENDRO_105 * DENDRO_129 + DENDRO_559) -
                          DENDRO_315 * (DENDRO_416 + DENDRO_626) -
                          DENDRO_315 * (DENDRO_210 * DENDRO_224 + DENDRO_599) -
                          DENDRO_315 * (DENDRO_216 * DENDRO_249 + DENDRO_419 +
                                        DENDRO_608) +
                          DENDRO_320 * (DENDRO_204 * DENDRO_246 + DENDRO_331) +
                          DENDRO_489 * (DENDRO_412 + DENDRO_601) +
                          DENDRO_489 * (DENDRO_486 + DENDRO_602) + DENDRO_583 +
                          DENDRO_584 + DENDRO_585 + DENDRO_587 + DENDRO_588 +
                          DENDRO_589 + DENDRO_590 + DENDRO_591 + DENDRO_592 +
                          DENDRO_593 + DENDRO_594 + DENDRO_595 + DENDRO_596 +
                          DENDRO_597 -
                          DENDRO_612 * (DENDRO_189 + DENDRO_493 + DENDRO_613) -
                          DENDRO_94 * (1.0 * DENDRO_321 + DENDRO_598) -
                          DENDRO_94 * (0.5 * DENDRO_349 + DENDRO_635) -
                          DENDRO_94 * (DENDRO_216 * DENDRO_246 + DENDRO_345 +
                                       DENDRO_616))) -
                DENDRO_17 *
                    (DENDRO_581 +
                     alpha[pp] *
                         (DENDRO_116 * (DENDRO_623 + DENDRO_633) +
                          DENDRO_116 * (DENDRO_624 - DENDRO_655) +
                          DENDRO_116 * (DENDRO_628 + DENDRO_647) +
                          DENDRO_116 * (-DENDRO_220 * DENDRO_232 + DENDRO_654) +
                          DENDRO_116 * (DENDRO_228 * DENDRO_246 + DENDRO_647) -
                          DENDRO_122 * (DENDRO_529 + DENDRO_609) -
                          DENDRO_122 * (DENDRO_529 + DENDRO_644) -
                          DENDRO_122 * (DENDRO_564 + DENDRO_634) +
                          DENDRO_127 * (DENDRO_410 + DENDRO_631) +
                          DENDRO_127 * (DENDRO_506 + DENDRO_653) +
                          DENDRO_127 * (DENDRO_607 + DENDRO_648) +
                          DENDRO_127 * (DENDRO_621 + DENDRO_643) +
                          DENDRO_127 * (DENDRO_642 + DENDRO_643) +
                          DENDRO_127 * (DENDRO_487 + DENDRO_649 + DENDRO_650) -
                          DENDRO_166 * (DENDRO_499 + DENDRO_651) -
                          DENDRO_166 * (DENDRO_185 * DENDRO_187 + DENDRO_637) -
                          DENDRO_315 * (DENDRO_608 + DENDRO_636) -
                          DENDRO_315 * (DENDRO_626 + DENDRO_656) -
                          DENDRO_315 * (DENDRO_228 * DENDRO_232 + DENDRO_400 +
                                        DENDRO_599) +
                          DENDRO_320 * (DENDRO_206 * DENDRO_232 + DENDRO_382) +
                          DENDRO_583 + DENDRO_584 + DENDRO_585 + DENDRO_587 +
                          DENDRO_588 + DENDRO_589 + DENDRO_590 + DENDRO_591 +
                          DENDRO_592 + DENDRO_593 + DENDRO_594 + DENDRO_595 +
                          DENDRO_596 + DENDRO_597 -
                          DENDRO_639 * (DENDRO_550 + DENDRO_641) + DENDRO_640 +
                          DENDRO_645 + DENDRO_646 + DENDRO_652 -
                          DENDRO_94 * (DENDRO_357 + DENDRO_635) -
                          DENDRO_94 * (DENDRO_616 + DENDRO_638))) +
                DENDRO_19 *
                    (DENDRO_299 * (DENDRO_268 + DENDRO_303 * DENDRO_374) +
                     4 * DENDRO_32 *
                         (DENDRO_33 * (DENDRO_375 + DENDRO_376 * DENDRO_377) +
                          DENDRO_378) -
                     4 * DENDRO_372 +
                     DENDRO_373 * (DENDRO_249 + DENDRO_61 * gt3[pp]) +
                     alpha[pp] *
                         (DENDRO_107 * (DENDRO_396 + DENDRO_397) +
                          DENDRO_107 * (DENDRO_108 * DENDRO_260 + DENDRO_399) +
                          DENDRO_107 * (DENDRO_187 * DENDRO_204 + DENDRO_408) +
                          DENDRO_108 * DENDRO_369 +
                          DENDRO_116 * (1.0 * DENDRO_395 + DENDRO_400) +
                          DENDRO_116 * (DENDRO_416 + DENDRO_417) +
                          DENDRO_116 * (DENDRO_419 - DENDRO_420) -
                          DENDRO_122 * (-1.0 * DENDRO_411 + DENDRO_413) -
                          DENDRO_122 * (DENDRO_240 * DENDRO_405 + DENDRO_402) -
                          DENDRO_122 * (DENDRO_258 * DENDRO_44 + DENDRO_410) -
                          DENDRO_122 * (DENDRO_260 * DENDRO_405 + DENDRO_406) +
                          DENDRO_127 * (1.0 * DENDRO_397 + DENDRO_401) +
                          DENDRO_127 * (DENDRO_421 + DENDRO_422) +
                          DENDRO_127 * (DENDRO_423 + DENDRO_424) -
                          DENDRO_187 * DENDRO_392 -
                          DENDRO_19 * DENDRO_232 * DENDRO_387 * DENDRO_93 -
                          DENDRO_195 * (DENDRO_198 * (DENDRO_229 * DENDRO_308 +
                                                      DENDRO_248 +
                                                      DENDRO_427 * DENDRO_74) +
                                        DENDRO_200 * (DENDRO_16 * DENDRO_427 +
                                                      DENDRO_224 * DENDRO_366 +
                                                      DENDRO_230) +
                                        DENDRO_201 * (DENDRO_229 * DENDRO_74 +
                                                      DENDRO_267 +
                                                      DENDRO_427 * DENDRO_76) +
                                        DENDRO_425) +
                          DENDRO_204 * DENDRO_368 + DENDRO_223 * DENDRO_367 -
                          DENDRO_258 * DENDRO_391 - DENDRO_260 * DENDRO_389 +
                          DENDRO_320 * (DENDRO_394 + DENDRO_395) +
                          DENDRO_320 * (DENDRO_108 * DENDRO_258 + DENDRO_407) +
                          DENDRO_320 * (DENDRO_204 * DENDRO_240 + DENDRO_398) +
                          DENDRO_371 * gt3[pp] + DENDRO_386 -
                          DENDRO_388 * (DENDRO_218 - DENDRO_327) -
                          DENDRO_393 * (DENDRO_129 + DENDRO_131))) +
                DENDRO_24 *
                    (DENDRO_459 +
                     alpha[pp] *
                         (DENDRO_116 * (DENDRO_326 + DENDRO_541) +
                          DENDRO_116 * (DENDRO_340 + DENDRO_513) +
                          DENDRO_116 * (DENDRO_347 + DENDRO_496) +
                          DENDRO_116 * (DENDRO_347 + DENDRO_514) +
                          DENDRO_116 * (DENDRO_502 + DENDRO_504) +
                          DENDRO_116 * (DENDRO_528 + DENDRO_530) -
                          DENDRO_121 * (DENDRO_246 * DENDRO_49 + DENDRO_333) -
                          DENDRO_122 * (DENDRO_509 + DENDRO_510) -
                          DENDRO_122 * (DENDRO_524 + DENDRO_525) -
                          DENDRO_122 * (DENDRO_185 * DENDRO_202 + DENDRO_509) -
                          DENDRO_122 * (DENDRO_221 * DENDRO_44 + DENDRO_536) -
                          DENDRO_122 * (DENDRO_246 * DENDRO_51 + DENDRO_537 +
                                        DENDRO_538) +
                          DENDRO_127 * (DENDRO_498 + DENDRO_500) +
                          DENDRO_127 * (DENDRO_498 + DENDRO_519) +
                          DENDRO_127 * (DENDRO_534 + DENDRO_535) +
                          DENDRO_127 * (DENDRO_531 + DENDRO_532 + DENDRO_533) -
                          DENDRO_166 * (DENDRO_126 + DENDRO_526) -
                          DENDRO_166 * (DENDRO_253 * DENDRO_46 + DENDRO_488) -
                          DENDRO_166 * (DENDRO_185 * DENDRO_199 + DENDRO_442 +
                                        DENDRO_520) -
                          DENDRO_315 * (DENDRO_505 + DENDRO_507) -
                          DENDRO_315 * (DENDRO_258 * DENDRO_43 + DENDRO_487) +
                          DENDRO_465 + DENDRO_467 + DENDRO_469 + DENDRO_471 +
                          DENDRO_472 + DENDRO_474 + DENDRO_476 + DENDRO_478 +
                          DENDRO_479 + DENDRO_480 + DENDRO_481 + DENDRO_482 +
                          DENDRO_483 + DENDRO_484 +
                          DENDRO_489 * (DENDRO_189 + DENDRO_494) +
                          DENDRO_489 * (DENDRO_490 + DENDRO_492) -
                          DENDRO_515 * (DENDRO_412 + DENDRO_516 + DENDRO_517) -
                          DENDRO_94 * (1.0 * DENDRO_343 + DENDRO_485) -
                          DENDRO_94 * (DENDRO_221 * DENDRO_405 + DENDRO_542) -
                          DENDRO_94 * (DENDRO_185 * DENDRO_246 + DENDRO_346 +
                                       DENDRO_522))) +
                DENDRO_24 *
                    (DENDRO_459 +
                     alpha[pp] *
                         (DENDRO_116 * (DENDRO_513 + DENDRO_529) +
                          DENDRO_116 * (DENDRO_514 + DENDRO_564) +
                          DENDRO_116 * (DENDRO_529 + DENDRO_553) +
                          DENDRO_116 * (DENDRO_541 + DENDRO_567) -
                          DENDRO_121 * (DENDRO_103 * DENDRO_202 + DENDRO_428) -
                          DENDRO_122 * (DENDRO_538 + DENDRO_556) -
                          DENDRO_122 * (DENDRO_246 * DENDRO_50 + DENDRO_556) -
                          DENDRO_122 * (-DENDRO_186 * DENDRO_202 + DENDRO_508 +
                                        DENDRO_510) +
                          DENDRO_127 * (DENDRO_430 + DENDRO_519) +
                          DENDRO_127 * (DENDRO_499 + DENDRO_565) +
                          DENDRO_127 * (DENDRO_533 + DENDRO_552) +
                          DENDRO_127 * (DENDRO_551 + DENDRO_552) -
                          DENDRO_166 * (1.0 * DENDRO_432 + DENDRO_520) -
                          DENDRO_166 * (DENDRO_202 * DENDRO_50 + DENDRO_437 +
                                        DENDRO_488) -
                          DENDRO_315 * (DENDRO_187 * DENDRO_216 + DENDRO_543) +
                          DENDRO_465 + DENDRO_467 + DENDRO_469 + DENDRO_471 +
                          DENDRO_472 + DENDRO_474 + DENDRO_476 + DENDRO_478 +
                          DENDRO_479 + DENDRO_480 + DENDRO_481 + DENDRO_482 +
                          DENDRO_483 + DENDRO_484 -
                          DENDRO_515 * (DENDRO_409 + DENDRO_486 + DENDRO_561) +
                          DENDRO_545 * (DENDRO_350 + DENDRO_550) + DENDRO_548 +
                          DENDRO_555 + DENDRO_560 + DENDRO_562 + DENDRO_563 +
                          DENDRO_568 + DENDRO_569 + DENDRO_570 -
                          DENDRO_94 * (DENDRO_522 + DENDRO_544) -
                          DENDRO_94 * (DENDRO_265 * DENDRO_50 + DENDRO_359 +
                                       DENDRO_485))) +
                DENDRO_27 *
                    (-4 * DENDRO_31 + DENDRO_373 * (DENDRO_199 + DENDRO_62) +
                     4 * DENDRO_54 +
                     4 * DENDRO_67 *
                         (-DENDRO_69 - DENDRO_71 + DENDRO_73 + DENDRO_78) +
                     alpha[pp] *
                         (-DENDRO_106 * DENDRO_167 * DENDRO_199 * DENDRO_27 +
                          DENDRO_107 * (DENDRO_435 + DENDRO_436) +
                          DENDRO_107 * (DENDRO_187 * DENDRO_49 + DENDRO_441) +
                          DENDRO_116 * (-1.0 * DENDRO_182 + DENDRO_184) +
                          DENDRO_116 * (-1.0 * DENDRO_188 + DENDRO_190) +
                          DENDRO_116 * (DENDRO_439 + 0.5 * DENDRO_440) +
                          DENDRO_116 * (DENDRO_253 * DENDRO_43 + DENDRO_430) +
                          DENDRO_116 * (DENDRO_260 * DENDRO_50 + DENDRO_431) -
                          DENDRO_121 * (DENDRO_433 + DENDRO_434) -
                          DENDRO_121 * (DENDRO_157 * DENDRO_49 + DENDRO_432) -
                          DENDRO_122 * (1.0 * DENDRO_434 + DENDRO_437) -
                          DENDRO_122 * (-DENDRO_194 * DENDRO_199 + DENDRO_442) +
                          DENDRO_127 * (1.0 * DENDRO_436 + DENDRO_438) +
                          DENDRO_127 * (DENDRO_119 * DENDRO_199 + DENDRO_443) +
                          DENDRO_134 - DENDRO_135 * DENDRO_428 -
                          DENDRO_141 * DENDRO_429 -
                          DENDRO_147 * (-DENDRO_117 + DENDRO_131) -
                          DENDRO_148 * DENDRO_270 -
                          DENDRO_158 * (-DENDRO_150 + DENDRO_152) -
                          DENDRO_165 * DENDRO_187 -
                          DENDRO_195 * (DENDRO_197 +
                                        DENDRO_198 * (DENDRO_308 * DENDRO_51 +
                                                      DENDRO_46 * DENDRO_74 +
                                                      DENDRO_65) +
                                        DENDRO_200 * (DENDRO_366 * DENDRO_44 +
                                                      DENDRO_52) +
                                        DENDRO_201 * (DENDRO_46 * DENDRO_76 +
                                                      DENDRO_51 * DENDRO_74 +
                                                      DENDRO_72)) +
                          DENDRO_367 * DENDRO_42 + DENDRO_368 * DENDRO_49 +
                          DENDRO_369 * DENDRO_45 + DENDRO_371 * gt0[pp])) +
                DENDRO_29 *
                    (-4 * DENDRO_297 + 4 * DENDRO_298 +
                     DENDRO_299 * (DENDRO_265 + DENDRO_304) +
                     4 * DENDRO_55 *
                         (-DENDRO_305 + DENDRO_306 - DENDRO_307 + DENDRO_310) +
                     alpha[pp] *
                         (DENDRO_103 * DENDRO_369 +
                          DENDRO_116 * (1.0 * DENDRO_341 + DENDRO_345) +
                          DENDRO_116 * (DENDRO_357 + DENDRO_358) -
                          DENDRO_121 * (DENDRO_103 * DENDRO_253 + DENDRO_343) -
                          DENDRO_121 * (DENDRO_157 * DENDRO_213 + DENDRO_342) -
                          DENDRO_122 * (1.0 * DENDRO_342 + DENDRO_346) -
                          DENDRO_122 * (DENDRO_353 + DENDRO_354) -
                          DENDRO_122 * (DENDRO_265 * DENDRO_360 + DENDRO_359) +
                          DENDRO_127 * (DENDRO_347 + 0.5 * DENDRO_348) +
                          DENDRO_127 * (DENDRO_157 * DENDRO_216 + DENDRO_340) +
                          DENDRO_127 * (DENDRO_185 * DENDRO_240 + DENDRO_339) +
                          DENDRO_127 * (DENDRO_258 * DENDRO_51 + DENDRO_351) -
                          DENDRO_141 * DENDRO_331 -
                          DENDRO_195 * (DENDRO_198 * (DENDRO_17 * DENDRO_363 +
                                                      DENDRO_214 * DENDRO_308 +
                                                      DENDRO_365 * DENDRO_74) +
                                        DENDRO_200 * (DENDRO_16 * DENDRO_365 +
                                                      DENDRO_215 +
                                                      DENDRO_363 * DENDRO_366) +
                                        DENDRO_201 * (DENDRO_16 * DENDRO_363 +
                                                      DENDRO_214 * DENDRO_74 +
                                                      DENDRO_365 * DENDRO_76) +
                                        DENDRO_361) +
                          DENDRO_206 * DENDRO_367 + DENDRO_213 * DENDRO_368 -
                          DENDRO_247 * DENDRO_334 - DENDRO_253 * DENDRO_335 -
                          DENDRO_258 * DENDRO_336 +
                          DENDRO_320 * (DENDRO_103 * DENDRO_258 + DENDRO_349) +
                          DENDRO_320 * (DENDRO_213 * DENDRO_240 + DENDRO_341) +
                          DENDRO_330 - DENDRO_332 * DENDRO_333 -
                          DENDRO_337 * (DENDRO_152 + DENDRO_185) -
                          DENDRO_338 * (DENDRO_216 + DENDRO_218) +
                          DENDRO_371 * gt5[pp]));
            double DENDRO_723 = DENDRO_14 * DENDRO_722;
            double DENDRO_724 = grad_1_beta0[pp];
            double DENDRO_725 = grad_1_beta2[pp];
            double DENDRO_726 = (1.0L / 3.0L) * At1[pp];
            double DENDRO_727 = (2.0L / 3.0L) * DENDRO_3;
            double DENDRO_728 = At4[pp] * DENDRO_17;
            double DENDRO_729 = -At3[pp] * DENDRO_19 + DENDRO_21 + DENDRO_728;
            double DENDRO_730 = -At1[pp] * DENDRO_27 + At3[pp] * DENDRO_16 -
                                At4[pp] * DENDRO_24;
            double DENDRO_731 = -At1[pp] * DENDRO_24 + At3[pp] * DENDRO_17 -
                                At4[pp] * DENDRO_29;
            double DENDRO_732 = 6.0 * DENDRO_67;
            double DENDRO_733 = 6.0 * DENDRO_32;
            double DENDRO_734 = 1.0 * DENDRO_17 * DENDRO_93;
            double DENDRO_735 = -DENDRO_237 + DENDRO_238 + DENDRO_239;
            double DENDRO_736 = DENDRO_119 * DENDRO_735;
            double DENDRO_737 = -DENDRO_255 + DENDRO_256 + DENDRO_257;
            double DENDRO_738 = DENDRO_119 * DENDRO_145;
            double DENDRO_739 = DENDRO_42 * DENDRO_737 + DENDRO_738;
            double DENDRO_740 = DENDRO_49 * DENDRO_737;
            double DENDRO_741 = DENDRO_119 * DENDRO_139;
            double DENDRO_742 = DENDRO_102 * DENDRO_139;
            double DENDRO_743 = DENDRO_225 - DENDRO_230 + DENDRO_231;
            double DENDRO_744 = DENDRO_145 * DENDRO_414;
            double DENDRO_745 = DENDRO_168 * DENDRO_735;
            double DENDRO_746 = DENDRO_216 * DENDRO_66;
            double DENDRO_747 = DENDRO_171 * DENDRO_610;
            double DENDRO_748 = DENDRO_169 + DENDRO_170;
            double DENDRO_749 = DENDRO_125 * DENDRO_139;
            double DENDRO_750 =
                DENDRO_149 * DENDRO_512 + DENDRO_497 * DENDRO_737;
            double DENDRO_751 = -DENDRO_162 * DENDRO_418;
            double DENDRO_752 = -0.25 * DENDRO_108 * DENDRO_17 +
                                0.25 * DENDRO_111 * DENDRO_29 +
                                0.25 * DENDRO_24 * DENDRO_42;
            double DENDRO_753 = DENDRO_295 * DENDRO_33;
            double DENDRO_754 = grad_2_beta0[pp];
            double DENDRO_755 = grad_2_beta1[pp];
            double DENDRO_756 = (1.0L / 3.0L) * At2[pp];
            double DENDRO_757 = (2.0L / 3.0L) * DENDRO_2;
            double DENDRO_758 =
                At2[pp] * DENDRO_16 - At4[pp] * DENDRO_19 + At5[pp] * DENDRO_17;
            double DENDRO_759 = -At2[pp] * DENDRO_27 + At4[pp] * DENDRO_16 -
                                At5[pp] * DENDRO_24;
            double DENDRO_760 = -At5[pp] * DENDRO_29 + DENDRO_25 + DENDRO_728;
            double DENDRO_761 = 6.0 * DENDRO_55;
            double DENDRO_762 = DENDRO_103 * DENDRO_145 + DENDRO_740;
            double DENDRO_763 = DENDRO_162 * DENDRO_344;
            double DENDRO_764 = -DENDRO_763;
            double DENDRO_765 = DENDRO_171 * DENDRO_327;
            double DENDRO_766 = DENDRO_179 + DENDRO_746;
            double DENDRO_767 = DENDRO_111 * DENDRO_145;
            double DENDRO_768 = DENDRO_150 * DENDRO_171;
            double DENDRO_769 = DENDRO_214 * DENDRO_66;
            double DENDRO_770 = DENDRO_243 - DENDRO_244 + DENDRO_245;
            double DENDRO_771 = DENDRO_139 * DENDRO_355;
            double DENDRO_772 = DENDRO_262 + DENDRO_263 - DENDRO_264;
            double DENDRO_773 = DENDRO_145 * DENDRO_355;
            double DENDRO_774 = DENDRO_111 * DENDRO_139;
            double DENDRO_775 = DENDRO_145 * DENDRO_193;
            double DENDRO_776 = (2.0L / 3.0L) * DENDRO_0;
            double DENDRO_777 = 2 * At4[pp];
            double DENDRO_778 = 2 * At3[pp] * DENDRO_14;
            double DENDRO_779 = 2 * At4[pp] * DENDRO_14;
            double DENDRO_780 = 12 * DENDRO_14 * DENDRO_67;
            double DENDRO_781 = DENDRO_204 * DENDRO_743;
            double DENDRO_782 = DENDRO_108 * DENDRO_743;
            double DENDRO_783 = DENDRO_162 * DENDRO_327;
            double DENDRO_784 = 1.0 * DENDRO_119 * DENDRO_27 -
                                1.0 * DENDRO_16 * DENDRO_204 +
                                1.0 * DENDRO_206 * DENDRO_24;
            double DENDRO_785 = 0.25 * DENDRO_738;
            double DENDRO_786 = (1.0L / 3.0L) * At4[pp];
            double DENDRO_787 = DENDRO_150 * DENDRO_735;
            double DENDRO_788 = DENDRO_619 - DENDRO_783;
            double DENDRO_789 = -DENDRO_327 * DENDRO_735 + DENDRO_627;
            double DENDRO_790 = 0.25 * DENDRO_119 * DENDRO_27 -
                                0.25 * DENDRO_16 * DENDRO_204 +
                                0.25 * DENDRO_206 * DENDRO_24;
            double DENDRO_791 = DENDRO_213 * DENDRO_735;
            double DENDRO_792 = DENDRO_206 * DENDRO_770;
            double DENDRO_793 = DENDRO_171 * DENDRO_213;
            double DENDRO_794 = DENDRO_103 * DENDRO_770;
            // Dendro: printing variables
            //--
            At_rhs00[pp] =
                (4.0L / 3.0L) * At0[pp] * DENDRO_0 - DENDRO_1 * DENDRO_2 -
                DENDRO_1 * DENDRO_3 +
                DENDRO_30 *
                    (-12 * DENDRO_31 + 12 * DENDRO_54 -
                     DENDRO_56 * (-DENDRO_62 + DENDRO_66) -
                     12 * DENDRO_67 *
                         (DENDRO_69 + DENDRO_71 - DENDRO_73 - DENDRO_78) +
                     DENDRO_723 * gt0[pp] +
                     DENDRO_79 *
                         (-DENDRO_107 * (DENDRO_175 + DENDRO_176) -
                          DENDRO_107 * (DENDRO_162 * DENDRO_49 + DENDRO_181) -
                          DENDRO_116 * (DENDRO_179 + 0.5 * DENDRO_180) -
                          DENDRO_116 * (DENDRO_182 - DENDRO_184) -
                          DENDRO_116 * (DENDRO_188 + DENDRO_191) -
                          DENDRO_116 * (DENDRO_139 * DENDRO_43 + DENDRO_169) -
                          DENDRO_116 * (DENDRO_145 * DENDRO_50 + DENDRO_170) +
                          DENDRO_121 * (DENDRO_173 + DENDRO_174) +
                          DENDRO_121 * (DENDRO_171 * DENDRO_49 + DENDRO_172) +
                          DENDRO_122 * (1.0 * DENDRO_174 + DENDRO_177) -
                          DENDRO_122 * (-DENDRO_171 * DENDRO_193 +
                                        DENDRO_194 * DENDRO_66) -
                          DENDRO_127 * (1.0 * DENDRO_176 + DENDRO_178) -
                          DENDRO_127 * (1.0 * DENDRO_119 * DENDRO_66 +
                                        DENDRO_162 * DENDRO_193) +
                          DENDRO_134 + DENDRO_135 * DENDRO_140 +
                          DENDRO_141 * DENDRO_146 +
                          DENDRO_147 * (DENDRO_117 + DENDRO_132) +
                          DENDRO_148 * DENDRO_149 * DENDRO_27 +
                          DENDRO_158 * (DENDRO_150 + DENDRO_153) +
                          DENDRO_162 * DENDRO_165 +
                          DENDRO_166 * DENDRO_167 * DENDRO_66 -
                          DENDRO_195 * (DENDRO_197 + DENDRO_198 * DENDRO_199 +
                                        DENDRO_200 * DENDRO_53 +
                                        DENDRO_201 * DENDRO_202) -
                          DENDRO_235 * DENDRO_42 - DENDRO_252 * DENDRO_49 -
                          DENDRO_271 * DENDRO_45 - DENDRO_296 * gt0[pp])) +
                DENDRO_4 * DENDRO_5 + DENDRO_6 * DENDRO_7 -
                alpha[pp] *
                    (-At0[pp] * K[pp] +
                     DENDRO_15 * (At0[pp] * DENDRO_16 - At1[pp] * DENDRO_19 +
                                  At2[pp] * DENDRO_17) +
                     DENDRO_20 *
                         (-At0[pp] * DENDRO_27 + DENDRO_21 + DENDRO_25) +
                     DENDRO_28 * (-At0[pp] * DENDRO_24 + At1[pp] * DENDRO_17 -
                                  At2[pp] * DENDRO_29)) +
                beta0[pp] * agrad_0_At0[pp] + beta1[pp] * agrad_1_At0[pp] +
                beta2[pp] * agrad_2_At0[pp];
            //--
            At_rhs01[pp] =
                At0[pp] * DENDRO_724 - At1[pp] * DENDRO_727 +
                At2[pp] * DENDRO_725 + At3[pp] * DENDRO_5 + At4[pp] * DENDRO_7 +
                DENDRO_0 * DENDRO_726 + DENDRO_2 * DENDRO_726 +
                DENDRO_30 *
                    (-12 * DENDRO_657 -
                     DENDRO_658 * (-6.0 * DENDRO_108 * DENDRO_17 +
                                   6.0 * DENDRO_111 * DENDRO_29 +
                                   6.0 * DENDRO_24 * DENDRO_42 -
                                   6.0 * DENDRO_33 * DENDRO_60 * gt1[pp]) +
                     DENDRO_659 * DENDRO_722 -
                     DENDRO_732 *
                         (-DENDRO_661 + DENDRO_662 + DENDRO_663 - DENDRO_664) +
                     DENDRO_733 *
                         (DENDRO_33 * (DENDRO_38 * DENDRO_659 + DENDRO_450) +
                          DENDRO_660) +
                     DENDRO_79 *
                         (-DENDRO_100 * DENDRO_252 -
                          DENDRO_107 * (DENDRO_108 * DENDRO_149 + DENDRO_146) +
                          DENDRO_116 * (-DENDRO_512 * DENDRO_743 + DENDRO_695) -
                          DENDRO_116 * (-DENDRO_620 + DENDRO_713 + DENDRO_714) +
                          DENDRO_116 * (DENDRO_139 * DENDRO_414 -
                                        DENDRO_145 * DENDRO_610 + DENDRO_650) +
                          DENDRO_122 * (DENDRO_749 + DENDRO_750) +
                          DENDRO_122 * (DENDRO_149 * DENDRO_405 + DENDRO_748) +
                          DENDRO_122 * (DENDRO_191 + DENDRO_716 + DENDRO_717) +
                          DENDRO_122 * (DENDRO_745 + DENDRO_746 + DENDRO_747) +
                          DENDRO_127 * (-DENDRO_43 * DENDRO_743 + DENDRO_692) +
                          DENDRO_127 * (-DENDRO_102 * DENDRO_752 -
                                        DENDRO_229 * DENDRO_66 + DENDRO_711) -
                          DENDRO_127 * (DENDRO_111 * DENDRO_752 +
                                        DENDRO_119 * DENDRO_752 +
                                        DENDRO_228 * DENDRO_66) +
                          DENDRO_127 * (-DENDRO_125 * DENDRO_145 +
                                        DENDRO_133 * DENDRO_149 + DENDRO_709) -
                          DENDRO_130 * DENDRO_235 +
                          DENDRO_166 * (DENDRO_149 * DENDRO_43 +
                                        DENDRO_149 * DENDRO_44 + DENDRO_178) +
                          DENDRO_166 *
                              (DENDRO_162 * DENDRO_168 + 0.5 * DENDRO_181 +
                               DENDRO_405 * DENDRO_66) -
                          DENDRO_195 * (DENDRO_114 * DENDRO_462 +
                                        DENDRO_187 * DENDRO_461 +
                                        DENDRO_260 * DENDRO_463 + DENDRO_666) -
                          DENDRO_271 * DENDRO_40 +
                          DENDRO_29 * DENDRO_93 *
                              (DENDRO_740 + DENDRO_741 + DENDRO_742) +
                          DENDRO_315 * (-DENDRO_684 + DENDRO_685) -
                          DENDRO_315 * (-DENDRO_117 * DENDRO_145 + DENDRO_700 +
                                        DENDRO_744) -
                          DENDRO_315 * (-DENDRO_162 * DENDRO_316 + DENDRO_704 +
                                        DENDRO_751) -
                          DENDRO_659 * DENDRO_753 + DENDRO_668 + DENDRO_669 +
                          DENDRO_670 + DENDRO_671 + DENDRO_672 + DENDRO_673 +
                          DENDRO_674 + DENDRO_675 + DENDRO_676 + DENDRO_677 +
                          DENDRO_678 + DENDRO_679 + DENDRO_680 + DENDRO_682 +
                          DENDRO_683 + DENDRO_687 + DENDRO_689 + DENDRO_693 +
                          DENDRO_697 -
                          DENDRO_734 * (DENDRO_108 * DENDRO_139 + DENDRO_739) -
                          DENDRO_734 * (DENDRO_162 * DENDRO_206 +
                                        DENDRO_171 * DENDRO_204 + DENDRO_736) +
                          DENDRO_94 *
                              (-DENDRO_339 + DENDRO_527 + DENDRO_617))) -
                alpha[pp] * (-At1[pp] * K[pp] + DENDRO_15 * DENDRO_729 +
                             DENDRO_20 * DENDRO_730 + DENDRO_28 * DENDRO_731) +
                beta0[pp] * agrad_0_At1[pp] + beta1[pp] * agrad_1_At1[pp] +
                beta2[pp] * agrad_2_At1[pp];
            //--
            At_rhs02[pp] =
                At0[pp] * DENDRO_754 + At1[pp] * DENDRO_755 -
                At2[pp] * DENDRO_757 + At4[pp] * DENDRO_5 + At5[pp] * DENDRO_7 +
                DENDRO_0 * DENDRO_756 + DENDRO_3 * DENDRO_756 +
                DENDRO_30 *
                    (-12 * DENDRO_444 + 6.0 * DENDRO_445 +
                     DENDRO_451 * DENDRO_722 -
                     DENDRO_732 *
                         (DENDRO_454 + DENDRO_455 - DENDRO_456 - DENDRO_458) -
                     DENDRO_761 *
                         (DENDRO_447 + DENDRO_448 - DENDRO_449 - DENDRO_452) +
                     DENDRO_79 *
                         (-DENDRO_116 *
                              (-DENDRO_539 + DENDRO_540 + DENDRO_566) +
                          DENDRO_116 * (DENDRO_355 * DENDRO_735 + DENDRO_764 -
                                        DENDRO_765) +
                          DENDRO_116 * (-DENDRO_43 * DENDRO_772 + DENDRO_773 -
                                        0.25 * DENDRO_774) -
                          DENDRO_116 * (DENDRO_512 * DENDRO_770 + DENDRO_763 +
                                        DENDRO_765) +
                          DENDRO_121 * (DENDRO_103 * DENDRO_149 + DENDRO_140) -
                          DENDRO_122 * (-DENDRO_139 * DENDRO_193 +
                                        DENDRO_149 * DENDRO_186 -
                                        DENDRO_46 * DENDRO_772) -
                          DENDRO_122 * (DENDRO_171 * DENDRO_355 - DENDRO_768 -
                                        DENDRO_769) +
                          DENDRO_122 * (DENDRO_50 * DENDRO_770 + DENDRO_768 +
                                        DENDRO_769) -
                          DENDRO_127 * (0.25 * DENDRO_180 + DENDRO_766) -
                          DENDRO_127 * (DENDRO_745 + DENDRO_766) -
                          DENDRO_127 * (DENDRO_750 + DENDRO_775) -
                          DENDRO_127 * (DENDRO_149 * DENDRO_518 + DENDRO_748) -
                          DENDRO_151 * DENDRO_252 +
                          DENDRO_166 *
                              (DENDRO_168 * DENDRO_171 + 1.0 * DENDRO_172) +
                          DENDRO_166 * (DENDRO_149 * DENDRO_50 +
                                        DENDRO_149 * DENDRO_51 + DENDRO_177) +
                          DENDRO_19 * DENDRO_93 * (DENDRO_739 + DENDRO_767) -
                          DENDRO_195 * (DENDRO_105 * DENDRO_462 +
                                        DENDRO_157 * DENDRO_461 +
                                        DENDRO_253 * DENDRO_463 + DENDRO_460) -
                          DENDRO_235 * DENDRO_98 - DENDRO_271 * DENDRO_47 +
                          DENDRO_315 *
                              (DENDRO_162 * DENDRO_216 + 0.25 * DENDRO_736) -
                          DENDRO_451 * DENDRO_753 + DENDRO_467 + DENDRO_469 +
                          DENDRO_471 + DENDRO_472 + DENDRO_474 + DENDRO_476 +
                          DENDRO_478 + DENDRO_479 + DENDRO_480 + DENDRO_481 +
                          DENDRO_482 + DENDRO_483 + DENDRO_484 + DENDRO_548 +
                          DENDRO_555 + DENDRO_560 + DENDRO_562 + DENDRO_563 +
                          DENDRO_568 + DENDRO_569 + DENDRO_570 -
                          DENDRO_734 * (DENDRO_741 + DENDRO_762) +
                          DENDRO_94 * (DENDRO_521 - DENDRO_544) -
                          DENDRO_94 * (-DENDRO_139 * DENDRO_150 -
                                       DENDRO_50 * DENDRO_772 + DENDRO_771))) -
                alpha[pp] * (-At2[pp] * K[pp] + DENDRO_15 * DENDRO_758 +
                             DENDRO_20 * DENDRO_759 + DENDRO_28 * DENDRO_760) +
                beta0[pp] * agrad_0_At2[pp] + beta1[pp] * agrad_1_At2[pp] +
                beta2[pp] * agrad_2_At2[pp];
            //--
            At_rhs11[pp] =
                (4.0L / 3.0L) * At3[pp] * DENDRO_2 - At3[pp] * DENDRO_727 -
                At3[pp] * DENDRO_776 +
                DENDRO_30 *
                    (12 * DENDRO_32 *
                         (DENDRO_33 * (DENDRO_375 + DENDRO_376 * DENDRO_38) +
                          DENDRO_378) -
                     12 * DENDRO_372 +
                     DENDRO_56 *
                         (DENDRO_249 -
                          DENDRO_374 * (-DENDRO_57 + DENDRO_58 + DENDRO_59)) +
                     DENDRO_723 * gt3[pp] +
                     DENDRO_780 *
                         (DENDRO_268 -
                          DENDRO_374 * (DENDRO_301 + DENDRO_302 - DENDRO_75)) +
                     DENDRO_79 *
                         (DENDRO_107 * (DENDRO_396 - DENDRO_782) +
                          DENDRO_107 * (-DENDRO_108 * DENDRO_145 + DENDRO_399) +
                          DENDRO_107 * (-DENDRO_162 * DENDRO_204 + DENDRO_408) -
                          DENDRO_108 * DENDRO_271 +
                          DENDRO_116 * (DENDRO_400 - 1.0 * DENDRO_781) -
                          DENDRO_116 * (DENDRO_415 - 1.0 * DENDRO_417) -
                          DENDRO_116 * (DENDRO_418 * DENDRO_735 + DENDRO_420) +
                          DENDRO_122 * (DENDRO_411 - DENDRO_413) +
                          DENDRO_122 * (DENDRO_123 * DENDRO_737 +
                                        DENDRO_145 * DENDRO_405) +
                          DENDRO_122 * (DENDRO_405 * DENDRO_735 + DENDRO_783) +
                          DENDRO_122 * (DENDRO_44 * DENDRO_784 + DENDRO_785) +
                          DENDRO_127 * (DENDRO_401 - 1.0 * DENDRO_782) +
                          DENDRO_127 * (DENDRO_422 + DENDRO_744) +
                          DENDRO_127 * (DENDRO_424 + DENDRO_751) +
                          DENDRO_145 * DENDRO_389 + DENDRO_162 * DENDRO_392 +
                          DENDRO_19 * DENDRO_387 * DENDRO_743 * DENDRO_93 -
                          DENDRO_195 * (DENDRO_198 * DENDRO_249 +
                                        DENDRO_200 * DENDRO_232 +
                                        DENDRO_201 * DENDRO_268 + DENDRO_425) -
                          DENDRO_204 * DENDRO_252 - DENDRO_223 * DENDRO_235 -
                          DENDRO_296 * gt3[pp] +
                          DENDRO_320 * (DENDRO_394 - DENDRO_781) +
                          DENDRO_320 * (-DENDRO_108 * DENDRO_737 + DENDRO_407) +
                          DENDRO_320 * (-DENDRO_204 * DENDRO_735 + DENDRO_398) +
                          DENDRO_386 + DENDRO_388 * (DENDRO_219 + DENDRO_327) +
                          DENDRO_391 * DENDRO_737 +
                          DENDRO_393 * (DENDRO_132 + DENDRO_426))) +
                DENDRO_4 * DENDRO_724 + DENDRO_725 * DENDRO_777 -
                alpha[pp] *
                    (-At3[pp] * K[pp] + DENDRO_15 * DENDRO_730 +
                     DENDRO_729 * DENDRO_778 + DENDRO_731 * DENDRO_779) +
                beta0[pp] * agrad_0_At3[pp] + beta1[pp] * agrad_1_At3[pp] +
                beta2[pp] * agrad_2_At3[pp];
            //--
            At_rhs12[pp] =
                At1[pp] * DENDRO_754 + At2[pp] * DENDRO_724 +
                At3[pp] * DENDRO_755 - At4[pp] * DENDRO_776 +
                At5[pp] * DENDRO_725 + DENDRO_2 * DENDRO_786 +
                DENDRO_3 * DENDRO_786 +
                DENDRO_30 *
                    (-12 * DENDRO_571 -
                     DENDRO_572 * (6.0 * DENDRO_119 * DENDRO_27 -
                                   6.0 * DENDRO_16 * DENDRO_204 +
                                   6.0 * DENDRO_206 * DENDRO_24 -
                                   6.0 * DENDRO_303 * DENDRO_33 * gt4[pp]) +
                     DENDRO_574 * DENDRO_722 +
                     DENDRO_733 *
                         (DENDRO_33 * (DENDRO_38 * DENDRO_574 + DENDRO_457) +
                          DENDRO_575) -
                     DENDRO_761 *
                         (-DENDRO_576 + DENDRO_577 + DENDRO_578 - DENDRO_580) +
                     DENDRO_79 *
                         (DENDRO_116 * (DENDRO_220 * DENDRO_743 + DENDRO_654) +
                          DENDRO_116 * (-DENDRO_228 * DENDRO_770 + DENDRO_789) +
                          DENDRO_116 * (DENDRO_328 * DENDRO_735 + DENDRO_789) -
                          DENDRO_116 * (DENDRO_111 * DENDRO_790 +
                                        DENDRO_129 * DENDRO_772 + DENDRO_655) +
                          DENDRO_116 * (-DENDRO_119 * DENDRO_790 -
                                        DENDRO_610 * DENDRO_737 + DENDRO_632) -
                          DENDRO_122 * (DENDRO_171 * DENDRO_328 + DENDRO_764 -
                                        DENDRO_787) -
                          DENDRO_122 * (-DENDRO_193 * DENDRO_737 -
                                        DENDRO_44 * DENDRO_772 + DENDRO_773) +
                          DENDRO_122 * (DENDRO_405 * DENDRO_770 + DENDRO_763 +
                                        DENDRO_787) +
                          DENDRO_127 * (-DENDRO_171 * DENDRO_316 + DENDRO_788) +
                          DENDRO_127 * (-DENDRO_518 * DENDRO_743 + DENDRO_648) +
                          DENDRO_127 * (-DENDRO_610 * DENDRO_735 + DENDRO_788) +
                          DENDRO_127 * (-DENDRO_117 * DENDRO_139 + DENDRO_650 -
                                        0.25 * DENDRO_767) +
                          DENDRO_127 * (-DENDRO_125 * DENDRO_737 + DENDRO_630 -
                                        DENDRO_785) +
                          DENDRO_127 * (-DENDRO_405 * DENDRO_743 + DENDRO_506 +
                                        DENDRO_600) +
                          DENDRO_166 * (DENDRO_162 * DENDRO_185 + DENDRO_747) +
                          DENDRO_166 * (DENDRO_169 + DENDRO_749 + DENDRO_775) -
                          DENDRO_195 * (DENDRO_210 * DENDRO_462 +
                                        DENDRO_240 * DENDRO_461 +
                                        DENDRO_258 * DENDRO_463 + DENDRO_582) -
                          DENDRO_217 * DENDRO_252 - DENDRO_226 * DENDRO_235 +
                          DENDRO_24 * DENDRO_93 * (DENDRO_742 + DENDRO_762) -
                          DENDRO_271 * DENDRO_95 -
                          DENDRO_315 * (-DENDRO_316 * DENDRO_735 + DENDRO_636) -
                          DENDRO_315 * (-DENDRO_117 * DENDRO_737 + DENDRO_625 +
                                        DENDRO_656) -
                          DENDRO_315 * (-DENDRO_228 * DENDRO_743 -
                                        DENDRO_229 * DENDRO_743 + DENDRO_400) +
                          DENDRO_320 * (-DENDRO_206 * DENDRO_743 + DENDRO_382) -
                          DENDRO_574 * DENDRO_753 + DENDRO_584 + DENDRO_585 +
                          DENDRO_587 + DENDRO_588 + DENDRO_589 + DENDRO_590 +
                          DENDRO_591 + DENDRO_592 + DENDRO_593 + DENDRO_594 +
                          DENDRO_595 + DENDRO_596 + DENDRO_597 + DENDRO_640 +
                          DENDRO_645 + DENDRO_646 + DENDRO_652 +
                          DENDRO_94 * (DENDRO_615 - DENDRO_638) -
                          DENDRO_94 * (-DENDRO_150 * DENDRO_737 +
                                       DENDRO_355 * DENDRO_737 -
                                       DENDRO_405 * DENDRO_772))) -
                alpha[pp] *
                    (-At4[pp] * K[pp] + DENDRO_15 * DENDRO_759 +
                     DENDRO_758 * DENDRO_778 + DENDRO_760 * DENDRO_779) +
                beta0[pp] * agrad_0_At4[pp] + beta1[pp] * agrad_1_At4[pp] +
                beta2[pp] * agrad_2_At4[pp];
            //--
            At_rhs22[pp] =
                (4.0L / 3.0L) * At5[pp] * DENDRO_3 - At5[pp] * DENDRO_757 -
                At5[pp] * DENDRO_776 +
                DENDRO_30 *
                    (-12 * DENDRO_297 + 12 * DENDRO_298 -
                     12 * DENDRO_55 *
                         (DENDRO_305 - DENDRO_306 + DENDRO_307 - DENDRO_310) +
                     DENDRO_723 * gt5[pp] -
                     DENDRO_780 * (-DENDRO_304 + DENDRO_772) +
                     DENDRO_79 *
                         (DENDRO_103 * DENDRO_171 * DENDRO_332 -
                          DENDRO_103 * DENDRO_271 -
                          DENDRO_116 * (DENDRO_356 - 1.0 * DENDRO_358) -
                          DENDRO_116 * (0.25 * DENDRO_791 + 1.0 * DENDRO_792) +
                          DENDRO_121 * (DENDRO_793 + DENDRO_794) +
                          DENDRO_121 * (DENDRO_103 * DENDRO_139 +
                                        DENDRO_49 * DENDRO_772) +
                          DENDRO_122 * (DENDRO_352 - 1.0 * DENDRO_354) +
                          DENDRO_122 * (0.25 * DENDRO_793 + 1.0 * DENDRO_794) -
                          DENDRO_122 * (-DENDRO_360 * DENDRO_772 + DENDRO_771) -
                          DENDRO_127 *
                              (DENDRO_168 * DENDRO_737 + 0.5 * DENDRO_774) -
                          DENDRO_127 * (DENDRO_171 * DENDRO_216 + DENDRO_787) -
                          DENDRO_127 * (DENDRO_185 * DENDRO_735 + DENDRO_765) -
                          DENDRO_127 *
                              (DENDRO_51 * DENDRO_784 + 0.25 * DENDRO_741) +
                          DENDRO_139 * DENDRO_335 +
                          DENDRO_141 * DENDRO_206 * DENDRO_735 -
                          DENDRO_195 * (DENDRO_198 * DENDRO_246 +
                                        DENDRO_200 * DENDRO_221 +
                                        DENDRO_201 * DENDRO_265 + DENDRO_361) -
                          DENDRO_206 * DENDRO_235 - DENDRO_213 * DENDRO_252 +
                          DENDRO_29 * DENDRO_334 * DENDRO_770 -
                          DENDRO_296 * gt5[pp] -
                          DENDRO_320 * (DENDRO_791 + DENDRO_792) -
                          DENDRO_320 * (DENDRO_103 * DENDRO_737 +
                                        DENDRO_119 * DENDRO_772) +
                          DENDRO_330 + DENDRO_336 * DENDRO_737 +
                          DENDRO_337 * (DENDRO_153 + DENDRO_364) +
                          DENDRO_338 * (DENDRO_219 + DENDRO_362))) +
                DENDRO_6 * DENDRO_754 + DENDRO_755 * DENDRO_777 -
                alpha[pp] *
                    (At5[pp] * DENDRO_288 * DENDRO_760 - At5[pp] * K[pp] +
                     DENDRO_28 * DENDRO_759 + DENDRO_758 * DENDRO_779) +
                beta0[pp] * agrad_0_At5[pp] + beta1[pp] * agrad_1_At5[pp] +
                beta2[pp] * agrad_2_At5[pp];
            // Dendro: reduced ops:  3541
            // Dendro: }}}
            /* debugging */
            /*unsigned int qi = 46 - 1;
            unsigned int qj = 10 - 1;
            unsigned int qk = 60 - 1;
            unsigned int qidx = qi + nx*(qj + ny*qk);
            if (0 && qidx == pp) {
            std::cout << ".... end OPTIMIZED debug stuff..." << std::endl;
            }*/
        }
    }
}
bssn::timer::t_rhs.stop();
