#![allow(unused_variables)]

use std::ops::{Index, IndexMut};

pub struct Matrix {
    value: f64,
}

impl Index<(usize, usize)> for Matrix {
    type Output = f64;
    fn index(&self, _: (usize, usize)) -> &f64 {
        &self.value
    }
}

impl IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, _: (usize, usize)) -> &mut f64 {
        &mut self.value
    }
}

impl Index<usize> for Matrix {
    type Output = f64;
    fn index(&self, _: usize) -> &f64 {
        &self.value
    }
}

impl IndexMut<usize> for Matrix {
    fn index_mut(&mut self, _: usize) -> &mut f64 {
        &mut self.value
    }
}

#[allow(non_snake_case)]
pub fn sympy_curvature_energy_hess_ii_raw(
    mu: f64,
    L_c: f64,
    a1: f64,
    a2: f64,
    a3: f64,
    dt: f64,
    phi_node_i: f64,
    phi_node_j: f64,
    grad_phi_node_i: &Matrix,
    grad_phi_node_j: &Matrix,
    omega_node_i: &Matrix,
    omega_node_j: &Matrix,
    omega_const: &Matrix,
    q_old: &Matrix,
    q_old_node_i: &Matrix,
    q_old_node_j: &Matrix,
    dRdx0_const: &Matrix,
    dRdx1_const: &Matrix,
    dRdx2_const: &Matrix,
) {
    let tmp0 = omega_const[1] + omega_node_i[1] * phi_node_i + omega_node_j[1] * phi_node_j;
    let tmp1 = omega_const[0] + omega_node_i[0] * phi_node_i + omega_node_j[0] * phi_node_j;
    let tmp2 = omega_const[2] + omega_node_i[2] * phi_node_i + omega_node_j[2] * phi_node_j;
    let tmp3 = 0.5 * dt;
    let tmp4 = q_old[0] + tmp3 * (-1.0 * q_old[1] * tmp2 + q_old[2] * tmp0 + q_old[3] * tmp1);
    let tmp5 = q_old[2] + tmp3 * (-1.0 * q_old[0] * tmp0 + q_old[1] * tmp1 + q_old[3] * tmp2);
    let tmp6 = tmp5.powi(2);
    let tmp7 = q_old[1] + tmp3 * (q_old[0] * tmp2 - 1.0 * q_old[2] * tmp1 + q_old[3] * tmp0);
    let tmp8 = tmp7.powi(2);
    let tmp9 = tmp4.powi(2);
    let tmp10 = q_old[3] + tmp3 * (-1.0 * q_old[0] * tmp1 - 1.0 * q_old[1] * tmp0 - 1.0 * q_old[2] * tmp2);
    let tmp11 = tmp10.powi(2);
    let tmp12 = tmp11 + tmp6 + tmp8 + tmp9;
    let tmp13 = tmp12.recip();
    let tmp14 = tmp13 * tmp5;
    let tmp15 = tmp14 * tmp4;
    let tmp16 = tmp13 * tmp7;
    let tmp17 = tmp10 * tmp16;
    let tmp18 = tmp15 + tmp17;
    let tmp19 = q_old_node_i[2]
        + tmp3
            * (omega_node_i[0] * q_old_node_i[1] - 1.0 * omega_node_i[1] * q_old_node_i[0]
                + omega_node_i[2] * q_old_node_i[3]);
    let tmp20 = q_old_node_i[3]
        + tmp3
            * (-1.0 * omega_node_i[0] * q_old_node_i[0]
                - 1.0 * omega_node_i[1] * q_old_node_i[1]
                - 1.0 * omega_node_i[2] * q_old_node_i[2]);
    let tmp21 = tmp19 * tmp20;
    let tmp22 = tmp19.powi(2);
    let tmp23 = q_old_node_i[1]
        + tmp3
            * (-1.0 * omega_node_i[0] * q_old_node_i[2]
                + omega_node_i[1] * q_old_node_i[3]
                + omega_node_i[2] * q_old_node_i[0]);
    let tmp24 = tmp23.powi(2);
    let tmp25 = q_old_node_i[0]
        + tmp3
            * (omega_node_i[0] * q_old_node_i[3] + omega_node_i[1] * q_old_node_i[2]
                - 1.0 * omega_node_i[2] * q_old_node_i[1]);
    let tmp26 = tmp25.powi(2);
    let tmp27 = tmp20.powi(2);
    let tmp28 = tmp22 + tmp24 + tmp26 + tmp27;
    let tmp29 = tmp28.powi(-2);
    let tmp30 = 1.0 * dt;
    let tmp31 = q_old_node_i[1] * tmp19;
    let tmp32 = tmp30 * tmp31;
    let tmp33 = q_old_node_i[2] * tmp23;
    let tmp34 = tmp30 * tmp33;
    let tmp35 = q_old_node_i[3] * tmp25;
    let tmp36 = tmp30 * tmp35;
    let tmp37 = q_old_node_i[0] * tmp20;
    let tmp38 = tmp30 * tmp37;
    let tmp39 = -1.0 * tmp32 + tmp34 - 1.0 * tmp36 + tmp38;
    let tmp40 = tmp29 * tmp39;
    let tmp41 = 2.0 * tmp40;
    let tmp42 = tmp21 * tmp41;
    let tmp43 = tmp25 * tmp41;
    let tmp44 = tmp28.recip();
    let tmp45 = q_old_node_i[3] * tmp23;
    let tmp46 = tmp30 * tmp45;
    let tmp47 = tmp44 * tmp46;
    let tmp48 = q_old_node_i[2] * tmp25;
    let tmp49 = tmp30 * tmp48;
    let tmp50 = tmp44 * tmp49;
    let tmp51 = -1.0 * tmp50;
    let tmp52 = tmp47 + tmp51;
    let tmp53 = tmp23 * tmp43 + tmp52;
    let tmp54 = q_old_node_i[0] * tmp19;
    let tmp55 = tmp30 * tmp54;
    let tmp56 = tmp44 * tmp55;
    let tmp57 = q_old_node_i[1] * tmp20;
    let tmp58 = tmp30 * tmp57;
    let tmp59 = tmp44 * tmp58;
    let tmp60 = -1.0 * tmp59;
    let tmp61 = tmp56 + tmp60;
    let tmp62 = -1.0 * tmp42 + tmp53 + tmp61;
    let tmp63 = grad_phi_node_i[0] * tmp62;
    let tmp64 = tmp18 * tmp63;
    let tmp65 = tmp16 * tmp4;
    let tmp66 = tmp10 * tmp14;
    let tmp67 = tmp65 + tmp66;
    let tmp68 = tmp20 * tmp43;
    let tmp69 = tmp30 * tmp44;
    let tmp70 = q_old_node_i[1] * tmp23;
    let tmp71 = tmp69 * tmp70;
    let tmp72 = tmp23 * tmp41;
    let tmp73 = q_old_node_i[2] * tmp19;
    let tmp74 = tmp69 * tmp73;
    let tmp75 = -1.0 * tmp74;
    let tmp76 = tmp19 * tmp72 + tmp71 + tmp75;
    let tmp77 = q_old_node_i[0] * tmp25;
    let tmp78 = tmp69 * tmp77;
    let tmp79 = q_old_node_i[3] * tmp20;
    let tmp80 = tmp69 * tmp79;
    let tmp81 = -1.0 * tmp80;
    let tmp82 = tmp78 + tmp81;
    let tmp83 = -1.0 * tmp68 + tmp76 + tmp82;
    let tmp84 = grad_phi_node_i[1] * tmp83;
    let tmp85 = tmp67 * tmp84;
    let tmp86 = tmp14 * tmp7;
    let tmp87 = tmp10 * tmp4;
    let tmp88 = tmp13 * tmp87;
    let tmp89 = tmp86 + tmp88;
    let tmp90 = tmp19 * tmp43;
    let tmp91 = tmp20 * tmp72;
    let tmp92 = q_old_node_i[2] * tmp20;
    let tmp93 = tmp30 * tmp92;
    let tmp94 = tmp44 * tmp93;
    let tmp95 = q_old_node_i[3] * tmp19;
    let tmp96 = tmp30 * tmp95;
    let tmp97 = tmp44 * tmp96;
    let tmp98 = tmp94 + tmp97;
    let tmp99 = q_old_node_i[0] * tmp23;
    let tmp100 = tmp30 * tmp99;
    let tmp101 = tmp100 * tmp44;
    let tmp102 = q_old_node_i[1] * tmp25;
    let tmp103 = tmp102 * tmp30;
    let tmp104 = tmp103 * tmp44;
    let tmp105 = tmp101 + tmp104;
    let tmp106 = tmp105 + tmp90 - 1.0 * tmp91 + tmp98;
    let tmp107 = grad_phi_node_i[2] * tmp106;
    let tmp108 = tmp107 * tmp89;
    let tmp109 = -1.0 * tmp65;
    let tmp110 = tmp109 + tmp66;
    let tmp111 = -1.0 * tmp101;
    let tmp112 = -1.0 * tmp94;
    let tmp113 = tmp111 + tmp112;
    let tmp114 = tmp104 + tmp113 + tmp97;
    let tmp115 = tmp114 + tmp90 + tmp91;
    let tmp116 = grad_phi_node_i[0] * tmp115;
    let tmp117 = tmp110 * tmp116;
    let tmp118 = -1.0 * tmp86;
    let tmp119 = tmp118 + tmp88;
    let tmp120 = -1.0 * tmp56;
    let tmp121 = tmp120 + tmp42 + tmp53 + tmp59;
    let tmp122 = grad_phi_node_i[1] * tmp121;
    let tmp123 = tmp119 * tmp122;
    let tmp124 = -1.0 * tmp15;
    let tmp125 = tmp124 + tmp17;
    let tmp126 = -1.0 * tmp78;
    let tmp127 = tmp126 + tmp80;
    let tmp128 = tmp127 + tmp68 + tmp76;
    let tmp129 = grad_phi_node_i[2] * tmp128;
    let tmp130 = tmp125 * tmp129;
    let tmp131 = -1.0 * tmp88;
    let tmp132 = tmp131 + tmp86;
    let tmp133 = tmp38 * tmp44;
    let tmp134 = -1.0 * tmp133;
    let tmp135 = tmp32 * tmp44;
    let tmp136 = -1.0 * tmp135;
    let tmp137 = tmp34 * tmp44;
    let tmp138 = -1.0 * tmp137;
    let tmp139 = tmp24 * tmp40;
    let tmp140 = tmp27 * tmp40;
    let tmp141 = tmp22 * tmp40;
    let tmp142 = tmp140 - 1.0 * tmp141;
    let tmp143 = tmp26 * tmp40;
    let tmp144 = tmp36 * tmp44;
    let tmp145 = -1.0 * tmp143 - 1.0 * tmp144;
    let tmp146 = tmp134 + tmp136 + tmp138 + tmp139 + tmp142 + tmp145;
    let tmp147 = grad_phi_node_i[0] * tmp146;
    let tmp148 = tmp132 * tmp147;
    let tmp149 = -1.0 * tmp17;
    let tmp150 = tmp149 + tmp15;
    let tmp151 = -1.0 * tmp139;
    let tmp152 = tmp134 + tmp137;
    let tmp153 = tmp135 + tmp140 + tmp141 + tmp145 + tmp151 + tmp152;
    let tmp154 = grad_phi_node_i[1] * tmp153;
    let tmp155 = tmp150 * tmp154;
    let tmp156 = -1.0 * tmp66;
    let tmp157 = tmp156 + tmp65;
    let tmp158 = tmp136 + tmp144;
    let tmp159 = tmp152 + tmp158;
    let tmp160 = tmp142 + tmp143 + tmp151 + tmp159;
    let tmp161 = grad_phi_node_i[2] * tmp160;
    let tmp162 = tmp157 * tmp161;
    let tmp163 = tmp118 + tmp131;
    let tmp164 = grad_phi_node_i[0] * tmp153;
    let tmp165 = tmp163 * tmp164;
    let tmp166 = tmp124 + tmp149;
    let tmp167 = grad_phi_node_i[1] * tmp160;
    let tmp168 = tmp166 * tmp167;
    let tmp169 = tmp109 + tmp156;
    let tmp170 = grad_phi_node_i[2] * tmp146;
    let tmp171 = tmp169 * tmp170;
    let tmp172 = tmp13 * tmp8;
    let tmp173 = tmp13 * tmp6;
    let tmp174 = -1.0 * tmp173;
    let tmp175 = tmp11 * tmp13;
    let tmp176 = tmp13 * tmp9;
    let tmp177 = tmp175 - 1.0 * tmp176;
    let tmp178 = tmp172 + tmp174 + tmp177;
    let tmp179 = (1_f64 / 2.0) * tmp178;
    let tmp180 = grad_phi_node_i[0] * tmp83;
    let tmp181 = -1.0 * tmp172;
    let tmp182 = tmp173 + tmp177 + tmp181;
    let tmp183 = (1_f64 / 2.0) * tmp182;
    let tmp184 = grad_phi_node_i[0] * tmp128;
    let tmp185 = tmp174 + tmp175 + tmp176 + tmp181;
    let tmp186 = (1_f64 / 2.0) * tmp185;
    let tmp187 = grad_phi_node_i[1] * tmp115;
    let tmp188 = grad_phi_node_i[1] * tmp106;
    let tmp189 = grad_phi_node_i[2] * tmp62;
    let tmp190 = grad_phi_node_i[2] * tmp121;
    let tmp191 = 2.0 * tmp44;
    let tmp192 = tmp191 * tmp25;
    let tmp193 = tmp192 * tmp23;
    let tmp194 = tmp191 * tmp21;
    let tmp195 = tmp193 - 1.0 * tmp194;
    let tmp196 = grad_phi_node_i[0] * tmp195;
    let tmp197 = q_old_node_j[1]
        + tmp3
            * (-1.0 * omega_node_j[0] * q_old_node_j[2]
                + omega_node_j[1] * q_old_node_j[3]
                + omega_node_j[2] * q_old_node_j[0]);
    let tmp198 = q_old_node_j[0]
        + tmp3
            * (omega_node_j[0] * q_old_node_j[3] + omega_node_j[1] * q_old_node_j[2]
                - 1.0 * omega_node_j[2] * q_old_node_j[1]);
    let tmp199 = q_old_node_j[2]
        + tmp3
            * (omega_node_j[0] * q_old_node_j[1] - 1.0 * omega_node_j[1] * q_old_node_j[0]
                + omega_node_j[2] * q_old_node_j[3]);
    let tmp200 = tmp199.powi(2);
    let tmp201 = tmp197.powi(2);
    let tmp202 = tmp198.powi(2);
    let tmp203 = q_old_node_j[3]
        + tmp3
            * (-1.0 * omega_node_j[0] * q_old_node_j[0]
                - 1.0 * omega_node_j[1] * q_old_node_j[1]
                - 1.0 * omega_node_j[2] * q_old_node_j[2]);
    let tmp204 = tmp203.powi(2);
    let tmp205 = (tmp200 + tmp201 + tmp202 + tmp204).recip();
    let tmp206 = 2.0 * tmp205;
    let tmp207 = tmp198 * tmp206;
    let tmp208 = tmp197 * tmp207;
    let tmp209 = tmp199 * tmp203 * tmp206;
    let tmp210 = tmp208 - 1.0 * tmp209;
    let tmp211 = grad_phi_node_j[0] * tmp210;
    let tmp212 = dRdx0_const[(0, 1)] + tmp196 + tmp211;
    let tmp213 = tmp12.powi(-2);
    let tmp214 = phi_node_i * tmp30;
    let tmp215 = q_old[1] * tmp5;
    let tmp216 = q_old[2] * tmp7;
    let tmp217 = q_old[3] * tmp4;
    let tmp218 = tmp214 * tmp217;
    let tmp219 = q_old[0] * tmp10;
    let tmp220 = tmp214 * tmp219;
    let tmp221 = -1.0 * tmp214 * tmp215 + tmp214 * tmp216 - 1.0 * tmp218 + tmp220;
    let tmp222 = tmp213 * tmp221;
    let tmp223 = tmp222 * tmp5;
    let tmp224 = tmp223 * tmp4;
    let tmp225 = tmp222 * tmp7;
    let tmp226 = tmp10 * tmp225;
    let tmp227 = q_old[1] * tmp4;
    let tmp228 = phi_node_i * tmp3;
    let tmp229 = tmp13 * tmp228;
    let tmp230 = tmp227 * tmp229;
    let tmp231 = q_old[2] * tmp10;
    let tmp232 = tmp229 * tmp231;
    let tmp233 = -1.0 * tmp232;
    let tmp234 = tmp230 + tmp233;
    let tmp235 = q_old[3] * tmp228;
    let tmp236 = tmp14 * tmp235;
    let tmp237 = q_old[0] * tmp16;
    let tmp238 = tmp228 * tmp237;
    let tmp239 = -1.0 * tmp238;
    let tmp240 = tmp236 + tmp239;
    let tmp241 = tmp234 + tmp240;
    let tmp242 = tmp224 + tmp226 + tmp241;
    let tmp243 = tmp212 * tmp242;
    let tmp244 = tmp191 * tmp23;
    let tmp245 = tmp19 * tmp244;
    let tmp246 = tmp192 * tmp20;
    let tmp247 = tmp245 - 1.0 * tmp246;
    let tmp248 = grad_phi_node_i[1] * tmp247;
    let tmp249 = tmp197 * tmp206;
    let tmp250 = tmp199 * tmp249;
    let tmp251 = tmp203 * tmp207;
    let tmp252 = tmp250 - 1.0 * tmp251;
    let tmp253 = grad_phi_node_j[1] * tmp252;
    let tmp254 = dRdx1_const[(1, 2)] + tmp248 + tmp253;
    let tmp255 = tmp10 * tmp223;
    let tmp256 = q_old[1] * tmp10;
    let tmp257 = tmp229 * tmp256;
    let tmp258 = q_old[0] * tmp14;
    let tmp259 = tmp228 * tmp258;
    let tmp260 = -1.0 * tmp259;
    let tmp261 = tmp225 * tmp4;
    let tmp262 = tmp16 * tmp235;
    let tmp263 = q_old[2] * tmp4;
    let tmp264 = tmp229 * tmp263;
    let tmp265 = -1.0 * tmp264;
    let tmp266 = tmp262 + tmp265;
    let tmp267 = tmp261 + tmp266;
    let tmp268 = tmp255 + tmp257 + tmp260 + tmp267;
    let tmp269 = tmp254 * tmp268;
    let tmp270 = tmp19 * tmp192;
    let tmp271 = tmp20 * tmp244;
    let tmp272 = tmp270 - 1.0 * tmp271;
    let tmp273 = grad_phi_node_i[2] * tmp272;
    let tmp274 = tmp199 * tmp207;
    let tmp275 = tmp203 * tmp249;
    let tmp276 = tmp274 - 1.0 * tmp275;
    let tmp277 = grad_phi_node_j[2] * tmp276;
    let tmp278 = dRdx2_const[(2, 0)] + tmp273 + tmp277;
    let tmp279 = tmp223 * tmp7;
    let tmp280 = tmp222 * tmp87;
    let tmp281 = q_old[3] * tmp10;
    let tmp282 = tmp229 * tmp281;
    let tmp283 = q_old[0] * tmp4;
    let tmp284 = tmp229 * tmp283;
    let tmp285 = -1.0 * tmp284;
    let tmp286 = tmp282 + tmp285;
    let tmp287 = tmp280 + tmp286;
    let tmp288 = q_old[1] * tmp16;
    let tmp289 = tmp228 * tmp288;
    let tmp290 = q_old[2] * tmp14;
    let tmp291 = tmp228 * tmp290;
    let tmp292 = -1.0 * tmp291;
    let tmp293 = tmp289 + tmp292;
    let tmp294 = tmp279 + tmp287 + tmp293;
    let tmp295 = tmp278 * tmp294;
    let tmp296 = tmp270 + tmp271;
    let tmp297 = grad_phi_node_i[0] * tmp296;
    let tmp298 = tmp274 + tmp275;
    let tmp299 = grad_phi_node_j[0] * tmp298;
    let tmp300 = dRdx0_const[(0, 2)] + tmp297 + tmp299;
    let tmp301 = -1.0 * tmp261;
    let tmp302 = tmp257 + tmp264;
    let tmp303 = -1.0 * tmp262;
    let tmp304 = tmp260 + tmp303;
    let tmp305 = tmp302 + tmp304;
    let tmp306 = tmp255 + tmp301 + tmp305;
    let tmp307 = tmp300 * tmp306;
    let tmp308 = tmp193 + tmp194;
    let tmp309 = grad_phi_node_i[1] * tmp308;
    let tmp310 = tmp208 + tmp209;
    let tmp311 = grad_phi_node_j[1] * tmp310;
    let tmp312 = dRdx1_const[(1, 0)] + tmp309 + tmp311;
    let tmp313 = -1.0 * tmp289;
    let tmp314 = tmp291 + tmp313;
    let tmp315 = -1.0 * tmp279 + tmp314;
    let tmp316 = tmp287 + tmp315;
    let tmp317 = tmp312 * tmp316;
    let tmp318 = tmp245 + tmp246;
    let tmp319 = grad_phi_node_i[2] * tmp318;
    let tmp320 = tmp250 + tmp251;
    let tmp321 = grad_phi_node_j[2] * tmp320;
    let tmp322 = dRdx2_const[(2, 1)] + tmp319 + tmp321;
    let tmp323 = -1.0 * tmp230;
    let tmp324 = -1.0 * tmp236;
    let tmp325 = tmp323 + tmp324;
    let tmp326 = -1.0 * tmp224 + tmp325;
    let tmp327 = tmp226 + tmp233 + tmp239 + tmp326;
    let tmp328 = tmp322 * tmp327;
    let tmp329 = tmp24 * tmp44;
    let tmp330 = tmp26 * tmp44;
    let tmp331 = -1.0 * tmp330;
    let tmp332 = tmp27 * tmp44;
    let tmp333 = tmp22 * tmp44;
    let tmp334 = tmp332 - 1.0 * tmp333;
    let tmp335 = tmp329 + tmp331 + tmp334;
    let tmp336 = grad_phi_node_i[0] * tmp335;
    let tmp337 = tmp201 * tmp205;
    let tmp338 = tmp202 * tmp205;
    let tmp339 = -1.0 * tmp338;
    let tmp340 = tmp204 * tmp205;
    let tmp341 = tmp200 * tmp205;
    let tmp342 = tmp340 - 1.0 * tmp341;
    let tmp343 = tmp337 + tmp339 + tmp342;
    let tmp344 = grad_phi_node_j[0] * tmp343;
    let tmp345 = dRdx0_const[(1, 1)] + tmp336 + tmp344;
    let tmp346 = -1.0 * tmp280;
    let tmp347 = -1.0 * tmp282;
    let tmp348 = tmp289 + tmp347;
    let tmp349 = tmp284 + tmp292;
    let tmp350 = tmp348 + tmp349;
    let tmp351 = tmp279 + tmp346 + tmp350;
    let tmp352 = tmp345 * tmp351;
    let tmp353 = -1.0 * tmp329;
    let tmp354 = tmp331 + tmp332 + tmp333 + tmp353;
    let tmp355 = grad_phi_node_i[1] * tmp354;
    let tmp356 = -1.0 * tmp337;
    let tmp357 = tmp339 + tmp340 + tmp341 + tmp356;
    let tmp358 = grad_phi_node_j[1] * tmp357;
    let tmp359 = dRdx1_const[(2, 2)] + tmp355 + tmp358;
    let tmp360 = tmp232 + tmp238;
    let tmp361 = -1.0 * tmp226 + tmp360;
    let tmp362 = tmp224 + tmp230 + tmp236 + tmp361;
    let tmp363 = tmp359 * tmp362;
    let tmp364 = tmp330 + tmp334 + tmp353;
    let tmp365 = grad_phi_node_i[2] * tmp364;
    let tmp366 = tmp338 + tmp342 + tmp356;
    let tmp367 = grad_phi_node_j[2] * tmp366;
    let tmp368 = dRdx2_const[(0, 0)] + tmp365 + tmp367;
    let tmp369 = -1.0 * tmp257;
    let tmp370 = tmp259 + tmp369;
    let tmp371 = -1.0 * tmp255 + tmp370;
    let tmp372 = tmp267 + tmp371;
    let tmp373 = tmp368 * tmp372;
    let tmp374 = grad_phi_node_i[0] * tmp354;
    let tmp375 = grad_phi_node_j[0] * tmp357;
    let tmp376 = dRdx0_const[(2, 2)] + tmp374 + tmp375;
    let tmp377 = tmp284 + tmp347;
    let tmp378 = tmp315 + tmp346 + tmp377;
    let tmp379 = tmp376 * tmp378;
    let tmp380 = grad_phi_node_i[1] * tmp364;
    let tmp381 = grad_phi_node_j[1] * tmp366;
    let tmp382 = dRdx1_const[(0, 0)] + tmp380 + tmp381;
    let tmp383 = tmp326 + tmp361;
    let tmp384 = tmp382 * tmp383;
    let tmp385 = grad_phi_node_i[2] * tmp335;
    let tmp386 = grad_phi_node_j[2] * tmp343;
    let tmp387 = dRdx2_const[(1, 1)] + tmp385 + tmp386;
    let tmp388 = tmp264 + tmp301 + tmp303 + tmp371;
    let tmp389 = tmp387 * tmp388;
    let tmp390 = grad_phi_node_i[0] * tmp247;
    let tmp391 = grad_phi_node_j[0] * tmp252;
    let tmp392 = -1_f64 / 2.0 * dRdx0_const[(1, 2)] - 1_f64 / 2.0 * tmp390 - 1_f64 / 2.0 * tmp391;
    let tmp393 = tmp222 * tmp8;
    let tmp394 = tmp16 * tmp214;
    let tmp395 = q_old[2] * tmp394;
    let tmp396 = tmp11 * tmp222;
    let tmp397 = tmp222 * tmp6;
    let tmp398 = -1.0 * tmp13 * tmp220;
    let tmp399 = q_old[1] * tmp14;
    let tmp400 = tmp214 * tmp399;
    let tmp401 = tmp396 - 1.0 * tmp397 + tmp398 - 1.0 * tmp400;
    let tmp402 = tmp222 * tmp9;
    let tmp403 = tmp13 * tmp218;
    let tmp404 = -1.0 * tmp402 - 1.0 * tmp403;
    let tmp405 = tmp393 - 1.0 * tmp395 + tmp401 + tmp404;
    let tmp406 = tmp392 * tmp405;
    let tmp407 = grad_phi_node_i[0] * tmp318;
    let tmp408 = grad_phi_node_j[0] * tmp320;
    let tmp409 = (1_f64 / 2.0) * dRdx0_const[(2, 1)] + (1_f64 / 2.0) * tmp407 + (1_f64 / 2.0) * tmp408;
    let tmp410 = -1.0 * tmp393 + tmp395;
    let tmp411 = tmp396 + tmp397 + tmp398 + tmp400 + tmp404 + tmp410;
    let tmp412 = tmp409 * tmp411;
    let tmp413 = grad_phi_node_i[1] * tmp296;
    let tmp414 = grad_phi_node_j[1] * tmp298;
    let tmp415 = (1_f64 / 2.0) * dRdx1_const[(0, 2)] + (1_f64 / 2.0) * tmp413 + (1_f64 / 2.0) * tmp414;
    let tmp416 = tmp401 + tmp402 + tmp403 + tmp410;
    let tmp417 = tmp415 * tmp416;
    let tmp418 = grad_phi_node_i[1] * tmp272;
    let tmp419 = grad_phi_node_j[1] * tmp276;
    let tmp420 = -1_f64 / 2.0 * dRdx1_const[(2, 0)] - 1_f64 / 2.0 * tmp418 - 1_f64 / 2.0 * tmp419;
    let tmp421 = tmp411 * tmp420;
    let tmp422 = grad_phi_node_i[2] * tmp195;
    let tmp423 = grad_phi_node_j[2] * tmp210;
    let tmp424 = -1_f64 / 2.0 * dRdx2_const[(0, 1)] - 1_f64 / 2.0 * tmp422 - 1_f64 / 2.0 * tmp423;
    let tmp425 = tmp416 * tmp424;
    let tmp426 = grad_phi_node_i[2] * tmp308;
    let tmp427 = grad_phi_node_j[2] * tmp310;
    let tmp428 = (1_f64 / 2.0) * dRdx2_const[(1, 0)] + (1_f64 / 2.0) * tmp426 + (1_f64 / 2.0) * tmp427;
    let tmp429 = tmp405 * tmp428;
    let tmp430 = tmp182 * tmp184;
    let tmp431 = tmp185 * tmp187;
    let tmp432 = tmp178 * tmp190;
    let tmp433 = tmp178 * tmp180;
    let tmp434 = tmp182 * tmp188;
    let tmp435 = tmp185 * tmp189;
    let tmp436 = 1.0 * a3;
    let tmp437 = tmp436
        * (2.0 * tmp108
            + 2.0 * tmp117
            + 2.0 * tmp123
            + 2.0 * tmp130
            + 2.0 * tmp148
            + 2.0 * tmp155
            + 2.0 * tmp162
            + 2.0 * tmp165
            + 2.0 * tmp168
            + 2.0 * tmp171
            + 2.0 * tmp243
            + 2.0 * tmp269
            + 2.0 * tmp295
            + 2.0 * tmp307
            + 2.0 * tmp317
            + 2.0 * tmp328
            + 2.0 * tmp352
            + 2.0 * tmp363
            + 2.0 * tmp373
            + 2.0 * tmp379
            + 2.0 * tmp384
            + 2.0 * tmp389
            + 2.0 * tmp406
            + 2.0 * tmp412
            + 2.0 * tmp417
            + 2.0 * tmp421
            + 2.0 * tmp425
            + 2.0 * tmp429
            + tmp430
            + tmp431
            + tmp432
            - 1.0 * tmp433
            - 1.0 * tmp434
            - 1.0 * tmp435
            + 2.0 * tmp64
            + 2.0 * tmp85);
    let tmp438 = 2.0 * dt;
    let tmp439 = tmp438 * tmp54;
    let tmp440 = tmp40 * tmp439;
    let tmp441 = tmp438 * tmp57;
    let tmp442 = tmp40 * tmp441;
    let tmp443 = tmp31 * tmp438;
    let tmp444 = tmp33 * tmp438;
    let tmp445 = tmp35 * tmp438;
    let tmp446 = tmp37 * tmp438;
    let tmp447 = tmp28.powi(-3);
    let tmp448 = tmp39 * tmp447;
    let tmp449 = tmp448 * (-1.0 * tmp443 + tmp444 - 1.0 * tmp445 + tmp446);
    let tmp450 = 2.0 * tmp449;
    let tmp451 = tmp21 * tmp450;
    let tmp452 = tmp438 * tmp45;
    let tmp453 = tmp25 * tmp450;
    let tmp454 = tmp438 * tmp48;
    let tmp455 = tmp23 * tmp453 + tmp40 * tmp452 - 1.0 * tmp40 * tmp454;
    let tmp456 = dt.powi(2);
    let tmp457 = 0.5 * tmp456;
    let tmp458 = q_old_node_i[0].powi(2) * tmp457;
    let tmp459 = q_old_node_i[1].powi(2) * tmp457;
    let tmp460 = q_old_node_i[2].powi(2) * tmp457;
    let tmp461 = q_old_node_i[3].powi(2) * tmp457;
    let tmp462 = tmp29 * (-1.0 * tmp458 - 1.0 * tmp459 - 1.0 * tmp460 - 1.0 * tmp461);
    let tmp463 = 2.0 * tmp462;
    let tmp464 = tmp25 * tmp463;
    let tmp465 = tmp23 * tmp464;
    let tmp466 = tmp21 * tmp463;
    let tmp467 = tmp465 - 1.0 * tmp466;
    let tmp468 = 1.0 * tmp456;
    let tmp469 = tmp44 * tmp468;
    let tmp470 = q_old_node_i[2] * q_old_node_i[3] * tmp469;
    let tmp471 = -1.0 * tmp470;
    let tmp472 = q_old_node_i[0] * tmp469;
    let tmp473 = q_old_node_i[1] * tmp472;
    let tmp474 = tmp471 + tmp473;
    let tmp475 = tmp440 - 1.0 * tmp442 - 1.0 * tmp451 + tmp455 + tmp467 + tmp474;
    let tmp476 = 2.0 * grad_phi_node_i[0];
    let tmp477 = tmp475 * tmp476;
    let tmp478 = tmp40 * tmp438;
    let tmp479 = tmp478 * tmp77;
    let tmp480 = tmp478 * tmp79;
    let tmp481 = tmp20 * tmp453;
    let tmp482 = tmp23 * tmp450;
    let tmp483 = tmp19 * tmp482 + tmp478 * tmp70 - 1.0 * tmp478 * tmp73;
    let tmp484 = q_old_node_i[1] * tmp469;
    let tmp485 = q_old_node_i[2] * tmp484;
    let tmp486 = -1.0 * tmp485;
    let tmp487 = tmp20 * tmp464;
    let tmp488 = -1.0 * tmp487;
    let tmp489 = tmp23 * tmp463;
    let tmp490 = tmp19 * tmp489;
    let tmp491 = q_old_node_i[3] * tmp472;
    let tmp492 = tmp490 + tmp491;
    let tmp493 = tmp486 + tmp488 + tmp492;
    let tmp494 = tmp479 - 1.0 * tmp480 - 1.0 * tmp481 + tmp483 + tmp493;
    let tmp495 = 2.0 * grad_phi_node_i[1];
    let tmp496 = tmp494 * tmp495;
    let tmp497 = tmp438 * tmp99;
    let tmp498 = tmp40 * tmp497;
    let tmp499 = tmp438 * tmp92;
    let tmp500 = tmp40 * tmp499;
    let tmp501 = tmp20 * tmp482;
    let tmp502 = tmp102 * tmp438;
    let tmp503 = tmp438 * tmp95;
    let tmp504 = tmp19 * tmp453 + tmp40 * tmp502 + tmp40 * tmp503;
    let tmp505 = q_old_node_i[3] * tmp484;
    let tmp506 = q_old_node_i[2] * tmp472;
    let tmp507 = -1.0 * tmp506;
    let tmp508 = tmp19 * tmp464;
    let tmp509 = tmp20 * tmp489;
    let tmp510 = tmp508 - 1.0 * tmp509;
    let tmp511 = tmp505 + tmp507 + tmp510;
    let tmp512 = tmp498 + tmp500 - 1.0 * tmp501 + tmp504 + tmp511;
    let tmp513 = 2.0 * grad_phi_node_i[2];
    let tmp514 = tmp512 * tmp513;
    let tmp515 = tmp508 + tmp509;
    let tmp516 = tmp505 + tmp506;
    let tmp517 = -1.0 * tmp498 - 1.0 * tmp500 + tmp501 + tmp504 + tmp515 + tmp516;
    let tmp518 = tmp476 * tmp517;
    let tmp519 = -1.0 * tmp473;
    let tmp520 = tmp465 + tmp466;
    let tmp521 = tmp471 + tmp519 + tmp520;
    let tmp522 = -1.0 * tmp440 + tmp442 + tmp451 + tmp455 + tmp521;
    let tmp523 = tmp495 * tmp522;
    let tmp524 = -1.0 * tmp491;
    let tmp525 = tmp490 + tmp524;
    let tmp526 = tmp486 + tmp487 + tmp525;
    let tmp527 = -1.0 * tmp479 + tmp480 + tmp481 + tmp483 + tmp526;
    let tmp528 = tmp513 * tmp527;
    let tmp529 = 4.0 * tmp63;
    let tmp530 = 4.0 * tmp84;
    let tmp531 = 4.0 * tmp107;
    let tmp532 = 4.0 * tmp116;
    let tmp533 = 4.0 * tmp122;
    let tmp534 = 4.0 * tmp129;
    let tmp535 = grad_phi_node_i[0] * tmp182;
    let tmp536 = tmp527 * tmp535;
    let tmp537 = grad_phi_node_i[1] * tmp185;
    let tmp538 = tmp517 * tmp537;
    let tmp539 = grad_phi_node_i[2] * tmp178;
    let tmp540 = tmp522 * tmp539;
    let tmp541 = grad_phi_node_i[0] * tmp178;
    let tmp542 = tmp494 * tmp541;
    let tmp543 = grad_phi_node_i[1] * tmp182;
    let tmp544 = tmp512 * tmp543;
    let tmp545 = grad_phi_node_i[2] * tmp185;
    let tmp546 = tmp475 * tmp545;
    let tmp547 = 4.0 * tmp147;
    let tmp548 = 4.0 * tmp154;
    let tmp549 = 4.0 * tmp161;
    let tmp550 = 4.0 * tmp164;
    let tmp551 = 4.0 * tmp167;
    let tmp552 = 4.0 * tmp170;
    let tmp553 = tmp128 * tmp411;
    let tmp554 = tmp180 * tmp405;
    let tmp555 = tmp115 * tmp416;
    let tmp556 = tmp106 * tmp411;
    let tmp557 = tmp190 * tmp405;
    let tmp558 = tmp416 * tmp62;
    let tmp559 = tmp27 * tmp449;
    let tmp560 = -1.0 * tmp40 * tmp446;
    let tmp561 = tmp24 * tmp449;
    let tmp562 = tmp40 * tmp444;
    let tmp563 = tmp26 * tmp449;
    let tmp564 = tmp40 * tmp445;
    let tmp565 = tmp44 * tmp461;
    let tmp566 = -1.0 * tmp565;
    let tmp567 = tmp27 * tmp462;
    let tmp568 = tmp26 * tmp462;
    let tmp569 = tmp567 - 1.0 * tmp568;
    let tmp570 = tmp566 + tmp569;
    let tmp571 = -1.0 * tmp563 - 1.0 * tmp564 + tmp570;
    let tmp572 = tmp22 * tmp462;
    let tmp573 = -1.0 * tmp572;
    let tmp574 = tmp22 * tmp449;
    let tmp575 = tmp40 * tmp443;
    let tmp576 = tmp573 - 1.0 * tmp574 - 1.0 * tmp575;
    let tmp577 = tmp44 * tmp458;
    let tmp578 = tmp44 * tmp459;
    let tmp579 = -1.0 * tmp578;
    let tmp580 = tmp44 * tmp460;
    let tmp581 = tmp24 * tmp462;
    let tmp582 = tmp577 + tmp579 + tmp580 + tmp581;
    let tmp583 = tmp559 + tmp560 + tmp561 - 1.0 * tmp562 + tmp571 + tmp576 + tmp582;
    let tmp584 = tmp476 * tmp583;
    let tmp585 = -1.0 * tmp580;
    let tmp586 = tmp577 + tmp585;
    let tmp587 = tmp559 + tmp560 - 1.0 * tmp561 + tmp562 + tmp586;
    let tmp588 = -1.0 * tmp581;
    let tmp589 = tmp578 + tmp588;
    let tmp590 = tmp572 + tmp589;
    let tmp591 = tmp571 + tmp574 + tmp575 + tmp587 + tmp590;
    let tmp592 = tmp495 * tmp591;
    let tmp593 = tmp565 + tmp579;
    let tmp594 = tmp588 + tmp593;
    let tmp595 = tmp567 + tmp568;
    let tmp596 = tmp563 + tmp564 + tmp576 + tmp587 + tmp594 + tmp595;
    let tmp597 = tmp513 * tmp596;
    let tmp598 = tmp476 * tmp591;
    let tmp599 = tmp495 * tmp596;
    let tmp600 = tmp513 * tmp583;
    let tmp601 = 2.0 * dRdx0_const[(0, 1)] + 2.0 * tmp196 + 2.0 * tmp211;
    let tmp602 = phi_node_i * tmp438;
    let tmp603 = tmp217 * tmp602;
    let tmp604 = tmp219 * tmp602;
    let tmp605 = tmp12.powi(-3);
    let tmp606 = tmp221 * tmp605;
    let tmp607 = tmp606 * (-1.0 * tmp215 * tmp602 + tmp216 * tmp602 - 1.0 * tmp603 + tmp604);
    let tmp608 = tmp5 * tmp607;
    let tmp609 = tmp4 * tmp608;
    let tmp610 = tmp607 * tmp7;
    let tmp611 = tmp10 * tmp610;
    let tmp612 = q_old[0].powi(2);
    let tmp613 = phi_node_i.powi(2);
    let tmp614 = tmp457 * tmp613;
    let tmp615 = tmp612 * tmp614;
    let tmp616 = q_old[1].powi(2);
    let tmp617 = tmp614 * tmp616;
    let tmp618 = q_old[2].powi(2);
    let tmp619 = tmp614 * tmp618;
    let tmp620 = q_old[3].powi(2);
    let tmp621 = tmp614 * tmp620;
    let tmp622 = tmp213 * (-1.0 * tmp615 - 1.0 * tmp617 - 1.0 * tmp619 - 1.0 * tmp621);
    let tmp623 = tmp5 * tmp622;
    let tmp624 = tmp4 * tmp623;
    let tmp625 = tmp622 * tmp7;
    let tmp626 = tmp10 * tmp625;
    let tmp627 = tmp624 + tmp626;
    let tmp628 = tmp13 * tmp614;
    let tmp629 = q_old[0] * q_old[2];
    let tmp630 = tmp628 * tmp629;
    let tmp631 = q_old[3] * tmp628;
    let tmp632 = q_old[1] * tmp631;
    let tmp633 = tmp630 + tmp632;
    let tmp634 = q_old[3] * tmp214;
    let tmp635 = tmp223 * tmp634;
    let tmp636 = q_old[0] * tmp225;
    let tmp637 = tmp214 * tmp636;
    let tmp638 = -1.0 * tmp637;
    let tmp639 = tmp214 * tmp227;
    let tmp640 = tmp222 * tmp639;
    let tmp641 = tmp214 * tmp231;
    let tmp642 = tmp222 * tmp641;
    let tmp643 = -1.0 * tmp642;
    let tmp644 = tmp640 + tmp643;
    let tmp645 = tmp635 + tmp638 + tmp644;
    let tmp646 = 2.0 * dRdx1_const[(1, 2)] + 2.0 * tmp248 + 2.0 * tmp253;
    let tmp647 = q_old[0] * tmp223;
    let tmp648 = tmp214 * tmp647;
    let tmp649 = -1.0 * tmp648;
    let tmp650 = tmp10 * tmp608;
    let tmp651 = tmp214 * tmp256;
    let tmp652 = tmp222 * tmp651;
    let tmp653 = tmp650 + tmp652;
    let tmp654 = tmp10 * tmp623;
    let tmp655 = tmp4 * tmp625;
    let tmp656 = tmp654 + tmp655;
    let tmp657 = q_old[0] * q_old[1];
    let tmp658 = tmp628 * tmp657;
    let tmp659 = -1.0 * tmp658;
    let tmp660 = q_old[2] * tmp631;
    let tmp661 = -1.0 * tmp660;
    let tmp662 = tmp659 + tmp661;
    let tmp663 = tmp656 + tmp662;
    let tmp664 = tmp4 * tmp610;
    let tmp665 = tmp225 * tmp634;
    let tmp666 = tmp214 * tmp263;
    let tmp667 = tmp222 * tmp666;
    let tmp668 = -1.0 * tmp667;
    let tmp669 = tmp665 + tmp668;
    let tmp670 = tmp664 + tmp669;
    let tmp671 = 2.0 * dRdx2_const[(2, 0)] + 2.0 * tmp273 + 2.0 * tmp277;
    let tmp672 = tmp608 * tmp7;
    let tmp673 = q_old[1] * tmp225;
    let tmp674 = tmp214 * tmp673;
    let tmp675 = q_old[2] * tmp223;
    let tmp676 = tmp214 * tmp675;
    let tmp677 = tmp672 + tmp674 - 1.0 * tmp676;
    let tmp678 = tmp607 * tmp87;
    let tmp679 = tmp214 * tmp222;
    let tmp680 = tmp281 * tmp679;
    let tmp681 = tmp283 * tmp679;
    let tmp682 = tmp678 + tmp680 - 1.0 * tmp681;
    let tmp683 = q_old[1] * q_old[2];
    let tmp684 = tmp628 * tmp683;
    let tmp685 = -1.0 * tmp684;
    let tmp686 = q_old[0] * tmp631;
    let tmp687 = -1.0 * tmp686;
    let tmp688 = tmp685 + tmp687;
    let tmp689 = tmp623 * tmp7;
    let tmp690 = tmp622 * tmp87;
    let tmp691 = tmp689 + tmp690;
    let tmp692 = tmp688 + tmp691;
    let tmp693 = 2.0 * dRdx0_const[(0, 2)] + 2.0 * tmp297 + 2.0 * tmp299;
    let tmp694 = -1.0 * tmp655;
    let tmp695 = tmp654 + tmp694;
    let tmp696 = tmp659 + tmp660;
    let tmp697 = -1.0 * tmp665;
    let tmp698 = tmp649 + tmp697;
    let tmp699 = -1.0 * tmp664 + tmp667;
    let tmp700 = 2.0 * dRdx1_const[(1, 0)] + 2.0 * tmp309 + 2.0 * tmp311;
    let tmp701 = -1.0 * tmp672 - 1.0 * tmp674 + tmp676;
    let tmp702 = tmp684 + tmp687;
    let tmp703 = -1.0 * tmp689;
    let tmp704 = tmp690 + tmp703;
    let tmp705 = tmp702 + tmp704;
    let tmp706 = 2.0 * dRdx2_const[(2, 1)] + 2.0 * tmp319 + 2.0 * tmp321;
    let tmp707 = -1.0 * tmp640;
    let tmp708 = -1.0 * tmp609 + tmp707;
    let tmp709 = -1.0 * tmp635;
    let tmp710 = tmp638 + tmp643 + tmp709;
    let tmp711 = -1.0 * tmp632;
    let tmp712 = tmp630 + tmp711;
    let tmp713 = -1.0 * tmp624;
    let tmp714 = tmp626 + tmp713;
    let tmp715 = tmp712 + tmp714;
    let tmp716 = 2.0 * dRdx0_const[(1, 1)] + 2.0 * tmp336 + 2.0 * tmp344;
    let tmp717 = -1.0 * tmp690;
    let tmp718 = tmp689 + tmp717;
    let tmp719 = tmp685 + tmp686;
    let tmp720 = tmp718 + tmp719;
    let tmp721 = -1.0 * tmp678 - 1.0 * tmp680 + tmp681;
    let tmp722 = 2.0 * dRdx1_const[(2, 2)] + 2.0 * tmp355 + 2.0 * tmp358;
    let tmp723 = -1.0 * tmp611 + tmp642;
    let tmp724 = -1.0 * tmp626;
    let tmp725 = tmp624 + tmp724;
    let tmp726 = -1.0 * tmp630;
    let tmp727 = tmp632 + tmp726;
    let tmp728 = tmp725 + tmp727;
    let tmp729 = 2.0 * dRdx2_const[(0, 0)] + 2.0 * tmp365 + 2.0 * tmp367;
    let tmp730 = tmp658 + tmp661;
    let tmp731 = -1.0 * tmp654;
    let tmp732 = tmp655 + tmp731;
    let tmp733 = -1.0 * tmp652;
    let tmp734 = tmp648 + tmp733;
    let tmp735 = -1.0 * tmp650 + tmp734;
    let tmp736 = 2.0 * dRdx0_const[(2, 2)] + 2.0 * tmp374 + 2.0 * tmp375;
    let tmp737 = tmp703 + tmp717;
    let tmp738 = tmp684 + tmp686;
    let tmp739 = tmp737 + tmp738;
    let tmp740 = 2.0 * dRdx1_const[(0, 0)] + 2.0 * tmp380 + 2.0 * tmp381;
    let tmp741 = tmp711 + tmp726;
    let tmp742 = tmp713 + tmp724;
    let tmp743 = tmp637 + tmp709;
    let tmp744 = 2.0 * dRdx2_const[(1, 1)] + 2.0 * tmp385 + 2.0 * tmp386;
    let tmp745 = tmp658 + tmp660;
    let tmp746 = tmp694 + tmp731;
    let tmp747 = tmp745 + tmp746;
    let tmp748 = dRdx0_const[(2, 1)] + tmp407 + tmp408;
    let tmp749 = tmp6 * tmp607;
    let tmp750 = q_old[1] * tmp223;
    let tmp751 = tmp602 * tmp750;
    let tmp752 = tmp11 * tmp607;
    let tmp753 = tmp607 * tmp8;
    let tmp754 = q_old[2] * tmp225;
    let tmp755 = tmp602 * tmp754;
    let tmp756 = -1.0 * tmp222 * tmp604;
    let tmp757 = tmp13 * tmp615;
    let tmp758 = tmp13 * tmp619;
    let tmp759 = -1.0 * tmp758;
    let tmp760 = tmp757 + tmp759;
    let tmp761 = tmp752 - 1.0 * tmp753 + tmp755 + tmp756 + tmp760;
    let tmp762 = tmp607 * tmp9;
    let tmp763 = tmp222 * tmp603;
    let tmp764 = -1.0 * tmp762 - 1.0 * tmp763;
    let tmp765 = tmp13 * tmp617;
    let tmp766 = tmp13 * tmp621;
    let tmp767 = -1.0 * tmp766;
    let tmp768 = tmp765 + tmp767;
    let tmp769 = tmp622 * tmp8;
    let tmp770 = -1.0 * tmp769;
    let tmp771 = tmp6 * tmp622;
    let tmp772 = tmp11 * tmp622;
    let tmp773 = tmp622 * tmp9;
    let tmp774 = tmp772 - 1.0 * tmp773;
    let tmp775 = tmp770 + tmp771 + tmp774;
    let tmp776 = tmp768 + tmp775;
    let tmp777 = tmp749 + tmp751 + tmp761 + tmp764 + tmp776;
    let tmp778 = dRdx1_const[(0, 2)] + tmp413 + tmp414;
    let tmp779 = -1.0 * tmp771;
    let tmp780 = tmp770 + tmp772 + tmp773 + tmp779;
    let tmp781 = -1.0 * tmp765;
    let tmp782 = tmp766 + tmp781;
    let tmp783 = -1.0 * tmp749 - 1.0 * tmp751;
    let tmp784 = tmp761 + tmp762 + tmp763 + tmp780 + tmp782 + tmp783;
    let tmp785 = dRdx2_const[(1, 0)] + tmp426 + tmp427;
    let tmp786 = tmp769 + tmp774 + tmp779;
    let tmp787 = tmp757 + tmp758 + tmp767 + tmp781 + tmp786;
    let tmp788 = tmp752 + tmp753 - 1.0 * tmp755 + tmp756 + tmp764 + tmp783 + tmp787;
    let tmp789 = -1.0 * dRdx0_const[(1, 2)] - 1.0 * tmp390 - 1.0 * tmp391;
    let tmp790 = -1.0 * dRdx1_const[(2, 0)] - 1.0 * tmp418 - 1.0 * tmp419;
    let tmp791 = -1.0 * dRdx2_const[(0, 1)] - 1.0 * tmp422 - 1.0 * tmp423;
    let tmp792 = 2.0 * tmp65;
    let tmp793 = 2.0 * tmp66;
    let tmp794 = tmp792 - 1.0 * tmp793;
    let tmp795 = (1_f64 / 2.0) * tmp794;
    let tmp796 = tmp792 + tmp793;
    let tmp797 = (1_f64 / 2.0) * tmp796;
    let tmp798 = 2.0 * tmp15;
    let tmp799 = 2.0 * tmp17;
    let tmp800 = tmp798 - 1.0 * tmp799;
    let tmp801 = (1_f64 / 2.0) * tmp800;
    let tmp802 = tmp798 + tmp799;
    let tmp803 = (1_f64 / 2.0) * tmp802;
    let tmp804 = 2.0 * tmp86;
    let tmp805 = 2.0 * tmp88;
    let tmp806 = tmp804 - 1.0 * tmp805;
    let tmp807 = (1_f64 / 2.0) * tmp806;
    let tmp808 = tmp804 + tmp805;
    let tmp809 = (1_f64 / 2.0) * tmp808;
    let tmp810 = dRdx0_const[(1, 2)] + tmp390 + tmp391;
    let tmp811 = dRdx1_const[(2, 0)] + tmp418 + tmp419;
    let tmp812 = dRdx2_const[(0, 1)] + tmp422 + tmp423;
    let tmp813 = tmp436
        * (tmp179 * tmp785 - 1.0 * tmp179 * tmp810 + tmp183 * tmp748 - 1.0 * tmp183 * tmp811 + tmp186 * tmp778
            - 1.0 * tmp186 * tmp812
            + tmp212 * tmp803
            + tmp254 * tmp797
            + tmp278 * tmp809
            - 1.0 * tmp300 * tmp795
            - 1.0 * tmp312 * tmp807
            - 1.0 * tmp322 * tmp801
            + tmp345 * tmp807
            + tmp359 * tmp801
            + tmp368 * tmp795
            - 1.0 * tmp376 * tmp809
            - 1.0 * tmp382 * tmp803
            - 1.0 * tmp387 * tmp797);
    let tmp814 = 0.5 * tmp65;
    let tmp815 = -1.0 * tmp814;
    let tmp816 = 0.5 * tmp66;
    let tmp817 = -1.0 * tmp816;
    let tmp818 = tmp815 + tmp817;
    let tmp819 = tmp180 * tmp818;
    let tmp820 = 0.5 * tmp86;
    let tmp821 = 0.5 * tmp88;
    let tmp822 = -1.0 * tmp821;
    let tmp823 = tmp820 + tmp822;
    let tmp824 = grad_phi_node_i[0] * tmp121;
    let tmp825 = tmp823 * tmp824;
    let tmp826 = tmp815 + tmp816;
    let tmp827 = tmp187 * tmp826;
    let tmp828 = 0.5 * tmp15;
    let tmp829 = 0.5 * tmp17;
    let tmp830 = tmp828 + tmp829;
    let tmp831 = grad_phi_node_i[1] * tmp62;
    let tmp832 = tmp830 * tmp831;
    let tmp833 = -1.0 * tmp828;
    let tmp834 = tmp829 + tmp833;
    let tmp835 = tmp164 * tmp834;
    let tmp836 = grad_phi_node_i[0] * tmp160;
    let tmp837 = tmp830 * tmp836;
    let tmp838 = -1.0 * tmp820;
    let tmp839 = tmp822 + tmp838;
    let tmp840 = tmp154 * tmp839;
    let tmp841 = grad_phi_node_i[1] * tmp146;
    let tmp842 = tmp823 * tmp841;
    let tmp843 = (1_f64 / 4.0) * tmp185;
    let tmp844 = tmp116 * tmp843;
    let tmp845 = (1_f64 / 4.0) * tmp182;
    let tmp846 = grad_phi_node_i[0] * tmp106;
    let tmp847 = tmp845 * tmp846;
    let tmp848 = (1_f64 / 4.0) * tmp178;
    let tmp849 = tmp84 * tmp848;
    let tmp850 = grad_phi_node_i[1] * tmp128;
    let tmp851 = tmp845 * tmp850;
    let tmp852 = grad_phi_node_i[0] * tmp308;
    let tmp853 = grad_phi_node_j[0] * tmp310;
    let tmp854 = dRdx0_const[(1, 0)] + tmp852 + tmp853;
    let tmp855 = 0.5 * tmp279;
    let tmp856 = 0.5 * tmp280;
    let tmp857 = dt * phi_node_i;
    let tmp858 = 0.25 * tmp857;
    let tmp859 = tmp13 * tmp858;
    let tmp860 = tmp283 * tmp859;
    let tmp861 = tmp281 * tmp859;
    let tmp862 = -1.0 * tmp861;
    let tmp863 = tmp860 + tmp862;
    let tmp864 = -1.0 * tmp856 + tmp863;
    let tmp865 = tmp16 * tmp858;
    let tmp866 = q_old[1] * tmp865;
    let tmp867 = tmp14 * tmp858;
    let tmp868 = q_old[2] * tmp867;
    let tmp869 = -1.0 * tmp868;
    let tmp870 = tmp866 + tmp869;
    let tmp871 = tmp855 + tmp864 + tmp870;
    let tmp872 = tmp854 * tmp871;
    let tmp873 = 0.5 * tmp255;
    let tmp874 = -1.0 * tmp873;
    let tmp875 = tmp258 * tmp858;
    let tmp876 = tmp256 * tmp859;
    let tmp877 = -1.0 * tmp876;
    let tmp878 = 0.5 * tmp261;
    let tmp879 = tmp263 * tmp859;
    let tmp880 = q_old[3] * tmp865;
    let tmp881 = -1.0 * tmp880;
    let tmp882 = tmp879 + tmp881;
    let tmp883 = -1.0 * tmp878 + tmp882;
    let tmp884 = tmp874 + tmp875 + tmp877 + tmp883;
    let tmp885 = tmp810 * tmp884;
    let tmp886 = grad_phi_node_i[1] * tmp195;
    let tmp887 = grad_phi_node_j[1] * tmp210;
    let tmp888 = dRdx1_const[(0, 1)] + tmp886 + tmp887;
    let tmp889 = 0.5 * tmp226;
    let tmp890 = q_old[0] * tmp865;
    let tmp891 = -1.0 * tmp890;
    let tmp892 = tmp231 * tmp859;
    let tmp893 = -1.0 * tmp892;
    let tmp894 = tmp891 + tmp893;
    let tmp895 = tmp889 + tmp894;
    let tmp896 = 0.5 * tmp224;
    let tmp897 = tmp227 * tmp859;
    let tmp898 = q_old[3] * tmp867;
    let tmp899 = tmp897 + tmp898;
    let tmp900 = tmp896 + tmp899;
    let tmp901 = tmp895 + tmp900;
    let tmp902 = tmp888 * tmp901;
    let tmp903 = -1.0 * tmp875;
    let tmp904 = tmp876 + tmp903;
    let tmp905 = tmp873 + tmp904;
    let tmp906 = tmp883 + tmp905;
    let tmp907 = tmp778 * tmp906;
    let tmp908 = grad_phi_node_i[0] * tmp364;
    let tmp909 = grad_phi_node_j[0] * tmp366;
    let tmp910 = dRdx0_const[(0, 0)] + tmp908 + tmp909;
    let tmp911 = tmp901 * tmp910;
    let tmp912 = -1.0 * tmp897;
    let tmp913 = -1.0 * tmp898;
    let tmp914 = tmp912 + tmp913;
    let tmp915 = -1.0 * tmp896 + tmp914;
    let tmp916 = tmp895 + tmp915;
    let tmp917 = tmp376 * tmp916;
    let tmp918 = grad_phi_node_i[1] * tmp335;
    let tmp919 = grad_phi_node_j[1] * tmp343;
    let tmp920 = dRdx1_const[(1, 1)] + tmp918 + tmp919;
    let tmp921 = tmp871 * tmp920;
    let tmp922 = -1.0 * tmp866;
    let tmp923 = tmp868 + tmp922;
    let tmp924 = -1.0 * tmp855 + tmp923;
    let tmp925 = tmp864 + tmp924;
    let tmp926 = tmp359 * tmp925;
    let tmp927 = (1_f64 / 4.0) * dRdx0_const[(0, 2)];
    let tmp928 = (1_f64 / 4.0) * tmp297;
    let tmp929 = (1_f64 / 4.0) * tmp299;
    let tmp930 = -1.0 * tmp927 - 1.0 * tmp928 - 1.0 * tmp929;
    let tmp931 = tmp416 * tmp930;
    let tmp932 = (1_f64 / 4.0) * dRdx0_const[(2, 0)];
    let tmp933 = grad_phi_node_i[0] * tmp272;
    let tmp934 = (1_f64 / 4.0) * tmp933;
    let tmp935 = grad_phi_node_j[0] * tmp276;
    let tmp936 = (1_f64 / 4.0) * tmp935;
    let tmp937 = tmp932 + tmp934 + tmp936;
    let tmp938 = tmp411 * tmp937;
    let tmp939 = (1_f64 / 4.0) * dRdx1_const[(1, 2)];
    let tmp940 = (1_f64 / 4.0) * tmp248;
    let tmp941 = (1_f64 / 4.0) * tmp253;
    let tmp942 = -1.0 * tmp939 - 1.0 * tmp940 - 1.0 * tmp941;
    let tmp943 = tmp405 * tmp942;
    let tmp944 = (1_f64 / 4.0) * dRdx1_const[(2, 1)];
    let tmp945 = grad_phi_node_i[1] * tmp318;
    let tmp946 = (1_f64 / 4.0) * tmp945;
    let tmp947 = grad_phi_node_j[1] * tmp320;
    let tmp948 = (1_f64 / 4.0) * tmp947;
    let tmp949 = tmp944 + tmp946 + tmp948;
    let tmp950 = tmp411 * tmp949;
    let tmp951 = tmp183 * tmp846;
    let tmp952 = tmp116 * tmp186;
    let tmp953 = tmp183 * tmp850;
    let tmp954 = tmp179 * tmp84;
    let tmp955 = 2.0 * tmp827
        + 2.0 * tmp832
        + 2.0 * tmp840
        + 2.0 * tmp842
        + 2.0 * tmp902
        + 2.0 * tmp907
        + 2.0 * tmp921
        + 2.0 * tmp926
        + 2.0 * tmp943
        + 2.0 * tmp950
        + tmp953
        - 1.0 * tmp954;
    let tmp956 = 2.0 * tmp819
        + 2.0 * tmp825
        + 2.0 * tmp835
        + 2.0 * tmp837
        + 2.0 * tmp872
        + 2.0 * tmp885
        + 2.0 * tmp911
        + 2.0 * tmp917
        + 2.0 * tmp931
        + 2.0 * tmp938
        + tmp951
        - 1.0 * tmp952
        + tmp955;
    let tmp957 = tmp184 * tmp834;
    let tmp958 = tmp820 + tmp821;
    let tmp959 = tmp846 * tmp958;
    let tmp960 = tmp814 + tmp817;
    let tmp961 = grad_phi_node_i[2] * tmp115;
    let tmp962 = tmp960 * tmp961;
    let tmp963 = -1.0 * tmp829;
    let tmp964 = tmp833 + tmp963;
    let tmp965 = tmp189 * tmp964;
    let tmp966 = tmp147 * tmp818;
    let tmp967 = tmp836 * tmp960;
    let tmp968 = tmp821 + tmp838;
    let tmp969 = tmp170 * tmp968;
    let tmp970 = grad_phi_node_i[2] * tmp153;
    let tmp971 = tmp958 * tmp970;
    let tmp972 = tmp63 * tmp843;
    let tmp973 = tmp824 * tmp848;
    let tmp974 = grad_phi_node_i[2] * tmp83;
    let tmp975 = tmp848 * tmp974;
    let tmp976 = tmp129 * tmp845;
    let tmp977 = dRdx0_const[(2, 0)] + tmp933 + tmp935;
    let tmp978 = -1.0 * tmp860;
    let tmp979 = tmp866 + tmp978;
    let tmp980 = tmp861 + tmp869;
    let tmp981 = tmp979 + tmp980;
    let tmp982 = tmp855 + tmp856 + tmp981;
    let tmp983 = tmp977 * tmp982;
    let tmp984 = tmp748 * tmp916;
    let tmp985 = tmp890 + tmp892;
    let tmp986 = -1.0 * tmp889 + tmp985;
    let tmp987 = tmp915 + tmp986;
    let tmp988 = tmp812 * tmp987;
    let tmp989 = grad_phi_node_i[2] * tmp296;
    let tmp990 = grad_phi_node_j[2] * tmp298;
    let tmp991 = dRdx2_const[(0, 2)] + tmp989 + tmp990;
    let tmp992 = -1.0 * tmp879;
    let tmp993 = tmp877 + tmp992;
    let tmp994 = tmp875 + tmp880;
    let tmp995 = tmp993 + tmp994;
    let tmp996 = tmp874 + tmp878 + tmp995;
    let tmp997 = tmp991 * tmp996;
    let tmp998 = tmp910 * tmp996;
    let tmp999 = tmp345 * tmp884;
    let tmp1000 = tmp861 + tmp978;
    let tmp1001 = tmp1000 + tmp856 + tmp924;
    let tmp1002 = tmp1001 * tmp387;
    let tmp1003 = grad_phi_node_i[2] * tmp354;
    let tmp1004 = grad_phi_node_j[2] * tmp357;
    let tmp1005 = dRdx2_const[(2, 2)] + tmp1003 + tmp1004;
    let tmp1006 = tmp1005 * tmp982;
    let tmp1007 = (1_f64 / 4.0) * dRdx0_const[(0, 1)];
    let tmp1008 = (1_f64 / 4.0) * tmp196;
    let tmp1009 = (1_f64 / 4.0) * tmp211;
    let tmp1010 = -1.0 * tmp1007 - 1.0 * tmp1008 - 1.0 * tmp1009;
    let tmp1011 = tmp1010 * tmp416;
    let tmp1012 = (1_f64 / 4.0) * dRdx0_const[(1, 0)];
    let tmp1013 = (1_f64 / 4.0) * tmp852;
    let tmp1014 = (1_f64 / 4.0) * tmp853;
    let tmp1015 = tmp1012 + tmp1013 + tmp1014;
    let tmp1016 = tmp1015 * tmp405;
    let tmp1017 = (1_f64 / 4.0) * dRdx2_const[(1, 2)];
    let tmp1018 = grad_phi_node_i[2] * tmp247;
    let tmp1019 = (1_f64 / 4.0) * tmp1018;
    let tmp1020 = grad_phi_node_j[2] * tmp252;
    let tmp1021 = (1_f64 / 4.0) * tmp1020;
    let tmp1022 = tmp1017 + tmp1019 + tmp1021;
    let tmp1023 = tmp1022 * tmp405;
    let tmp1024 = (1_f64 / 4.0) * dRdx2_const[(2, 1)];
    let tmp1025 = (1_f64 / 4.0) * tmp319;
    let tmp1026 = (1_f64 / 4.0) * tmp321;
    let tmp1027 = -1.0 * tmp1024 - 1.0 * tmp1025 - 1.0 * tmp1026;
    let tmp1028 = tmp1027 * tmp411;
    let tmp1029 = tmp179 * tmp974;
    let tmp1030 = tmp129 * tmp183;
    let tmp1031 = tmp179 * tmp824;
    let tmp1032 = tmp186 * tmp63;
    let tmp1033 = 2.0 * tmp1011 + 2.0 * tmp1016 + tmp1031 - 1.0 * tmp1032
        + 2.0 * tmp957
        + 2.0 * tmp959
        + 2.0 * tmp966
        + 2.0 * tmp967
        + 2.0 * tmp983
        + 2.0 * tmp984
        + 2.0 * tmp998
        + 2.0 * tmp999;
    let tmp1034 = 2.0 * tmp1002 + 2.0 * tmp1006 + 2.0 * tmp1023 + 2.0 * tmp1028 + tmp1029 - 1.0 * tmp1030
        + tmp1033
        + 2.0 * tmp962
        + 2.0 * tmp965
        + 2.0 * tmp969
        + 2.0 * tmp971
        + 2.0 * tmp988
        + 2.0 * tmp997;
    let tmp1035 = tmp828 + tmp963;
    let tmp1036 = tmp1035 * tmp184;
    let tmp1037 = tmp839 * tmp846;
    let tmp1038 = tmp826 * tmp961;
    let tmp1039 = tmp189 * tmp830;
    let tmp1040 = tmp826 * tmp836;
    let tmp1041 = tmp814 + tmp816;
    let tmp1042 = tmp1041 * tmp147;
    let tmp1043 = tmp839 * tmp970;
    let tmp1044 = tmp170 * tmp823;
    let tmp1045 = tmp925 * tmp977;
    let tmp1046 = tmp900 + tmp986;
    let tmp1047 = tmp1046 * tmp748;
    let tmp1048 = tmp812 * tmp901;
    let tmp1049 = tmp906 * tmp991;
    let tmp1050 = tmp906 * tmp910;
    let tmp1051 = tmp878 + tmp880 + tmp905 + tmp992;
    let tmp1052 = tmp1051 * tmp345;
    let tmp1053 = tmp387 * tmp871;
    let tmp1054 = tmp1005 * tmp925;
    let tmp1055 = tmp1007 + tmp1008 + tmp1009;
    let tmp1056 = tmp1055 * tmp416;
    let tmp1057 = -1.0 * tmp1012 - 1.0 * tmp1013 - 1.0 * tmp1014;
    let tmp1058 = tmp1057 * tmp405;
    let tmp1059 = -1.0 * tmp1017 - 1.0 * tmp1019 - 1.0 * tmp1021;
    let tmp1060 = tmp1059 * tmp405;
    let tmp1061 = tmp1024 + tmp1025 + tmp1026;
    let tmp1062 = tmp1061 * tmp411;
    let tmp1063 = -1.0 * tmp1029
        + tmp1030
        + 2.0 * tmp1038
        + 2.0 * tmp1039
        + 2.0 * tmp1043
        + 2.0 * tmp1044
        + 2.0 * tmp1048
        + 2.0 * tmp1049
        + 2.0 * tmp1053
        + 2.0 * tmp1054
        + 2.0 * tmp1060
        + 2.0 * tmp1062;
    let tmp1064 = -1.0 * tmp1031
        + tmp1032
        + 2.0 * tmp1036
        + 2.0 * tmp1037
        + 2.0 * tmp1040
        + 2.0 * tmp1042
        + 2.0 * tmp1045
        + 2.0 * tmp1047
        + 2.0 * tmp1050
        + 2.0 * tmp1052
        + 2.0 * tmp1056
        + 2.0 * tmp1058
        + tmp1063;
    let tmp1065 = tmp1041 * tmp180;
    let tmp1066 = tmp824 * tmp968;
    let tmp1067 = tmp187 * tmp960;
    let tmp1068 = tmp831 * tmp964;
    let tmp1069 = tmp836 * tmp964;
    let tmp1070 = tmp1035 * tmp164;
    let tmp1071 = tmp841 * tmp968;
    let tmp1072 = tmp154 * tmp958;
    let tmp1073 = tmp1001 * tmp854;
    let tmp1074 = tmp1051 * tmp810;
    let tmp1075 = tmp888 * tmp987;
    let tmp1076 = tmp778 * tmp996;
    let tmp1077 = tmp910 * tmp987;
    let tmp1078 = tmp1046 * tmp376;
    let tmp1079 = tmp1001 * tmp920;
    let tmp1080 = tmp359 * tmp982;
    let tmp1081 = tmp927 + tmp928 + tmp929;
    let tmp1082 = tmp1081 * tmp416;
    let tmp1083 = -1.0 * tmp932 - 1.0 * tmp934 - 1.0 * tmp936;
    let tmp1084 = tmp1083 * tmp411;
    let tmp1085 = tmp939 + tmp940 + tmp941;
    let tmp1086 = tmp1085 * tmp405;
    let tmp1087 = -1.0 * tmp944 - 1.0 * tmp946 - 1.0 * tmp948;
    let tmp1088 = tmp1087 * tmp411;
    let tmp1089 = 2.0 * tmp1065
        + 2.0 * tmp1066
        + 2.0 * tmp1069
        + 2.0 * tmp1070
        + 2.0 * tmp1073
        + 2.0 * tmp1074
        + 2.0 * tmp1077
        + 2.0 * tmp1078
        + 2.0 * tmp1082
        + 2.0 * tmp1084
        - 1.0 * tmp951
        + tmp952;
    let tmp1090 = 2.0 * tmp1067
        + 2.0 * tmp1068
        + 2.0 * tmp1071
        + 2.0 * tmp1072
        + 2.0 * tmp1075
        + 2.0 * tmp1076
        + 2.0 * tmp1079
        + 2.0 * tmp1080
        + 2.0 * tmp1086
        + 2.0 * tmp1088
        + tmp1089
        - 1.0 * tmp953
        + tmp954;
    let tmp1091 = tmp834 * tmp850;
    let tmp1092 = tmp188 * tmp958;
    let tmp1093 = tmp818 * tmp974;
    let tmp1094 = tmp190 * tmp823;
    let tmp1095 = tmp818 * tmp841;
    let tmp1096 = tmp167 * tmp960;
    let tmp1097 = tmp834 * tmp970;
    let tmp1098 = tmp161 * tmp830;
    let tmp1099 = tmp831 * tmp843;
    let tmp1100 = tmp122 * tmp848;
    let tmp1101 = tmp843 * tmp961;
    let tmp1102 = tmp107 * tmp845;
    let tmp1103 = tmp811 * tmp982;
    let tmp1104 = dRdx1_const[(2, 1)] + tmp945 + tmp947;
    let tmp1105 = tmp1104 * tmp916;
    let tmp1106 = tmp785 * tmp871;
    let tmp1107 = dRdx2_const[(1, 2)] + tmp1018 + tmp1020;
    let tmp1108 = tmp1107 * tmp884;
    let tmp1109 = tmp382 * tmp996;
    let tmp1110 = tmp884 * tmp920;
    let tmp1111 = tmp368 * tmp901;
    let tmp1112 = tmp1005 * tmp916;
    let tmp1113 = (1_f64 / 4.0) * dRdx1_const[(0, 1)];
    let tmp1114 = (1_f64 / 4.0) * tmp886;
    let tmp1115 = (1_f64 / 4.0) * tmp887;
    let tmp1116 = -1.0 * tmp1113 - 1.0 * tmp1114 - 1.0 * tmp1115;
    let tmp1117 = tmp1116 * tmp416;
    let tmp1118 = (1_f64 / 4.0) * dRdx1_const[(1, 0)];
    let tmp1119 = (1_f64 / 4.0) * tmp309;
    let tmp1120 = (1_f64 / 4.0) * tmp311;
    let tmp1121 = tmp1118 + tmp1119 + tmp1120;
    let tmp1122 = tmp1121 * tmp405;
    let tmp1123 = (1_f64 / 4.0) * dRdx2_const[(0, 2)];
    let tmp1124 = (1_f64 / 4.0) * tmp989;
    let tmp1125 = (1_f64 / 4.0) * tmp990;
    let tmp1126 = -1.0 * tmp1123 - 1.0 * tmp1124 - 1.0 * tmp1125;
    let tmp1127 = tmp1126 * tmp416;
    let tmp1128 = (1_f64 / 4.0) * dRdx2_const[(2, 0)];
    let tmp1129 = (1_f64 / 4.0) * tmp273;
    let tmp1130 = (1_f64 / 4.0) * tmp277;
    let tmp1131 = tmp1128 + tmp1129 + tmp1130;
    let tmp1132 = tmp1131 * tmp411;
    let tmp1133 = tmp107 * tmp183;
    let tmp1134 = tmp186 * tmp961;
    let tmp1135 = tmp122 * tmp179;
    let tmp1136 = tmp186 * tmp831;
    let tmp1137 = 2.0 * tmp1091
        + 2.0 * tmp1092
        + 2.0 * tmp1095
        + 2.0 * tmp1096
        + 2.0 * tmp1103
        + 2.0 * tmp1105
        + 2.0 * tmp1109
        + 2.0 * tmp1110
        + 2.0 * tmp1117
        + 2.0 * tmp1122
        + tmp1135
        - 1.0 * tmp1136;
    let tmp1138 = 2.0 * tmp1093
        + 2.0 * tmp1094
        + 2.0 * tmp1097
        + 2.0 * tmp1098
        + 2.0 * tmp1106
        + 2.0 * tmp1108
        + 2.0 * tmp1111
        + 2.0 * tmp1112
        + 2.0 * tmp1127
        + 2.0 * tmp1132
        + tmp1133
        - 1.0 * tmp1134
        + tmp1137;
    let tmp1139 = tmp1035 * tmp850;
    let tmp1140 = tmp188 * tmp839;
    let tmp1141 = tmp1041 * tmp974;
    let tmp1142 = tmp190 * tmp968;
    let tmp1143 = tmp167 * tmp826;
    let tmp1144 = tmp1041 * tmp841;
    let tmp1145 = tmp161 * tmp964;
    let tmp1146 = tmp1035 * tmp970;
    let tmp1147 = tmp811 * tmp925;
    let tmp1148 = tmp1046 * tmp1104;
    let tmp1149 = tmp1001 * tmp785;
    let tmp1150 = tmp1051 * tmp1107;
    let tmp1151 = tmp382 * tmp906;
    let tmp1152 = tmp1051 * tmp920;
    let tmp1153 = tmp368 * tmp987;
    let tmp1154 = tmp1005 * tmp1046;
    let tmp1155 = tmp1113 + tmp1114 + tmp1115;
    let tmp1156 = tmp1155 * tmp416;
    let tmp1157 = -1.0 * tmp1118 - 1.0 * tmp1119 - 1.0 * tmp1120;
    let tmp1158 = tmp1157 * tmp405;
    let tmp1159 = tmp1123 + tmp1124 + tmp1125;
    let tmp1160 = tmp1159 * tmp416;
    let tmp1161 = -1.0 * tmp1128 - 1.0 * tmp1129 - 1.0 * tmp1130;
    let tmp1162 = tmp1161 * tmp411;
    let tmp1163 = -1.0 * tmp1133
        + tmp1134
        + 2.0 * tmp1141
        + 2.0 * tmp1142
        + 2.0 * tmp1145
        + 2.0 * tmp1146
        + 2.0 * tmp1149
        + 2.0 * tmp1150
        + 2.0 * tmp1153
        + 2.0 * tmp1154
        + 2.0 * tmp1160
        + 2.0 * tmp1162;
    let tmp1164 = -1.0 * tmp1135
        + tmp1136
        + 2.0 * tmp1139
        + 2.0 * tmp1140
        + 2.0 * tmp1143
        + 2.0 * tmp1144
        + 2.0 * tmp1147
        + 2.0 * tmp1148
        + 2.0 * tmp1151
        + 2.0 * tmp1152
        + 2.0 * tmp1156
        + 2.0 * tmp1158
        + tmp1163;
    let tmp1165 = (1_f64 / 4.0) * tmp794;
    let tmp1166 = tmp1165 * tmp991;
    let tmp1167 = (1_f64 / 4.0) * tmp800;
    let tmp1168 = tmp1167 * tmp748;
    let tmp1169 = (1_f64 / 4.0) * tmp802;
    let tmp1170 = tmp1169 * tmp812;
    let tmp1171 = (1_f64 / 4.0) * tmp808;
    let tmp1172 = tmp1171 * tmp977;
    let tmp1173 = tmp1165 * tmp910;
    let tmp1174 = (1_f64 / 4.0) * tmp796;
    let tmp1175 = tmp1174 * tmp345;
    let tmp1176 = (1_f64 / 4.0) * tmp806;
    let tmp1177 = tmp1176 * tmp387;
    let tmp1178 = tmp1005 * tmp1171;
    let tmp1179 = tmp212 * tmp843;
    let tmp1180 = tmp848 * tmp854;
    let tmp1181 = tmp1107 * tmp848;
    let tmp1182 = tmp322 * tmp845;
    let tmp1183 = -1.0 * tmp1166 + tmp1168 + tmp1170 - 1.0 * tmp1172 - 1.0 * tmp1173 + tmp1175 + tmp1177
        - 1.0 * tmp1178
        + tmp1179
        - 1.0 * tmp1180
        - 1.0 * tmp1181
        + tmp1182;
    let tmp1184 = (1_f64 / 2.0) * dRdx0_const[(0, 1)];
    let tmp1185 = (1_f64 / 2.0) * tmp196;
    let tmp1186 = (1_f64 / 2.0) * tmp211;
    let tmp1187 = tmp1184 + tmp1185 + tmp1186;
    let tmp1188 = 2.0 * dRdx0_const[(0, 0)] + 2.0 * tmp908 + 2.0 * tmp909;
    let tmp1189 = 0.5 * tmp650;
    let tmp1190 = 0.5 * tmp664;
    let tmp1191 = -1.0 * tmp1190;
    let tmp1192 = 0.5 * tmp654;
    let tmp1193 = 0.5 * tmp655;
    let tmp1194 = -1.0 * tmp1193;
    let tmp1195 = tmp1192 + tmp1194;
    let tmp1196 = tmp13 * tmp613;
    let tmp1197 = tmp1196 * tmp456;
    let tmp1198 = 0.25 * tmp1197;
    let tmp1199 = tmp1198 * tmp657;
    let tmp1200 = -1.0 * tmp1199;
    let tmp1201 = q_old[2] * tmp1198;
    let tmp1202 = q_old[3] * tmp1201;
    let tmp1203 = tmp1200 + tmp1202;
    let tmp1204 = tmp222 * tmp228;
    let tmp1205 = tmp1204 * tmp256;
    let tmp1206 = tmp223 * tmp228;
    let tmp1207 = q_old[0] * tmp1206;
    let tmp1208 = -1.0 * tmp1207;
    let tmp1209 = tmp1204 * tmp263;
    let tmp1210 = tmp225 * tmp235;
    let tmp1211 = -1.0 * tmp1210;
    let tmp1212 = tmp1205 + tmp1208 + tmp1209 + tmp1211;
    let tmp1213 = tmp1189 + tmp1191 + tmp1195 + tmp1203 + tmp1212;
    let tmp1214 = tmp1205 + tmp1210;
    let tmp1215 = -1.0 * tmp1209;
    let tmp1216 = tmp1208 + tmp1215;
    let tmp1217 = tmp1192 + tmp1193;
    let tmp1218 = -1.0 * tmp1202;
    let tmp1219 = tmp1200 + tmp1218;
    let tmp1220 = tmp1217 + tmp1219;
    let tmp1221 = tmp1189 + tmp1190 + tmp1214 + tmp1216 + tmp1220;
    let tmp1222 = 2.0 * dRdx0_const[(2, 0)] + 2.0 * tmp933 + 2.0 * tmp935;
    let tmp1223 = 0.5 * tmp672;
    let tmp1224 = -1.0 * tmp1223;
    let tmp1225 = 0.5 * tmp678;
    let tmp1226 = -1.0 * tmp1225;
    let tmp1227 = 0.5 * tmp689;
    let tmp1228 = -1.0 * tmp1227;
    let tmp1229 = 0.5 * tmp690;
    let tmp1230 = -1.0 * tmp1229;
    let tmp1231 = tmp1228 + tmp1230;
    let tmp1232 = q_old[1] * tmp1201;
    let tmp1233 = q_old[3] * tmp1198;
    let tmp1234 = q_old[0] * tmp1233;
    let tmp1235 = tmp1232 + tmp1234;
    let tmp1236 = tmp1231 + tmp1235;
    let tmp1237 = q_old[2] * tmp1206;
    let tmp1238 = tmp225 * tmp228;
    let tmp1239 = q_old[1] * tmp1238;
    let tmp1240 = -1.0 * tmp1239;
    let tmp1241 = tmp1204 * tmp283;
    let tmp1242 = tmp1204 * tmp281;
    let tmp1243 = -1.0 * tmp1242;
    let tmp1244 = tmp1237 + tmp1240 + tmp1241 + tmp1243;
    let tmp1245 = tmp1224 + tmp1226 + tmp1236 + tmp1244;
    let tmp1246 = 2.0 * dRdx0_const[(2, 1)] + 2.0 * tmp407 + 2.0 * tmp408;
    let tmp1247 = 0.5 * tmp611;
    let tmp1248 = -1.0 * tmp1247;
    let tmp1249 = q_old[0] * tmp1238;
    let tmp1250 = tmp1204 * tmp231;
    let tmp1251 = 0.5 * tmp626;
    let tmp1252 = -1.0 * tmp1251;
    let tmp1253 = 0.5 * tmp624;
    let tmp1254 = tmp1252 + tmp1253;
    let tmp1255 = q_old[0] * tmp1201;
    let tmp1256 = -1.0 * tmp1255;
    let tmp1257 = q_old[1] * tmp1233;
    let tmp1258 = tmp1256 + tmp1257;
    let tmp1259 = tmp1254 + tmp1258;
    let tmp1260 = 0.5 * tmp609;
    let tmp1261 = tmp1204 * tmp227;
    let tmp1262 = tmp223 * tmp235;
    let tmp1263 = tmp1261 + tmp1262;
    let tmp1264 = tmp1260 + tmp1263;
    let tmp1265 = tmp1248 + tmp1249 + tmp1250 + tmp1259 + tmp1264;
    let tmp1266 = (1_f64 / 2.0) * dRdx0_const[(1, 0)];
    let tmp1267 = (1_f64 / 2.0) * tmp852;
    let tmp1268 = (1_f64 / 2.0) * tmp853;
    let tmp1269 = -1.0 * tmp1266 - 1.0 * tmp1267 - 1.0 * tmp1268;
    let tmp1270 = grad_phi_node_i[0] * tmp558;
    let tmp1271 = tmp186 * tmp475;
    let tmp1272 = grad_phi_node_i[0] * tmp1271;
    let tmp1273 = tmp405 * tmp824;
    let tmp1274 = tmp1035 * tmp476;
    let tmp1275 = tmp476 * tmp596;
    let tmp1276 = tmp512 * tmp839;
    let tmp1277 = 4.0 * tmp1051;
    let tmp1278 = 4.0 * tmp1046;
    let tmp1279 = 4.0 * tmp906;
    let tmp1280 = 4.0 * tmp925;
    let tmp1281 = grad_phi_node_i[0] * tmp179;
    let tmp1282 = tmp1281 * tmp522;
    let tmp1283 = (1_f64 / 2.0) * dRdx2_const[(2, 1)];
    let tmp1284 = (1_f64 / 2.0) * tmp319;
    let tmp1285 = (1_f64 / 2.0) * tmp321;
    let tmp1286 = tmp1283 + tmp1284 + tmp1285;
    let tmp1287 = 2.0 * dRdx2_const[(0, 1)] + 2.0 * tmp422 + 2.0 * tmp423;
    let tmp1288 = tmp1255 + tmp1257;
    let tmp1289 = tmp1251 + tmp1253;
    let tmp1290 = -1.0 * tmp1249;
    let tmp1291 = -1.0 * tmp1250;
    let tmp1292 = tmp1290 + tmp1291;
    let tmp1293 = tmp1247 + tmp1292;
    let tmp1294 = tmp1264 + tmp1288 + tmp1289 + tmp1293;
    let tmp1295 = 2.0 * dRdx2_const[(0, 2)] + 2.0 * tmp989 + 2.0 * tmp990;
    let tmp1296 = -1.0 * tmp1232;
    let tmp1297 = tmp1234 + tmp1296;
    let tmp1298 = tmp1227 + tmp1230;
    let tmp1299 = tmp1297 + tmp1298;
    let tmp1300 = tmp1239 + tmp1243;
    let tmp1301 = -1.0 * tmp1237;
    let tmp1302 = tmp1241 + tmp1301;
    let tmp1303 = tmp1300 + tmp1302;
    let tmp1304 = tmp1223 + tmp1226 + tmp1299 + tmp1303;
    let tmp1305 = 2.0 * dRdx2_const[(2, 2)] + 2.0 * tmp1003 + 2.0 * tmp1004;
    let tmp1306 = (1_f64 / 2.0) * dRdx2_const[(1, 2)];
    let tmp1307 = (1_f64 / 2.0) * tmp1018;
    let tmp1308 = (1_f64 / 2.0) * tmp1020;
    let tmp1309 = -1.0 * tmp1306 - 1.0 * tmp1307 - 1.0 * tmp1308;
    let tmp1310 = grad_phi_node_i[2] * tmp553;
    let tmp1311 = grad_phi_node_i[2] * tmp183;
    let tmp1312 = tmp1311 * tmp527;
    let tmp1313 = tmp405 * tmp974;
    let tmp1314 = tmp475 * tmp830;
    let tmp1315 = tmp517 * tmp826;
    let tmp1316 = tmp591 * tmp839;
    let tmp1317 = 4.0 * tmp901;
    let tmp1318 = grad_phi_node_i[2] * tmp179;
    let tmp1319 = tmp1318 * tmp494;
    let tmp1320 = tmp1213 * tmp1295
        + tmp1245 * tmp1305
        + tmp1279 * tmp961
        + tmp1280 * tmp970
        + tmp1286 * tmp777
        + tmp1287 * tmp1294
        + tmp1304 * tmp744
        + tmp1309 * tmp788
        + tmp1310
        + tmp1312
        - 1.0 * tmp1313
        + tmp1314 * tmp513
        + tmp1315 * tmp513
        + tmp1316 * tmp513
        + tmp1317 * tmp189
        - 1.0 * tmp1319
        + tmp552 * tmp871
        + tmp600 * tmp823;
    let tmp1321 = tmp1166 - 1.0 * tmp1168 - 1.0 * tmp1170 + tmp1172 + tmp1173 - 1.0 * tmp1175 - 1.0 * tmp1177 + tmp1178
        - 1.0 * tmp1179
        + tmp1180
        + tmp1181
        - 1.0 * tmp1182;
    let tmp1322 = tmp1306 + tmp1307 + tmp1308;
    let tmp1323 = -1.0 * tmp1260;
    let tmp1324 = -1.0 * tmp1253;
    let tmp1325 = tmp1252 + tmp1324;
    let tmp1326 = -1.0 * tmp1257;
    let tmp1327 = tmp1256 + tmp1326;
    let tmp1328 = -1.0 * tmp1262;
    let tmp1329 = tmp1249 + tmp1328;
    let tmp1330 = -1.0 * tmp1261;
    let tmp1331 = tmp1250 + tmp1330;
    let tmp1332 = tmp1329 + tmp1331;
    let tmp1333 = tmp1248 + tmp1323 + tmp1325 + tmp1327 + tmp1332;
    let tmp1334 = -1.0 * tmp1189;
    let tmp1335 = tmp1199 + tmp1218;
    let tmp1336 = -1.0 * tmp1192;
    let tmp1337 = tmp1193 + tmp1336;
    let tmp1338 = -1.0 * tmp1205;
    let tmp1339 = tmp1207 + tmp1210 + tmp1215 + tmp1338;
    let tmp1340 = tmp1190 + tmp1334 + tmp1335 + tmp1337 + tmp1339;
    let tmp1341 = tmp1228 + tmp1229;
    let tmp1342 = -1.0 * tmp1234;
    let tmp1343 = tmp1232 + tmp1342;
    let tmp1344 = tmp1341 + tmp1343;
    let tmp1345 = -1.0 * tmp1241;
    let tmp1346 = tmp1237 + tmp1345;
    let tmp1347 = tmp1240 + tmp1242;
    let tmp1348 = tmp1346 + tmp1347;
    let tmp1349 = tmp1224 + tmp1225 + tmp1344 + tmp1348;
    let tmp1350 = tmp1296 + tmp1342;
    let tmp1351 = tmp1227 + tmp1229;
    let tmp1352 = tmp1350 + tmp1351;
    let tmp1353 = tmp1239 + tmp1242 + tmp1301 + tmp1345;
    let tmp1354 = tmp1223 + tmp1225 + tmp1352 + tmp1353;
    let tmp1355 = -1.0 * tmp1283 - 1.0 * tmp1284 - 1.0 * tmp1285;
    let tmp1356 = tmp517 * tmp960;
    let tmp1357 = tmp513 * tmp958;
    let tmp1358 = tmp475 * tmp964;
    let tmp1359 = 4.0 * tmp996;
    let tmp1360 = 4.0 * tmp982;
    let tmp1361 = 4.0 * tmp1001;
    let tmp1362 = 4.0 * tmp987;
    let tmp1363 = tmp1266 + tmp1267 + tmp1268;
    let tmp1364 = tmp1207 + tmp1209;
    let tmp1365 = tmp1211 + tmp1338;
    let tmp1366 = tmp1199 + tmp1202;
    let tmp1367 = tmp1194 + tmp1336;
    let tmp1368 = tmp1366 + tmp1367;
    let tmp1369 = tmp1191 + tmp1334 + tmp1364 + tmp1365 + tmp1368;
    let tmp1370 = tmp1255 + tmp1326;
    let tmp1371 = tmp1251 + tmp1324;
    let tmp1372 = tmp1370 + tmp1371;
    let tmp1373 = tmp1293 + tmp1323 + tmp1328 + tmp1330 + tmp1372;
    let tmp1374 = -1.0 * tmp1184 - 1.0 * tmp1185 - 1.0 * tmp1186;
    let tmp1375 = tmp596 * tmp960;
    let tmp1376 = tmp512 * tmp958;
    let tmp1377 = tmp476 * tmp834;
    let tmp1378 = 4.0 * tmp916;
    let tmp1379 = tmp1188 * tmp1340 + tmp1222 * tmp1354 + tmp1246 * tmp1373 - 1.0 * tmp1270 - 1.0 * tmp1272
        + tmp1273
        + tmp1282
        + tmp1359 * tmp836
        + tmp1360 * tmp846
        + tmp1363 * tmp788
        + tmp1369 * tmp716
        + tmp1374 * tmp784
        + tmp1375 * tmp476
        + tmp1376 * tmp476
        + tmp1377 * tmp527
        + tmp1378 * tmp184
        + tmp547 * tmp884
        + tmp584 * tmp818;
    let tmp1380 = tmp1107 * tmp1174;
    let tmp1381 = tmp1104 * tmp1167;
    let tmp1382 = tmp1176 * tmp785;
    let tmp1383 = tmp1171 * tmp811;
    let tmp1384 = tmp1165 * tmp382;
    let tmp1385 = tmp1174 * tmp920;
    let tmp1386 = tmp1005 * tmp1167;
    let tmp1387 = tmp1169 * tmp368;
    let tmp1388 = tmp843 * tmp888;
    let tmp1389 = tmp312 * tmp848;
    let tmp1390 = tmp843 * tmp991;
    let tmp1391 = tmp278 * tmp845;
    let tmp1392 = tmp1380 + tmp1381 - 1.0 * tmp1382 - 1.0 * tmp1383 - 1.0 * tmp1384 + tmp1385 + tmp1386 - 1.0 * tmp1387
        + tmp1388
        - 1.0 * tmp1389
        + tmp1390
        - 1.0 * tmp1391;
    let tmp1393 = (1_f64 / 2.0) * dRdx1_const[(0, 1)];
    let tmp1394 = (1_f64 / 2.0) * tmp886;
    let tmp1395 = (1_f64 / 2.0) * tmp887;
    let tmp1396 = tmp1393 + tmp1394 + tmp1395;
    let tmp1397 = 2.0 * dRdx1_const[(1, 1)] + 2.0 * tmp918 + 2.0 * tmp919;
    let tmp1398 = 2.0 * dRdx1_const[(2, 0)] + 2.0 * tmp418 + 2.0 * tmp419;
    let tmp1399 = 2.0 * dRdx1_const[(2, 1)] + 2.0 * tmp945 + 2.0 * tmp947;
    let tmp1400 = (1_f64 / 2.0) * dRdx1_const[(1, 0)];
    let tmp1401 = (1_f64 / 2.0) * tmp309;
    let tmp1402 = (1_f64 / 2.0) * tmp311;
    let tmp1403 = -1.0 * tmp1400 - 1.0 * tmp1401 - 1.0 * tmp1402;
    let tmp1404 = grad_phi_node_i[1] * tmp558;
    let tmp1405 = grad_phi_node_i[1] * tmp1271;
    let tmp1406 = grad_phi_node_i[1] * tmp405;
    let tmp1407 = tmp121 * tmp1406;
    let tmp1408 = tmp495 * tmp583;
    let tmp1409 = grad_phi_node_i[1] * tmp527;
    let tmp1410 = grad_phi_node_i[1] * tmp179;
    let tmp1411 = tmp1410 * tmp522;
    let tmp1412 = (1_f64 / 2.0) * dRdx2_const[(0, 2)];
    let tmp1413 = (1_f64 / 2.0) * tmp989;
    let tmp1414 = (1_f64 / 2.0) * tmp990;
    let tmp1415 = tmp1412 + tmp1413 + tmp1414;
    let tmp1416 = 2.0 * dRdx2_const[(1, 0)] + 2.0 * tmp426 + 2.0 * tmp427;
    let tmp1417 = 2.0 * dRdx2_const[(1, 2)] + 2.0 * tmp1018 + 2.0 * tmp1020;
    let tmp1418 = (1_f64 / 2.0) * dRdx2_const[(2, 0)];
    let tmp1419 = (1_f64 / 2.0) * tmp273;
    let tmp1420 = (1_f64 / 2.0) * tmp277;
    let tmp1421 = -1.0 * tmp1418 - 1.0 * tmp1419 - 1.0 * tmp1420;
    let tmp1422 = grad_phi_node_i[2] * tmp555;
    let tmp1423 = tmp186 * tmp517;
    let tmp1424 = grad_phi_node_i[2] * tmp1423;
    let tmp1425 = grad_phi_node_i[2] * tmp556;
    let tmp1426 = tmp1041 * tmp494;
    let tmp1427 = tmp1035 * tmp591;
    let tmp1428 = tmp522 * tmp968;
    let tmp1429 = tmp596 * tmp964;
    let tmp1430 = tmp183 * tmp512;
    let tmp1431 = grad_phi_node_i[2] * tmp1430;
    let tmp1432 = tmp1221 * tmp1417
        + tmp1265 * tmp1305
        + tmp1277 * tmp974
        + tmp1278 * tmp970
        + tmp1333 * tmp729
        + tmp1349 * tmp1416
        + tmp1361 * tmp190
        + tmp1362 * tmp161
        + tmp1415 * tmp784
        + tmp1421 * tmp777
        + tmp1422
        + tmp1424
        - 1.0 * tmp1425
        + tmp1426 * tmp513
        + tmp1427 * tmp513
        + tmp1428 * tmp513
        + tmp1429 * tmp513
        - 1.0 * tmp1431;
    let tmp1433 = -1.0 * tmp1380 - 1.0 * tmp1381 + tmp1382 + tmp1383 + tmp1384 - 1.0 * tmp1385 - 1.0 * tmp1386
        + tmp1387
        - 1.0 * tmp1388
        + tmp1389
        - 1.0 * tmp1390
        + tmp1391;
    let tmp1434 = tmp1418 + tmp1419 + tmp1420;
    let tmp1435 = -1.0 * tmp1412 - 1.0 * tmp1413 - 1.0 * tmp1414;
    let tmp1436 = tmp522 * tmp823;
    let tmp1437 = 2.0 * tmp834;
    let tmp1438 = grad_phi_node_i[2] * tmp1437;
    let tmp1439 = tmp494 * tmp818;
    let tmp1440 = 4.0 * tmp871;
    let tmp1441 = 4.0 * tmp884;
    let tmp1442 = tmp1400 + tmp1401 + tmp1402;
    let tmp1443 = -1.0 * tmp1393 - 1.0 * tmp1394 - 1.0 * tmp1395;
    let tmp1444 = 4.0 * tmp841;
    let tmp1445 = tmp1340 * tmp740
        + tmp1354 * tmp1398
        + tmp1359 * tmp167
        + tmp1360 * tmp188
        + tmp1369 * tmp1397
        + tmp1373 * tmp1399
        + tmp1375 * tmp495
        + tmp1376 * tmp495
        + tmp1378 * tmp850
        - 1.0 * tmp1404
        - 1.0 * tmp1405
        + tmp1407
        + tmp1408 * tmp818
        + tmp1409 * tmp1437
        + tmp1411
        + tmp1442 * tmp788
        + tmp1443 * tmp784
        + tmp1444 * tmp884;
    let tmp1446 = tmp1165 * tmp778;
    let tmp1447 = tmp1174 * tmp810;
    let tmp1448 = tmp1169 * tmp888;
    let tmp1449 = tmp1176 * tmp854;
    let tmp1450 = tmp1167 * tmp376;
    let tmp1451 = tmp1169 * tmp910;
    let tmp1452 = tmp1176 * tmp920;
    let tmp1453 = tmp1171 * tmp359;
    let tmp1454 = tmp300 * tmp843;
    let tmp1455 = tmp845 * tmp977;
    let tmp1456 = tmp254 * tmp848;
    let tmp1457 = tmp1104 * tmp845;
    let tmp1458 = -1.0 * tmp1446 - 1.0 * tmp1447 + tmp1448 + tmp1449 - 1.0 * tmp1450 + tmp1451 + tmp1452
        - 1.0 * tmp1453
        - 1.0 * tmp1454
        + tmp1455
        - 1.0 * tmp1456
        + tmp1457;
    let tmp1459 = (1_f64 / 2.0) * dRdx0_const[(2, 0)];
    let tmp1460 = (1_f64 / 2.0) * tmp933;
    let tmp1461 = (1_f64 / 2.0) * tmp935;
    let tmp1462 = tmp1459 + tmp1460 + tmp1461;
    let tmp1463 = 2.0 * dRdx0_const[(1, 0)] + 2.0 * tmp852 + 2.0 * tmp853;
    let tmp1464 = 2.0 * dRdx0_const[(1, 2)] + 2.0 * tmp390 + 2.0 * tmp391;
    let tmp1465 = (1_f64 / 2.0) * dRdx0_const[(0, 2)];
    let tmp1466 = (1_f64 / 2.0) * tmp297;
    let tmp1467 = (1_f64 / 2.0) * tmp299;
    let tmp1468 = -1.0 * tmp1465 - 1.0 * tmp1466 - 1.0 * tmp1467;
    let tmp1469 = grad_phi_node_i[0] * tmp556;
    let tmp1470 = grad_phi_node_i[0] * tmp1430;
    let tmp1471 = grad_phi_node_i[0] * tmp555;
    let tmp1472 = grad_phi_node_i[0] * tmp1423;
    let tmp1473 = (1_f64 / 2.0) * dRdx1_const[(2, 1)];
    let tmp1474 = (1_f64 / 2.0) * tmp945;
    let tmp1475 = (1_f64 / 2.0) * tmp947;
    let tmp1476 = tmp1473 + tmp1474 + tmp1475;
    let tmp1477 = 2.0 * dRdx1_const[(0, 1)] + 2.0 * tmp886 + 2.0 * tmp887;
    let tmp1478 = 2.0 * dRdx1_const[(0, 2)] + 2.0 * tmp413 + 2.0 * tmp414;
    let tmp1479 = (1_f64 / 2.0) * dRdx1_const[(1, 2)];
    let tmp1480 = (1_f64 / 2.0) * tmp248;
    let tmp1481 = (1_f64 / 2.0) * tmp253;
    let tmp1482 = -1.0 * tmp1479 - 1.0 * tmp1480 - 1.0 * tmp1481;
    let tmp1483 = grad_phi_node_i[1] * tmp553;
    let tmp1484 = tmp1409 * tmp183;
    let tmp1485 = tmp1406 * tmp83;
    let tmp1486 = tmp1410 * tmp494;
    let tmp1487 = tmp1213 * tmp1478
        + tmp1245 * tmp722
        + tmp1279 * tmp187
        + tmp1280 * tmp154
        + tmp1294 * tmp1477
        + tmp1304 * tmp1397
        + tmp1314 * tmp495
        + tmp1315 * tmp495
        + tmp1316 * tmp495
        + tmp1317 * tmp831
        + tmp1408 * tmp823
        + tmp1444 * tmp871
        + tmp1476 * tmp777
        + tmp1482 * tmp788
        + tmp1483
        + tmp1484
        - 1.0 * tmp1485
        - 1.0 * tmp1486;
    let tmp1488 =
        tmp1446 + tmp1447 - 1.0 * tmp1448 - 1.0 * tmp1449 + tmp1450 - 1.0 * tmp1451 - 1.0 * tmp1452 + tmp1453 + tmp1454
            - 1.0 * tmp1455
            + tmp1456
            - 1.0 * tmp1457;
    let tmp1489 = tmp1479 + tmp1480 + tmp1481;
    let tmp1490 = -1.0 * tmp1473 - 1.0 * tmp1474 - 1.0 * tmp1475;
    let tmp1491 = tmp1465 + tmp1466 + tmp1467;
    let tmp1492 = -1.0 * tmp1459 - 1.0 * tmp1460 - 1.0 * tmp1461;
    let tmp1493 = tmp1188 * tmp1333
        + tmp1221 * tmp1464
        + tmp1265 * tmp736
        + tmp1277 * tmp180
        + tmp1278 * tmp164
        + tmp1349 * tmp1463
        + tmp1361 * tmp824
        + tmp1362 * tmp836
        + tmp1426 * tmp476
        + tmp1427 * tmp476
        + tmp1428 * tmp476
        + tmp1429 * tmp476
        - 1.0 * tmp1469
        - 1.0 * tmp1470
        + tmp1471
        + tmp1472
        + tmp1491 * tmp784
        + tmp1492 * tmp777;
    let tmp1494 = 1.0 * a2;
    let tmp1495 = tmp1033 + tmp1063;
    let tmp1496 = tmp1089 + tmp955;
    let tmp1497 = tmp1137 + tmp1163;
    let tmp1498 = -1.0 * tmp1005 * tmp809 - 1.0 * tmp1107 * tmp179 + tmp179 * tmp854 + tmp183 * tmp322
        - 1.0 * tmp186 * tmp212
        - 1.0 * tmp345 * tmp797
        + tmp387 * tmp807
        - 1.0 * tmp748 * tmp801
        + tmp795 * tmp910
        - 1.0 * tmp795 * tmp991
        + tmp803 * tmp812
        + tmp809 * tmp977;
    let tmp1499 = tmp1005 * tmp801 - 1.0 * tmp1104 * tmp801 + tmp1107 * tmp797 + tmp179 * tmp312
        - 1.0 * tmp183 * tmp278
        - 1.0 * tmp186 * tmp888
        + tmp186 * tmp991
        - 1.0 * tmp368 * tmp803
        + tmp382 * tmp795
        - 1.0 * tmp785 * tmp807
        - 1.0 * tmp797 * tmp920
        + tmp809 * tmp811;
    let tmp1500 = tmp1104 * tmp183 - 1.0 * tmp179 * tmp254 - 1.0 * tmp183 * tmp977 + tmp186 * tmp300
        - 1.0 * tmp359 * tmp809
        + tmp376 * tmp801
        - 1.0 * tmp778 * tmp795
        + tmp797 * tmp810
        + tmp803 * tmp888
        - 1.0 * tmp803 * tmp910
        - 1.0 * tmp807 * tmp854
        + tmp807 * tmp920;
    let tmp1501 = 0.666666666666667 * tmp261;
    let tmp1502 = 0.666666666666667 * tmp255;
    let tmp1503 = -1.0 * tmp1502;
    let tmp1504 = 0.333333333333333 * tmp857;
    let tmp1505 = tmp1504 * tmp258;
    let tmp1506 = q_old[3] * tmp1504;
    let tmp1507 = tmp1506 * tmp16;
    let tmp1508 = tmp1505 + tmp1507;
    let tmp1509 = tmp13 * tmp1504;
    let tmp1510 = tmp1509 * tmp256;
    let tmp1511 = -1.0 * tmp1510;
    let tmp1512 = tmp1509 * tmp263;
    let tmp1513 = -1.0 * tmp1512;
    let tmp1514 = tmp1511 + tmp1513;
    let tmp1515 = tmp1508 + tmp1514;
    let tmp1516 = tmp1501 + tmp1503 + tmp1515;
    let tmp1517 = tmp1516 * tmp368;
    let tmp1518 = -1.0 * tmp1507;
    let tmp1519 = tmp1512 + tmp1518;
    let tmp1520 = -1.0 * tmp1501 + tmp1519;
    let tmp1521 = tmp1503 + tmp1505 + tmp1511 + tmp1520;
    let tmp1522 = tmp1521 * tmp387;
    let tmp1523 = 0.666666666666667 * tmp280;
    let tmp1524 = 0.666666666666667 * tmp279;
    let tmp1525 = tmp1509 * tmp281;
    let tmp1526 = tmp1504 * tmp290;
    let tmp1527 = -1.0 * tmp1526;
    let tmp1528 = tmp1525 + tmp1527;
    let tmp1529 = tmp1504 * tmp288;
    let tmp1530 = tmp1509 * tmp283;
    let tmp1531 = -1.0 * tmp1530;
    let tmp1532 = tmp1529 + tmp1531;
    let tmp1533 = tmp1528 + tmp1532;
    let tmp1534 = tmp1523 + tmp1524 + tmp1533;
    let tmp1535 = tmp1534 * tmp278;
    let tmp1536 = 0.666666666666667 * tmp226;
    let tmp1537 = 0.666666666666667 * tmp224;
    let tmp1538 = -1.0 * tmp1537;
    let tmp1539 = tmp1504 * tmp237;
    let tmp1540 = -1.0 * tmp1539;
    let tmp1541 = tmp1509 * tmp227;
    let tmp1542 = -1.0 * tmp1541;
    let tmp1543 = tmp1540 + tmp1542;
    let tmp1544 = tmp1509 * tmp231;
    let tmp1545 = -1.0 * tmp1544;
    let tmp1546 = tmp14 * tmp1506;
    let tmp1547 = -1.0 * tmp1546;
    let tmp1548 = tmp1545 + tmp1547;
    let tmp1549 = tmp1536 + tmp1538 + tmp1543 + tmp1548;
    let tmp1550 = tmp1549 * tmp322;
    let tmp1551 = (1_f64 / 3.0) * dRdx2_const[(0, 1)];
    let tmp1552 = (1_f64 / 3.0) * tmp422;
    let tmp1553 = (1_f64 / 3.0) * tmp423;
    let tmp1554 = -1.0 * tmp1551 - 1.0 * tmp1552 - 1.0 * tmp1553;
    let tmp1555 = tmp1554 * tmp416;
    let tmp1556 = (1_f64 / 3.0) * dRdx2_const[(1, 0)];
    let tmp1557 = (1_f64 / 3.0) * tmp426;
    let tmp1558 = (1_f64 / 3.0) * tmp427;
    let tmp1559 = tmp1556 + tmp1557 + tmp1558;
    let tmp1560 = tmp1559 * tmp405;
    let tmp1561 = 0.666666666666667 * tmp65;
    let tmp1562 = 0.666666666666667 * tmp66;
    let tmp1563 = -1.0 * tmp1562;
    let tmp1564 = tmp1561 + tmp1563;
    let tmp1565 = tmp1564 * tmp161;
    let tmp1566 = 0.666666666666667 * tmp86;
    let tmp1567 = 0.666666666666667 * tmp88;
    let tmp1568 = tmp1566 + tmp1567;
    let tmp1569 = tmp107 * tmp1568;
    let tmp1570 = 0.666666666666667 * tmp15;
    let tmp1571 = -1.0 * tmp1570;
    let tmp1572 = 0.666666666666667 * tmp17;
    let tmp1573 = tmp1571 + tmp1572;
    let tmp1574 = tmp129 * tmp1573;
    let tmp1575 = -1.0 * tmp1561;
    let tmp1576 = tmp1563 + tmp1575;
    let tmp1577 = tmp1576 * tmp170;
    let tmp1578 = (1_f64 / 3.0) * tmp435;
    let tmp1579 = (1_f64 / 3.0) * tmp178;
    let tmp1580 = tmp1579 * tmp190;
    let tmp1581 = 0.333333333333333 * tmp224;
    let tmp1582 = 0.333333333333333 * tmp226;
    let tmp1583 = 0.166666666666667 * tmp857;
    let tmp1584 = q_old[3] * tmp1583;
    let tmp1585 = tmp14 * tmp1584;
    let tmp1586 = tmp1583 * tmp237;
    let tmp1587 = -1.0 * tmp1586;
    let tmp1588 = tmp1585 + tmp1587;
    let tmp1589 = tmp13 * tmp1583;
    let tmp1590 = tmp1589 * tmp227;
    let tmp1591 = tmp1589 * tmp231;
    let tmp1592 = -1.0 * tmp1591;
    let tmp1593 = tmp1590 + tmp1592;
    let tmp1594 = tmp1588 + tmp1593;
    let tmp1595 = tmp1581 + tmp1582 + tmp1594;
    let tmp1596 = tmp1595 * tmp382;
    let tmp1597 = 0.333333333333333 * tmp279;
    let tmp1598 = tmp1583 * tmp288;
    let tmp1599 = tmp1583 * tmp290;
    let tmp1600 = -1.0 * tmp1599;
    let tmp1601 = tmp1598 + tmp1600;
    let tmp1602 = tmp1597 + tmp1601;
    let tmp1603 = 0.333333333333333 * tmp280;
    let tmp1604 = tmp1589 * tmp283;
    let tmp1605 = tmp1589 * tmp281;
    let tmp1606 = -1.0 * tmp1605;
    let tmp1607 = tmp1604 + tmp1606;
    let tmp1608 = -1.0 * tmp1603 + tmp1607;
    let tmp1609 = tmp1602 + tmp1608;
    let tmp1610 = tmp1609 * tmp312;
    let tmp1611 = 0.333333333333333 * tmp255;
    let tmp1612 = -1.0 * tmp1611;
    let tmp1613 = tmp1583 * tmp258;
    let tmp1614 = tmp1589 * tmp256;
    let tmp1615 = -1.0 * tmp1614;
    let tmp1616 = 0.333333333333333 * tmp261;
    let tmp1617 = tmp1589 * tmp263;
    let tmp1618 = tmp1584 * tmp16;
    let tmp1619 = -1.0 * tmp1618;
    let tmp1620 = tmp1617 + tmp1619;
    let tmp1621 = -1.0 * tmp1616 + tmp1620;
    let tmp1622 = tmp1612 + tmp1613 + tmp1615 + tmp1621;
    let tmp1623 = tmp1622 * tmp254;
    let tmp1624 = -1.0 * tmp1590;
    let tmp1625 = -1.0 * tmp1585;
    let tmp1626 = -1.0 * tmp1581;
    let tmp1627 = tmp1582 + tmp1587 + tmp1592 + tmp1624 + tmp1625 + tmp1626;
    let tmp1628 = tmp1627 * tmp359;
    let tmp1629 = -1_f64 / 6.0 * dRdx1_const[(0, 2)] - 1_f64 / 6.0 * tmp413 - 1_f64 / 6.0 * tmp414;
    let tmp1630 = tmp1629 * tmp416;
    let tmp1631 = (1_f64 / 6.0) * dRdx1_const[(2, 0)] + (1_f64 / 6.0) * tmp418 + (1_f64 / 6.0) * tmp419;
    let tmp1632 = tmp1631 * tmp411;
    let tmp1633 = 0.333333333333333 * tmp15;
    let tmp1634 = 0.333333333333333 * tmp17;
    let tmp1635 = tmp1633 + tmp1634;
    let tmp1636 = tmp1635 * tmp167;
    let tmp1637 = 0.333333333333333 * tmp86;
    let tmp1638 = 0.333333333333333 * tmp88;
    let tmp1639 = -1.0 * tmp1638;
    let tmp1640 = tmp1637 + tmp1639;
    let tmp1641 = tmp122 * tmp1640;
    let tmp1642 = -1.0 * tmp1633;
    let tmp1643 = tmp1634 + tmp1642;
    let tmp1644 = tmp154 * tmp1643;
    let tmp1645 = 0.333333333333333 * tmp65;
    let tmp1646 = -1.0 * tmp1645;
    let tmp1647 = 0.333333333333333 * tmp66;
    let tmp1648 = -1.0 * tmp1647;
    let tmp1649 = tmp1646 + tmp1648;
    let tmp1650 = tmp1649 * tmp84;
    let tmp1651 = (1_f64 / 6.0) * tmp182;
    let tmp1652 = tmp1596
        + tmp1610
        + tmp1623
        + tmp1628
        + tmp1630
        + tmp1632
        + tmp1636
        + tmp1641
        + tmp1644
        + tmp1650
        + tmp1651 * tmp188
        - 1_f64 / 6.0 * tmp431;
    let tmp1653 = -1.0 * tmp1582;
    let tmp1654 = tmp1591 + tmp1624;
    let tmp1655 = tmp1586 + tmp1625;
    let tmp1656 = tmp1654 + tmp1655;
    let tmp1657 = tmp1626 + tmp1653 + tmp1656;
    let tmp1658 = tmp1657 * tmp212;
    let tmp1659 = -1.0 * tmp1617;
    let tmp1660 = tmp1615 + tmp1659;
    let tmp1661 = tmp1613 + tmp1618;
    let tmp1662 = tmp1660 + tmp1661;
    let tmp1663 = tmp1612 + tmp1616 + tmp1662;
    let tmp1664 = tmp1663 * tmp300;
    let tmp1665 = -1.0 * tmp1604;
    let tmp1666 = tmp1605 + tmp1665;
    let tmp1667 = tmp1603 + tmp1666;
    let tmp1668 = -1.0 * tmp1598;
    let tmp1669 = tmp1599 + tmp1668;
    let tmp1670 = -1.0 * tmp1597 + tmp1669;
    let tmp1671 = tmp1667 + tmp1670;
    let tmp1672 = tmp1671 * tmp345;
    let tmp1673 = tmp1602 + tmp1667;
    let tmp1674 = tmp1673 * tmp376;
    let tmp1675 = -1_f64 / 6.0 * dRdx0_const[(2, 1)] - 1_f64 / 6.0 * tmp407 - 1_f64 / 6.0 * tmp408;
    let tmp1676 = tmp1675 * tmp411;
    let tmp1677 = (1_f64 / 6.0) * dRdx0_const[(1, 2)] + (1_f64 / 6.0) * tmp390 + (1_f64 / 6.0) * tmp391;
    let tmp1678 = tmp1677 * tmp405;
    let tmp1679 = tmp1645 + tmp1648;
    let tmp1680 = tmp116 * tmp1679;
    let tmp1681 = tmp1637 + tmp1638;
    let tmp1682 = tmp164 * tmp1681;
    let tmp1683 = -1.0 * tmp1637;
    let tmp1684 = tmp1638 + tmp1683;
    let tmp1685 = tmp147 * tmp1684;
    let tmp1686 = -1.0 * tmp1634;
    let tmp1687 = tmp1642 + tmp1686;
    let tmp1688 = tmp1687 * tmp63;
    let tmp1689 = (1_f64 / 6.0) * tmp178;
    let tmp1690 = -1.0 * tmp1651 * tmp184
        + tmp1658
        + tmp1664
        + tmp1672
        + tmp1674
        + tmp1676
        + tmp1678
        + tmp1680
        + tmp1682
        + tmp1685
        + tmp1688
        + tmp1689 * tmp180;
    let tmp1691 = (1_f64 / 3.0) * tmp431;
    let tmp1692 = (1_f64 / 3.0) * tmp182;
    let tmp1693 = tmp1692 * tmp188;
    let tmp1694 = 2.0 * tmp1596
        + 2.0 * tmp1610
        + 2.0 * tmp1623
        + 2.0 * tmp1628
        + 2.0 * tmp1630
        + 2.0 * tmp1632
        + 2.0 * tmp1636
        + 2.0 * tmp1641
        + 2.0 * tmp1644
        + 2.0 * tmp1650
        - 1.0 * tmp1691
        + tmp1693;
    let tmp1695 = tmp1692 * tmp184;
    let tmp1696 = tmp1579 * tmp180;
    let tmp1697 = 2.0 * tmp1658
        + 2.0 * tmp1664
        + 2.0 * tmp1672
        + 2.0 * tmp1674
        + 2.0 * tmp1676
        + 2.0 * tmp1678
        + 2.0 * tmp1680
        + 2.0 * tmp1682
        + 2.0 * tmp1685
        + 2.0 * tmp1688
        - 1.0 * tmp1695
        + tmp1696;
    let tmp1698 = 2.0 * tmp1517
        + 2.0 * tmp1522
        + 2.0 * tmp1535
        + 2.0 * tmp1550
        + 2.0 * tmp1555
        + 2.0 * tmp1560
        + 2.0 * tmp1565
        + 2.0 * tmp1569
        + 2.0 * tmp1574
        + 2.0 * tmp1577
        + tmp1694
        + tmp1697
        + (2_f64 / 3.0) * tmp432
        - 2_f64 / 3.0 * tmp435;
    let tmp1699 = tmp1539 + tmp1544;
    let tmp1700 = -1.0 * tmp1536 + tmp1699;
    let tmp1701 = tmp1542 + tmp1547;
    let tmp1702 = tmp1538 + tmp1700 + tmp1701;
    let tmp1703 = tmp1702 * tmp382;
    let tmp1704 = -1.0 * tmp1529;
    let tmp1705 = tmp1526 + tmp1704;
    let tmp1706 = -1.0 * tmp1524 + tmp1705;
    let tmp1707 = tmp1525 + tmp1531;
    let tmp1708 = tmp1523 + tmp1706 + tmp1707;
    let tmp1709 = tmp1708 * tmp312;
    let tmp1710 = -1.0 * tmp1505;
    let tmp1711 = tmp1510 + tmp1710;
    let tmp1712 = tmp1502 + tmp1711;
    let tmp1713 = tmp1501 + tmp1507 + tmp1513 + tmp1712;
    let tmp1714 = tmp1713 * tmp254;
    let tmp1715 = tmp1537 + tmp1541 + tmp1546;
    let tmp1716 = tmp1700 + tmp1715;
    let tmp1717 = tmp1716 * tmp359;
    let tmp1718 = (1_f64 / 3.0) * dRdx1_const[(2, 0)];
    let tmp1719 = (1_f64 / 3.0) * tmp418;
    let tmp1720 = (1_f64 / 3.0) * tmp419;
    let tmp1721 = -1.0 * tmp1718 - 1.0 * tmp1719 - 1.0 * tmp1720;
    let tmp1722 = tmp1721 * tmp411;
    let tmp1723 = (1_f64 / 3.0) * dRdx1_const[(0, 2)];
    let tmp1724 = (1_f64 / 3.0) * tmp413;
    let tmp1725 = (1_f64 / 3.0) * tmp414;
    let tmp1726 = tmp1723 + tmp1724 + tmp1725;
    let tmp1727 = tmp1726 * tmp416;
    let tmp1728 = tmp1561 + tmp1562;
    let tmp1729 = tmp1728 * tmp84;
    let tmp1730 = -1.0 * tmp1572;
    let tmp1731 = tmp1570 + tmp1730;
    let tmp1732 = tmp154 * tmp1731;
    let tmp1733 = -1.0 * tmp1566;
    let tmp1734 = tmp1567 + tmp1733;
    let tmp1735 = tmp122 * tmp1734;
    let tmp1736 = tmp1571 + tmp1730;
    let tmp1737 = tmp167 * tmp1736;
    let tmp1738 = -1.0 * tmp1613;
    let tmp1739 = tmp1614 + tmp1738;
    let tmp1740 = tmp1611 + tmp1739;
    let tmp1741 = tmp1621 + tmp1740;
    let tmp1742 = tmp1741 * tmp368;
    let tmp1743 = tmp1616 + tmp1618 + tmp1659 + tmp1740;
    let tmp1744 = tmp1743 * tmp387;
    let tmp1745 = tmp1608 + tmp1670;
    let tmp1746 = tmp1745 * tmp278;
    let tmp1747 = tmp1581 + tmp1585 + tmp1586 + tmp1590 + tmp1591 + tmp1653;
    let tmp1748 = tmp1747 * tmp322;
    let tmp1749 = -1_f64 / 6.0 * dRdx2_const[(1, 0)] - 1_f64 / 6.0 * tmp426 - 1_f64 / 6.0 * tmp427;
    let tmp1750 = tmp1749 * tmp405;
    let tmp1751 = (1_f64 / 6.0) * dRdx2_const[(0, 1)] + (1_f64 / 6.0) * tmp422 + (1_f64 / 6.0) * tmp423;
    let tmp1752 = tmp1751 * tmp416;
    let tmp1753 = tmp1645 + tmp1647;
    let tmp1754 = tmp170 * tmp1753;
    let tmp1755 = tmp1633 + tmp1686;
    let tmp1756 = tmp129 * tmp1755;
    let tmp1757 = tmp1646 + tmp1647;
    let tmp1758 = tmp161 * tmp1757;
    let tmp1759 = tmp1639 + tmp1683;
    let tmp1760 = tmp107 * tmp1759;
    let tmp1761 = -1.0 * tmp1689 * tmp190
        + tmp1742
        + tmp1744
        + tmp1746
        + tmp1748
        + tmp1750
        + tmp1752
        + tmp1754
        + tmp1756
        + tmp1758
        + tmp1760
        + (1_f64 / 6.0) * tmp435;
    let tmp1762 = tmp1578 - 1.0 * tmp1580
        + 2.0 * tmp1742
        + 2.0 * tmp1744
        + 2.0 * tmp1746
        + 2.0 * tmp1748
        + 2.0 * tmp1750
        + 2.0 * tmp1752
        + 2.0 * tmp1754
        + 2.0 * tmp1756
        + 2.0 * tmp1758
        + 2.0 * tmp1760;
    let tmp1763 = tmp1697
        + 2.0 * tmp1703
        + 2.0 * tmp1709
        + 2.0 * tmp1714
        + 2.0 * tmp1717
        + 2.0 * tmp1722
        + 2.0 * tmp1727
        + 2.0 * tmp1729
        + 2.0 * tmp1732
        + 2.0 * tmp1735
        + 2.0 * tmp1737
        + tmp1762
        + (2_f64 / 3.0) * tmp431
        - 2_f64 / 3.0 * tmp434;
    let tmp1764 = tmp1540 + tmp1545;
    let tmp1765 = tmp1536 + tmp1715 + tmp1764;
    let tmp1766 = tmp1765 * tmp212;
    let tmp1767 = tmp1520 + tmp1712;
    let tmp1768 = tmp1767 * tmp300;
    let tmp1769 = -1.0 * tmp1525;
    let tmp1770 = tmp1530 + tmp1769;
    let tmp1771 = -1.0 * tmp1523 + tmp1770;
    let tmp1772 = tmp1527 + tmp1529;
    let tmp1773 = tmp1524 + tmp1771 + tmp1772;
    let tmp1774 = tmp1773 * tmp345;
    let tmp1775 = tmp1706 + tmp1771;
    let tmp1776 = tmp1775 * tmp376;
    let tmp1777 = (1_f64 / 3.0) * dRdx0_const[(1, 2)];
    let tmp1778 = (1_f64 / 3.0) * tmp390;
    let tmp1779 = (1_f64 / 3.0) * tmp391;
    let tmp1780 = -1.0 * tmp1777 - 1.0 * tmp1778 - 1.0 * tmp1779;
    let tmp1781 = tmp1780 * tmp405;
    let tmp1782 = (1_f64 / 3.0) * dRdx0_const[(2, 1)];
    let tmp1783 = (1_f64 / 3.0) * tmp407;
    let tmp1784 = (1_f64 / 3.0) * tmp408;
    let tmp1785 = tmp1782 + tmp1783 + tmp1784;
    let tmp1786 = tmp1785 * tmp411;
    let tmp1787 = tmp1570 + tmp1572;
    let tmp1788 = tmp1787 * tmp63;
    let tmp1789 = -1.0 * tmp1567;
    let tmp1790 = tmp1566 + tmp1789;
    let tmp1791 = tmp147 * tmp1790;
    let tmp1792 = tmp1562 + tmp1575;
    let tmp1793 = tmp116 * tmp1792;
    let tmp1794 = tmp1733 + tmp1789;
    let tmp1795 = tmp164 * tmp1794;
    let tmp1796 = tmp1694
        + tmp1762
        + 2.0 * tmp1766
        + 2.0 * tmp1768
        + 2.0 * tmp1774
        + 2.0 * tmp1776
        + 2.0 * tmp1781
        + 2.0 * tmp1786
        + 2.0 * tmp1788
        + 2.0 * tmp1791
        + 2.0 * tmp1793
        + 2.0 * tmp1795
        + (2_f64 / 3.0) * tmp430
        - 2_f64 / 3.0 * tmp433;
    let tmp1797 = (1_f64 / 3.0) * tmp794;
    let tmp1798 = (1_f64 / 3.0) * tmp808;
    let tmp1799 = (1_f64 / 3.0) * tmp802;
    let tmp1800 = (1_f64 / 3.0) * tmp806;
    let tmp1801 = tmp254 * tmp796;
    let tmp1802 = tmp359 * tmp800;
    let tmp1803 = tmp185 * tmp778;
    let tmp1804 = (1_f64 / 6.0) * tmp802;
    let tmp1805 = (1_f64 / 6.0) * tmp806;
    let tmp1806 = tmp1651 * tmp811 - 1_f64 / 6.0 * tmp1801 - 1_f64 / 6.0 * tmp1802 - 1_f64 / 6.0 * tmp1803
        + tmp1804 * tmp382
        + tmp1805 * tmp312;
    let tmp1807 = (1_f64 / 6.0) * tmp794;
    let tmp1808 = (1_f64 / 6.0) * tmp808;
    let tmp1809 = tmp387 * tmp796;
    let tmp1810 = tmp322 * tmp800;
    let tmp1811 = tmp185 * tmp812;
    let tmp1812 = -1.0 * tmp1689 * tmp785 - 1.0 * tmp1807 * tmp368 - 1.0 * tmp1808 * tmp278
        + (1_f64 / 6.0) * tmp1809
        + (1_f64 / 6.0) * tmp1810
        + (1_f64 / 6.0) * tmp1811;
    let tmp1813 = -1.0 * tmp1579 * tmp810 + tmp1692 * tmp748 - 1.0 * tmp1797 * tmp300 - 1.0 * tmp1798 * tmp376
        + tmp1799 * tmp212
        + tmp1800 * tmp345
        + tmp1806
        + tmp1812;
    let tmp1814 = 0.666666666666667 * tmp626;
    let tmp1815 = 0.333333333333333 * tmp1197;
    let tmp1816 = q_old[2] * tmp1815;
    let tmp1817 = q_old[0] * tmp1816;
    let tmp1818 = 0.666666666666667 * tmp611;
    let tmp1819 = 0.666666666666667 * tmp857;
    let tmp1820 = tmp1819 * tmp225;
    let tmp1821 = q_old[0] * tmp1820;
    let tmp1822 = tmp1819 * tmp222;
    let tmp1823 = tmp1822 * tmp231;
    let tmp1824 = tmp1818 - 1.0 * tmp1821 - 1.0 * tmp1823;
    let tmp1825 = 0.666666666666667 * tmp609;
    let tmp1826 = tmp1822 * tmp227;
    let tmp1827 = q_old[3] * tmp223;
    let tmp1828 = tmp1819 * tmp1827;
    let tmp1829 = 0.666666666666667 * tmp624;
    let tmp1830 = q_old[3] * tmp1815;
    let tmp1831 = q_old[1] * tmp1830;
    let tmp1832 = tmp1829 + tmp1831;
    let tmp1833 = tmp1825 + tmp1826 + tmp1828 + tmp1832;
    let tmp1834 = 0.666666666666667 * tmp664;
    let tmp1835 = tmp1822 * tmp263;
    let tmp1836 = q_old[3] * tmp1820;
    let tmp1837 = -1.0 * tmp1834 + tmp1835 - 1.0 * tmp1836;
    let tmp1838 = q_old[2] * tmp1830;
    let tmp1839 = tmp1815 * tmp657;
    let tmp1840 = -1.0 * tmp1839;
    let tmp1841 = tmp1838 + tmp1840;
    let tmp1842 = 0.666666666666667 * tmp655;
    let tmp1843 = -1.0 * tmp1842;
    let tmp1844 = 0.666666666666667 * tmp654;
    let tmp1845 = tmp1843 + tmp1844;
    let tmp1846 = 0.666666666666667 * tmp650;
    let tmp1847 = tmp1822 * tmp256;
    let tmp1848 = tmp1819 * tmp647;
    let tmp1849 = tmp1846 + tmp1847 - 1.0 * tmp1848;
    let tmp1850 = 0.666666666666667 * tmp672;
    let tmp1851 = tmp1819 * tmp673;
    let tmp1852 = tmp1819 * tmp675;
    let tmp1853 = tmp1850 + tmp1851 - 1.0 * tmp1852;
    let tmp1854 = 0.666666666666667 * tmp689;
    let tmp1855 = q_old[1] * tmp1816;
    let tmp1856 = -1.0 * tmp1855;
    let tmp1857 = 0.666666666666667 * tmp690;
    let tmp1858 = -1.0 * tmp1857;
    let tmp1859 = q_old[0] * tmp1830;
    let tmp1860 = tmp1858 + tmp1859;
    let tmp1861 = tmp1854 + tmp1856 + tmp1860;
    let tmp1862 = 0.666666666666667 * tmp678;
    let tmp1863 = tmp1822 * tmp283;
    let tmp1864 = tmp1822 * tmp281;
    let tmp1865 = -1.0 * tmp1862 + tmp1863 - 1.0 * tmp1864;
    let tmp1866 = -1.0 * tmp1854;
    let tmp1867 = tmp1855 + tmp1866;
    let tmp1868 = tmp1860 + tmp1867;
    let tmp1869 = -1.0 * tmp1850 - 1.0 * tmp1851 + tmp1852;
    let tmp1870 = -2_f64 / 3.0 * dRdx0_const[(1, 2)] - 2_f64 / 3.0 * tmp390 - 2_f64 / 3.0 * tmp391;
    let tmp1871 = (2_f64 / 3.0) * dRdx0_const[(2, 1)] + (2_f64 / 3.0) * tmp407 + (2_f64 / 3.0) * tmp408;
    let tmp1872 = grad_phi_node_i[0] * tmp553;
    let tmp1873 = 0.333333333333333 * tmp609;
    let tmp1874 = 0.333333333333333 * tmp611;
    let tmp1875 = 0.166666666666667 * tmp1197;
    let tmp1876 = q_old[2] * tmp1875;
    let tmp1877 = q_old[0] * tmp1876;
    let tmp1878 = q_old[3] * tmp1875;
    let tmp1879 = q_old[1] * tmp1878;
    let tmp1880 = tmp1877 + tmp1879;
    let tmp1881 = 0.333333333333333 * tmp624;
    let tmp1882 = 0.333333333333333 * tmp626;
    let tmp1883 = tmp1881 + tmp1882;
    let tmp1884 = tmp1506 * tmp223;
    let tmp1885 = q_old[0] * tmp1504;
    let tmp1886 = tmp1885 * tmp225;
    let tmp1887 = -1.0 * tmp1886;
    let tmp1888 = tmp1884 + tmp1887;
    let tmp1889 = tmp1504 * tmp222;
    let tmp1890 = tmp1889 * tmp227;
    let tmp1891 = tmp1889 * tmp231;
    let tmp1892 = -1.0 * tmp1891;
    let tmp1893 = tmp1890 + tmp1892;
    let tmp1894 = tmp1888 + tmp1893;
    let tmp1895 = 0.333333333333333 * tmp672;
    let tmp1896 = 0.333333333333333 * tmp678;
    let tmp1897 = -1.0 * tmp1896;
    let tmp1898 = tmp1504 * tmp673;
    let tmp1899 = tmp1504 * tmp675;
    let tmp1900 = -1.0 * tmp1899;
    let tmp1901 = tmp1898 + tmp1900;
    let tmp1902 = tmp1889 * tmp283;
    let tmp1903 = tmp1889 * tmp281;
    let tmp1904 = -1.0 * tmp1903;
    let tmp1905 = tmp1902 + tmp1904;
    let tmp1906 = tmp1901 + tmp1905;
    let tmp1907 = q_old[1] * tmp1876;
    let tmp1908 = -1.0 * tmp1907;
    let tmp1909 = q_old[0] * tmp1878;
    let tmp1910 = tmp1908 + tmp1909;
    let tmp1911 = 0.333333333333333 * tmp689;
    let tmp1912 = 0.333333333333333 * tmp690;
    let tmp1913 = -1.0 * tmp1912;
    let tmp1914 = tmp1911 + tmp1913;
    let tmp1915 = tmp1910 + tmp1914;
    let tmp1916 = 0.333333333333333 * tmp650;
    let tmp1917 = -1.0 * tmp1916;
    let tmp1918 = tmp1885 * tmp223;
    let tmp1919 = tmp1889 * tmp256;
    let tmp1920 = -1.0 * tmp1919;
    let tmp1921 = 0.333333333333333 * tmp654;
    let tmp1922 = -1.0 * tmp1921;
    let tmp1923 = 0.333333333333333 * tmp655;
    let tmp1924 = -1.0 * tmp1923;
    let tmp1925 = tmp1922 + tmp1924;
    let tmp1926 = tmp1875 * tmp657;
    let tmp1927 = q_old[3] * tmp1876;
    let tmp1928 = tmp1926 + tmp1927;
    let tmp1929 = tmp1925 + tmp1928;
    let tmp1930 = 0.333333333333333 * tmp664;
    let tmp1931 = tmp1889 * tmp263;
    let tmp1932 = tmp1506 * tmp225;
    let tmp1933 = -1.0 * tmp1932;
    let tmp1934 = tmp1931 + tmp1933;
    let tmp1935 = -1.0 * tmp1930 + tmp1934;
    let tmp1936 = -1.0 * tmp1873;
    let tmp1937 = -1.0 * tmp1890;
    let tmp1938 = -1.0 * tmp1884;
    let tmp1939 = -1.0 * tmp1881;
    let tmp1940 = tmp1882 + tmp1939;
    let tmp1941 = -1.0 * tmp1879;
    let tmp1942 = tmp1877 + tmp1941;
    let tmp1943 = tmp1940 + tmp1942;
    let tmp1944 = -1.0 * tmp1723 - 1.0 * tmp1724 - 1.0 * tmp1725;
    let tmp1945 = tmp1718 + tmp1719 + tmp1720;
    let tmp1946 = (2_f64 / 3.0) * grad_phi_node_i[1];
    let tmp1947 = grad_phi_node_i[1] * tmp1692;
    let tmp1948 = tmp1595 * tmp551
        + tmp1609 * tmp533
        + tmp1622 * tmp530
        + tmp1627 * tmp548
        + tmp1635 * tmp599
        + tmp1640 * tmp523
        + tmp1643 * tmp592
        + tmp1649 * tmp496
        + tmp1944 * tmp784
        + tmp1945 * tmp777
        - 1.0 * tmp1946 * tmp555
        + tmp1946 * tmp556
        + tmp1947 * tmp512
        - 1_f64 / 3.0 * tmp538
        + tmp646 * (tmp1917 + tmp1918 + tmp1920 + tmp1929 + tmp1935)
        + tmp700 * (tmp1895 + tmp1897 + tmp1906 + tmp1915)
        + tmp722 * (tmp1874 + tmp1887 + tmp1892 + tmp1936 + tmp1937 + tmp1938 + tmp1943)
        + tmp740 * (tmp1873 + tmp1874 + tmp1880 + tmp1883 + tmp1894);
    let tmp1949 = -1.0 * tmp1926;
    let tmp1950 = tmp1927 + tmp1949;
    let tmp1951 = tmp1921 + tmp1924;
    let tmp1952 = -1.0 * tmp1918;
    let tmp1953 = tmp1919 + tmp1952;
    let tmp1954 = tmp1916 + tmp1953;
    let tmp1955 = -1.0 * tmp1931;
    let tmp1956 = -1.0 * tmp1927;
    let tmp1957 = tmp1949 + tmp1956;
    let tmp1958 = tmp1921 + tmp1923;
    let tmp1959 = tmp1957 + tmp1958;
    let tmp1960 = -1.0 * tmp1895;
    let tmp1961 = -1.0 * tmp1898;
    let tmp1962 = tmp1899 + tmp1905 + tmp1961;
    let tmp1963 = tmp1907 + tmp1909;
    let tmp1964 = -1.0 * tmp1911;
    let tmp1965 = tmp1913 + tmp1964;
    let tmp1966 = tmp1963 + tmp1965;
    let tmp1967 = -1.0 * tmp1874;
    let tmp1968 = -1.0 * tmp1877;
    let tmp1969 = tmp1879 + tmp1968;
    let tmp1970 = -1.0 * tmp1882;
    let tmp1971 = tmp1881 + tmp1970;
    let tmp1972 = tmp1969 + tmp1971;
    let tmp1973 = -1.0 * tmp1556 - 1.0 * tmp1557 - 1.0 * tmp1558;
    let tmp1974 = tmp1551 + tmp1552 + tmp1553;
    let tmp1975 = grad_phi_node_i[2] * tmp1579;
    let tmp1976 = grad_phi_node_i[2] * tmp558;
    let tmp1977 = tmp1741 * tmp549
        + tmp1743 * tmp552
        + tmp1745 * tmp531
        + tmp1747 * tmp534
        + tmp1753 * tmp600
        + tmp1755 * tmp528
        + tmp1757 * tmp597
        + tmp1759 * tmp514
        + tmp1973 * tmp788
        + tmp1974 * tmp784
        - 1.0 * tmp1975 * tmp522
        + (2_f64 / 3.0) * tmp1976
        + (1_f64 / 3.0) * tmp546
        - 2_f64 / 3.0 * tmp557
        + tmp671 * (tmp1897 + tmp1960 + tmp1962 + tmp1966)
        + tmp706 * (tmp1873 + tmp1884 + tmp1886 + tmp1890 + tmp1891 + tmp1967 + tmp1972)
        + tmp729 * (tmp1935 + tmp1950 + tmp1951 + tmp1954)
        + tmp744 * (tmp1930 + tmp1932 + tmp1954 + tmp1955 + tmp1959);
    let tmp1978 = -1.0 * tmp1651 * tmp748 + tmp1689 * tmp810 - 1.0 * tmp1804 * tmp212 - 1.0 * tmp1805 * tmp345
        + tmp1807 * tmp300
        + tmp1808 * tmp376;
    let tmp1979 = -1.0 * tmp1692 * tmp811 - 1.0 * tmp1799 * tmp382 - 1.0 * tmp1800 * tmp312
        + (1_f64 / 3.0) * tmp1801
        + (1_f64 / 3.0) * tmp1802
        + (1_f64 / 3.0) * tmp1803
        + tmp1812
        + tmp1978;
    let tmp1980 = -1.0 * tmp1829;
    let tmp1981 = -1.0 * tmp1831;
    let tmp1982 = -1.0 * tmp1825 - 1.0 * tmp1826 - 1.0 * tmp1828;
    let tmp1983 = -1.0 * tmp1814;
    let tmp1984 = -1.0 * tmp1817;
    let tmp1985 = tmp1983 + tmp1984;
    let tmp1986 = -1.0 * tmp1818 + tmp1821 + tmp1823 + tmp1985;
    let tmp1987 = tmp1862 - 1.0 * tmp1863 + tmp1864;
    let tmp1988 = -1.0 * tmp1859;
    let tmp1989 = tmp1857 + tmp1867 + tmp1988;
    let tmp1990 = tmp1834 - 1.0 * tmp1835 + tmp1836;
    let tmp1991 = -1.0 * tmp1838;
    let tmp1992 = tmp1840 + tmp1991;
    let tmp1993 = tmp1842 + tmp1844;
    let tmp1994 = tmp1992 + tmp1993;
    let tmp1995 = -2_f64 / 3.0 * dRdx1_const[(2, 0)] - 2_f64 / 3.0 * tmp418 - 2_f64 / 3.0 * tmp419;
    let tmp1996 = (2_f64 / 3.0) * dRdx1_const[(0, 2)] + (2_f64 / 3.0) * tmp413 + (2_f64 / 3.0) * tmp414;
    let tmp1997 = (4_f64 / 3.0) * grad_phi_node_i[1];
    let tmp1998 = tmp1939 + tmp1970;
    let tmp1999 = tmp1941 + tmp1968;
    let tmp2000 = tmp1891 + tmp1937;
    let tmp2001 = tmp1886 + tmp1938;
    let tmp2002 = tmp2000 + tmp2001;
    let tmp2003 = tmp1922 + tmp1923;
    let tmp2004 = tmp1926 + tmp1956;
    let tmp2005 = tmp1918 + tmp1932;
    let tmp2006 = tmp1920 + tmp1955;
    let tmp2007 = tmp2005 + tmp2006;
    let tmp2008 = tmp1912 + tmp1964;
    let tmp2009 = -1.0 * tmp1909;
    let tmp2010 = tmp1907 + tmp2009;
    let tmp2011 = tmp2008 + tmp2010;
    let tmp2012 = tmp1903 + tmp1961;
    let tmp2013 = -1.0 * tmp1902;
    let tmp2014 = tmp1899 + tmp2013;
    let tmp2015 = tmp2012 + tmp2014;
    let tmp2016 = tmp1911 + tmp1912;
    let tmp2017 = tmp1908 + tmp2009;
    let tmp2018 = tmp2016 + tmp2017;
    let tmp2019 = tmp1901 + tmp1903 + tmp2013;
    let tmp2020 = -1.0 * tmp1782 - 1.0 * tmp1783 - 1.0 * tmp1784;
    let tmp2021 = tmp1777 + tmp1778 + tmp1779;
    let tmp2022 = grad_phi_node_i[0] * tmp1692;
    let tmp2023 = grad_phi_node_i[0] * tmp1579;
    let tmp2024 = tmp1657 * tmp529
        + tmp1663 * tmp532
        + tmp1671 * tmp547
        + tmp1673 * tmp550
        + tmp1679 * tmp518
        + tmp1681 * tmp598
        + tmp1684 * tmp584
        + tmp1687 * tmp477
        - 2_f64 / 3.0 * tmp1872
        + tmp2020 * tmp777
        + tmp2021 * tmp788
        - 1.0 * tmp2022 * tmp527
        + tmp2023 * tmp494
        + (2_f64 / 3.0) * tmp554
        + tmp601 * (tmp1936 + tmp1967 + tmp1998 + tmp1999 + tmp2002)
        + tmp693 * (tmp1917 + tmp1930 + tmp2003 + tmp2004 + tmp2007)
        + tmp716 * (tmp1896 + tmp1960 + tmp2011 + tmp2015)
        + tmp736 * (tmp1895 + tmp1896 + tmp2018 + tmp2019);
    let tmp2025 = tmp1579 * tmp785 + tmp1797 * tmp368 + tmp1798 * tmp278 + tmp1806
        - 1_f64 / 3.0 * tmp1809
        - 1_f64 / 3.0 * tmp1810
        - 1_f64 / 3.0 * tmp1811
        + tmp1978;
    let tmp2026 = -1.0 * tmp1844;
    let tmp2027 = tmp1842 + tmp2026;
    let tmp2028 = tmp1839 + tmp1991;
    let tmp2029 = -1.0 * tmp1846 - 1.0 * tmp1847 + tmp1848;
    let tmp2030 = tmp1843 + tmp2026;
    let tmp2031 = tmp1838 + tmp1839;
    let tmp2032 = tmp2030 + tmp2031;
    let tmp2033 = tmp1854 + tmp1988;
    let tmp2034 = tmp1856 + tmp1857;
    let tmp2035 = tmp2033 + tmp2034;
    let tmp2036 = tmp1814 + tmp1981;
    let tmp2037 = tmp1817 + tmp1980;
    let tmp2038 = tmp2036 + tmp2037;
    let tmp2039 = -2_f64 / 3.0 * dRdx2_const[(0, 1)] - 2_f64 / 3.0 * tmp422 - 2_f64 / 3.0 * tmp423;
    let tmp2040 = (2_f64 / 3.0) * dRdx2_const[(1, 0)] + (2_f64 / 3.0) * tmp426 + (2_f64 / 3.0) * tmp427;
    let tmp2041 = 1.0 * a1;
    let tmp2042 = L_c.powi(2) * mu;
    let tmp2043 = -1.0 * tmp46 - 1.0 * tmp49 + tmp55 + tmp58;
    let tmp2044 = tmp2043 * tmp29;
    let tmp2045 = 2.0 * tmp2044;
    let tmp2046 = tmp2045 * tmp25;
    let tmp2047 = tmp2046 * tmp23;
    let tmp2048 = tmp2045 * tmp21;
    let tmp2049 = tmp133 + tmp135 + tmp137 + tmp144 + tmp2047 - 1.0 * tmp2048;
    let tmp2050 = grad_phi_node_i[0] * tmp2049;
    let tmp2051 = tmp18 * tmp2050;
    let tmp2052 = tmp2045 * tmp23;
    let tmp2053 = tmp19 * tmp2052;
    let tmp2054 = tmp20 * tmp2046;
    let tmp2055 = tmp114 + tmp2053 - 1.0 * tmp2054;
    let tmp2056 = grad_phi_node_i[1] * tmp2055;
    let tmp2057 = tmp2056 * tmp67;
    let tmp2058 = tmp20 * tmp2052;
    let tmp2059 = tmp19 * tmp2046 + tmp74;
    let tmp2060 = tmp126 - 1.0 * tmp2058 + tmp2059 + tmp71 + tmp81;
    let tmp2061 = grad_phi_node_i[2] * tmp2060;
    let tmp2062 = tmp2061 * tmp89;
    let tmp2063 = -1.0 * tmp71;
    let tmp2064 = tmp127 + tmp2058 + tmp2059 + tmp2063;
    let tmp2065 = grad_phi_node_i[0] * tmp2064;
    let tmp2066 = tmp110 * tmp2065;
    let tmp2067 = tmp159 + tmp2047 + tmp2048;
    let tmp2068 = grad_phi_node_i[1] * tmp2067;
    let tmp2069 = tmp119 * tmp2068;
    let tmp2070 = -1.0 * tmp104;
    let tmp2071 = tmp111 + tmp2053 + tmp2054 + tmp2070 + tmp98;
    let tmp2072 = grad_phi_node_i[2] * tmp2071;
    let tmp2073 = tmp125 * tmp2072;
    let tmp2074 = tmp2044 * tmp22;
    let tmp2075 = -1.0 * tmp2074;
    let tmp2076 = tmp2044 * tmp24;
    let tmp2077 = tmp52 + tmp61;
    let tmp2078 = tmp2044 * tmp27;
    let tmp2079 = tmp2044 * tmp26;
    let tmp2080 = tmp2078 - 1.0 * tmp2079;
    let tmp2081 = tmp2075 + tmp2076 + tmp2077 + tmp2080;
    let tmp2082 = grad_phi_node_i[0] * tmp2081;
    let tmp2083 = tmp132 * tmp2082;
    let tmp2084 = -1.0 * tmp2076 - 1.0 * tmp47;
    let tmp2085 = tmp120 + tmp2074 + tmp2080 + tmp2084 + tmp51 + tmp60;
    let tmp2086 = grad_phi_node_i[1] * tmp2085;
    let tmp2087 = tmp150 * tmp2086;
    let tmp2088 = tmp2075 + tmp2078 + tmp2079 + tmp2084 + tmp50 + tmp61;
    let tmp2089 = grad_phi_node_i[2] * tmp2088;
    let tmp2090 = tmp157 * tmp2089;
    let tmp2091 = grad_phi_node_i[0] * tmp2085;
    let tmp2092 = tmp163 * tmp2091;
    let tmp2093 = grad_phi_node_i[1] * tmp2088;
    let tmp2094 = tmp166 * tmp2093;
    let tmp2095 = grad_phi_node_i[2] * tmp2081;
    let tmp2096 = tmp169 * tmp2095;
    let tmp2097 = grad_phi_node_i[0] * tmp2055;
    let tmp2098 = grad_phi_node_i[0] * tmp2071;
    let tmp2099 = grad_phi_node_i[1] * tmp2064;
    let tmp2100 = grad_phi_node_i[1] * tmp2060;
    let tmp2101 = grad_phi_node_i[2] * tmp2049;
    let tmp2102 = grad_phi_node_i[2] * tmp2067;
    let tmp2103 = q_old[0] * tmp5;
    let tmp2104 = q_old[3] * tmp7;
    let tmp2105 = tmp2103 * tmp214 - 1.0 * tmp2104 * tmp214 + tmp651 - 1.0 * tmp666;
    let tmp2106 = tmp2105 * tmp213;
    let tmp2107 = tmp2106 * tmp5;
    let tmp2108 = tmp2107 * tmp4;
    let tmp2109 = tmp2106 * tmp7;
    let tmp2110 = tmp10 * tmp2109;
    let tmp2111 = tmp2108 + tmp2110 + tmp286 + tmp314;
    let tmp2112 = tmp2111 * tmp212;
    let tmp2113 = tmp10 * tmp2107;
    let tmp2114 = tmp2109 * tmp4;
    let tmp2115 = q_old[2] * tmp16;
    let tmp2116 = tmp2115 * tmp228;
    let tmp2117 = tmp217 * tmp229;
    let tmp2118 = tmp2114 + tmp2116 + tmp2117;
    let tmp2119 = tmp219 * tmp229;
    let tmp2120 = -1.0 * tmp2119;
    let tmp2121 = tmp228 * tmp399;
    let tmp2122 = -1.0 * tmp2121;
    let tmp2123 = tmp2120 + tmp2122;
    let tmp2124 = tmp2113 + tmp2118 + tmp2123;
    let tmp2125 = tmp2124 * tmp254;
    let tmp2126 = tmp2107 * tmp7;
    let tmp2127 = tmp2106 * tmp87;
    let tmp2128 = tmp2126 + tmp2127 + tmp232 + tmp240 + tmp323;
    let tmp2129 = tmp2128 * tmp278;
    let tmp2130 = -1.0 * tmp2114;
    let tmp2131 = -1.0 * tmp2116;
    let tmp2132 = tmp2122 + tmp2131;
    let tmp2133 = -1.0 * tmp2117;
    let tmp2134 = tmp2120 + tmp2133;
    let tmp2135 = tmp2113 + tmp2130 + tmp2132 + tmp2134;
    let tmp2136 = tmp2135 * tmp300;
    let tmp2137 = -1.0 * tmp2126;
    let tmp2138 = tmp2127 + tmp2137 + tmp325 + tmp360;
    let tmp2139 = tmp2138 * tmp312;
    let tmp2140 = -1.0 * tmp2108;
    let tmp2141 = tmp282 + tmp313 + tmp349;
    let tmp2142 = tmp2110 + tmp2140 + tmp2141;
    let tmp2143 = tmp2142 * tmp322;
    let tmp2144 = -1.0 * tmp2127;
    let tmp2145 = tmp2126 + tmp2144 + tmp241;
    let tmp2146 = tmp2145 * tmp345;
    let tmp2147 = -1.0 * tmp2110;
    let tmp2148 = tmp285 + tmp291 + tmp348;
    let tmp2149 = tmp2108 + tmp2147 + tmp2148;
    let tmp2150 = tmp2149 * tmp359;
    let tmp2151 = tmp2119 + tmp2121;
    let tmp2152 = -1.0 * tmp2113 + tmp2151;
    let tmp2153 = tmp2118 + tmp2152;
    let tmp2154 = tmp2153 * tmp368;
    let tmp2155 = tmp2137 + tmp2144 + tmp234 + tmp238 + tmp324;
    let tmp2156 = tmp2155 * tmp376;
    let tmp2157 = tmp2140 + tmp2147 + tmp350;
    let tmp2158 = tmp2157 * tmp382;
    let tmp2159 = tmp2131 + tmp2133;
    let tmp2160 = tmp2130 + tmp2152 + tmp2159;
    let tmp2161 = tmp2160 * tmp387;
    let tmp2162 = tmp11 * tmp2106;
    let tmp2163 = -1.0 * tmp13 * tmp651;
    let tmp2164 = tmp2106 * tmp8;
    let tmp2165 = q_old[3] * tmp394;
    let tmp2166 = tmp2106 * tmp6;
    let tmp2167 = tmp214 * tmp258;
    let tmp2168 = -1.0 * tmp2166 + tmp2167;
    let tmp2169 = tmp2106 * tmp9;
    let tmp2170 = tmp13 * tmp666;
    let tmp2171 = -1.0 * tmp2169 - 1.0 * tmp2170;
    let tmp2172 = tmp2162 + tmp2163 + tmp2164 + tmp2165 + tmp2168 + tmp2171;
    let tmp2173 = tmp2172 * tmp392;
    let tmp2174 = tmp2162 + tmp2163 - 1.0 * tmp2164 - 1.0 * tmp2165;
    let tmp2175 = tmp2166 - 1.0 * tmp2167 + tmp2171 + tmp2174;
    let tmp2176 = tmp2175 * tmp409;
    let tmp2177 = tmp2168 + tmp2169 + tmp2170 + tmp2174;
    let tmp2178 = tmp2177 * tmp415;
    let tmp2179 = tmp2175 * tmp420;
    let tmp2180 = tmp2177 * tmp424;
    let tmp2181 = tmp2172 * tmp428;
    let tmp2182 = -1.0 * tmp179 * tmp2097 + tmp179 * tmp2102 + tmp183 * tmp2098 - 1.0 * tmp183 * tmp2100
        + tmp186 * tmp2099
        - 1.0 * tmp186 * tmp2101
        + tmp2051
        + tmp2057
        + tmp2062
        + tmp2066
        + tmp2069
        + tmp2073
        + tmp2083
        + tmp2087
        + tmp2090
        + tmp2092
        + tmp2094
        + tmp2096
        + tmp2112
        + tmp2125
        + tmp2129
        + tmp2136
        + tmp2139
        + tmp2143
        + tmp2146
        + tmp2150
        + tmp2154
        + tmp2156
        + tmp2158
        + tmp2161
        + tmp2173
        + tmp2176
        + tmp2178
        + tmp2179
        + tmp2180
        + tmp2181;
    let tmp2183 = 2.0 * tmp2050;
    let tmp2184 = 2.0 * tmp63;
    let tmp2185 = 2.0 * tmp2056;
    let tmp2186 = 2.0 * tmp84;
    let tmp2187 = 2.0 * tmp107;
    let tmp2188 = 2.0 * tmp2061;
    let tmp2189 = 2.0 * tmp116;
    let tmp2190 = 2.0 * tmp2065;
    let tmp2191 = 2.0 * tmp2068;
    let tmp2192 = 2.0 * tmp122;
    let tmp2193 = 2.0 * tmp2072;
    let tmp2194 = 2.0 * tmp129;
    let tmp2195 = tmp100 * tmp2044;
    let tmp2196 = tmp30 * tmp40;
    let tmp2197 = tmp2196 * tmp70;
    let tmp2198 = tmp2044 * tmp93;
    let tmp2199 = tmp439 + tmp441 - 1.0 * tmp452 - 1.0 * tmp454;
    let tmp2200 = tmp2199 * tmp448;
    let tmp2201 = 2.0 * tmp2200;
    let tmp2202 = tmp2201 * tmp23;
    let tmp2203 = tmp20 * tmp2202;
    let tmp2204 = tmp2196 * tmp73;
    let tmp2205 = tmp2196 * tmp79;
    let tmp2206 = tmp2204 - 1.0 * tmp2205;
    let tmp2207 = tmp2201 * tmp25;
    let tmp2208 = tmp2196 * tmp77;
    let tmp2209 = tmp103 * tmp2044;
    let tmp2210 = tmp2044 * tmp96;
    let tmp2211 = tmp2209 + tmp2210;
    let tmp2212 = tmp19 * tmp2207 - 1.0 * tmp2208 + tmp2211;
    let tmp2213 = tmp2195 + tmp2197 + tmp2198 - 1.0 * tmp2203 + tmp2206 + tmp2212;
    let tmp2214 = tmp2213 * tmp513;
    let tmp2215 = tmp40 * tmp93;
    let tmp2216 = tmp20 * tmp2207;
    let tmp2217 = tmp2044 * tmp30;
    let tmp2218 = tmp2217 * tmp77;
    let tmp2219 = tmp103 * tmp40;
    let tmp2220 = -1.0 * tmp2219;
    let tmp2221 = tmp2217 * tmp79;
    let tmp2222 = tmp2217 * tmp73;
    let tmp2223 = -1.0 * tmp2222;
    let tmp2224 = tmp2221 + tmp2223;
    let tmp2225 = tmp2217 * tmp70;
    let tmp2226 = tmp40 * tmp96;
    let tmp2227 = tmp100 * tmp40;
    let tmp2228 = -1.0 * tmp2227;
    let tmp2229 = tmp2226 + tmp2228;
    let tmp2230 = tmp19 * tmp2202 + tmp2225 + tmp2229;
    let tmp2231 = tmp2215 + tmp2216 - 1.0 * tmp2218 + tmp2220 + tmp2224 + tmp2230;
    let tmp2232 = tmp2231 * tmp513;
    let tmp2233 = tmp2218 - 1.0 * tmp2221;
    let tmp2234 = -1.0 * tmp2215;
    let tmp2235 = tmp2219 + tmp2234;
    let tmp2236 = tmp2235 + tmp516;
    let tmp2237 = -1.0 * tmp2216 + tmp2223 + tmp2230 + tmp2233 + tmp2236;
    let tmp2238 = tmp2237 * tmp495;
    let tmp2239 = -1.0 * tmp2197;
    let tmp2240 = tmp2205 + tmp2239;
    let tmp2241 = -1.0 * tmp2198;
    let tmp2242 = -1.0 * tmp2195 + tmp2241;
    let tmp2243 = tmp2242 + tmp485 + tmp524;
    let tmp2244 = tmp2203 + tmp2204 + tmp2212 + tmp2240 + tmp2243;
    let tmp2245 = tmp2244 * tmp476;
    let tmp2246 = 2.0 * tmp147;
    let tmp2247 = 2.0 * tmp2082;
    let tmp2248 = 2.0 * tmp2086;
    let tmp2249 = 2.0 * tmp154;
    let tmp2250 = 2.0 * tmp161;
    let tmp2251 = 2.0 * tmp2089;
    let tmp2252 = 2.0 * tmp164;
    let tmp2253 = 2.0 * tmp2091;
    let tmp2254 = 2.0 * tmp2093;
    let tmp2255 = 2.0 * tmp167;
    let tmp2256 = 2.0 * tmp170;
    let tmp2257 = 2.0 * tmp2095;
    let tmp2258 = tmp2231 * tmp535;
    let tmp2259 = tmp2213 * tmp543;
    let tmp2260 = tmp36 * tmp40;
    let tmp2261 = tmp32 * tmp40;
    let tmp2262 = tmp21 * tmp2201;
    let tmp2263 = tmp34 * tmp40;
    let tmp2264 = tmp2044 * tmp46;
    let tmp2265 = tmp2044 * tmp49;
    let tmp2266 = tmp2264 - 1.0 * tmp2265;
    let tmp2267 = tmp2207 * tmp23 + tmp2263 + tmp2266;
    let tmp2268 = -1.0 * tmp577;
    let tmp2269 = tmp38 * tmp40;
    let tmp2270 = tmp2268 + tmp2269;
    let tmp2271 = tmp2044 * tmp55;
    let tmp2272 = tmp2044 * tmp58;
    let tmp2273 = tmp565 + tmp578 + tmp585;
    let tmp2274 = tmp2271 - 1.0 * tmp2272 + tmp2273;
    let tmp2275 = tmp2260 + tmp2261 - 1.0 * tmp2262 + tmp2267 + tmp2270 + tmp2274;
    let tmp2276 = tmp2275 * tmp476;
    let tmp2277 = tmp2260 - 1.0 * tmp2261 + tmp593;
    let tmp2278 = -1.0 * tmp2269 + tmp2277 + tmp586;
    let tmp2279 = tmp2262 + tmp2267 - 1.0 * tmp2271 + tmp2272 + tmp2278;
    let tmp2280 = tmp2279 * tmp495;
    let tmp2281 = tmp22 * tmp2200;
    let tmp2282 = tmp2044 * tmp32;
    let tmp2283 = tmp40 * tmp55;
    let tmp2284 = tmp2200 * tmp27;
    let tmp2285 = tmp2200 * tmp24;
    let tmp2286 = tmp40 * tmp46;
    let tmp2287 = tmp2044 * tmp34;
    let tmp2288 = tmp2044 * tmp38;
    let tmp2289 = -1.0 * tmp2288;
    let tmp2290 = tmp2287 + tmp2289;
    let tmp2291 = tmp2284 - 1.0 * tmp2285 - 1.0 * tmp2286 + tmp2290;
    let tmp2292 = tmp2200 * tmp26;
    let tmp2293 = tmp2044 * tmp36;
    let tmp2294 = tmp40 * tmp58;
    let tmp2295 = -1.0 * tmp2294;
    let tmp2296 = tmp40 * tmp49;
    let tmp2297 = tmp2295 - 1.0 * tmp2296;
    let tmp2298 = -1.0 * tmp2292 - 1.0 * tmp2293 + tmp2297;
    let tmp2299 = tmp2281 + tmp2282 - 1.0 * tmp2283 + tmp2291 + tmp2298;
    let tmp2300 = tmp2299 * tmp495;
    let tmp2301 = tmp2299 * tmp476;
    let tmp2302 = tmp2244 * tmp537;
    let tmp2303 = tmp2237 * tmp541;
    let tmp2304 = tmp2098 * tmp411;
    let tmp2305 = tmp184 * tmp2175;
    let tmp2306 = tmp187 * tmp2177;
    let tmp2307 = tmp2099 * tmp416;
    let tmp2308 = tmp2102 * tmp405;
    let tmp2309 = tmp190 * tmp2172;
    let tmp2310 = tmp2097 * tmp405;
    let tmp2311 = tmp180 * tmp2172;
    let tmp2312 = tmp188 * tmp2175;
    let tmp2313 = tmp2100 * tmp411;
    let tmp2314 = tmp2101 * tmp416;
    let tmp2315 = tmp189 * tmp2177;
    let tmp2316 = -1.0 * tmp2282;
    let tmp2317 = -1.0 * tmp2287;
    let tmp2318 = -1.0 * tmp2281;
    let tmp2319 = tmp2283 + tmp474;
    let tmp2320 = tmp2284 + tmp2285 + tmp2286 + tmp2289 + tmp2298 + tmp2316 + tmp2317 + tmp2318 + tmp2319;
    let tmp2321 = tmp2320 * tmp476;
    let tmp2322 = tmp2283 + tmp2296;
    let tmp2323 = tmp470 + tmp473;
    let tmp2324 = tmp2293 + tmp2316;
    let tmp2325 = tmp2291 + tmp2292 + tmp2295 + tmp2318 + tmp2322 + tmp2323 + tmp2324;
    let tmp2326 = tmp2325 * tmp513;
    let tmp2327 = tmp2325 * tmp495;
    let tmp2328 = tmp2320 * tmp513;
    let tmp2329 = tmp2279 * tmp539;
    let tmp2330 = tmp2275 * tmp545;
    let tmp2331 = tmp263 * tmp602;
    let tmp2332 = tmp256 * tmp602;
    let tmp2333 = tmp2103 * tmp602 - 1.0 * tmp2104 * tmp602 - 1.0 * tmp2331 + tmp2332;
    let tmp2334 = tmp2333 * tmp606;
    let tmp2335 = tmp2334 * tmp5;
    let tmp2336 = tmp2335 * tmp7;
    let tmp2337 = q_old[1] * tmp2109;
    let tmp2338 = tmp228 * tmp2337;
    let tmp2339 = q_old[2] * tmp2107;
    let tmp2340 = tmp228 * tmp2339;
    let tmp2341 = -1.0 * tmp2340;
    let tmp2342 = tmp2338 + tmp2341;
    let tmp2343 = tmp2336 + tmp2342;
    let tmp2344 = tmp2334 * tmp87;
    let tmp2345 = tmp2106 * tmp228;
    let tmp2346 = tmp2345 * tmp281;
    let tmp2347 = tmp2345 * tmp283;
    let tmp2348 = -1.0 * tmp2347;
    let tmp2349 = tmp2346 + tmp2348;
    let tmp2350 = tmp2344 + tmp2349;
    let tmp2351 = tmp1302 + tmp1347;
    let tmp2352 = tmp2334 * tmp7;
    let tmp2353 = tmp10 * tmp2352;
    let tmp2354 = q_old[0] * tmp228;
    let tmp2355 = tmp2109 * tmp2354;
    let tmp2356 = -1.0 * tmp2355;
    let tmp2357 = tmp231 * tmp2345;
    let tmp2358 = -1.0 * tmp2357;
    let tmp2359 = tmp2356 + tmp2358;
    let tmp2360 = tmp2353 + tmp2359;
    let tmp2361 = tmp2335 * tmp4;
    let tmp2362 = tmp227 * tmp2345;
    let tmp2363 = -1.0 * tmp2362;
    let tmp2364 = tmp2107 * tmp235;
    let tmp2365 = -1.0 * tmp2364;
    let tmp2366 = tmp2363 + tmp2365;
    let tmp2367 = -1.0 * tmp2361 + tmp2366;
    let tmp2368 = tmp1300 + tmp1346;
    let tmp2369 = tmp2362 + tmp2364;
    let tmp2370 = tmp2361 + tmp2369;
    let tmp2371 = tmp2355 + tmp2357;
    let tmp2372 = -1.0 * tmp2353 + tmp2371;
    let tmp2373 = -1.0 * tmp2346;
    let tmp2374 = tmp2347 + tmp2373;
    let tmp2375 = -1.0 * tmp2344 + tmp2374;
    let tmp2376 = -1.0 * tmp2338;
    let tmp2377 = tmp2340 + tmp2376;
    let tmp2378 = -1.0 * tmp2336 + tmp2377;
    let tmp2379 = tmp10 * tmp2335;
    let tmp2380 = tmp2345 * tmp256;
    let tmp2381 = tmp2107 * tmp2354;
    let tmp2382 = -1.0 * tmp2381;
    let tmp2383 = tmp2352 * tmp4;
    let tmp2384 = tmp1198 * tmp620;
    let tmp2385 = tmp1204 * tmp217;
    let tmp2386 = tmp2384 + tmp2385;
    let tmp2387 = tmp1198 * tmp618;
    let tmp2388 = -1.0 * tmp2387;
    let tmp2389 = tmp228 * tmp754;
    let tmp2390 = tmp2388 + tmp2389;
    let tmp2391 = tmp2109 * tmp235;
    let tmp2392 = tmp2345 * tmp263;
    let tmp2393 = -1.0 * tmp2392;
    let tmp2394 = tmp2391 + tmp2393;
    let tmp2395 = tmp2383 + tmp2386 + tmp2390 + tmp2394;
    let tmp2396 = tmp1204 * tmp219;
    let tmp2397 = -1.0 * tmp2396;
    let tmp2398 = tmp228 * tmp750;
    let tmp2399 = -1.0 * tmp2398;
    let tmp2400 = tmp1198 * tmp612;
    let tmp2401 = tmp1198 * tmp616;
    let tmp2402 = -1.0 * tmp2401;
    let tmp2403 = tmp2400 + tmp2402;
    let tmp2404 = tmp2397 + tmp2399 + tmp2403;
    let tmp2405 = -1.0 * tmp2383;
    let tmp2406 = -1.0 * tmp2389;
    let tmp2407 = tmp2387 + tmp2399 + tmp2402 + tmp2406;
    let tmp2408 = -1.0 * tmp2384;
    let tmp2409 = -1.0 * tmp2385;
    let tmp2410 = tmp2397 + tmp2400 + tmp2408 + tmp2409;
    let tmp2411 = tmp2380 + tmp2392;
    let tmp2412 = -1.0 * tmp2391;
    let tmp2413 = tmp2382 + tmp2412;
    let tmp2414 = tmp2411 + tmp2413;
    let tmp2415 = q_old[1] * tmp214;
    let tmp2416 = tmp2107 * tmp2415;
    let tmp2417 = tmp2334 * tmp6;
    let tmp2418 = q_old[2] * tmp214;
    let tmp2419 = tmp2109 * tmp2418;
    let tmp2420 = tmp11 * tmp2334;
    let tmp2421 = tmp2334 * tmp8;
    let tmp2422 = tmp2419 + tmp2420 - 1.0 * tmp2421;
    let tmp2423 = tmp2334 * tmp9;
    let tmp2424 = tmp2106 * tmp220;
    let tmp2425 = -1.0 * tmp2424;
    let tmp2426 = tmp2106 * tmp218;
    let tmp2427 = -1.0 * tmp2426;
    let tmp2428 = tmp2425 + tmp2427;
    let tmp2429 = -1.0 * tmp2423 + tmp2428;
    let tmp2430 = tmp2416 + tmp2417 + tmp2422 + tmp2429 + tmp668 + tmp698 + tmp733;
    let tmp2431 = -1.0 * tmp2400;
    let tmp2432 = tmp2396 + tmp2431;
    let tmp2433 = tmp2398 + tmp2401;
    let tmp2434 = tmp2432 + tmp2433;
    let tmp2435 = -1.0 * tmp2380;
    let tmp2436 = tmp2381 + tmp2435;
    let tmp2437 = -1.0 * tmp2379 + tmp2434 + tmp2436;
    let tmp2438 = tmp2387 + tmp2408;
    let tmp2439 = tmp2406 + tmp2409 + tmp2438;
    let tmp2440 = tmp1196 * tmp468;
    let tmp2441 = q_old[3] * tmp2440;
    let tmp2442 = q_old[2] * tmp2441;
    let tmp2443 = -1.0 * tmp2416;
    let tmp2444 = -1.0 * tmp2417 + tmp2440 * tmp657 + tmp2443 + tmp734;
    let tmp2445 = tmp2422 + tmp2423 + tmp2425 + tmp2426 + tmp2442 + tmp2444 + tmp667 + tmp697;
    let tmp2446 = -1.0 * tmp2419;
    let tmp2447 = tmp2420 + tmp2421 + tmp2429 - 1.0 * tmp2442 + tmp2444 + tmp2446 + tmp669;
    let tmp2448 = tmp2097 * tmp818;
    let tmp2449 = grad_phi_node_i[0] * tmp2067;
    let tmp2450 = tmp2449 * tmp823;
    let tmp2451 = tmp2099 * tmp826;
    let tmp2452 = grad_phi_node_i[1] * tmp2049;
    let tmp2453 = tmp2452 * tmp830;
    let tmp2454 = tmp2091 * tmp834;
    let tmp2455 = grad_phi_node_i[0] * tmp2088;
    let tmp2456 = tmp2455 * tmp830;
    let tmp2457 = tmp2086 * tmp839;
    let tmp2458 = grad_phi_node_i[1] * tmp2081;
    let tmp2459 = tmp2458 * tmp823;
    let tmp2460 = tmp2065 * tmp843;
    let tmp2461 = grad_phi_node_i[0] * tmp2060;
    let tmp2462 = tmp2461 * tmp845;
    let tmp2463 = tmp2056 * tmp848;
    let tmp2464 = grad_phi_node_i[1] * tmp2071;
    let tmp2465 = tmp2464 * tmp845;
    let tmp2466 = 0.5 * tmp2126;
    let tmp2467 = 0.5 * tmp2127;
    let tmp2468 = -1.0 * tmp2467;
    let tmp2469 = tmp2466 + tmp2468 + tmp894 + tmp899;
    let tmp2470 = tmp2469 * tmp854;
    let tmp2471 = 0.5 * tmp2114;
    let tmp2472 = q_old[2] * tmp865;
    let tmp2473 = -1.0 * tmp2472;
    let tmp2474 = tmp217 * tmp859;
    let tmp2475 = -1.0 * tmp2474;
    let tmp2476 = tmp2473 + tmp2475;
    let tmp2477 = -1.0 * tmp2471 + tmp2476;
    let tmp2478 = 0.5 * tmp2113;
    let tmp2479 = tmp219 * tmp859;
    let tmp2480 = q_old[1] * tmp867;
    let tmp2481 = tmp2479 + tmp2480;
    let tmp2482 = -1.0 * tmp2478 + tmp2481;
    let tmp2483 = tmp2477 + tmp2482;
    let tmp2484 = tmp2483 * tmp810;
    let tmp2485 = 0.5 * tmp2108;
    let tmp2486 = 0.5 * tmp2110;
    let tmp2487 = tmp1000 + tmp2485 + tmp2486 + tmp923;
    let tmp2488 = tmp2487 * tmp888;
    let tmp2489 = -1.0 * tmp2479;
    let tmp2490 = -1.0 * tmp2480;
    let tmp2491 = tmp2477 + tmp2478 + tmp2489 + tmp2490;
    let tmp2492 = tmp2491 * tmp778;
    let tmp2493 = tmp2487 * tmp910;
    let tmp2494 = -1.0 * tmp2485;
    let tmp2495 = tmp860 + tmp922 + tmp980;
    let tmp2496 = tmp2486 + tmp2494 + tmp2495;
    let tmp2497 = tmp2496 * tmp376;
    let tmp2498 = tmp2469 * tmp920;
    let tmp2499 = -1.0 * tmp2466;
    let tmp2500 = tmp2468 + tmp2499 + tmp890 + tmp893 + tmp897 + tmp913;
    let tmp2501 = tmp2500 * tmp359;
    let tmp2502 = tmp2177 * tmp930;
    let tmp2503 = tmp2175 * tmp937;
    let tmp2504 = tmp2172 * tmp942;
    let tmp2505 = tmp2175 * tmp949;
    let tmp2506 = tmp2448 + tmp2450 + tmp2451 + tmp2453 + tmp2454 + tmp2456 + tmp2457 + tmp2459 - 1.0 * tmp2460
        + tmp2462
        - 1.0 * tmp2463
        + tmp2465
        + tmp2470
        + tmp2484
        + tmp2488
        + tmp2492
        + tmp2493
        + tmp2497
        + tmp2498
        + tmp2501
        + tmp2502
        + tmp2503
        + tmp2504
        + tmp2505;
    let tmp2507 = tmp2098 * tmp834;
    let tmp2508 = tmp2461 * tmp958;
    let tmp2509 = grad_phi_node_i[2] * tmp2064;
    let tmp2510 = tmp2509 * tmp960;
    let tmp2511 = tmp2101 * tmp964;
    let tmp2512 = tmp2082 * tmp818;
    let tmp2513 = tmp2455 * tmp960;
    let tmp2514 = tmp2095 * tmp968;
    let tmp2515 = grad_phi_node_i[2] * tmp2085;
    let tmp2516 = tmp2515 * tmp958;
    let tmp2517 = tmp2050 * tmp843;
    let tmp2518 = tmp2449 * tmp848;
    let tmp2519 = grad_phi_node_i[2] * tmp2055;
    let tmp2520 = tmp2519 * tmp848;
    let tmp2521 = tmp2072 * tmp845;
    let tmp2522 = tmp2466 + tmp2467 + tmp891 + tmp892 + tmp898 + tmp912;
    let tmp2523 = tmp2522 * tmp977;
    let tmp2524 = tmp2496 * tmp748;
    let tmp2525 = -1.0 * tmp2486;
    let tmp2526 = tmp2494 + tmp2525 + tmp863 + tmp870;
    let tmp2527 = tmp2526 * tmp812;
    let tmp2528 = tmp2471 + tmp2472 + tmp2474 + tmp2482;
    let tmp2529 = tmp2528 * tmp991;
    let tmp2530 = tmp2528 * tmp910;
    let tmp2531 = tmp2483 * tmp345;
    let tmp2532 = tmp2467 + tmp2499 + tmp914 + tmp985;
    let tmp2533 = tmp2532 * tmp387;
    let tmp2534 = tmp1005 * tmp2522;
    let tmp2535 = tmp1010 * tmp2177;
    let tmp2536 = tmp1015 * tmp2172;
    let tmp2537 = tmp1022 * tmp2172;
    let tmp2538 = tmp1027 * tmp2175;
    let tmp2539 = tmp2507 + tmp2508 + tmp2510 + tmp2511 + tmp2512 + tmp2513 + tmp2514 + tmp2516 - 1.0 * tmp2517
        + tmp2518
        + tmp2520
        - 1.0 * tmp2521
        + tmp2523
        + tmp2524
        + tmp2527
        + tmp2529
        + tmp2530
        + tmp2531
        + tmp2533
        + tmp2534
        + tmp2535
        + tmp2536
        + tmp2537
        + tmp2538;
    let tmp2540 = tmp1035 * tmp2098;
    let tmp2541 = tmp2461 * tmp839;
    let tmp2542 = tmp2509 * tmp826;
    let tmp2543 = tmp2101 * tmp830;
    let tmp2544 = tmp2455 * tmp826;
    let tmp2545 = tmp1041 * tmp2082;
    let tmp2546 = tmp2515 * tmp839;
    let tmp2547 = tmp2095 * tmp823;
    let tmp2548 = tmp2500 * tmp977;
    let tmp2549 = tmp862 + tmp868 + tmp979;
    let tmp2550 = tmp2485 + tmp2525 + tmp2549;
    let tmp2551 = tmp2550 * tmp748;
    let tmp2552 = tmp2487 * tmp812;
    let tmp2553 = tmp2491 * tmp991;
    let tmp2554 = tmp2491 * tmp910;
    let tmp2555 = tmp2472 + tmp2489;
    let tmp2556 = tmp2474 + tmp2490;
    let tmp2557 = tmp2555 + tmp2556;
    let tmp2558 = tmp2471 + tmp2478 + tmp2557;
    let tmp2559 = tmp2558 * tmp345;
    let tmp2560 = tmp2469 * tmp387;
    let tmp2561 = tmp1005 * tmp2500;
    let tmp2562 = tmp1055 * tmp2177;
    let tmp2563 = tmp1057 * tmp2172;
    let tmp2564 = tmp1059 * tmp2172;
    let tmp2565 = tmp1061 * tmp2175;
    let tmp2566 = tmp2517 - 1.0 * tmp2518 - 1.0 * tmp2520
        + tmp2521
        + tmp2540
        + tmp2541
        + tmp2542
        + tmp2543
        + tmp2544
        + tmp2545
        + tmp2546
        + tmp2547
        + tmp2548
        + tmp2551
        + tmp2552
        + tmp2553
        + tmp2554
        + tmp2559
        + tmp2560
        + tmp2561
        + tmp2562
        + tmp2563
        + tmp2564
        + tmp2565;
    let tmp2567 = tmp1041 * tmp2097;
    let tmp2568 = tmp2449 * tmp968;
    let tmp2569 = tmp2099 * tmp960;
    let tmp2570 = tmp2452 * tmp964;
    let tmp2571 = tmp2455 * tmp964;
    let tmp2572 = tmp1035 * tmp2091;
    let tmp2573 = tmp2458 * tmp968;
    let tmp2574 = tmp2086 * tmp958;
    let tmp2575 = tmp2532 * tmp854;
    let tmp2576 = tmp2558 * tmp810;
    let tmp2577 = tmp2526 * tmp888;
    let tmp2578 = tmp2528 * tmp778;
    let tmp2579 = tmp2526 * tmp910;
    let tmp2580 = tmp2550 * tmp376;
    let tmp2581 = tmp2532 * tmp920;
    let tmp2582 = tmp2522 * tmp359;
    let tmp2583 = tmp1081 * tmp2177;
    let tmp2584 = tmp1083 * tmp2175;
    let tmp2585 = tmp1085 * tmp2172;
    let tmp2586 = tmp1087 * tmp2175;
    let tmp2587 = tmp2460 - 1.0 * tmp2462 + tmp2463 - 1.0 * tmp2465
        + tmp2567
        + tmp2568
        + tmp2569
        + tmp2570
        + tmp2571
        + tmp2572
        + tmp2573
        + tmp2574
        + tmp2575
        + tmp2576
        + tmp2577
        + tmp2578
        + tmp2579
        + tmp2580
        + tmp2581
        + tmp2582
        + tmp2583
        + tmp2584
        + tmp2585
        + tmp2586;
    let tmp2588 = tmp2464 * tmp834;
    let tmp2589 = tmp2100 * tmp958;
    let tmp2590 = tmp2519 * tmp818;
    let tmp2591 = tmp2102 * tmp823;
    let tmp2592 = tmp2458 * tmp818;
    let tmp2593 = tmp2093 * tmp960;
    let tmp2594 = tmp2515 * tmp834;
    let tmp2595 = tmp2089 * tmp830;
    let tmp2596 = tmp2452 * tmp843;
    let tmp2597 = tmp2068 * tmp848;
    let tmp2598 = tmp2509 * tmp843;
    let tmp2599 = tmp2061 * tmp845;
    let tmp2600 = tmp2522 * tmp811;
    let tmp2601 = tmp1104 * tmp2496;
    let tmp2602 = tmp2469 * tmp785;
    let tmp2603 = tmp1107 * tmp2483;
    let tmp2604 = tmp2528 * tmp382;
    let tmp2605 = tmp2483 * tmp920;
    let tmp2606 = tmp2487 * tmp368;
    let tmp2607 = tmp1005 * tmp2496;
    let tmp2608 = tmp1116 * tmp2177;
    let tmp2609 = tmp1121 * tmp2172;
    let tmp2610 = tmp1126 * tmp2177;
    let tmp2611 = tmp1131 * tmp2175;
    let tmp2612 = tmp2588 + tmp2589 + tmp2590 + tmp2591 + tmp2592 + tmp2593 + tmp2594 + tmp2595 - 1.0 * tmp2596
        + tmp2597
        - 1.0 * tmp2598
        + tmp2599
        + tmp2600
        + tmp2601
        + tmp2602
        + tmp2603
        + tmp2604
        + tmp2605
        + tmp2606
        + tmp2607
        + tmp2608
        + tmp2609
        + tmp2610
        + tmp2611;
    let tmp2613 = tmp1035 * tmp2464;
    let tmp2614 = tmp2100 * tmp839;
    let tmp2615 = tmp1041 * tmp2519;
    let tmp2616 = tmp2102 * tmp968;
    let tmp2617 = tmp2093 * tmp826;
    let tmp2618 = tmp1041 * tmp2458;
    let tmp2619 = tmp2089 * tmp964;
    let tmp2620 = tmp1035 * tmp2515;
    let tmp2621 = tmp2500 * tmp811;
    let tmp2622 = tmp1104 * tmp2550;
    let tmp2623 = tmp2532 * tmp785;
    let tmp2624 = tmp1107 * tmp2558;
    let tmp2625 = tmp2491 * tmp382;
    let tmp2626 = tmp2558 * tmp920;
    let tmp2627 = tmp2526 * tmp368;
    let tmp2628 = tmp1005 * tmp2550;
    let tmp2629 = tmp1155 * tmp2177;
    let tmp2630 = tmp1157 * tmp2172;
    let tmp2631 = tmp1159 * tmp2177;
    let tmp2632 = tmp1161 * tmp2175;
    let tmp2633 = tmp2596 - 1.0 * tmp2597 + tmp2598 - 1.0 * tmp2599
        + tmp2613
        + tmp2614
        + tmp2615
        + tmp2616
        + tmp2617
        + tmp2618
        + tmp2619
        + tmp2620
        + tmp2621
        + tmp2622
        + tmp2623
        + tmp2624
        + tmp2625
        + tmp2626
        + tmp2627
        + tmp2628
        + tmp2629
        + tmp2630
        + tmp2631
        + tmp2632;
    let tmp2634 = tmp675 * tmp858;
    let tmp2635 = tmp222 * tmp858;
    let tmp2636 = tmp2635 * tmp283;
    let tmp2637 = -1.0 * tmp2636;
    let tmp2638 = tmp2635 * tmp281;
    let tmp2639 = tmp673 * tmp858;
    let tmp2640 = -1.0 * tmp2639;
    let tmp2641 = 0.5 * tmp2361;
    let tmp2642 = tmp2106 * tmp858;
    let tmp2643 = tmp227 * tmp2642;
    let tmp2644 = tmp2107 * tmp858;
    let tmp2645 = q_old[3] * tmp2644;
    let tmp2646 = tmp2641 + tmp2643 + tmp2645;
    let tmp2647 = 0.5 * tmp2353;
    let tmp2648 = tmp2109 * tmp858;
    let tmp2649 = q_old[0] * tmp2648;
    let tmp2650 = tmp231 * tmp2642;
    let tmp2651 = tmp2647 - 1.0 * tmp2649 - 1.0 * tmp2650;
    let tmp2652 = tmp1343 + tmp2634 + tmp2637 + tmp2638 + tmp2640 + tmp2646 + tmp2651;
    let tmp2653 = tmp2642 * tmp283;
    let tmp2654 = tmp2642 * tmp281;
    let tmp2655 = -1.0 * tmp2654;
    let tmp2656 = q_old[1] * tmp2648;
    let tmp2657 = q_old[2] * tmp2644;
    let tmp2658 = -1.0 * tmp2657;
    let tmp2659 = 0.5 * tmp2344;
    let tmp2660 = tmp227 * tmp2635;
    let tmp2661 = tmp231 * tmp2635;
    let tmp2662 = -1.0 * tmp2659 + tmp2660 - 1.0 * tmp2661;
    let tmp2663 = 0.5 * tmp2336;
    let tmp2664 = tmp1827 * tmp858;
    let tmp2665 = tmp636 * tmp858;
    let tmp2666 = tmp2663 + tmp2664 - 1.0 * tmp2665;
    let tmp2667 = tmp1288 + tmp2653 + tmp2655 + tmp2656 + tmp2658 + tmp2662 + tmp2666;
    let tmp2668 = tmp263 * tmp2642;
    let tmp2669 = q_old[3] * tmp2648;
    let tmp2670 = -1.0 * tmp2669;
    let tmp2671 = q_old[0] * tmp2644;
    let tmp2672 = tmp256 * tmp2642;
    let tmp2673 = -1.0 * tmp2672;
    let tmp2674 = 0.5 * tmp2383;
    let tmp2675 = 0.125 * tmp1197;
    let tmp2676 = tmp2675 * tmp620;
    let tmp2677 = -1.0 * tmp2676;
    let tmp2678 = tmp217 * tmp2635;
    let tmp2679 = tmp2677 - 1.0 * tmp2678;
    let tmp2680 = tmp2675 * tmp618;
    let tmp2681 = tmp754 * tmp858;
    let tmp2682 = tmp2680 - 1.0 * tmp2681;
    let tmp2683 = tmp2679 + tmp2682;
    let tmp2684 = -1.0 * tmp2674 + tmp2683;
    let tmp2685 = 0.5 * tmp2379;
    let tmp2686 = tmp2675 * tmp616;
    let tmp2687 = tmp750 * tmp858;
    let tmp2688 = tmp2686 + tmp2687;
    let tmp2689 = tmp2675 * tmp612;
    let tmp2690 = -1.0 * tmp2689;
    let tmp2691 = tmp219 * tmp2635;
    let tmp2692 = tmp2690 + tmp2691;
    let tmp2693 = tmp2688 + tmp2692;
    let tmp2694 = -1.0 * tmp2685 + tmp2693;
    let tmp2695 = tmp2668 + tmp2670 + tmp2671 + tmp2673 + tmp2684 + tmp2694;
    let tmp2696 = -1.0 * tmp2634;
    let tmp2697 = tmp2638 + tmp2696;
    let tmp2698 = tmp2636 + tmp2640;
    let tmp2699 = tmp2697 + tmp2698;
    let tmp2700 = -1.0 * tmp2641 - 1.0 * tmp2643 - 1.0 * tmp2645;
    let tmp2701 = tmp2651 + tmp2699 + tmp2700;
    let tmp2702 = grad_phi_node_i[0] * tmp183;
    let tmp2703 = tmp2213 * tmp2702;
    let tmp2704 = (1_f64 / 2.0) * tmp2175;
    let tmp2705 = tmp2704 * tmp846;
    let tmp2706 = (1_f64 / 2.0) * tmp411;
    let tmp2707 = tmp2461 * tmp2706;
    let tmp2708 = tmp2325 * tmp476;
    let tmp2709 = tmp2279 * tmp823;
    let tmp2710 = tmp2237 * tmp818;
    let tmp2711 = 2.0 * tmp2487;
    let tmp2712 = 2.0 * tmp901;
    let tmp2713 = 2.0 * tmp2469;
    let tmp2714 = 2.0 * tmp871;
    let tmp2715 = 2.0 * tmp2496;
    let tmp2716 = 2.0 * tmp916;
    let tmp2717 = 2.0 * tmp2483;
    let tmp2718 = 2.0 * tmp884;
    let tmp2719 = tmp186 * tmp2244;
    let tmp2720 = grad_phi_node_i[0] * tmp2719;
    let tmp2721 = (1_f64 / 2.0) * tmp2177;
    let tmp2722 = tmp116 * tmp2721;
    let tmp2723 = (1_f64 / 2.0) * tmp416;
    let tmp2724 = tmp2065 * tmp2723;
    let tmp2725 = -1.0 * tmp2691;
    let tmp2726 = -1.0 * tmp2687;
    let tmp2727 = -1.0 * tmp2671;
    let tmp2728 = tmp2670 + tmp2689;
    let tmp2729 = -1.0 * tmp2686;
    let tmp2730 = tmp2668 + tmp2729;
    let tmp2731 = tmp2672 + tmp2727 + tmp2728 + tmp2730;
    let tmp2732 = tmp2684 + tmp2685 + tmp2725 + tmp2726 + tmp2731;
    let tmp2733 = -1.0 * tmp2656;
    let tmp2734 = tmp2653 + tmp2733;
    let tmp2735 = tmp2655 + tmp2657;
    let tmp2736 = tmp2734 + tmp2735;
    let tmp2737 = -1.0 * tmp2663 - 1.0 * tmp2664 + tmp2665;
    let tmp2738 = tmp2662 + tmp2736 + tmp2737;
    let tmp2739 = grad_phi_node_i[1] * tmp183;
    let tmp2740 = tmp2231 * tmp2739;
    let tmp2741 = tmp2704 * tmp850;
    let tmp2742 = tmp2464 * tmp2706;
    let tmp2743 = tmp2275 * tmp830;
    let tmp2744 = tmp2320 * tmp495;
    let tmp2745 = tmp2244 * tmp826;
    let tmp2746 = tmp2299 * tmp839;
    let tmp2747 = 2.0 * tmp841;
    let tmp2748 = 2.0 * tmp2458;
    let tmp2749 = 2.0 * tmp2491;
    let tmp2750 = 2.0 * tmp906;
    let tmp2751 = 2.0 * tmp2500;
    let tmp2752 = 2.0 * tmp925;
    let tmp2753 = tmp1410 * tmp2237;
    let tmp2754 = (1_f64 / 2.0) * tmp2172;
    let tmp2755 = tmp2754 * tmp84;
    let tmp2756 = (1_f64 / 2.0) * tmp1406;
    let tmp2757 = tmp2055 * tmp2756;
    let tmp2758 = tmp1397 * tmp2667
        + tmp1476 * tmp2430
        + tmp1477 * tmp2652
        + tmp1478 * tmp2732
        + tmp1482 * tmp2447
        + tmp154 * tmp2751
        + tmp187 * tmp2749
        + tmp2086 * tmp2752
        + tmp2099 * tmp2750
        + tmp2452 * tmp2712
        + tmp2469 * tmp2747
        + tmp2711 * tmp831
        + tmp2738 * tmp722
        + tmp2740
        + tmp2741
        + tmp2742
        + tmp2743 * tmp495
        + tmp2744 * tmp823
        + tmp2745 * tmp495
        + tmp2746 * tmp495
        + tmp2748 * tmp871
        - 1.0 * tmp2753
        - 1.0 * tmp2755
        - 1.0 * tmp2757;
    let tmp2759 = -1.0 * tmp2638;
    let tmp2760 = -1.0 * tmp2647 + tmp2649 + tmp2650;
    let tmp2761 = tmp1297 + tmp2636 + tmp2639 + tmp2696 + tmp2700 + tmp2759 + tmp2760;
    let tmp2762 = -1.0 * tmp2668;
    let tmp2763 = -1.0 * tmp2680;
    let tmp2764 = tmp2673 + tmp2763;
    let tmp2765 = tmp2671 + tmp2676;
    let tmp2766 = tmp2669 + tmp2762 + tmp2764 + tmp2765;
    let tmp2767 = tmp2674 + tmp2678 + tmp2681 + tmp2694 + tmp2766;
    let tmp2768 = -1.0 * tmp2653;
    let tmp2769 = tmp2659 - 1.0 * tmp2660 + tmp2661;
    let tmp2770 = tmp1327 + tmp2654 + tmp2657 + tmp2733 + tmp2737 + tmp2768 + tmp2769;
    let tmp2771 = tmp2654 + tmp2658;
    let tmp2772 = tmp2656 + tmp2768;
    let tmp2773 = tmp2771 + tmp2772;
    let tmp2774 = tmp2666 + tmp2769 + tmp2773;
    let tmp2775 = tmp2244 * tmp960;
    let tmp2776 = tmp2275 * tmp964;
    let tmp2777 = 2.0 * tmp2528;
    let tmp2778 = 2.0 * tmp996;
    let tmp2779 = 2.0 * tmp2522;
    let tmp2780 = 2.0 * tmp2532;
    let tmp2781 = 2.0 * tmp982;
    let tmp2782 = 2.0 * tmp1001;
    let tmp2783 = 2.0 * tmp2526;
    let tmp2784 = 2.0 * tmp987;
    let tmp2785 = tmp2727 + tmp2762;
    let tmp2786 = tmp2669 + tmp2672;
    let tmp2787 = tmp2681 + tmp2689 + tmp2725 + tmp2763;
    let tmp2788 = tmp2676 + tmp2678 + tmp2726 + tmp2729;
    let tmp2789 = tmp2787 + tmp2788;
    let tmp2790 = tmp2674 + tmp2685 + tmp2785 + tmp2786 + tmp2789;
    let tmp2791 = tmp2637 + tmp2639;
    let tmp2792 = tmp2634 + tmp2759;
    let tmp2793 = tmp2791 + tmp2792;
    let tmp2794 = tmp2646 + tmp2760 + tmp2793;
    let tmp2795 = tmp1041 * tmp2237;
    let tmp2796 = tmp2279 * tmp968;
    let tmp2797 = tmp2325 * tmp964;
    let tmp2798 = 2.0 * tmp2558;
    let tmp2799 = 2.0 * tmp1051;
    let tmp2800 = 2.0 * tmp2550;
    let tmp2801 = 2.0 * tmp1046;
    let tmp2802 = tmp1188 * tmp2761
        + tmp1274 * tmp2299
        + tmp1463 * tmp2770
        + tmp1464 * tmp2790
        + tmp1491 * tmp2445
        + tmp1492 * tmp2430
        + tmp164 * tmp2800
        + tmp180 * tmp2798
        + tmp2091 * tmp2801
        + tmp2097 * tmp2799
        + tmp2449 * tmp2782
        + tmp2455 * tmp2784
        - 1.0 * tmp2703
        - 1.0 * tmp2705
        - 1.0 * tmp2707
        + tmp2720
        + tmp2722
        + tmp2724
        + tmp2780 * tmp824
        + tmp2783 * tmp836
        + tmp2794 * tmp736
        + tmp2795 * tmp476
        + tmp2796 * tmp476
        + tmp2797 * tmp476;
    let tmp2803 = tmp186 * tmp2275;
    let tmp2804 = grad_phi_node_i[0] * tmp2803;
    let tmp2805 = tmp2721 * tmp63;
    let tmp2806 = tmp2050 * tmp2723;
    let tmp2807 = tmp2213 * tmp839;
    let tmp2808 = tmp1281 * tmp2279;
    let tmp2809 = tmp2754 * tmp824;
    let tmp2810 = (1_f64 / 2.0) * tmp405;
    let tmp2811 = tmp2449 * tmp2810;
    let tmp2812 = tmp1311 * tmp2231;
    let tmp2813 = tmp129 * tmp2704;
    let tmp2814 = tmp2072 * tmp2706;
    let tmp2815 = tmp1318 * tmp2237;
    let tmp2816 = tmp2754 * tmp974;
    let tmp2817 = tmp2519 * tmp2810;
    let tmp2818 = tmp1286 * tmp2430
        + tmp1287 * tmp2652
        + tmp1295 * tmp2732
        + tmp1305 * tmp2738
        + tmp1309 * tmp2447
        + tmp189 * tmp2711
        + tmp2101 * tmp2712
        + tmp2256 * tmp2469
        + tmp2257 * tmp871
        + tmp2328 * tmp823
        + tmp2509 * tmp2750
        + tmp2515 * tmp2752
        + tmp2667 * tmp744
        + tmp2743 * tmp513
        + tmp2745 * tmp513
        + tmp2746 * tmp513
        + tmp2749 * tmp961
        + tmp2751 * tmp970
        + tmp2812
        + tmp2813
        + tmp2814
        - 1.0 * tmp2815
        - 1.0 * tmp2816
        - 1.0 * tmp2817;
    let tmp2819 = tmp2325 * tmp960;
    let tmp2820 = tmp2213 * tmp958;
    let tmp2821 = tmp1188 * tmp2767
        + tmp1222 * tmp2774
        + tmp1246 * tmp2701
        + tmp1363 * tmp2447
        + tmp1374 * tmp2445
        + tmp1377 * tmp2231
        + tmp184 * tmp2715
        + tmp2098 * tmp2716
        + tmp2246 * tmp2483
        + tmp2247 * tmp884
        + tmp2321 * tmp818
        + tmp2455 * tmp2778
        + tmp2461 * tmp2781
        + tmp2695 * tmp716
        + tmp2777 * tmp836
        + tmp2779 * tmp846
        - 1.0 * tmp2804
        - 1.0 * tmp2805
        - 1.0 * tmp2806
        + tmp2808
        + tmp2809
        + tmp2811
        + tmp2819 * tmp476
        + tmp2820 * tmp476;
    let tmp2822 = grad_phi_node_i[1] * tmp2803;
    let tmp2823 = tmp2721 * tmp831;
    let tmp2824 = tmp2452 * tmp2723;
    let tmp2825 = tmp2231 * tmp495;
    let tmp2826 = tmp1410 * tmp2279;
    let tmp2827 = tmp122 * tmp2754;
    let tmp2828 = tmp2067 * tmp2756;
    let tmp2829 = grad_phi_node_i[2] * tmp2719;
    let tmp2830 = tmp2721 * tmp961;
    let tmp2831 = tmp2509 * tmp2723;
    let tmp2832 = tmp1035 * tmp513;
    let tmp2833 = tmp1311 * tmp2213;
    let tmp2834 = tmp107 * tmp2704;
    let tmp2835 = tmp2061 * tmp2706;
    let tmp2836 = tmp1305 * tmp2794
        + tmp1415 * tmp2445
        + tmp1416 * tmp2770
        + tmp1417 * tmp2790
        + tmp1421 * tmp2430
        + tmp161 * tmp2783
        + tmp190 * tmp2780
        + tmp2089 * tmp2784
        + tmp2102 * tmp2782
        + tmp2299 * tmp2832
        + tmp2515 * tmp2801
        + tmp2519 * tmp2799
        + tmp2761 * tmp729
        + tmp2795 * tmp513
        + tmp2796 * tmp513
        + tmp2797 * tmp513
        + tmp2798 * tmp974
        + tmp2800 * tmp970
        + tmp2829
        + tmp2830
        + tmp2831
        - 1.0 * tmp2833
        - 1.0 * tmp2834
        - 1.0 * tmp2835;
    let tmp2837 = tmp1397 * tmp2695
        + tmp1398 * tmp2774
        + tmp1399 * tmp2701
        + tmp1442 * tmp2447
        + tmp1443 * tmp2445
        + tmp167 * tmp2777
        + tmp188 * tmp2779
        + tmp2093 * tmp2778
        + tmp2100 * tmp2781
        + tmp2464 * tmp2716
        + tmp2483 * tmp2747
        + tmp2715 * tmp850
        + tmp2744 * tmp818
        + tmp2748 * tmp884
        + tmp2767 * tmp740
        + tmp2819 * tmp495
        + tmp2820 * tmp495
        - 1.0 * tmp2822
        - 1.0 * tmp2823
        - 1.0 * tmp2824
        + tmp2825 * tmp834
        + tmp2826
        + tmp2827
        + tmp2828;
    let tmp2838 = tmp179 * tmp2449;
    let tmp2839 = tmp186 * tmp2050;
    let tmp2840 = 2.0 * tmp2507
        + 2.0 * tmp2508
        + 2.0 * tmp2512
        + 2.0 * tmp2513
        + 2.0 * tmp2523
        + 2.0 * tmp2524
        + 2.0 * tmp2530
        + 2.0 * tmp2531
        + 2.0 * tmp2535
        + 2.0 * tmp2536
        + tmp2838
        - 1.0 * tmp2839;
    let tmp2841 = tmp183 * tmp2072;
    let tmp2842 = tmp179 * tmp2519;
    let tmp2843 = 2.0 * tmp2542
        + 2.0 * tmp2543
        + 2.0 * tmp2546
        + 2.0 * tmp2547
        + 2.0 * tmp2552
        + 2.0 * tmp2553
        + 2.0 * tmp2560
        + 2.0 * tmp2561
        + 2.0 * tmp2564
        + 2.0 * tmp2565
        + tmp2841
        - 1.0 * tmp2842;
    let tmp2844 = tmp2840 + tmp2843;
    let tmp2845 = tmp183 * tmp2464;
    let tmp2846 = tmp179 * tmp2056;
    let tmp2847 = 2.0 * tmp2451
        + 2.0 * tmp2453
        + 2.0 * tmp2457
        + 2.0 * tmp2459
        + 2.0 * tmp2488
        + 2.0 * tmp2492
        + 2.0 * tmp2498
        + 2.0 * tmp2501
        + 2.0 * tmp2504
        + 2.0 * tmp2505
        + tmp2845
        - 1.0 * tmp2846;
    let tmp2848 = tmp186 * tmp2065;
    let tmp2849 = tmp183 * tmp2461;
    let tmp2850 = 2.0 * tmp2567
        + 2.0 * tmp2568
        + 2.0 * tmp2571
        + 2.0 * tmp2572
        + 2.0 * tmp2575
        + 2.0 * tmp2576
        + 2.0 * tmp2579
        + 2.0 * tmp2580
        + 2.0 * tmp2583
        + 2.0 * tmp2584
        + tmp2848
        - 1.0 * tmp2849;
    let tmp2851 = tmp2847 + tmp2850;
    let tmp2852 = tmp179 * tmp2068;
    let tmp2853 = tmp186 * tmp2452;
    let tmp2854 = 2.0 * tmp2588
        + 2.0 * tmp2589
        + 2.0 * tmp2592
        + 2.0 * tmp2593
        + 2.0 * tmp2600
        + 2.0 * tmp2601
        + 2.0 * tmp2604
        + 2.0 * tmp2605
        + 2.0 * tmp2608
        + 2.0 * tmp2609
        + tmp2852
        - 1.0 * tmp2853;
    let tmp2855 = tmp186 * tmp2509;
    let tmp2856 = tmp183 * tmp2061;
    let tmp2857 = 2.0 * tmp2615
        + 2.0 * tmp2616
        + 2.0 * tmp2619
        + 2.0 * tmp2620
        + 2.0 * tmp2623
        + 2.0 * tmp2624
        + 2.0 * tmp2627
        + 2.0 * tmp2628
        + 2.0 * tmp2631
        + 2.0 * tmp2632
        + tmp2855
        - 1.0 * tmp2856;
    let tmp2858 = tmp2854 + tmp2857;
    let tmp2859 = 0.666666666666667 * tmp2113;
    let tmp2860 = tmp1509 * tmp219;
    let tmp2861 = tmp1504 * tmp399;
    let tmp2862 = tmp2860 + tmp2861;
    let tmp2863 = -1.0 * tmp2859 + tmp2862;
    let tmp2864 = 0.666666666666667 * tmp2114;
    let tmp2865 = tmp1504 * tmp2115;
    let tmp2866 = tmp1509 * tmp217;
    let tmp2867 = tmp2865 + tmp2866;
    let tmp2868 = tmp2864 + tmp2867;
    let tmp2869 = tmp2863 + tmp2868;
    let tmp2870 = tmp2869 * tmp368;
    let tmp2871 = -1.0 * tmp2865;
    let tmp2872 = -1.0 * tmp2866;
    let tmp2873 = tmp2871 + tmp2872;
    let tmp2874 = -1.0 * tmp2864 + tmp2873;
    let tmp2875 = tmp2863 + tmp2874;
    let tmp2876 = tmp2875 * tmp387;
    let tmp2877 = 0.666666666666667 * tmp2127;
    let tmp2878 = 0.666666666666667 * tmp2126;
    let tmp2879 = tmp1546 + tmp2878;
    let tmp2880 = tmp1543 + tmp1544 + tmp2877 + tmp2879;
    let tmp2881 = tmp278 * tmp2880;
    let tmp2882 = 0.666666666666667 * tmp2110;
    let tmp2883 = 0.666666666666667 * tmp2108;
    let tmp2884 = -1.0 * tmp2883;
    let tmp2885 = tmp1528 + tmp1530 + tmp1704;
    let tmp2886 = tmp2882 + tmp2884 + tmp2885;
    let tmp2887 = tmp2886 * tmp322;
    let tmp2888 = tmp1554 * tmp2177;
    let tmp2889 = tmp1559 * tmp2172;
    let tmp2890 = tmp1564 * tmp2089;
    let tmp2891 = tmp1568 * tmp2061;
    let tmp2892 = tmp1573 * tmp2072;
    let tmp2893 = tmp1576 * tmp2095;
    let tmp2894 = tmp2049 * tmp545;
    let tmp2895 = (1_f64 / 3.0) * tmp2894;
    let tmp2896 = tmp1579 * tmp2102;
    let tmp2897 = 0.333333333333333 * tmp2108;
    let tmp2898 = 0.333333333333333 * tmp2110;
    let tmp2899 = tmp1666 + tmp1669 + tmp2897 + tmp2898;
    let tmp2900 = tmp2899 * tmp382;
    let tmp2901 = 0.333333333333333 * tmp2126;
    let tmp2902 = 0.333333333333333 * tmp2127;
    let tmp2903 = -1.0 * tmp2902;
    let tmp2904 = tmp1594 + tmp2901 + tmp2903;
    let tmp2905 = tmp2904 * tmp312;
    let tmp2906 = 0.333333333333333 * tmp2113;
    let tmp2907 = -1.0 * tmp2906;
    let tmp2908 = 0.333333333333333 * tmp2114;
    let tmp2909 = tmp1583 * tmp2115;
    let tmp2910 = -1.0 * tmp2909;
    let tmp2911 = tmp1589 * tmp217;
    let tmp2912 = -1.0 * tmp2911;
    let tmp2913 = -1.0 * tmp2908 + tmp2910 + tmp2912;
    let tmp2914 = tmp1589 * tmp219;
    let tmp2915 = tmp1583 * tmp399;
    let tmp2916 = tmp2914 + tmp2915;
    let tmp2917 = tmp2907 + tmp2913 + tmp2916;
    let tmp2918 = tmp254 * tmp2917;
    let tmp2919 = -1.0 * tmp2897;
    let tmp2920 = tmp1600 + tmp1604 + tmp1605 + tmp1668;
    let tmp2921 = tmp2898 + tmp2919 + tmp2920;
    let tmp2922 = tmp2921 * tmp359;
    let tmp2923 = tmp1629 * tmp2177;
    let tmp2924 = tmp1631 * tmp2175;
    let tmp2925 = tmp1635 * tmp2093;
    let tmp2926 = tmp1640 * tmp2068;
    let tmp2927 = tmp1643 * tmp2086;
    let tmp2928 = tmp1649 * tmp2056;
    let tmp2929 = tmp2064 * tmp537;
    let tmp2930 = tmp1651 * tmp2100
        + tmp2900
        + tmp2905
        + tmp2918
        + tmp2922
        + tmp2923
        + tmp2924
        + tmp2925
        + tmp2926
        + tmp2927
        + tmp2928
        - 1_f64 / 6.0 * tmp2929;
    let tmp2931 = -1.0 * tmp2898;
    let tmp2932 = tmp1601 + tmp1607 + tmp2919 + tmp2931;
    let tmp2933 = tmp212 * tmp2932;
    let tmp2934 = tmp2909 + tmp2915;
    let tmp2935 = tmp2911 + tmp2914;
    let tmp2936 = tmp2907 + tmp2908 + tmp2934 + tmp2935;
    let tmp2937 = tmp2936 * tmp300;
    let tmp2938 = -1.0 * tmp2901;
    let tmp2939 = tmp1656 + tmp2902 + tmp2938;
    let tmp2940 = tmp2939 * tmp345;
    let tmp2941 = tmp1588 + tmp1654 + tmp2901 + tmp2902;
    let tmp2942 = tmp2941 * tmp376;
    let tmp2943 = tmp1675 * tmp2175;
    let tmp2944 = tmp1677 * tmp2172;
    let tmp2945 = tmp1679 * tmp2065;
    let tmp2946 = tmp1681 * tmp2091;
    let tmp2947 = tmp1684 * tmp2082;
    let tmp2948 = tmp1687 * tmp2050;
    let tmp2949 = -1.0 * tmp1651 * tmp2098
        + tmp1689 * tmp2097
        + tmp2933
        + tmp2937
        + tmp2940
        + tmp2942
        + tmp2943
        + tmp2944
        + tmp2945
        + tmp2946
        + tmp2947
        + tmp2948;
    let tmp2950 = tmp2870 + tmp2876 + tmp2881 + tmp2887 + tmp2888 + tmp2889 + tmp2890 + tmp2891 + tmp2892 + tmp2893
        - 1.0 * tmp2895
        + tmp2896
        + tmp2930
        + tmp2949;
    let tmp2951 = -1.0 * tmp2882;
    let tmp2952 = tmp1770 + tmp1772 + tmp2884 + tmp2951;
    let tmp2953 = tmp2952 * tmp382;
    let tmp2954 = -1.0 * tmp2878;
    let tmp2955 = tmp1699 + tmp1701 + tmp2877 + tmp2954;
    let tmp2956 = tmp2955 * tmp312;
    let tmp2957 = -1.0 * tmp2860;
    let tmp2958 = -1.0 * tmp2861;
    let tmp2959 = tmp2957 + tmp2958;
    let tmp2960 = tmp2859 + tmp2959;
    let tmp2961 = tmp2868 + tmp2960;
    let tmp2962 = tmp254 * tmp2961;
    let tmp2963 = tmp1526 + tmp1532 + tmp1769;
    let tmp2964 = tmp2883 + tmp2951 + tmp2963;
    let tmp2965 = tmp2964 * tmp359;
    let tmp2966 = tmp1721 * tmp2175;
    let tmp2967 = tmp1726 * tmp2177;
    let tmp2968 = tmp1728 * tmp2056;
    let tmp2969 = tmp1731 * tmp2086;
    let tmp2970 = tmp1734 * tmp2068;
    let tmp2971 = tmp1736 * tmp2093;
    let tmp2972 = tmp1692 * tmp2100;
    let tmp2973 = (1_f64 / 3.0) * tmp2929;
    let tmp2974 = -1.0 * tmp2914;
    let tmp2975 = -1.0 * tmp2915;
    let tmp2976 = tmp2906 + tmp2974 + tmp2975;
    let tmp2977 = tmp2913 + tmp2976;
    let tmp2978 = tmp2977 * tmp368;
    let tmp2979 = tmp2909 + tmp2911;
    let tmp2980 = tmp2908 + tmp2976 + tmp2979;
    let tmp2981 = tmp2980 * tmp387;
    let tmp2982 = tmp1593 + tmp1655 + tmp2903 + tmp2938;
    let tmp2983 = tmp278 * tmp2982;
    let tmp2984 = tmp1598 + tmp1599 + tmp1606 + tmp1665;
    let tmp2985 = tmp2897 + tmp2931 + tmp2984;
    let tmp2986 = tmp2985 * tmp322;
    let tmp2987 = tmp1749 * tmp2172;
    let tmp2988 = tmp1751 * tmp2177;
    let tmp2989 = tmp1753 * tmp2095;
    let tmp2990 = tmp1755 * tmp2072;
    let tmp2991 = tmp1757 * tmp2089;
    let tmp2992 = tmp1759 * tmp2061;
    let tmp2993 = -1.0 * tmp1689 * tmp2102
        + (1_f64 / 6.0) * tmp2894
        + tmp2978
        + tmp2981
        + tmp2983
        + tmp2986
        + tmp2987
        + tmp2988
        + tmp2989
        + tmp2990
        + tmp2991
        + tmp2992;
    let tmp2994 =
        tmp2949 + tmp2953 + tmp2956 + tmp2962 + tmp2965 + tmp2966 + tmp2967 + tmp2968 + tmp2969 + tmp2970 + tmp2971
            - 1.0 * tmp2972
            + tmp2973
            + tmp2993;
    let tmp2995 = tmp1705 + tmp1707 + tmp2882 + tmp2883;
    let tmp2996 = tmp212 * tmp2995;
    let tmp2997 = tmp2874 + tmp2960;
    let tmp2998 = tmp2997 * tmp300;
    let tmp2999 = tmp1541 - 1.0 * tmp2877;
    let tmp3000 = tmp1764 + tmp2879 + tmp2999;
    let tmp3001 = tmp3000 * tmp345;
    let tmp3002 = tmp1539 + tmp1548 + tmp2954 + tmp2999;
    let tmp3003 = tmp3002 * tmp376;
    let tmp3004 = tmp1780 * tmp2172;
    let tmp3005 = tmp1785 * tmp2175;
    let tmp3006 = tmp1787 * tmp2050;
    let tmp3007 = tmp1790 * tmp2082;
    let tmp3008 = tmp1792 * tmp2065;
    let tmp3009 = tmp1794 * tmp2091;
    let tmp3010 = tmp1579 * tmp2097;
    let tmp3011 = tmp1692 * tmp2098;
    let tmp3012 = tmp2930
        + tmp2993
        + tmp2996
        + tmp2998
        + tmp3001
        + tmp3003
        + tmp3004
        + tmp3005
        + tmp3006
        + tmp3007
        + tmp3008
        + tmp3009
        - 1.0 * tmp3010
        + tmp3011;
    let tmp3013 = 0.666666666666667 * tmp2353;
    let tmp3014 = tmp1885 * tmp2109;
    let tmp3015 = -1.0 * tmp3014;
    let tmp3016 = tmp1504 * tmp2106;
    let tmp3017 = tmp231 * tmp3016;
    let tmp3018 = -1.0 * tmp3017;
    let tmp3019 = tmp3015 + tmp3018;
    let tmp3020 = tmp3013 + tmp3019;
    let tmp3021 = 0.666666666666667 * tmp2361;
    let tmp3022 = tmp227 * tmp3016;
    let tmp3023 = tmp1506 * tmp2107;
    let tmp3024 = tmp3022 + tmp3023;
    let tmp3025 = tmp3021 + tmp3024;
    let tmp3026 = 0.666666666666667 * tmp2383;
    let tmp3027 = tmp1875 * tmp618;
    let tmp3028 = tmp1504 * tmp754;
    let tmp3029 = tmp3027 - 1.0 * tmp3028;
    let tmp3030 = tmp1875 * tmp620;
    let tmp3031 = -1.0 * tmp3030;
    let tmp3032 = tmp1889 * tmp217;
    let tmp3033 = tmp3031 - 1.0 * tmp3032;
    let tmp3034 = tmp3029 + tmp3033;
    let tmp3035 = tmp263 * tmp3016;
    let tmp3036 = tmp1506 * tmp2109;
    let tmp3037 = -1.0 * tmp3036;
    let tmp3038 = tmp3035 + tmp3037;
    let tmp3039 = -1.0 * tmp3026 + tmp3034 + tmp3038;
    let tmp3040 = 0.666666666666667 * tmp2379;
    let tmp3041 = tmp1875 * tmp616;
    let tmp3042 = -1.0 * tmp3041;
    let tmp3043 = tmp1504 * tmp750;
    let tmp3044 = tmp3042 - 1.0 * tmp3043;
    let tmp3045 = tmp1875 * tmp612;
    let tmp3046 = tmp1889 * tmp219;
    let tmp3047 = tmp3045 - 1.0 * tmp3046;
    let tmp3048 = tmp3044 + tmp3047;
    let tmp3049 = tmp256 * tmp3016;
    let tmp3050 = tmp1885 * tmp2107;
    let tmp3051 = -1.0 * tmp3050;
    let tmp3052 = tmp3049 + tmp3051;
    let tmp3053 = tmp3040 + tmp3048 + tmp3052;
    let tmp3054 = 0.666666666666667 * tmp2336;
    let tmp3055 = tmp1504 * tmp2337;
    let tmp3056 = tmp1504 * tmp2339;
    let tmp3057 = -1.0 * tmp3056;
    let tmp3058 = tmp3055 + tmp3057;
    let tmp3059 = tmp3054 + tmp3058;
    let tmp3060 = 0.666666666666667 * tmp2344;
    let tmp3061 = tmp283 * tmp3016;
    let tmp3062 = tmp281 * tmp3016;
    let tmp3063 = -1.0 * tmp3062;
    let tmp3064 = tmp3061 + tmp3063;
    let tmp3065 = -1.0 * tmp3060 + tmp3064;
    let tmp3066 = -1.0 * tmp3055;
    let tmp3067 = tmp3056 + tmp3066;
    let tmp3068 = -1.0 * tmp3054 + tmp3067;
    let tmp3069 = tmp1583 * tmp222;
    let tmp3070 = tmp281 * tmp3069;
    let tmp3071 = tmp1583 * tmp673;
    let tmp3072 = -1.0 * tmp3071;
    let tmp3073 = tmp1583 * tmp675;
    let tmp3074 = tmp283 * tmp3069;
    let tmp3075 = -1.0 * tmp3074;
    let tmp3076 = 0.333333333333333 * tmp2353;
    let tmp3077 = q_old[0] * tmp1583;
    let tmp3078 = tmp2109 * tmp3077;
    let tmp3079 = tmp1583 * tmp2106;
    let tmp3080 = tmp231 * tmp3079;
    let tmp3081 = tmp3076 - 1.0 * tmp3078 - 1.0 * tmp3080;
    let tmp3082 = 0.333333333333333 * tmp2361;
    let tmp3083 = tmp227 * tmp3079;
    let tmp3084 = tmp1584 * tmp2107;
    let tmp3085 = tmp3082 + tmp3083 + tmp3084;
    let tmp3086 = tmp1583 * tmp2337;
    let tmp3087 = tmp1583 * tmp2339;
    let tmp3088 = -1.0 * tmp3087;
    let tmp3089 = tmp283 * tmp3079;
    let tmp3090 = tmp281 * tmp3079;
    let tmp3091 = -1.0 * tmp3090;
    let tmp3092 = 0.333333333333333 * tmp2336;
    let tmp3093 = tmp1584 * tmp223;
    let tmp3094 = tmp1583 * tmp636;
    let tmp3095 = tmp3092 + tmp3093 - 1.0 * tmp3094;
    let tmp3096 = 0.333333333333333 * tmp2344;
    let tmp3097 = tmp227 * tmp3069;
    let tmp3098 = tmp231 * tmp3069;
    let tmp3099 = -1.0 * tmp3096 + tmp3097 - 1.0 * tmp3098;
    let tmp3100 = 0.333333333333333 * tmp2379;
    let tmp3101 = -1.0 * tmp3100;
    let tmp3102 = tmp2107 * tmp3077;
    let tmp3103 = tmp256 * tmp3079;
    let tmp3104 = -1.0 * tmp3103;
    let tmp3105 = 0.333333333333333 * tmp2383;
    let tmp3106 = 0.0833333333333333 * tmp1197;
    let tmp3107 = tmp3106 * tmp620;
    let tmp3108 = -1.0 * tmp3107;
    let tmp3109 = tmp217 * tmp3069;
    let tmp3110 = tmp3108 - 1.0 * tmp3109;
    let tmp3111 = tmp3106 * tmp618;
    let tmp3112 = tmp1583 * tmp754;
    let tmp3113 = tmp3111 - 1.0 * tmp3112;
    let tmp3114 = tmp263 * tmp3079;
    let tmp3115 = tmp1584 * tmp2109;
    let tmp3116 = -1.0 * tmp3115;
    let tmp3117 = tmp3114 + tmp3116;
    let tmp3118 = -1.0 * tmp3105 + tmp3110 + tmp3113 + tmp3117;
    let tmp3119 = tmp219 * tmp3069;
    let tmp3120 = tmp1583 * tmp750;
    let tmp3121 = tmp3106 * tmp616;
    let tmp3122 = tmp3106 * tmp612;
    let tmp3123 = -1.0 * tmp3122;
    let tmp3124 = tmp3121 + tmp3123;
    let tmp3125 = tmp3119 + tmp3120 + tmp3124;
    let tmp3126 = tmp3072 + tmp3074;
    let tmp3127 = -1.0 * tmp3073;
    let tmp3128 = tmp3070 + tmp3127;
    let tmp3129 = tmp3126 + tmp3128;
    let tmp3130 = -1.0 * tmp3082 - 1.0 * tmp3083 - 1.0 * tmp3084;
    let tmp3131 = tmp1595 * tmp2254
        + tmp1609 * tmp2191
        + tmp1622 * tmp2185
        + tmp1627 * tmp2248
        + tmp1635 * tmp2327
        + tmp1640 * tmp2280
        + tmp1643 * tmp2300
        + tmp1649 * tmp2238
        + tmp1944 * tmp2445
        + tmp1945 * tmp2430
        + tmp1947 * tmp2213
        + tmp2186 * tmp2917
        + tmp2192 * tmp2904
        + tmp2249 * tmp2921
        + tmp2255 * tmp2899
        - 1_f64 / 3.0 * tmp2302
        - 1_f64 / 3.0 * tmp2306
        - 1_f64 / 3.0 * tmp2307
        + (1_f64 / 3.0) * tmp2312
        + (1_f64 / 3.0) * tmp2313
        + tmp646 * (tmp3101 + tmp3102 + tmp3104 + tmp3118 + tmp3125)
        + tmp700 * (tmp1880 + tmp3086 + tmp3088 + tmp3089 + tmp3091 + tmp3095 + tmp3099)
        + tmp722 * (tmp3081 + tmp3129 + tmp3130)
        + tmp740 * (tmp2010 + tmp3070 + tmp3072 + tmp3073 + tmp3075 + tmp3081 + tmp3085);
    let tmp3132 = -1.0 * tmp3119 + tmp3122;
    let tmp3133 = -1.0 * tmp3121;
    let tmp3134 = -1.0 * tmp3120 + tmp3133;
    let tmp3135 = tmp3132 + tmp3134;
    let tmp3136 = -1.0 * tmp3102;
    let tmp3137 = tmp3103 + tmp3136;
    let tmp3138 = tmp3100 + tmp3135 + tmp3137;
    let tmp3139 = -1.0 * tmp3114;
    let tmp3140 = -1.0 * tmp3111;
    let tmp3141 = tmp3107 + tmp3140;
    let tmp3142 = tmp3109 + tmp3112 + tmp3141;
    let tmp3143 = tmp3087 + tmp3091;
    let tmp3144 = -1.0 * tmp3086;
    let tmp3145 = tmp3089 + tmp3144;
    let tmp3146 = tmp3143 + tmp3145;
    let tmp3147 = -1.0 * tmp3092 - 1.0 * tmp3093 + tmp3094;
    let tmp3148 = -1.0 * tmp3070;
    let tmp3149 = tmp3073 + tmp3148;
    let tmp3150 = tmp3071 + tmp3075;
    let tmp3151 = tmp3149 + tmp3150;
    let tmp3152 = -1.0 * tmp3076 + tmp3078 + tmp3080;
    let tmp3153 = tmp1741 * tmp2251
        + tmp1743 * tmp2257
        + tmp1745 * tmp2188
        + tmp1747 * tmp2193
        + tmp1753 * tmp2328
        + tmp1755 * tmp2232
        + tmp1757 * tmp2326
        + tmp1759 * tmp2214
        + tmp1973 * tmp2447
        + tmp1974 * tmp2445
        - 1.0 * tmp1975 * tmp2279
        + tmp2187 * tmp2982
        + tmp2194 * tmp2985
        + tmp2250 * tmp2977
        + tmp2256 * tmp2980
        - 1_f64 / 3.0 * tmp2308
        - 1_f64 / 3.0 * tmp2309
        + (1_f64 / 3.0) * tmp2314
        + (1_f64 / 3.0) * tmp2315
        + (1_f64 / 3.0) * tmp2330
        + tmp671 * (tmp3099 + tmp3146 + tmp3147)
        + tmp706 * (tmp3085 + tmp3151 + tmp3152)
        + tmp729 * (tmp3118 + tmp3138)
        + tmp744 * (tmp3105 + tmp3115 + tmp3138 + tmp3139 + tmp3142);
    let tmp3154 = -1.0 * tmp3022;
    let tmp3155 = -1.0 * tmp3023;
    let tmp3156 = tmp3154 + tmp3155;
    let tmp3157 = -1.0 * tmp3021 + tmp3156;
    let tmp3158 = tmp3014 + tmp3017;
    let tmp3159 = -1.0 * tmp3013 + tmp3158;
    let tmp3160 = -1.0 * tmp3061;
    let tmp3161 = tmp3062 + tmp3160;
    let tmp3162 = tmp3060 + tmp3161;
    let tmp3163 = -1.0 * tmp3035;
    let tmp3164 = tmp3030 + tmp3032;
    let tmp3165 = -1.0 * tmp3027;
    let tmp3166 = tmp3028 + tmp3165;
    let tmp3167 = tmp3164 + tmp3166;
    let tmp3168 = tmp3026 + tmp3167;
    let tmp3169 = tmp1898 + tmp1904 + tmp2014;
    let tmp3170 = tmp3112 + tmp3120 + tmp3121 + tmp3140;
    let tmp3171 = tmp3107 + tmp3109 + tmp3119 + tmp3123;
    let tmp3172 = tmp3104 + tmp3139;
    let tmp3173 = tmp3102 + tmp3115;
    let tmp3174 = tmp3172 + tmp3173;
    let tmp3175 = -1.0 * tmp3089;
    let tmp3176 = tmp3096 - 1.0 * tmp3097 + tmp3098;
    let tmp3177 = tmp3086 + tmp3175;
    let tmp3178 = tmp3088 + tmp3090;
    let tmp3179 = tmp3177 + tmp3178;
    let tmp3180 = tmp1657 * tmp2183
        + tmp1663 * tmp2190
        + tmp1671 * tmp2247
        + tmp1673 * tmp2253
        + tmp1679 * tmp2245
        + tmp1681 * tmp2301
        + tmp1684 * tmp2321
        + tmp1687 * tmp2276
        + tmp2020 * tmp2430
        + tmp2021 * tmp2447
        - 1.0 * tmp2022 * tmp2231
        + tmp2023 * tmp2237
        + tmp2184 * tmp2932
        + tmp2189 * tmp2936
        + tmp2246 * tmp2939
        + tmp2252 * tmp2941
        - 1_f64 / 3.0 * tmp2304
        - 1_f64 / 3.0 * tmp2305
        + (1_f64 / 3.0) * tmp2310
        + (1_f64 / 3.0) * tmp2311
        + tmp601 * (tmp1910 + tmp3071 + tmp3074 + tmp3127 + tmp3130 + tmp3148 + tmp3152)
        + tmp693 * (tmp3101 + tmp3105 + tmp3170 + tmp3171 + tmp3174)
        + tmp716 * (tmp1999 + tmp3087 + tmp3090 + tmp3144 + tmp3147 + tmp3175 + tmp3176)
        + tmp736 * (tmp3095 + tmp3176 + tmp3179);
    let tmp3181 = -1.0 * tmp3045;
    let tmp3182 = tmp3041 + tmp3043 + tmp3046 + tmp3181;
    let tmp3183 = -1.0 * tmp3040 + tmp3182;
    let tmp3184 = tmp3036 + tmp3050;
    let tmp3185 = -1.0 * tmp3049;
    let tmp3186 = tmp3163 + tmp3185;
    let tmp3187 = tmp3184 + tmp3186;
    let tmp3188 = tmp1900 + tmp1902 + tmp2012;
    let tmp3189 = -1.0 * tmp100 + tmp103 + tmp93 - 1.0 * tmp96;
    let tmp3190 = tmp29 * tmp3189;
    let tmp3191 = 2.0 * tmp3190;
    let tmp3192 = tmp21 * tmp3191;
    let tmp3193 = tmp25 * tmp3191;
    let tmp3194 = tmp2063 + tmp23 * tmp3193;
    let tmp3195 = -1.0 * tmp3192 + tmp3194 + tmp74 + tmp82;
    let tmp3196 = grad_phi_node_i[0] * tmp3195;
    let tmp3197 = tmp18 * tmp3196;
    let tmp3198 = tmp23 * tmp3191;
    let tmp3199 = tmp19 * tmp3198;
    let tmp3200 = tmp20 * tmp3193;
    let tmp3201 = tmp3199 - 1.0 * tmp3200 + tmp47 + tmp50 + tmp56 + tmp59;
    let tmp3202 = grad_phi_node_i[1] * tmp3201;
    let tmp3203 = tmp3202 * tmp67;
    let tmp3204 = tmp19 * tmp3193;
    let tmp3205 = tmp20 * tmp3198;
    let tmp3206 = tmp159 + tmp3204 - 1.0 * tmp3205;
    let tmp3207 = grad_phi_node_i[2] * tmp3206;
    let tmp3208 = tmp3207 * tmp89;
    let tmp3209 = tmp133 + tmp138 + tmp158 + tmp3204 + tmp3205;
    let tmp3210 = grad_phi_node_i[0] * tmp3209;
    let tmp3211 = tmp110 * tmp3210;
    let tmp3212 = tmp3192 + tmp3194 + tmp75 + tmp78 + tmp80;
    let tmp3213 = grad_phi_node_i[1] * tmp3212;
    let tmp3214 = tmp119 * tmp3213;
    let tmp3215 = tmp2077 + tmp3199 + tmp3200;
    let tmp3216 = grad_phi_node_i[2] * tmp3215;
    let tmp3217 = tmp125 * tmp3216;
    let tmp3218 = tmp27 * tmp3190;
    let tmp3219 = tmp26 * tmp3190;
    let tmp3220 = -1.0 * tmp3219;
    let tmp3221 = tmp24 * tmp3190;
    let tmp3222 = tmp22 * tmp3190;
    let tmp3223 = -1.0 * tmp3222 - 1.0 * tmp97;
    let tmp3224 = tmp105 + tmp112 + tmp3218 + tmp3220 + tmp3221 + tmp3223;
    let tmp3225 = grad_phi_node_i[0] * tmp3224;
    let tmp3226 = tmp132 * tmp3225;
    let tmp3227 = tmp3218 - 1.0 * tmp3221;
    let tmp3228 = tmp114 + tmp3220 + tmp3222 + tmp3227;
    let tmp3229 = grad_phi_node_i[1] * tmp3228;
    let tmp3230 = tmp150 * tmp3229;
    let tmp3231 = tmp113 + tmp2070 + tmp3219 + tmp3223 + tmp3227;
    let tmp3232 = grad_phi_node_i[2] * tmp3231;
    let tmp3233 = tmp157 * tmp3232;
    let tmp3234 = grad_phi_node_i[0] * tmp3228;
    let tmp3235 = tmp163 * tmp3234;
    let tmp3236 = grad_phi_node_i[1] * tmp3231;
    let tmp3237 = tmp166 * tmp3236;
    let tmp3238 = grad_phi_node_i[2] * tmp3224;
    let tmp3239 = tmp169 * tmp3238;
    let tmp3240 = grad_phi_node_i[0] * tmp3201;
    let tmp3241 = grad_phi_node_i[0] * tmp3215;
    let tmp3242 = grad_phi_node_i[1] * tmp3209;
    let tmp3243 = grad_phi_node_i[1] * tmp3206;
    let tmp3244 = grad_phi_node_i[2] * tmp3195;
    let tmp3245 = grad_phi_node_i[2] * tmp3212;
    let tmp3246 = q_old[0] * tmp7;
    let tmp3247 = q_old[3] * tmp5;
    let tmp3248 = -1.0 * tmp214 * tmp3246 - 1.0 * tmp214 * tmp3247 + tmp639 + tmp641;
    let tmp3249 = tmp213 * tmp3248;
    let tmp3250 = tmp3249 * tmp7;
    let tmp3251 = tmp10 * tmp3250;
    let tmp3252 = tmp3249 * tmp5;
    let tmp3253 = tmp3252 * tmp4;
    let tmp3254 = tmp2117 + tmp3253;
    let tmp3255 = tmp2119 + tmp2132 + tmp3251 + tmp3254;
    let tmp3256 = tmp212 * tmp3255;
    let tmp3257 = tmp10 * tmp3252;
    let tmp3258 = tmp3250 * tmp4;
    let tmp3259 = tmp2141 + tmp3257 + tmp3258;
    let tmp3260 = tmp254 * tmp3259;
    let tmp3261 = tmp3252 * tmp7;
    let tmp3262 = tmp3249 * tmp87;
    let tmp3263 = tmp266 + tmp3261 + tmp3262 + tmp370;
    let tmp3264 = tmp278 * tmp3263;
    let tmp3265 = -1.0 * tmp3258;
    let tmp3266 = tmp286 + tmp293 + tmp3257 + tmp3265;
    let tmp3267 = tmp300 * tmp3266;
    let tmp3268 = -1.0 * tmp3261;
    let tmp3269 = tmp265 + tmp304 + tmp3262 + tmp3268 + tmp369;
    let tmp3270 = tmp312 * tmp3269;
    let tmp3271 = -1.0 * tmp3253;
    let tmp3272 = tmp2151 + tmp2159 + tmp3251 + tmp3271;
    let tmp3273 = tmp322 * tmp3272;
    let tmp3274 = -1.0 * tmp3262;
    let tmp3275 = tmp259 + tmp262 + tmp302 + tmp3261 + tmp3274;
    let tmp3276 = tmp3275 * tmp345;
    let tmp3277 = tmp2116 - 1.0 * tmp3251;
    let tmp3278 = tmp2123 + tmp3254 + tmp3277;
    let tmp3279 = tmp3278 * tmp359;
    let tmp3280 = -1.0 * tmp3257;
    let tmp3281 = tmp314 + tmp3258 + tmp3280 + tmp377;
    let tmp3282 = tmp3281 * tmp368;
    let tmp3283 = tmp305 + tmp3268 + tmp3274;
    let tmp3284 = tmp3283 * tmp376;
    let tmp3285 = tmp2121 + tmp2134 + tmp3271 + tmp3277;
    let tmp3286 = tmp3285 * tmp382;
    let tmp3287 = tmp2148 + tmp3265 + tmp3280;
    let tmp3288 = tmp3287 * tmp387;
    let tmp3289 = tmp11 * tmp3249;
    let tmp3290 = -1.0 * tmp13 * tmp641;
    let tmp3291 = tmp3249 * tmp8;
    let tmp3292 = q_old[0] * tmp394;
    let tmp3293 = tmp3249 * tmp6;
    let tmp3294 = tmp14 * tmp634;
    let tmp3295 = -1.0 * tmp3293 - 1.0 * tmp3294;
    let tmp3296 = tmp3249 * tmp9;
    let tmp3297 = tmp13 * tmp639;
    let tmp3298 = -1.0 * tmp3296 + tmp3297;
    let tmp3299 = tmp3289 + tmp3290 + tmp3291 + tmp3292 + tmp3295 + tmp3298;
    let tmp3300 = tmp3299 * tmp392;
    let tmp3301 = tmp3289 + tmp3290 - 1.0 * tmp3291 - 1.0 * tmp3292;
    let tmp3302 = tmp3293 + tmp3294 + tmp3298 + tmp3301;
    let tmp3303 = tmp3302 * tmp409;
    let tmp3304 = tmp3295 + tmp3296 - 1.0 * tmp3297 + tmp3301;
    let tmp3305 = tmp3304 * tmp415;
    let tmp3306 = tmp3302 * tmp420;
    let tmp3307 = tmp3304 * tmp424;
    let tmp3308 = tmp3299 * tmp428;
    let tmp3309 = -1.0 * tmp179 * tmp3240 + tmp179 * tmp3245 + tmp183 * tmp3241 - 1.0 * tmp183 * tmp3243
        + tmp186 * tmp3242
        - 1.0 * tmp186 * tmp3244
        + tmp3197
        + tmp3203
        + tmp3208
        + tmp3211
        + tmp3214
        + tmp3217
        + tmp3226
        + tmp3230
        + tmp3233
        + tmp3235
        + tmp3237
        + tmp3239
        + tmp3256
        + tmp3260
        + tmp3264
        + tmp3267
        + tmp3270
        + tmp3273
        + tmp3276
        + tmp3279
        + tmp3282
        + tmp3284
        + tmp3286
        + tmp3288
        + tmp3300
        + tmp3303
        + tmp3305
        + tmp3306
        + tmp3307
        + tmp3308;
    let tmp3310 = 2.0 * tmp3196;
    let tmp3311 = 2.0 * tmp3202;
    let tmp3312 = 2.0 * tmp3207;
    let tmp3313 = 2.0 * tmp3210;
    let tmp3314 = 2.0 * tmp3213;
    let tmp3315 = 2.0 * tmp3216;
    let tmp3316 = tmp30 * tmp3190;
    let tmp3317 = tmp3316 * tmp77;
    let tmp3318 = -1.0 * tmp497 + tmp499 + tmp502 - 1.0 * tmp503;
    let tmp3319 = tmp3318 * tmp448;
    let tmp3320 = 2.0 * tmp3319;
    let tmp3321 = tmp25 * tmp3320;
    let tmp3322 = tmp20 * tmp3321;
    let tmp3323 = tmp3316 * tmp70;
    let tmp3324 = tmp3316 * tmp79;
    let tmp3325 = tmp3323 - 1.0 * tmp3324;
    let tmp3326 = tmp23 * tmp3320;
    let tmp3327 = tmp3316 * tmp73;
    let tmp3328 = tmp19 * tmp3326 + tmp2286 - 1.0 * tmp3327;
    let tmp3329 = tmp2294 + tmp2322 + tmp3317 - 1.0 * tmp3322 + tmp3325 + tmp3328;
    let tmp3330 = tmp3329 * tmp495;
    let tmp3331 = tmp3190 * tmp58;
    let tmp3332 = tmp21 * tmp3320;
    let tmp3333 = tmp3190 * tmp55;
    let tmp3334 = tmp3190 * tmp49;
    let tmp3335 = -1.0 * tmp3334;
    let tmp3336 = -1.0 * tmp3333 + tmp3335;
    let tmp3337 = tmp3190 * tmp46;
    let tmp3338 = tmp2208 + tmp23 * tmp3321 + tmp3337;
    let tmp3339 = -1.0 * tmp2204 + tmp2240 + tmp3331 + tmp3332 + tmp3336 + tmp3338;
    let tmp3340 = tmp3339 * tmp495;
    let tmp3341 = -1.0 * tmp3331;
    let tmp3342 = tmp3333 + tmp3341;
    let tmp3343 = tmp3335 + tmp3342 + tmp485 + tmp491;
    let tmp3344 = tmp2206 + tmp2239 - 1.0 * tmp3332 + tmp3338 + tmp3343;
    let tmp3345 = tmp3344 * tmp476;
    let tmp3346 = -1.0 * tmp3317;
    let tmp3347 = tmp3324 + tmp3346;
    let tmp3348 = tmp2297 + tmp2319 + tmp3322 + tmp3323 + tmp3328 + tmp3347;
    let tmp3349 = tmp3348 * tmp513;
    let tmp3350 = 2.0 * tmp3225;
    let tmp3351 = 2.0 * tmp3229;
    let tmp3352 = 2.0 * tmp3232;
    let tmp3353 = 2.0 * tmp3234;
    let tmp3354 = 2.0 * tmp3236;
    let tmp3355 = 2.0 * tmp3238;
    let tmp3356 = tmp3339 * tmp539;
    let tmp3357 = tmp3329 * tmp541;
    let tmp3358 = tmp3190 * tmp96;
    let tmp3359 = tmp100 * tmp3190;
    let tmp3360 = tmp3190 * tmp93;
    let tmp3361 = tmp20 * tmp3326;
    let tmp3362 = tmp103 * tmp3190;
    let tmp3363 = tmp19 * tmp3321 + tmp3362;
    let tmp3364 = tmp2263 + tmp2278 + tmp3358 + tmp3359 + tmp3360 - 1.0 * tmp3361 + tmp3363;
    let tmp3365 = tmp3364 * tmp513;
    let tmp3366 = tmp3358 - 1.0 * tmp3359;
    let tmp3367 = -1.0 * tmp3360 + tmp3366 + tmp580;
    let tmp3368 = -1.0 * tmp2263 + tmp2270 + tmp2277 + tmp3361 + tmp3363 + tmp3367;
    let tmp3369 = tmp3368 * tmp476;
    let tmp3370 = tmp24 * tmp3319;
    let tmp3371 = tmp3190 * tmp34;
    let tmp3372 = tmp27 * tmp3319;
    let tmp3373 = tmp26 * tmp3319;
    let tmp3374 = tmp3190 * tmp36;
    let tmp3375 = tmp3372 - 1.0 * tmp3373 - 1.0 * tmp3374;
    let tmp3376 = tmp22 * tmp3319;
    let tmp3377 = tmp3190 * tmp38;
    let tmp3378 = -1.0 * tmp3377;
    let tmp3379 = tmp3190 * tmp32;
    let tmp3380 = tmp3378 - 1.0 * tmp3379;
    let tmp3381 = -1.0 * tmp2226 - 1.0 * tmp3376 + tmp3380;
    let tmp3382 = tmp2227 + tmp2235 + tmp3370 - 1.0 * tmp3371 + tmp3375 + tmp3381;
    let tmp3383 = tmp3382 * tmp476;
    let tmp3384 = tmp3382 * tmp513;
    let tmp3385 = tmp3348 * tmp535;
    let tmp3386 = tmp3344 * tmp545;
    let tmp3387 = tmp3241 * tmp411;
    let tmp3388 = tmp184 * tmp3302;
    let tmp3389 = tmp3242 * tmp416;
    let tmp3390 = tmp187 * tmp3304;
    let tmp3391 = tmp190 * tmp3299;
    let tmp3392 = tmp3245 * tmp405;
    let tmp3393 = tmp3240 * tmp405;
    let tmp3394 = tmp180 * tmp3299;
    let tmp3395 = tmp3243 * tmp411;
    let tmp3396 = tmp188 * tmp3302;
    let tmp3397 = tmp189 * tmp3304;
    let tmp3398 = tmp3244 * tmp416;
    let tmp3399 = -1.0 * tmp3370;
    let tmp3400 = tmp3371 + tmp3379;
    let tmp3401 = tmp2229 + tmp2236 + tmp3375 + tmp3376 + tmp3378 + tmp3399 + tmp3400;
    let tmp3402 = tmp3401 * tmp495;
    let tmp3403 = -1.0 * tmp505;
    let tmp3404 = tmp3403 + tmp506;
    let tmp3405 = tmp3371 + tmp3404;
    let tmp3406 = tmp2220 + tmp2228 + tmp2234 + tmp3372 + tmp3373 + tmp3374 + tmp3381 + tmp3399 + tmp3405;
    let tmp3407 = tmp3406 * tmp513;
    let tmp3408 = tmp3401 * tmp476;
    let tmp3409 = tmp3406 * tmp495;
    let tmp3410 = tmp3368 * tmp537;
    let tmp3411 = tmp3364 * tmp543;
    let tmp3412 = tmp227 * tmp602;
    let tmp3413 = tmp231 * tmp602;
    let tmp3414 = -1.0 * tmp3246 * tmp602 - 1.0 * tmp3247 * tmp602 + tmp3412 + tmp3413;
    let tmp3415 = tmp3414 * tmp606;
    let tmp3416 = tmp3415 * tmp7;
    let tmp3417 = tmp3416 * tmp4;
    let tmp3418 = tmp3415 * tmp5;
    let tmp3419 = tmp10 * tmp3418;
    let tmp3420 = tmp2354 * tmp3252;
    let tmp3421 = -1.0 * tmp3420;
    let tmp3422 = tmp228 * tmp3249;
    let tmp3423 = tmp263 * tmp3422;
    let tmp3424 = -1.0 * tmp3423;
    let tmp3425 = tmp3421 + tmp3424;
    let tmp3426 = tmp256 * tmp3422;
    let tmp3427 = tmp235 * tmp3250;
    let tmp3428 = tmp3426 + tmp3427;
    let tmp3429 = tmp3415 * tmp87;
    let tmp3430 = tmp3418 * tmp7;
    let tmp3431 = -1.0 * tmp3430;
    let tmp3432 = tmp281 * tmp3422;
    let tmp3433 = q_old[1] * tmp3250;
    let tmp3434 = tmp228 * tmp3433;
    let tmp3435 = -1.0 * tmp3434;
    let tmp3436 = tmp3432 + tmp3435;
    let tmp3437 = q_old[2] * tmp3252;
    let tmp3438 = tmp228 * tmp3437;
    let tmp3439 = tmp283 * tmp3422;
    let tmp3440 = -1.0 * tmp3439;
    let tmp3441 = tmp3438 + tmp3440;
    let tmp3442 = tmp3436 + tmp3441;
    let tmp3443 = -1.0 * tmp3429;
    let tmp3444 = -1.0 * tmp3438;
    let tmp3445 = tmp3439 + tmp3444;
    let tmp3446 = -1.0 * tmp3432;
    let tmp3447 = tmp3434 + tmp3446;
    let tmp3448 = tmp3445 + tmp3447;
    let tmp3449 = -1.0 * tmp3417;
    let tmp3450 = -1.0 * tmp3419;
    let tmp3451 = -1.0 * tmp3426;
    let tmp3452 = -1.0 * tmp3427;
    let tmp3453 = tmp3451 + tmp3452;
    let tmp3454 = tmp3420 + tmp3423;
    let tmp3455 = tmp3432 + tmp3434 + tmp3440 + tmp3444;
    let tmp3456 = tmp3421 + tmp3423 + tmp3426 + tmp3452;
    let tmp3457 = tmp3420 + tmp3424 + tmp3427 + tmp3451;
    let tmp3458 = tmp3435 + tmp3438 + tmp3439 + tmp3446;
    let tmp3459 = tmp10 * tmp3416;
    let tmp3460 = tmp3418 * tmp4;
    let tmp3461 = tmp2386 + tmp3460;
    let tmp3462 = tmp227 * tmp3422;
    let tmp3463 = tmp231 * tmp3422;
    let tmp3464 = -1.0 * tmp3463;
    let tmp3465 = tmp3462 + tmp3464;
    let tmp3466 = tmp235 * tmp3252;
    let tmp3467 = tmp2354 * tmp3250;
    let tmp3468 = -1.0 * tmp3467;
    let tmp3469 = tmp3466 + tmp3468;
    let tmp3470 = tmp3465 + tmp3469;
    let tmp3471 = -1.0 * tmp3460;
    let tmp3472 = -1.0 * tmp3462;
    let tmp3473 = -1.0 * tmp3466;
    let tmp3474 = tmp3415 * tmp8;
    let tmp3475 = tmp2415 * tmp3252;
    let tmp3476 = -1.0 * tmp3475;
    let tmp3477 = tmp2418 * tmp3250;
    let tmp3478 = -1.0 * tmp3477;
    let tmp3479 = tmp3476 + tmp3478;
    let tmp3480 = tmp11 * tmp3415;
    let tmp3481 = tmp3415 * tmp6;
    let tmp3482 = tmp3480 - 1.0 * tmp3481;
    let tmp3483 = tmp3415 * tmp9;
    let tmp3484 = tmp220 * tmp3249;
    let tmp3485 = -1.0 * tmp3484;
    let tmp3486 = tmp218 * tmp3249;
    let tmp3487 = -1.0 * tmp3486;
    let tmp3488 = tmp3485 + tmp3487;
    let tmp3489 = -1.0 * tmp3483 + tmp3488;
    let tmp3490 = tmp3474 + tmp3479 + tmp3482 + tmp3489 + tmp644 + tmp743;
    let tmp3491 = tmp2390 - 1.0 * tmp3459;
    let tmp3492 = tmp3467 + tmp3473;
    let tmp3493 = tmp3463 + tmp3472;
    let tmp3494 = tmp3492 + tmp3493;
    let tmp3495 = q_old[1] * tmp2441;
    let tmp3496 = tmp2440 * tmp629 - 1.0 * tmp3474 + tmp3477;
    let tmp3497 = tmp3475 + tmp3480 + tmp3481 + tmp3489 + tmp3495 + tmp3496 + tmp645;
    let tmp3498 = tmp3476 + tmp3485;
    let tmp3499 = tmp3482 + tmp3483 + tmp3486 - 1.0 * tmp3495 + tmp3496 + tmp3498 + tmp707 + tmp710;
    let tmp3500 = tmp3240 * tmp818;
    let tmp3501 = grad_phi_node_i[0] * tmp3212;
    let tmp3502 = tmp3501 * tmp823;
    let tmp3503 = tmp3242 * tmp826;
    let tmp3504 = grad_phi_node_i[1] * tmp3195;
    let tmp3505 = tmp3504 * tmp830;
    let tmp3506 = tmp3234 * tmp834;
    let tmp3507 = grad_phi_node_i[0] * tmp3231;
    let tmp3508 = tmp3507 * tmp830;
    let tmp3509 = tmp3229 * tmp839;
    let tmp3510 = grad_phi_node_i[1] * tmp3224;
    let tmp3511 = tmp3510 * tmp823;
    let tmp3512 = tmp3210 * tmp843;
    let tmp3513 = grad_phi_node_i[0] * tmp3206;
    let tmp3514 = tmp3513 * tmp845;
    let tmp3515 = tmp3202 * tmp848;
    let tmp3516 = grad_phi_node_i[1] * tmp3215;
    let tmp3517 = tmp3516 * tmp845;
    let tmp3518 = 0.5 * tmp3261;
    let tmp3519 = 0.5 * tmp3262;
    let tmp3520 = -1.0 * tmp3519;
    let tmp3521 = tmp3518 + tmp3520 + tmp876 + tmp879 + tmp994;
    let tmp3522 = tmp3521 * tmp854;
    let tmp3523 = 0.5 * tmp3258;
    let tmp3524 = -1.0 * tmp3523;
    let tmp3525 = 0.5 * tmp3257;
    let tmp3526 = -1.0 * tmp3525;
    let tmp3527 = tmp2549 + tmp3524 + tmp3526;
    let tmp3528 = tmp3527 * tmp810;
    let tmp3529 = 0.5 * tmp3253;
    let tmp3530 = 0.5 * tmp3251;
    let tmp3531 = tmp2473 + tmp2479 + tmp2556 + tmp3529 + tmp3530;
    let tmp3532 = tmp3531 * tmp888;
    let tmp3533 = tmp3524 + tmp3525 + tmp981;
    let tmp3534 = tmp3533 * tmp778;
    let tmp3535 = tmp3531 * tmp910;
    let tmp3536 = -1.0 * tmp3529;
    let tmp3537 = tmp2476 + tmp2481 + tmp3530 + tmp3536;
    let tmp3538 = tmp3537 * tmp376;
    let tmp3539 = tmp3521 * tmp920;
    let tmp3540 = -1.0 * tmp3518;
    let tmp3541 = tmp3520 + tmp3540 + tmp882 + tmp904;
    let tmp3542 = tmp3541 * tmp359;
    let tmp3543 = tmp3304 * tmp930;
    let tmp3544 = tmp3302 * tmp937;
    let tmp3545 = tmp3299 * tmp942;
    let tmp3546 = tmp3302 * tmp949;
    let tmp3547 = tmp3500 + tmp3502 + tmp3503 + tmp3505 + tmp3506 + tmp3508 + tmp3509 + tmp3511 - 1.0 * tmp3512
        + tmp3514
        - 1.0 * tmp3515
        + tmp3517
        + tmp3522
        + tmp3528
        + tmp3532
        + tmp3534
        + tmp3535
        + tmp3538
        + tmp3539
        + tmp3542
        + tmp3543
        + tmp3544
        + tmp3545
        + tmp3546;
    let tmp3548 = tmp3241 * tmp834;
    let tmp3549 = tmp3513 * tmp958;
    let tmp3550 = grad_phi_node_i[2] * tmp3209;
    let tmp3551 = tmp3550 * tmp960;
    let tmp3552 = tmp3244 * tmp964;
    let tmp3553 = tmp3225 * tmp818;
    let tmp3554 = tmp3507 * tmp960;
    let tmp3555 = tmp3238 * tmp968;
    let tmp3556 = grad_phi_node_i[2] * tmp3228;
    let tmp3557 = tmp3556 * tmp958;
    let tmp3558 = tmp3196 * tmp843;
    let tmp3559 = tmp3501 * tmp848;
    let tmp3560 = grad_phi_node_i[2] * tmp3201;
    let tmp3561 = tmp3560 * tmp848;
    let tmp3562 = tmp3216 * tmp845;
    let tmp3563 = tmp3518 + tmp3519 + tmp995;
    let tmp3564 = tmp3563 * tmp977;
    let tmp3565 = tmp3537 * tmp748;
    let tmp3566 = -1.0 * tmp3530;
    let tmp3567 = tmp2475 + tmp2480 + tmp2555 + tmp3536 + tmp3566;
    let tmp3568 = tmp3567 * tmp812;
    let tmp3569 = tmp3523 + tmp3526 + tmp863 + tmp923;
    let tmp3570 = tmp3569 * tmp991;
    let tmp3571 = tmp3569 * tmp910;
    let tmp3572 = tmp345 * tmp3527;
    let tmp3573 = tmp3519 + tmp3540 + tmp881 + tmp903 + tmp993;
    let tmp3574 = tmp3573 * tmp387;
    let tmp3575 = tmp1005 * tmp3563;
    let tmp3576 = tmp1010 * tmp3304;
    let tmp3577 = tmp1015 * tmp3299;
    let tmp3578 = tmp1022 * tmp3299;
    let tmp3579 = tmp1027 * tmp3302;
    let tmp3580 = tmp3548 + tmp3549 + tmp3551 + tmp3552 + tmp3553 + tmp3554 + tmp3555 + tmp3557 - 1.0 * tmp3558
        + tmp3559
        + tmp3561
        - 1.0 * tmp3562
        + tmp3564
        + tmp3565
        + tmp3568
        + tmp3570
        + tmp3571
        + tmp3572
        + tmp3574
        + tmp3575
        + tmp3576
        + tmp3577
        + tmp3578
        + tmp3579;
    let tmp3581 = tmp1035 * tmp3241;
    let tmp3582 = tmp3513 * tmp839;
    let tmp3583 = tmp3550 * tmp826;
    let tmp3584 = tmp3244 * tmp830;
    let tmp3585 = tmp3507 * tmp826;
    let tmp3586 = tmp1041 * tmp3225;
    let tmp3587 = tmp3556 * tmp839;
    let tmp3588 = tmp3238 * tmp823;
    let tmp3589 = tmp3541 * tmp977;
    let tmp3590 = tmp2557 + tmp3529 + tmp3566;
    let tmp3591 = tmp3590 * tmp748;
    let tmp3592 = tmp3531 * tmp812;
    let tmp3593 = tmp3533 * tmp991;
    let tmp3594 = tmp3533 * tmp910;
    let tmp3595 = tmp2495 + tmp3523 + tmp3525;
    let tmp3596 = tmp345 * tmp3595;
    let tmp3597 = tmp3521 * tmp387;
    let tmp3598 = tmp1005 * tmp3541;
    let tmp3599 = tmp1055 * tmp3304;
    let tmp3600 = tmp1057 * tmp3299;
    let tmp3601 = tmp1059 * tmp3299;
    let tmp3602 = tmp1061 * tmp3302;
    let tmp3603 = tmp3558 - 1.0 * tmp3559 - 1.0 * tmp3561
        + tmp3562
        + tmp3581
        + tmp3582
        + tmp3583
        + tmp3584
        + tmp3585
        + tmp3586
        + tmp3587
        + tmp3588
        + tmp3589
        + tmp3591
        + tmp3592
        + tmp3593
        + tmp3594
        + tmp3596
        + tmp3597
        + tmp3598
        + tmp3599
        + tmp3600
        + tmp3601
        + tmp3602;
    let tmp3604 = tmp1041 * tmp3240;
    let tmp3605 = tmp3501 * tmp968;
    let tmp3606 = tmp3242 * tmp960;
    let tmp3607 = tmp3504 * tmp964;
    let tmp3608 = tmp3507 * tmp964;
    let tmp3609 = tmp1035 * tmp3234;
    let tmp3610 = tmp3510 * tmp968;
    let tmp3611 = tmp3229 * tmp958;
    let tmp3612 = tmp3573 * tmp854;
    let tmp3613 = tmp3595 * tmp810;
    let tmp3614 = tmp3567 * tmp888;
    let tmp3615 = tmp3569 * tmp778;
    let tmp3616 = tmp3567 * tmp910;
    let tmp3617 = tmp3590 * tmp376;
    let tmp3618 = tmp3573 * tmp920;
    let tmp3619 = tmp3563 * tmp359;
    let tmp3620 = tmp1081 * tmp3304;
    let tmp3621 = tmp1083 * tmp3302;
    let tmp3622 = tmp1085 * tmp3299;
    let tmp3623 = tmp1087 * tmp3302;
    let tmp3624 = tmp3512 - 1.0 * tmp3514 + tmp3515 - 1.0 * tmp3517
        + tmp3604
        + tmp3605
        + tmp3606
        + tmp3607
        + tmp3608
        + tmp3609
        + tmp3610
        + tmp3611
        + tmp3612
        + tmp3613
        + tmp3614
        + tmp3615
        + tmp3616
        + tmp3617
        + tmp3618
        + tmp3619
        + tmp3620
        + tmp3621
        + tmp3622
        + tmp3623;
    let tmp3625 = tmp3516 * tmp834;
    let tmp3626 = tmp3243 * tmp958;
    let tmp3627 = tmp3560 * tmp818;
    let tmp3628 = tmp3245 * tmp823;
    let tmp3629 = tmp3510 * tmp818;
    let tmp3630 = tmp3236 * tmp960;
    let tmp3631 = tmp3556 * tmp834;
    let tmp3632 = tmp3232 * tmp830;
    let tmp3633 = tmp3504 * tmp843;
    let tmp3634 = tmp3213 * tmp848;
    let tmp3635 = tmp3550 * tmp843;
    let tmp3636 = tmp3207 * tmp845;
    let tmp3637 = tmp3563 * tmp811;
    let tmp3638 = tmp1104 * tmp3537;
    let tmp3639 = tmp3521 * tmp785;
    let tmp3640 = tmp1107 * tmp3527;
    let tmp3641 = tmp3569 * tmp382;
    let tmp3642 = tmp3527 * tmp920;
    let tmp3643 = tmp3531 * tmp368;
    let tmp3644 = tmp1005 * tmp3537;
    let tmp3645 = tmp1116 * tmp3304;
    let tmp3646 = tmp1121 * tmp3299;
    let tmp3647 = tmp1126 * tmp3304;
    let tmp3648 = tmp1131 * tmp3302;
    let tmp3649 = tmp3625 + tmp3626 + tmp3627 + tmp3628 + tmp3629 + tmp3630 + tmp3631 + tmp3632 - 1.0 * tmp3633
        + tmp3634
        - 1.0 * tmp3635
        + tmp3636
        + tmp3637
        + tmp3638
        + tmp3639
        + tmp3640
        + tmp3641
        + tmp3642
        + tmp3643
        + tmp3644
        + tmp3645
        + tmp3646
        + tmp3647
        + tmp3648;
    let tmp3650 = tmp1035 * tmp3516;
    let tmp3651 = tmp3243 * tmp839;
    let tmp3652 = tmp1041 * tmp3560;
    let tmp3653 = tmp3245 * tmp968;
    let tmp3654 = tmp3236 * tmp826;
    let tmp3655 = tmp1041 * tmp3510;
    let tmp3656 = tmp3232 * tmp964;
    let tmp3657 = tmp1035 * tmp3556;
    let tmp3658 = tmp3541 * tmp811;
    let tmp3659 = tmp1104 * tmp3590;
    let tmp3660 = tmp3573 * tmp785;
    let tmp3661 = tmp1107 * tmp3595;
    let tmp3662 = tmp3533 * tmp382;
    let tmp3663 = tmp3595 * tmp920;
    let tmp3664 = tmp3567 * tmp368;
    let tmp3665 = tmp1005 * tmp3590;
    let tmp3666 = tmp1155 * tmp3304;
    let tmp3667 = tmp1157 * tmp3299;
    let tmp3668 = tmp1159 * tmp3304;
    let tmp3669 = tmp1161 * tmp3302;
    let tmp3670 = tmp3633 - 1.0 * tmp3634 + tmp3635 - 1.0 * tmp3636
        + tmp3650
        + tmp3651
        + tmp3652
        + tmp3653
        + tmp3654
        + tmp3655
        + tmp3656
        + tmp3657
        + tmp3658
        + tmp3659
        + tmp3660
        + tmp3661
        + tmp3662
        + tmp3663
        + tmp3664
        + tmp3665
        + tmp3666
        + tmp3667
        + tmp3668
        + tmp3669;
    let tmp3671 = 0.5 * tmp3419;
    let tmp3672 = tmp3249 * tmp858;
    let tmp3673 = tmp256 * tmp3672;
    let tmp3674 = tmp3252 * tmp858;
    let tmp3675 = q_old[0] * tmp3674;
    let tmp3676 = tmp3671 + tmp3673 - 1.0 * tmp3675;
    let tmp3677 = 0.5 * tmp3417;
    let tmp3678 = tmp263 * tmp3672;
    let tmp3679 = tmp3250 * tmp858;
    let tmp3680 = q_old[3] * tmp3679;
    let tmp3681 = -1.0 * tmp3677 + tmp3678 - 1.0 * tmp3680;
    let tmp3682 = tmp1350 + tmp2697 + tmp2791 + tmp3676 + tmp3681;
    let tmp3683 = tmp3677 - 1.0 * tmp3678 + tmp3680;
    let tmp3684 = tmp2699 + tmp3676 + tmp3683;
    let tmp3685 = q_old[2] * tmp3674;
    let tmp3686 = q_old[1] * tmp3679;
    let tmp3687 = -1.0 * tmp3686;
    let tmp3688 = tmp283 * tmp3672;
    let tmp3689 = tmp281 * tmp3672;
    let tmp3690 = -1.0 * tmp3689;
    let tmp3691 = 0.5 * tmp3430;
    let tmp3692 = tmp647 * tmp858;
    let tmp3693 = q_old[3] * tmp225 * tmp858;
    let tmp3694 = -1.0 * tmp3691 - 1.0 * tmp3692 - 1.0 * tmp3693;
    let tmp3695 = 0.5 * tmp3429;
    let tmp3696 = tmp256 * tmp2635;
    let tmp3697 = tmp263 * tmp2635;
    let tmp3698 = -1.0 * tmp3695 + tmp3696 + tmp3697;
    let tmp3699 = tmp1203 + tmp3685 + tmp3687 + tmp3688 + tmp3690 + tmp3694 + tmp3698;
    let tmp3700 = 0.5 * tmp3459;
    let tmp3701 = -1.0 * tmp3700;
    let tmp3702 = q_old[0] * tmp3679;
    let tmp3703 = tmp231 * tmp3672;
    let tmp3704 = 0.5 * tmp3460;
    let tmp3705 = tmp227 * tmp3672;
    let tmp3706 = q_old[3] * tmp3674;
    let tmp3707 = tmp3705 + tmp3706;
    let tmp3708 = tmp3704 + tmp3707;
    let tmp3709 = tmp2789 + tmp3701 + tmp3702 + tmp3703 + tmp3708;
    let tmp3710 = tmp186 * tmp3344;
    let tmp3711 = grad_phi_node_i[0] * tmp3710;
    let tmp3712 = tmp2723 * tmp3196;
    let tmp3713 = (1_f64 / 2.0) * tmp3304;
    let tmp3714 = tmp3713 * tmp63;
    let tmp3715 = tmp3406 * tmp476;
    let tmp3716 = tmp3364 * tmp839;
    let tmp3717 = 2.0 * tmp3590;
    let tmp3718 = 2.0 * tmp836;
    let tmp3719 = 2.0 * tmp3541;
    let tmp3720 = tmp1281 * tmp3339;
    let tmp3721 = tmp2810 * tmp3501;
    let tmp3722 = (1_f64 / 2.0) * tmp3299;
    let tmp3723 = tmp3722 * tmp824;
    let tmp3724 = -1.0 * tmp3702;
    let tmp3725 = -1.0 * tmp3703;
    let tmp3726 = tmp3724 + tmp3725;
    let tmp3727 = tmp3700 + tmp3726;
    let tmp3728 = tmp2682 + tmp2692 + tmp2788 + tmp3708 + tmp3727;
    let tmp3729 = tmp3686 + tmp3690;
    let tmp3730 = -1.0 * tmp3685;
    let tmp3731 = tmp3688 + tmp3730;
    let tmp3732 = tmp3729 + tmp3731;
    let tmp3733 = tmp3691 + tmp3692 + tmp3693;
    let tmp3734 = tmp3698 + tmp3732 + tmp3733;
    let tmp3735 = tmp1311 * tmp3348;
    let tmp3736 = tmp2706 * tmp3216;
    let tmp3737 = (1_f64 / 2.0) * tmp3302;
    let tmp3738 = tmp129 * tmp3737;
    let tmp3739 = tmp3344 * tmp830;
    let tmp3740 = tmp3368 * tmp826;
    let tmp3741 = tmp3401 * tmp839;
    let tmp3742 = 2.0 * tmp3531;
    let tmp3743 = 2.0 * tmp3533;
    let tmp3744 = 2.0 * tmp970;
    let tmp3745 = tmp1318 * tmp3329;
    let tmp3746 = tmp3722 * tmp974;
    let tmp3747 = tmp2810 * tmp3560;
    let tmp3748 = tmp1286 * tmp3497
        + tmp1287 * tmp3728
        + tmp1295 * tmp3682
        + tmp1305 * tmp3699
        + tmp1309 * tmp3490
        + tmp189 * tmp3742
        + tmp2256 * tmp3521
        + tmp2712 * tmp3244
        + tmp2714 * tmp3238
        + tmp2750 * tmp3550
        + tmp2752 * tmp3556
        + tmp3384 * tmp823
        + tmp3541 * tmp3744
        + tmp3734 * tmp744
        + tmp3735
        + tmp3736
        + tmp3738
        + tmp3739 * tmp513
        + tmp3740 * tmp513
        + tmp3741 * tmp513
        + tmp3743 * tmp961
        - 1.0 * tmp3745
        - 1.0 * tmp3746
        - 1.0 * tmp3747;
    let tmp3749 = -1.0 * tmp3704;
    let tmp3750 = -1.0 * tmp3706;
    let tmp3751 = tmp3702 + tmp3750;
    let tmp3752 = -1.0 * tmp3705;
    let tmp3753 = tmp3703 + tmp3752;
    let tmp3754 = tmp3751 + tmp3753;
    let tmp3755 = tmp2679 + tmp2688 + tmp2787 + tmp3701 + tmp3749 + tmp3754;
    let tmp3756 = -1.0 * tmp3671 - 1.0 * tmp3673 + tmp3675;
    let tmp3757 = tmp1235 + tmp2698 + tmp2792 + tmp3683 + tmp3756;
    let tmp3758 = -1.0 * tmp3688;
    let tmp3759 = tmp3685 + tmp3758;
    let tmp3760 = tmp3687 + tmp3689;
    let tmp3761 = tmp3759 + tmp3760;
    let tmp3762 = tmp3695 - 1.0 * tmp3696 - 1.0 * tmp3697;
    let tmp3763 = tmp3694 + tmp3761 + tmp3762;
    let tmp3764 = tmp1335 + tmp3686 + tmp3689 + tmp3730 + tmp3733 + tmp3758 + tmp3762;
    let tmp3765 = tmp3368 * tmp960;
    let tmp3766 = tmp3344 * tmp964;
    let tmp3767 = 2.0 * tmp3569;
    let tmp3768 = 2.0 * tmp3563;
    let tmp3769 = 2.0 * tmp3567;
    let tmp3770 = tmp2793 + tmp3681 + tmp3756;
    let tmp3771 = tmp2683 + tmp2693 + tmp3727 + tmp3749 + tmp3750 + tmp3752;
    let tmp3772 = tmp3406 * tmp960;
    let tmp3773 = tmp3364 * tmp958;
    let tmp3774 = 2.0 * tmp3537;
    let tmp3775 = tmp1188 * tmp3757
        + tmp1222 * tmp3764
        + tmp1246 * tmp3771
        + tmp1363 * tmp3490
        + tmp1374 * tmp3499
        + tmp1377 * tmp3348
        + tmp184 * tmp3774
        + tmp2246 * tmp3527
        + tmp2716 * tmp3241
        + tmp2718 * tmp3225
        + tmp2778 * tmp3507
        + tmp2781 * tmp3513
        + tmp3383 * tmp818
        + tmp3569 * tmp3718
        - 1.0 * tmp3711
        - 1.0 * tmp3712
        - 1.0 * tmp3714
        + tmp3720
        + tmp3721
        + tmp3723
        + tmp3768 * tmp846
        + tmp3770 * tmp716
        + tmp3772 * tmp476
        + tmp3773 * tmp476;
    let tmp3776 = grad_phi_node_i[1] * tmp3710;
    let tmp3777 = tmp2723 * tmp3504;
    let tmp3778 = tmp3713 * tmp831;
    let tmp3779 = tmp3382 * tmp495;
    let tmp3780 = tmp3348 * tmp495;
    let tmp3781 = tmp1410 * tmp3339;
    let tmp3782 = tmp2756 * tmp3212;
    let tmp3783 = tmp122 * tmp3722;
    let tmp3784 = tmp186 * tmp3368;
    let tmp3785 = grad_phi_node_i[2] * tmp3784;
    let tmp3786 = tmp2723 * tmp3550;
    let tmp3787 = tmp3713 * tmp961;
    let tmp3788 = tmp1041 * tmp3329;
    let tmp3789 = tmp3339 * tmp968;
    let tmp3790 = tmp3406 * tmp964;
    let tmp3791 = 2.0 * tmp3595;
    let tmp3792 = 2.0 * tmp3573;
    let tmp3793 = tmp1311 * tmp3364;
    let tmp3794 = tmp107 * tmp3737;
    let tmp3795 = tmp2706 * tmp3207;
    let tmp3796 = tmp1305 * tmp3709
        + tmp1415 * tmp3499
        + tmp1416 * tmp3763
        + tmp1417 * tmp3684
        + tmp1421 * tmp3497
        + tmp190 * tmp3792
        + tmp2250 * tmp3567
        + tmp2782 * tmp3245
        + tmp2784 * tmp3232
        + tmp2799 * tmp3560
        + tmp2801 * tmp3556
        + tmp2832 * tmp3401
        + tmp3590 * tmp3744
        + tmp3755 * tmp729
        + tmp3785
        + tmp3786
        + tmp3787
        + tmp3788 * tmp513
        + tmp3789 * tmp513
        + tmp3790 * tmp513
        + tmp3791 * tmp974
        - 1.0 * tmp3793
        - 1.0 * tmp3794
        - 1.0 * tmp3795;
    let tmp3797 = tmp3339 * tmp823;
    let tmp3798 = tmp3329 * tmp818;
    let tmp3799 = 2.0 * tmp3521;
    let tmp3800 = 2.0 * tmp3527;
    let tmp3801 = tmp1397 * tmp3770
        + tmp1398 * tmp3764
        + tmp1399 * tmp3771
        + tmp1442 * tmp3490
        + tmp1443 * tmp3499
        + tmp188 * tmp3768
        + tmp2255 * tmp3569
        + tmp2716 * tmp3516
        + tmp2718 * tmp3510
        + tmp2747 * tmp3527
        + tmp2778 * tmp3236
        + tmp2781 * tmp3243
        + tmp3757 * tmp740
        + tmp3772 * tmp495
        + tmp3773 * tmp495
        + tmp3774 * tmp850
        - 1.0 * tmp3776
        - 1.0 * tmp3777
        - 1.0 * tmp3778
        + tmp3779 * tmp818
        + tmp3780 * tmp834
        + tmp3781
        + tmp3782
        + tmp3783;
    let tmp3802 = tmp2702 * tmp3364;
    let tmp3803 = tmp3737 * tmp846;
    let tmp3804 = tmp2706 * tmp3513;
    let tmp3805 = grad_phi_node_i[0] * tmp3784;
    let tmp3806 = tmp2723 * tmp3210;
    let tmp3807 = tmp116 * tmp3713;
    let tmp3808 = tmp2739 * tmp3348;
    let tmp3809 = tmp2706 * tmp3516;
    let tmp3810 = tmp3737 * tmp850;
    let tmp3811 = tmp1410 * tmp3329;
    let tmp3812 = tmp3722 * tmp84;
    let tmp3813 = tmp2756 * tmp3201;
    let tmp3814 = tmp1397 * tmp3734
        + tmp1476 * tmp3497
        + tmp1477 * tmp3728
        + tmp1478 * tmp3682
        + tmp1482 * tmp3490
        + tmp187 * tmp3743
        + tmp2249 * tmp3541
        + tmp2712 * tmp3504
        + tmp2714 * tmp3510
        + tmp2747 * tmp3521
        + tmp2750 * tmp3242
        + tmp2752 * tmp3229
        + tmp3699 * tmp722
        + tmp3739 * tmp495
        + tmp3740 * tmp495
        + tmp3741 * tmp495
        + tmp3742 * tmp831
        + tmp3779 * tmp823
        + tmp3808
        + tmp3809
        + tmp3810
        - 1.0 * tmp3811
        - 1.0 * tmp3812
        - 1.0 * tmp3813;
    let tmp3815 = tmp1188 * tmp3755
        + tmp1274 * tmp3401
        + tmp1463 * tmp3763
        + tmp1464 * tmp3684
        + tmp1491 * tmp3499
        + tmp1492 * tmp3497
        + tmp180 * tmp3791
        + tmp2252 * tmp3590
        + tmp2782 * tmp3501
        + tmp2784 * tmp3507
        + tmp2799 * tmp3240
        + tmp2801 * tmp3234
        + tmp3567 * tmp3718
        + tmp3709 * tmp736
        + tmp3788 * tmp476
        + tmp3789 * tmp476
        + tmp3790 * tmp476
        + tmp3792 * tmp824
        - 1.0 * tmp3802
        - 1.0 * tmp3803
        - 1.0 * tmp3804
        + tmp3805
        + tmp3806
        + tmp3807;
    let tmp3816 = tmp179 * tmp3501;
    let tmp3817 = tmp186 * tmp3196;
    let tmp3818 = 2.0 * tmp3548
        + 2.0 * tmp3549
        + 2.0 * tmp3553
        + 2.0 * tmp3554
        + 2.0 * tmp3564
        + 2.0 * tmp3565
        + 2.0 * tmp3571
        + 2.0 * tmp3572
        + 2.0 * tmp3576
        + 2.0 * tmp3577
        + tmp3816
        - 1.0 * tmp3817;
    let tmp3819 = tmp183 * tmp3216;
    let tmp3820 = tmp179 * tmp3560;
    let tmp3821 = 2.0 * tmp3583
        + 2.0 * tmp3584
        + 2.0 * tmp3587
        + 2.0 * tmp3588
        + 2.0 * tmp3592
        + 2.0 * tmp3593
        + 2.0 * tmp3597
        + 2.0 * tmp3598
        + 2.0 * tmp3601
        + 2.0 * tmp3602
        + tmp3819
        - 1.0 * tmp3820;
    let tmp3822 = tmp3818 + tmp3821;
    let tmp3823 = tmp183 * tmp3516;
    let tmp3824 = tmp179 * tmp3202;
    let tmp3825 = 2.0 * tmp3503
        + 2.0 * tmp3505
        + 2.0 * tmp3509
        + 2.0 * tmp3511
        + 2.0 * tmp3532
        + 2.0 * tmp3534
        + 2.0 * tmp3539
        + 2.0 * tmp3542
        + 2.0 * tmp3545
        + 2.0 * tmp3546
        + tmp3823
        - 1.0 * tmp3824;
    let tmp3826 = tmp186 * tmp3210;
    let tmp3827 = tmp183 * tmp3513;
    let tmp3828 = 2.0 * tmp3604
        + 2.0 * tmp3605
        + 2.0 * tmp3608
        + 2.0 * tmp3609
        + 2.0 * tmp3612
        + 2.0 * tmp3613
        + 2.0 * tmp3616
        + 2.0 * tmp3617
        + 2.0 * tmp3620
        + 2.0 * tmp3621
        + tmp3826
        - 1.0 * tmp3827;
    let tmp3829 = tmp3825 + tmp3828;
    let tmp3830 = tmp179 * tmp3213;
    let tmp3831 = tmp186 * tmp3504;
    let tmp3832 = 2.0 * tmp3625
        + 2.0 * tmp3626
        + 2.0 * tmp3629
        + 2.0 * tmp3630
        + 2.0 * tmp3637
        + 2.0 * tmp3638
        + 2.0 * tmp3641
        + 2.0 * tmp3642
        + 2.0 * tmp3645
        + 2.0 * tmp3646
        + tmp3830
        - 1.0 * tmp3831;
    let tmp3833 = tmp186 * tmp3550;
    let tmp3834 = tmp183 * tmp3207;
    let tmp3835 = 2.0 * tmp3652
        + 2.0 * tmp3653
        + 2.0 * tmp3656
        + 2.0 * tmp3657
        + 2.0 * tmp3660
        + 2.0 * tmp3661
        + 2.0 * tmp3664
        + 2.0 * tmp3665
        + 2.0 * tmp3668
        + 2.0 * tmp3669
        + tmp3833
        - 1.0 * tmp3834;
    let tmp3836 = tmp3832 + tmp3835;
    let tmp3837 = 0.666666666666667 * tmp3258;
    let tmp3838 = 0.666666666666667 * tmp3257;
    let tmp3839 = -1.0 * tmp3838;
    let tmp3840 = tmp1705 + tmp1770 + tmp3837 + tmp3839;
    let tmp3841 = tmp368 * tmp3840;
    let tmp3842 = -1.0 * tmp3837;
    let tmp3843 = tmp2963 + tmp3839 + tmp3842;
    let tmp3844 = tmp3843 * tmp387;
    let tmp3845 = 0.666666666666667 * tmp3262;
    let tmp3846 = 0.666666666666667 * tmp3261;
    let tmp3847 = tmp1515 + tmp3845 + tmp3846;
    let tmp3848 = tmp278 * tmp3847;
    let tmp3849 = 0.666666666666667 * tmp3251;
    let tmp3850 = 0.666666666666667 * tmp3253;
    let tmp3851 = -1.0 * tmp3850;
    let tmp3852 = tmp2862 + tmp2873 + tmp3849 + tmp3851;
    let tmp3853 = tmp322 * tmp3852;
    let tmp3854 = tmp1554 * tmp3304;
    let tmp3855 = tmp1559 * tmp3299;
    let tmp3856 = tmp1564 * tmp3232;
    let tmp3857 = tmp1568 * tmp3207;
    let tmp3858 = tmp1573 * tmp3216;
    let tmp3859 = tmp1576 * tmp3238;
    let tmp3860 = tmp3195 * tmp545;
    let tmp3861 = (1_f64 / 3.0) * tmp3860;
    let tmp3862 = tmp1579 * tmp3245;
    let tmp3863 = 0.333333333333333 * tmp3251;
    let tmp3864 = tmp2910 + tmp3863;
    let tmp3865 = 0.333333333333333 * tmp3253;
    let tmp3866 = tmp2975 + tmp3865;
    let tmp3867 = tmp2935 + tmp3864 + tmp3866;
    let tmp3868 = tmp382 * tmp3867;
    let tmp3869 = 0.333333333333333 * tmp3261;
    let tmp3870 = 0.333333333333333 * tmp3262;
    let tmp3871 = -1.0 * tmp3870;
    let tmp3872 = tmp1614 + tmp1617 + tmp1661 + tmp3869 + tmp3871;
    let tmp3873 = tmp312 * tmp3872;
    let tmp3874 = 0.333333333333333 * tmp3257;
    let tmp3875 = -1.0 * tmp3874;
    let tmp3876 = 0.333333333333333 * tmp3258;
    let tmp3877 = -1.0 * tmp3876;
    let tmp3878 = tmp2984 + tmp3875 + tmp3877;
    let tmp3879 = tmp254 * tmp3878;
    let tmp3880 = tmp2912 - 1.0 * tmp3865;
    let tmp3881 = tmp2916 + tmp3864 + tmp3880;
    let tmp3882 = tmp359 * tmp3881;
    let tmp3883 = tmp1629 * tmp3304;
    let tmp3884 = tmp1631 * tmp3302;
    let tmp3885 = tmp1635 * tmp3236;
    let tmp3886 = tmp1640 * tmp3213;
    let tmp3887 = tmp1643 * tmp3229;
    let tmp3888 = tmp1649 * tmp3202;
    let tmp3889 = tmp3209 * tmp537;
    let tmp3890 = tmp1651 * tmp3243
        + tmp3868
        + tmp3873
        + tmp3879
        + tmp3882
        + tmp3883
        + tmp3884
        + tmp3885
        + tmp3886
        + tmp3887
        + tmp3888
        - 1_f64 / 6.0 * tmp3889;
    let tmp3891 = tmp2974 - 1.0 * tmp3863;
    let tmp3892 = tmp2934 + tmp3880 + tmp3891;
    let tmp3893 = tmp212 * tmp3892;
    let tmp3894 = tmp1607 + tmp1669 + tmp3875 + tmp3876;
    let tmp3895 = tmp300 * tmp3894;
    let tmp3896 = -1.0 * tmp3869;
    let tmp3897 = tmp1619 + tmp1660 + tmp1738 + tmp3870 + tmp3896;
    let tmp3898 = tmp345 * tmp3897;
    let tmp3899 = tmp1662 + tmp3869 + tmp3870;
    let tmp3900 = tmp376 * tmp3899;
    let tmp3901 = tmp1675 * tmp3302;
    let tmp3902 = tmp1677 * tmp3299;
    let tmp3903 = tmp1679 * tmp3210;
    let tmp3904 = tmp1681 * tmp3234;
    let tmp3905 = tmp1684 * tmp3225;
    let tmp3906 = tmp1687 * tmp3196;
    let tmp3907 = -1.0 * tmp1651 * tmp3241
        + tmp1689 * tmp3240
        + tmp3893
        + tmp3895
        + tmp3898
        + tmp3900
        + tmp3901
        + tmp3902
        + tmp3903
        + tmp3904
        + tmp3905
        + tmp3906;
    let tmp3908 = tmp3841 + tmp3844 + tmp3848 + tmp3853 + tmp3854 + tmp3855 + tmp3856 + tmp3857 + tmp3858 + tmp3859
        - 1.0 * tmp3861
        + tmp3862
        + tmp3890
        + tmp3907;
    let tmp3909 = -1.0 * tmp3849;
    let tmp3910 = tmp2861 + tmp2865 + tmp2872 + tmp2957 + tmp3851 + tmp3909;
    let tmp3911 = tmp382 * tmp3910;
    let tmp3912 = -1.0 * tmp3846;
    let tmp3913 = tmp1514 + tmp1518 + tmp1710 + tmp3845 + tmp3912;
    let tmp3914 = tmp312 * tmp3913;
    let tmp3915 = tmp2885 + tmp3837 + tmp3838;
    let tmp3916 = tmp254 * tmp3915;
    let tmp3917 = tmp2867 + tmp2959 + tmp3850 + tmp3909;
    let tmp3918 = tmp359 * tmp3917;
    let tmp3919 = tmp1721 * tmp3302;
    let tmp3920 = tmp1726 * tmp3304;
    let tmp3921 = tmp1728 * tmp3202;
    let tmp3922 = tmp1731 * tmp3229;
    let tmp3923 = tmp1734 * tmp3213;
    let tmp3924 = tmp1736 * tmp3236;
    let tmp3925 = tmp1692 * tmp3243;
    let tmp3926 = (1_f64 / 3.0) * tmp3889;
    let tmp3927 = tmp1601 + tmp1666 + tmp3874 + tmp3877;
    let tmp3928 = tmp368 * tmp3927;
    let tmp3929 = tmp2920 + tmp3874 + tmp3876;
    let tmp3930 = tmp387 * tmp3929;
    let tmp3931 = tmp1620 + tmp1739 + tmp3871 + tmp3896;
    let tmp3932 = tmp278 * tmp3931;
    let tmp3933 = tmp2979 + tmp3866 + tmp3891;
    let tmp3934 = tmp322 * tmp3933;
    let tmp3935 = tmp1749 * tmp3299;
    let tmp3936 = tmp1751 * tmp3304;
    let tmp3937 = tmp1753 * tmp3238;
    let tmp3938 = tmp1755 * tmp3216;
    let tmp3939 = tmp1757 * tmp3232;
    let tmp3940 = tmp1759 * tmp3207;
    let tmp3941 = -1.0 * tmp1689 * tmp3245
        + (1_f64 / 6.0) * tmp3860
        + tmp3928
        + tmp3930
        + tmp3932
        + tmp3934
        + tmp3935
        + tmp3936
        + tmp3937
        + tmp3938
        + tmp3939
        + tmp3940;
    let tmp3942 =
        tmp3907 + tmp3911 + tmp3914 + tmp3916 + tmp3918 + tmp3919 + tmp3920 + tmp3921 + tmp3922 + tmp3923 + tmp3924
            - 1.0 * tmp3925
            + tmp3926
            + tmp3941;
    let tmp3943 = tmp2860 + tmp2866 + tmp2871 + tmp2958 + tmp3849 + tmp3850;
    let tmp3944 = tmp212 * tmp3943;
    let tmp3945 = tmp1533 + tmp3838 + tmp3842;
    let tmp3946 = tmp300 * tmp3945;
    let tmp3947 = -1.0 * tmp3845;
    let tmp3948 = tmp1508 + tmp1510 + tmp1512 + tmp3846 + tmp3947;
    let tmp3949 = tmp345 * tmp3948;
    let tmp3950 = tmp1519 + tmp1711 + tmp3912 + tmp3947;
    let tmp3951 = tmp376 * tmp3950;
    let tmp3952 = tmp1780 * tmp3299;
    let tmp3953 = tmp1785 * tmp3302;
    let tmp3954 = tmp1787 * tmp3196;
    let tmp3955 = tmp1790 * tmp3225;
    let tmp3956 = tmp1792 * tmp3210;
    let tmp3957 = tmp1794 * tmp3234;
    let tmp3958 = tmp1579 * tmp3240;
    let tmp3959 = tmp1692 * tmp3241;
    let tmp3960 = tmp3890
        + tmp3941
        + tmp3944
        + tmp3946
        + tmp3949
        + tmp3951
        + tmp3952
        + tmp3953
        + tmp3954
        + tmp3955
        + tmp3956
        + tmp3957
        - 1.0 * tmp3958
        + tmp3959;
    let tmp3961 = 0.666666666666667 * tmp3459;
    let tmp3962 = 0.666666666666667 * tmp3460;
    let tmp3963 = tmp1504 * tmp3249;
    let tmp3964 = tmp227 * tmp3963;
    let tmp3965 = tmp3962 + tmp3964;
    let tmp3966 = tmp1885 * tmp3250;
    let tmp3967 = -1.0 * tmp3966;
    let tmp3968 = tmp231 * tmp3963;
    let tmp3969 = -1.0 * tmp3968;
    let tmp3970 = tmp1506 * tmp3252;
    let tmp3971 = tmp3181 + tmp3970;
    let tmp3972 = tmp3967 + tmp3969 + tmp3971;
    let tmp3973 = 0.666666666666667 * tmp3417;
    let tmp3974 = -1.0 * tmp3973;
    let tmp3975 = 0.666666666666667 * tmp3419;
    let tmp3976 = tmp1856 + tmp1988;
    let tmp3977 = tmp263 * tmp3963;
    let tmp3978 = tmp1506 * tmp3250;
    let tmp3979 = -1.0 * tmp3978;
    let tmp3980 = tmp256 * tmp3963;
    let tmp3981 = tmp1885 * tmp3252;
    let tmp3982 = -1.0 * tmp3981;
    let tmp3983 = tmp3977 + tmp3979 + tmp3980 + tmp3982;
    let tmp3984 = 0.666666666666667 * tmp3430;
    let tmp3985 = tmp1504 * tmp3433;
    let tmp3986 = tmp1504 * tmp3437;
    let tmp3987 = -1.0 * tmp3986;
    let tmp3988 = tmp3985 + tmp3987;
    let tmp3989 = tmp3984 + tmp3988;
    let tmp3990 = 0.666666666666667 * tmp3429;
    let tmp3991 = tmp283 * tmp3963;
    let tmp3992 = tmp281 * tmp3963;
    let tmp3993 = -1.0 * tmp3992;
    let tmp3994 = tmp3991 + tmp3993;
    let tmp3995 = -1.0 * tmp3990 + tmp3994;
    let tmp3996 = -1.0 * tmp3985;
    let tmp3997 = tmp3986 + tmp3996;
    let tmp3998 = -1.0 * tmp3984 + tmp3997;
    let tmp3999 = 0.333333333333333 * tmp3460;
    let tmp4000 = 0.333333333333333 * tmp3459;
    let tmp4001 = tmp3113 + tmp4000;
    let tmp4002 = tmp1584 * tmp3252;
    let tmp4003 = tmp3077 * tmp3250;
    let tmp4004 = -1.0 * tmp4003;
    let tmp4005 = tmp4002 + tmp4004;
    let tmp4006 = tmp1583 * tmp3249;
    let tmp4007 = tmp227 * tmp4006;
    let tmp4008 = tmp231 * tmp4006;
    let tmp4009 = -1.0 * tmp4008;
    let tmp4010 = tmp4007 + tmp4009;
    let tmp4011 = tmp4005 + tmp4010;
    let tmp4012 = tmp283 * tmp4006;
    let tmp4013 = tmp1583 * tmp3437;
    let tmp4014 = -1.0 * tmp4013;
    let tmp4015 = tmp4012 + tmp4014;
    let tmp4016 = tmp1583 * tmp3433;
    let tmp4017 = tmp281 * tmp4006;
    let tmp4018 = -1.0 * tmp4017;
    let tmp4019 = tmp4016 + tmp4018;
    let tmp4020 = tmp4015 + tmp4019;
    let tmp4021 = 0.333333333333333 * tmp3430;
    let tmp4022 = tmp1583 * tmp647;
    let tmp4023 = tmp1584 * tmp225;
    let tmp4024 = tmp4021 + tmp4022 + tmp4023;
    let tmp4025 = 0.333333333333333 * tmp3429;
    let tmp4026 = tmp256 * tmp3069;
    let tmp4027 = tmp263 * tmp3069;
    let tmp4028 = -1.0 * tmp4025 + tmp4026 + tmp4027;
    let tmp4029 = 0.333333333333333 * tmp3419;
    let tmp4030 = tmp3077 * tmp3252;
    let tmp4031 = tmp256 * tmp4006;
    let tmp4032 = -1.0 * tmp4029 + tmp4030 - 1.0 * tmp4031;
    let tmp4033 = 0.333333333333333 * tmp3417;
    let tmp4034 = tmp263 * tmp4006;
    let tmp4035 = tmp1584 * tmp3250;
    let tmp4036 = -1.0 * tmp4033 + tmp4034 - 1.0 * tmp4035;
    let tmp4037 = -1.0 * tmp4007;
    let tmp4038 = -1.0 * tmp4002;
    let tmp4039 = tmp3110 - 1.0 * tmp3999;
    let tmp4040 = tmp1595 * tmp3354
        + tmp1609 * tmp3314
        + tmp1622 * tmp3311
        + tmp1627 * tmp3351
        + tmp1635 * tmp3409
        + tmp1640 * tmp3340
        + tmp1643 * tmp3402
        + tmp1649 * tmp3330
        + tmp1944 * tmp3499
        + tmp1945 * tmp3497
        + tmp1947 * tmp3364
        + tmp2186 * tmp3878
        + tmp2192 * tmp3872
        + tmp2249 * tmp3881
        + tmp2255 * tmp3867
        - 1_f64 / 3.0 * tmp3389
        - 1_f64 / 3.0 * tmp3390
        + (1_f64 / 3.0) * tmp3395
        + (1_f64 / 3.0) * tmp3396
        - 1_f64 / 3.0 * tmp3410
        + tmp646 * (tmp3151 + tmp4032 + tmp4036)
        + tmp700 * (tmp4020 + tmp4024 + tmp4028)
        + tmp722 * (tmp3125 + tmp4001 + tmp4004 + tmp4009 + tmp4037 + tmp4038 + tmp4039)
        + tmp740 * (tmp3134 + tmp3171 + tmp3999 + tmp4001 + tmp4011);
    let tmp4041 = tmp4029 - 1.0 * tmp4030 + tmp4031;
    let tmp4042 = tmp4033 - 1.0 * tmp4034 + tmp4035;
    let tmp4043 = -1.0 * tmp4016;
    let tmp4044 = -1.0 * tmp4021 - 1.0 * tmp4022 - 1.0 * tmp4023;
    let tmp4045 = -1.0 * tmp4000;
    let tmp4046 = tmp1741 * tmp3352
        + tmp1743 * tmp3355
        + tmp1745 * tmp3312
        + tmp1747 * tmp3315
        + tmp1753 * tmp3384
        + tmp1755 * tmp3349
        + tmp1757 * tmp3407
        + tmp1759 * tmp3365
        + tmp1973 * tmp3490
        + tmp1974 * tmp3499
        - 1.0 * tmp1975 * tmp3339
        + tmp2187 * tmp3931
        + tmp2194 * tmp3933
        + tmp2250 * tmp3927
        + tmp2256 * tmp3929
        + (1_f64 / 3.0) * tmp3386
        - 1_f64 / 3.0 * tmp3391
        - 1_f64 / 3.0 * tmp3392
        + (1_f64 / 3.0) * tmp3397
        + (1_f64 / 3.0) * tmp3398
        + tmp671 * (tmp1950 + tmp4012 + tmp4013 + tmp4018 + tmp4028 + tmp4043 + tmp4044)
        + tmp706 * (tmp3135 + tmp3142 + tmp3999 + tmp4002 + tmp4003 + tmp4007 + tmp4008 + tmp4045)
        + tmp729 * (tmp2017 + tmp3128 + tmp3150 + tmp4036 + tmp4041)
        + tmp744 * (tmp3129 + tmp4041 + tmp4042);
    let tmp4047 = -1.0 * tmp3962;
    let tmp4048 = -1.0 * tmp3961 + tmp3966;
    let tmp4049 = -1.0 * tmp3964;
    let tmp4050 = -1.0 * tmp3970;
    let tmp4051 = tmp3041 + tmp3968;
    let tmp4052 = tmp4049 + tmp4050 + tmp4051;
    let tmp4053 = -1.0 * tmp3991;
    let tmp4054 = tmp3992 + tmp4053;
    let tmp4055 = tmp3990 + tmp4054;
    let tmp4056 = -1.0 * tmp3977;
    let tmp4057 = tmp3982 + tmp4056;
    let tmp4058 = tmp3978 + tmp3980;
    let tmp4059 = tmp4008 + tmp4037;
    let tmp4060 = tmp4003 + tmp4038;
    let tmp4061 = tmp4059 + tmp4060;
    let tmp4062 = tmp4017 + tmp4043;
    let tmp4063 = -1.0 * tmp4012;
    let tmp4064 = tmp4013 + tmp4063;
    let tmp4065 = tmp4062 + tmp4064;
    let tmp4066 = tmp4025 - 1.0 * tmp4026 - 1.0 * tmp4027;
    let tmp4067 = tmp1657 * tmp3310
        + tmp1663 * tmp3313
        + tmp1671 * tmp3350
        + tmp1673 * tmp3353
        + tmp1679 * tmp3369
        + tmp1681 * tmp3408
        + tmp1684 * tmp3383
        + tmp1687 * tmp3345
        + tmp2020 * tmp3497
        + tmp2021 * tmp3490
        - 1.0 * tmp2022 * tmp3348
        + tmp2023 * tmp3329
        + tmp2184 * tmp3892
        + tmp2189 * tmp3894
        + tmp2246 * tmp3897
        + tmp2252 * tmp3899
        - 1_f64 / 3.0 * tmp3387
        - 1_f64 / 3.0 * tmp3388
        + (1_f64 / 3.0) * tmp3393
        + (1_f64 / 3.0) * tmp3394
        + tmp601 * (tmp3132 + tmp3170 + tmp4039 + tmp4045 + tmp4061)
        + tmp693 * (tmp1963 + tmp3126 + tmp3149 + tmp4032 + tmp4042)
        + tmp716 * (tmp4044 + tmp4065 + tmp4066)
        + tmp736 * (tmp2004 + tmp4014 + tmp4016 + tmp4017 + tmp4024 + tmp4063 + tmp4066);
    let tmp4068 = -1.0 * tmp3975;
    let tmp4069 = -1.0 * tmp3980;
    let tmp4070 = tmp3978 + tmp3981 + tmp4056 + tmp4069;
    let tmp4071 = tmp1855 + tmp1859;
    let tmp4072 = tmp3979 + tmp4069;
    let tmp4073 = tmp3977 + tmp3981;
    let tmp4074 = tmp3967 + tmp4049;
    let tmp4075 = tmp3969 + tmp4050;
    let tmp4076 = tmp2071 * tmp535;
    let tmp4077 = tmp2067 * tmp539;
    let tmp4078 = tmp2055 * tmp541;
    let tmp4079 = tmp2060 * tmp543;
    let tmp4080 = tmp436
        * (2.0 * tmp2051
            + 2.0 * tmp2057
            + 2.0 * tmp2062
            + 2.0 * tmp2066
            + 2.0 * tmp2069
            + 2.0 * tmp2073
            + 2.0 * tmp2083
            + 2.0 * tmp2087
            + 2.0 * tmp2090
            + 2.0 * tmp2092
            + 2.0 * tmp2094
            + 2.0 * tmp2096
            + 2.0 * tmp2112
            + 2.0 * tmp2125
            + 2.0 * tmp2129
            + 2.0 * tmp2136
            + 2.0 * tmp2139
            + 2.0 * tmp2143
            + 2.0 * tmp2146
            + 2.0 * tmp2150
            + 2.0 * tmp2154
            + 2.0 * tmp2156
            + 2.0 * tmp2158
            + 2.0 * tmp2161
            + 2.0 * tmp2173
            + 2.0 * tmp2176
            + 2.0 * tmp2178
            + 2.0 * tmp2179
            + 2.0 * tmp2180
            + 2.0 * tmp2181
            - 1.0 * tmp2894
            + tmp2929
            + tmp4076
            + tmp4077
            - 1.0 * tmp4078
            - 1.0 * tmp4079);
    let tmp4081 = tmp2044 * tmp446;
    let tmp4082 = tmp2044 * tmp443;
    let tmp4083 = tmp2043 * tmp447;
    let tmp4084 = tmp2199 * tmp4083;
    let tmp4085 = 2.0 * tmp4084;
    let tmp4086 = tmp21 * tmp4085;
    let tmp4087 = tmp25 * tmp4085;
    let tmp4088 = tmp2044 * tmp444 + tmp2044 * tmp445 + tmp23 * tmp4087;
    let tmp4089 = tmp467 + tmp470 + tmp519;
    let tmp4090 = tmp4081 + tmp4082 - 1.0 * tmp4086 + tmp4088 + tmp4089;
    let tmp4091 = tmp4090 * tmp476;
    let tmp4092 = tmp2044 * tmp502;
    let tmp4093 = tmp2044 * tmp499;
    let tmp4094 = tmp20 * tmp4087;
    let tmp4095 = tmp23 * tmp4085;
    let tmp4096 = tmp19 * tmp4095 - 1.0 * tmp2044 * tmp497 + tmp2044 * tmp503;
    let tmp4097 = tmp4092 - 1.0 * tmp4093 - 1.0 * tmp4094 + tmp4096 + tmp485 + tmp488 + tmp525;
    let tmp4098 = tmp4097 * tmp495;
    let tmp4099 = tmp2044 * tmp438;
    let tmp4100 = tmp4099 * tmp70;
    let tmp4101 = tmp4099 * tmp79;
    let tmp4102 = tmp20 * tmp4095;
    let tmp4103 = tmp19 * tmp4087 + tmp4099 * tmp73 - 1.0 * tmp4099 * tmp77;
    let tmp4104 = tmp4100 - 1.0 * tmp4101 - 1.0 * tmp4102 + tmp4103 + tmp511;
    let tmp4105 = tmp4104 * tmp513;
    let tmp4106 = tmp3403 + tmp507 + tmp515;
    let tmp4107 = -1.0 * tmp4100 + tmp4101 + tmp4102 + tmp4103 + tmp4106;
    let tmp4108 = tmp4107 * tmp476;
    let tmp4109 = tmp2323 - 1.0 * tmp4081 - 1.0 * tmp4082 + tmp4086 + tmp4088 + tmp520;
    let tmp4110 = tmp4109 * tmp495;
    let tmp4111 = -1.0 * tmp4092 + tmp4093 + tmp4094 + tmp4096 + tmp526;
    let tmp4112 = tmp4111 * tmp513;
    let tmp4113 = 4.0 * tmp2050;
    let tmp4114 = 4.0 * tmp2056;
    let tmp4115 = 4.0 * tmp2061;
    let tmp4116 = 4.0 * tmp2065;
    let tmp4117 = 4.0 * tmp2068;
    let tmp4118 = 4.0 * tmp2072;
    let tmp4119 = tmp4111 * tmp535;
    let tmp4120 = tmp4107 * tmp537;
    let tmp4121 = tmp4109 * tmp539;
    let tmp4122 = tmp4097 * tmp541;
    let tmp4123 = tmp4104 * tmp543;
    let tmp4124 = tmp4090 * tmp545;
    let tmp4125 = 4.0 * tmp2082;
    let tmp4126 = 4.0 * tmp2086;
    let tmp4127 = 4.0 * tmp2089;
    let tmp4128 = 4.0 * tmp2091;
    let tmp4129 = 4.0 * tmp2093;
    let tmp4130 = 4.0 * tmp2095;
    let tmp4131 = 2.0 * tmp2175;
    let tmp4132 = tmp2097 * tmp2172;
    let tmp4133 = tmp2099 * tmp2177;
    let tmp4134 = tmp2102 * tmp2172;
    let tmp4135 = tmp2101 * tmp2177;
    let tmp4136 = tmp27 * tmp4084;
    let tmp4137 = -1.0 * tmp2044 * tmp441;
    let tmp4138 = tmp24 * tmp4084;
    let tmp4139 = tmp2044 * tmp452;
    let tmp4140 = tmp26 * tmp4084;
    let tmp4141 = tmp2044 * tmp454;
    let tmp4142 = -1.0 * tmp4140 - 1.0 * tmp4141;
    let tmp4143 = tmp22 * tmp4084;
    let tmp4144 = tmp2044 * tmp439;
    let tmp4145 = -1.0 * tmp4143 + tmp4144 + tmp573;
    let tmp4146 = tmp2268 + tmp2273 + tmp4136 + tmp4137 + tmp4138 + tmp4139 + tmp4142 + tmp4145 + tmp569 + tmp581;
    let tmp4147 = tmp4146 * tmp476;
    let tmp4148 = tmp4136 + tmp4137 - 1.0 * tmp4138 - 1.0 * tmp4139;
    let tmp4149 = tmp4142 + tmp4143 - 1.0 * tmp4144 + tmp4148 + tmp570 + tmp586 + tmp590;
    let tmp4150 = tmp4149 * tmp495;
    let tmp4151 = tmp2268 + tmp580;
    let tmp4152 = tmp4151 + tmp566 + tmp589 + tmp595;
    let tmp4153 = tmp4140 + tmp4141 + tmp4145 + tmp4148 + tmp4152;
    let tmp4154 = tmp4153 * tmp513;
    let tmp4155 = tmp4149 * tmp476;
    let tmp4156 = tmp4153 * tmp495;
    let tmp4157 = tmp4146 * tmp513;
    let tmp4158 = tmp627 + tmp741;
    let tmp4159 = tmp2105 * tmp605;
    let tmp4160 = tmp2333 * tmp4159;
    let tmp4161 = tmp4160 * tmp5;
    let tmp4162 = tmp4 * tmp4161;
    let tmp4163 = tmp2107 * tmp2418;
    let tmp4164 = tmp2106 * tmp214;
    let tmp4165 = tmp283 * tmp4164;
    let tmp4166 = tmp4162 + tmp4163 - 1.0 * tmp4165;
    let tmp4167 = tmp4160 * tmp7;
    let tmp4168 = tmp10 * tmp4167;
    let tmp4169 = tmp281 * tmp4164;
    let tmp4170 = tmp2109 * tmp2415;
    let tmp4171 = tmp4168 + tmp4169 - 1.0 * tmp4170;
    let tmp4172 = tmp10 * tmp4161;
    let tmp4173 = tmp2443 + tmp4172;
    let tmp4174 = tmp4 * tmp4167;
    let tmp4175 = tmp2419 + tmp2426 + tmp4174;
    let tmp4176 = q_old[0] * tmp2109;
    let tmp4177 = tmp214 * tmp4176;
    let tmp4178 = -1.0 * tmp4177;
    let tmp4179 = tmp4161 * tmp7;
    let tmp4180 = tmp2107 * tmp634;
    let tmp4181 = tmp4179 + tmp4180;
    let tmp4182 = tmp4160 * tmp87;
    let tmp4183 = tmp2106 * tmp641;
    let tmp4184 = tmp2106 * tmp639;
    let tmp4185 = -1.0 * tmp4184;
    let tmp4186 = tmp4182 + tmp4183 + tmp4185;
    let tmp4187 = tmp695 + tmp730;
    let tmp4188 = tmp2446 - 1.0 * tmp4174;
    let tmp4189 = -1.0 * tmp4180;
    let tmp4190 = tmp4177 - 1.0 * tmp4179 + tmp4189;
    let tmp4191 = -1.0 * tmp4162 - 1.0 * tmp4163 + tmp4165;
    let tmp4192 = -1.0 * tmp4183;
    let tmp4193 = tmp4178 + tmp4192;
    let tmp4194 = -1.0 * tmp4182 + tmp4184;
    let tmp4195 = -1.0 * tmp4168 - 1.0 * tmp4169 + tmp4170;
    let tmp4196 = tmp696 + tmp732;
    let tmp4197 = tmp2416 + tmp2424 - 1.0 * tmp4172;
    let tmp4198 = tmp633 + tmp742;
    let tmp4199 = tmp4160 * tmp6;
    let tmp4200 = q_old[0] * tmp602;
    let tmp4201 = tmp2107 * tmp4200;
    let tmp4202 = tmp11 * tmp4160;
    let tmp4203 = tmp4160 * tmp8;
    let tmp4204 = -1.0 * tmp2106 * tmp2332;
    let tmp4205 = q_old[3] * tmp602;
    let tmp4206 = tmp2109 * tmp4205;
    let tmp4207 = tmp4202 - 1.0 * tmp4203 + tmp4204 - 1.0 * tmp4206;
    let tmp4208 = tmp4160 * tmp9;
    let tmp4209 = tmp2106 * tmp2331;
    let tmp4210 = -1.0 * tmp4208 - 1.0 * tmp4209;
    let tmp4211 = tmp4199 - 1.0 * tmp4201 + tmp4207 + tmp4210 + tmp760 + tmp776;
    let tmp4212 = -1.0 * tmp757;
    let tmp4213 = tmp4212 + tmp758;
    let tmp4214 = tmp4213 + tmp768 + tmp780;
    let tmp4215 = -1.0 * tmp4199 + tmp4201;
    let tmp4216 = tmp4207 + tmp4208 + tmp4209 + tmp4214 + tmp4215;
    let tmp4217 =
        tmp4202 + tmp4203 + tmp4204 + tmp4206 + tmp4210 + tmp4212 + tmp4215 + tmp759 + tmp765 + tmp766 + tmp786;
    let tmp4218 = 2.0 * tmp2448
        + 2.0 * tmp2450
        + 2.0 * tmp2454
        + 2.0 * tmp2456
        + 2.0 * tmp2470
        + 2.0 * tmp2484
        + 2.0 * tmp2493
        + 2.0 * tmp2497
        + 2.0 * tmp2502
        + 2.0 * tmp2503
        + tmp2847
        - 1.0 * tmp2848
        + tmp2849;
    let tmp4219 = 2.0 * tmp2510
        + 2.0 * tmp2511
        + 2.0 * tmp2514
        + 2.0 * tmp2516
        + 2.0 * tmp2527
        + 2.0 * tmp2529
        + 2.0 * tmp2533
        + 2.0 * tmp2534
        + 2.0 * tmp2537
        + 2.0 * tmp2538
        + tmp2840
        - 1.0 * tmp2841
        + tmp2842;
    let tmp4220 = 2.0 * tmp2540
        + 2.0 * tmp2541
        + 2.0 * tmp2544
        + 2.0 * tmp2545
        + 2.0 * tmp2548
        + 2.0 * tmp2551
        + 2.0 * tmp2554
        + 2.0 * tmp2559
        + 2.0 * tmp2562
        + 2.0 * tmp2563
        - 1.0 * tmp2838
        + tmp2839
        + tmp2843;
    let tmp4221 = 2.0 * tmp2569
        + 2.0 * tmp2570
        + 2.0 * tmp2573
        + 2.0 * tmp2574
        + 2.0 * tmp2577
        + 2.0 * tmp2578
        + 2.0 * tmp2581
        + 2.0 * tmp2582
        + 2.0 * tmp2585
        + 2.0 * tmp2586
        - 1.0 * tmp2845
        + tmp2846
        + tmp2850;
    let tmp4222 = 2.0 * tmp2590
        + 2.0 * tmp2591
        + 2.0 * tmp2594
        + 2.0 * tmp2595
        + 2.0 * tmp2602
        + 2.0 * tmp2603
        + 2.0 * tmp2606
        + 2.0 * tmp2607
        + 2.0 * tmp2610
        + 2.0 * tmp2611
        + tmp2854
        - 1.0 * tmp2855
        + tmp2856;
    let tmp4223 = 2.0 * tmp2613
        + 2.0 * tmp2614
        + 2.0 * tmp2617
        + 2.0 * tmp2618
        + 2.0 * tmp2621
        + 2.0 * tmp2622
        + 2.0 * tmp2625
        + 2.0 * tmp2626
        + 2.0 * tmp2629
        + 2.0 * tmp2630
        - 1.0 * tmp2852
        + tmp2853
        + tmp2857;
    let tmp4224 = 0.5 * tmp4172;
    let tmp4225 = tmp219 * tmp2345;
    let tmp4226 = -1.0 * tmp4225;
    let tmp4227 = q_old[1] * tmp2107;
    let tmp4228 = tmp228 * tmp4227;
    let tmp4229 = -1.0 * tmp4228;
    let tmp4230 = tmp1195 + tmp1335;
    let tmp4231 = 0.5 * tmp4174;
    let tmp4232 = q_old[2] * tmp2109;
    let tmp4233 = tmp228 * tmp4232;
    let tmp4234 = -1.0 * tmp4233;
    let tmp4235 = tmp217 * tmp2345;
    let tmp4236 = -1.0 * tmp4235;
    let tmp4237 = tmp4234 + tmp4236;
    let tmp4238 = -1.0 * tmp4231 + tmp4237;
    let tmp4239 = tmp4224 + tmp4226 + tmp4229 + tmp4230 + tmp4238;
    let tmp4240 = tmp4229 + tmp4235;
    let tmp4241 = tmp4226 + tmp4233;
    let tmp4242 = tmp4240 + tmp4241;
    let tmp4243 = tmp1217 + tmp1366 + tmp4224 + tmp4231 + tmp4242;
    let tmp4244 = 0.5 * tmp4179;
    let tmp4245 = -1.0 * tmp4244;
    let tmp4246 = 0.5 * tmp4182;
    let tmp4247 = -1.0 * tmp4246;
    let tmp4248 = tmp1236 + tmp2355 + tmp2358 + tmp2362 + tmp2365 + tmp4245 + tmp4247;
    let tmp4249 = 0.5 * tmp4168;
    let tmp4250 = -1.0 * tmp4249;
    let tmp4251 = 0.5 * tmp4162;
    let tmp4252 = tmp2338 + tmp2340 + tmp2348 + tmp2373;
    let tmp4253 = tmp1259 + tmp4250 + tmp4251 + tmp4252;
    let tmp4254 = tmp2050 * tmp2177;
    let tmp4255 = tmp186 * tmp4090;
    let tmp4256 = grad_phi_node_i[0] * tmp4255;
    let tmp4257 = tmp2172 * tmp2449;
    let tmp4258 = tmp4153 * tmp476;
    let tmp4259 = tmp4104 * tmp839;
    let tmp4260 = 4.0 * tmp2558;
    let tmp4261 = 4.0 * tmp2550;
    let tmp4262 = 4.0 * tmp2491;
    let tmp4263 = 4.0 * tmp2500;
    let tmp4264 = tmp1281 * tmp4109;
    let tmp4265 = tmp1289 + tmp1327;
    let tmp4266 = tmp2349 + tmp2377 + tmp4249 + tmp4251 + tmp4265;
    let tmp4267 = tmp1298 + tmp1343 + tmp2359 + tmp2369 + tmp4244 + tmp4247;
    let tmp4268 = tmp2072 * tmp2175;
    let tmp4269 = tmp1311 * tmp4111;
    let tmp4270 = tmp2172 * tmp2519;
    let tmp4271 = tmp4090 * tmp830;
    let tmp4272 = tmp4107 * tmp826;
    let tmp4273 = tmp4149 * tmp839;
    let tmp4274 = 4.0 * tmp2487;
    let tmp4275 = tmp1318 * tmp4097;
    let tmp4276 = tmp1286 * tmp4211
        + tmp1287 * tmp4266
        + tmp1295 * tmp4239
        + tmp1305 * tmp4248
        + tmp1309 * tmp4217
        + tmp2101 * tmp4274
        + tmp2469 * tmp4130
        + tmp2509 * tmp4262
        + tmp2515 * tmp4263
        + tmp4157 * tmp823
        + tmp4267 * tmp744
        + tmp4268
        + tmp4269
        - 1.0 * tmp4270
        + tmp4271 * tmp513
        + tmp4272 * tmp513
        + tmp4273 * tmp513
        - 1.0 * tmp4275;
    let tmp4277 = -1.0 * tmp4251;
    let tmp4278 = tmp1288 + tmp1325;
    let tmp4279 = tmp2342 + tmp2374 + tmp4250 + tmp4277 + tmp4278;
    let tmp4280 = tmp1203 + tmp1337;
    let tmp4281 = tmp4225 + tmp4228;
    let tmp4282 = -1.0 * tmp4224 + tmp4281;
    let tmp4283 = tmp4231 + tmp4233 + tmp4235 + tmp4280 + tmp4282;
    let tmp4284 = tmp1297 + tmp1341 + tmp2366 + tmp2371 + tmp4245 + tmp4246;
    let tmp4285 = tmp1352 + tmp2356 + tmp2357 + tmp2363 + tmp2364 + tmp4244 + tmp4246;
    let tmp4286 = tmp4107 * tmp960;
    let tmp4287 = tmp4090 * tmp964;
    let tmp4288 = 4.0 * tmp2528;
    let tmp4289 = 4.0 * tmp2522;
    let tmp4290 = 4.0 * tmp2532;
    let tmp4291 = 4.0 * tmp2526;
    let tmp4292 = tmp1219 + tmp1367 + tmp4238 + tmp4282;
    let tmp4293 = tmp2341 + tmp2346 + tmp2347 + tmp2376;
    let tmp4294 = tmp1372 + tmp4249 + tmp4277 + tmp4293;
    let tmp4295 = tmp4153 * tmp960;
    let tmp4296 = tmp4104 * tmp958;
    let tmp4297 = 4.0 * tmp2496;
    let tmp4298 = tmp1188 * tmp4283
        + tmp1222 * tmp4285
        + tmp1246 * tmp4294
        + tmp1363 * tmp4217
        + tmp1374 * tmp4216
        + tmp1377 * tmp4111
        + tmp2098 * tmp4297
        + tmp2455 * tmp4288
        + tmp2461 * tmp4289
        + tmp2483 * tmp4125
        + tmp4147 * tmp818
        - 1.0 * tmp4254
        - 1.0 * tmp4256
        + tmp4257
        + tmp4264
        + tmp4292 * tmp716
        + tmp4295 * tmp476
        + tmp4296 * tmp476;
    let tmp4299 = tmp2177 * tmp2452;
    let tmp4300 = grad_phi_node_i[1] * tmp4255;
    let tmp4301 = tmp2068 * tmp2172;
    let tmp4302 = tmp4146 * tmp495;
    let tmp4303 = tmp4111 * tmp495;
    let tmp4304 = tmp1410 * tmp4109;
    let tmp4305 = tmp2177 * tmp2509;
    let tmp4306 = tmp186 * tmp4107;
    let tmp4307 = grad_phi_node_i[2] * tmp4306;
    let tmp4308 = tmp2061 * tmp2175;
    let tmp4309 = tmp1041 * tmp4097;
    let tmp4310 = tmp4109 * tmp968;
    let tmp4311 = tmp4153 * tmp964;
    let tmp4312 = tmp1311 * tmp4104;
    let tmp4313 = tmp1305 * tmp4253
        + tmp1415 * tmp4216
        + tmp1416 * tmp4284
        + tmp1417 * tmp4243
        + tmp1421 * tmp4211
        + tmp2089 * tmp4291
        + tmp2102 * tmp4290
        + tmp2515 * tmp4261
        + tmp2519 * tmp4260
        + tmp2832 * tmp4149
        + tmp4279 * tmp729
        + tmp4305
        + tmp4307
        - 1.0 * tmp4308
        + tmp4309 * tmp513
        + tmp4310 * tmp513
        + tmp4311 * tmp513
        - 1.0 * tmp4312;
    let tmp4314 = tmp4109 * tmp823;
    let tmp4315 = tmp4097 * tmp818;
    let tmp4316 = 4.0 * tmp2469;
    let tmp4317 = 4.0 * tmp2483;
    let tmp4318 = 4.0 * tmp2458;
    let tmp4319 = tmp1397 * tmp4292
        + tmp1398 * tmp4285
        + tmp1399 * tmp4294
        + tmp1442 * tmp4217
        + tmp1443 * tmp4216
        + tmp2093 * tmp4288
        + tmp2100 * tmp4289
        + tmp2464 * tmp4297
        + tmp2483 * tmp4318
        + tmp4283 * tmp740
        + tmp4295 * tmp495
        + tmp4296 * tmp495
        - 1.0 * tmp4299
        - 1.0 * tmp4300
        + tmp4301
        + tmp4302 * tmp818
        + tmp4303 * tmp834
        + tmp4304;
    let tmp4320 = tmp2175 * tmp2461;
    let tmp4321 = tmp2702 * tmp4104;
    let tmp4322 = tmp2065 * tmp2177;
    let tmp4323 = grad_phi_node_i[0] * tmp4306;
    let tmp4324 = tmp2175 * tmp2464;
    let tmp4325 = tmp2739 * tmp4111;
    let tmp4326 = tmp2056 * tmp2172;
    let tmp4327 = tmp1410 * tmp4097;
    let tmp4328 = tmp1397 * tmp4267
        + tmp1476 * tmp4211
        + tmp1477 * tmp4266
        + tmp1478 * tmp4239
        + tmp1482 * tmp4217
        + tmp2086 * tmp4263
        + tmp2099 * tmp4262
        + tmp2452 * tmp4274
        + tmp2469 * tmp4318
        + tmp4248 * tmp722
        + tmp4271 * tmp495
        + tmp4272 * tmp495
        + tmp4273 * tmp495
        + tmp4302 * tmp823
        + tmp4324
        + tmp4325
        - 1.0 * tmp4326
        - 1.0 * tmp4327;
    let tmp4329 = tmp1188 * tmp4279
        + tmp1274 * tmp4149
        + tmp1463 * tmp4284
        + tmp1464 * tmp4243
        + tmp1491 * tmp4216
        + tmp1492 * tmp4211
        + tmp2091 * tmp4261
        + tmp2097 * tmp4260
        + tmp2449 * tmp4290
        + tmp2455 * tmp4291
        + tmp4253 * tmp736
        + tmp4309 * tmp476
        + tmp4310 * tmp476
        + tmp4311 * tmp476
        - 1.0 * tmp4320
        - 1.0 * tmp4321
        + tmp4322
        + tmp4323;
    let tmp4330 = 2.0 * tmp2900
        + 2.0 * tmp2905
        + 2.0 * tmp2918
        + 2.0 * tmp2922
        + 2.0 * tmp2923
        + 2.0 * tmp2924
        + 2.0 * tmp2925
        + 2.0 * tmp2926
        + 2.0 * tmp2927
        + 2.0 * tmp2928
        + tmp2972
        - 1.0 * tmp2973;
    let tmp4331 = 2.0 * tmp2933
        + 2.0 * tmp2937
        + 2.0 * tmp2940
        + 2.0 * tmp2942
        + 2.0 * tmp2943
        + 2.0 * tmp2944
        + 2.0 * tmp2945
        + 2.0 * tmp2946
        + 2.0 * tmp2947
        + 2.0 * tmp2948
        + tmp3010
        - 1.0 * tmp3011;
    let tmp4332 = 2.0 * tmp2870
        + 2.0 * tmp2876
        + 2.0 * tmp2881
        + 2.0 * tmp2887
        + 2.0 * tmp2888
        + 2.0 * tmp2889
        + 2.0 * tmp2890
        + 2.0 * tmp2891
        + 2.0 * tmp2892
        + 2.0 * tmp2893
        - 2_f64 / 3.0 * tmp2894
        + (2_f64 / 3.0) * tmp4077
        + tmp4330
        + tmp4331;
    let tmp4333 = tmp2895 - 1.0 * tmp2896
        + 2.0 * tmp2978
        + 2.0 * tmp2981
        + 2.0 * tmp2983
        + 2.0 * tmp2986
        + 2.0 * tmp2987
        + 2.0 * tmp2988
        + 2.0 * tmp2989
        + 2.0 * tmp2990
        + 2.0 * tmp2991
        + 2.0 * tmp2992;
    let tmp4334 = (2_f64 / 3.0) * tmp2929
        + 2.0 * tmp2953
        + 2.0 * tmp2956
        + 2.0 * tmp2962
        + 2.0 * tmp2965
        + 2.0 * tmp2966
        + 2.0 * tmp2967
        + 2.0 * tmp2968
        + 2.0 * tmp2969
        + 2.0 * tmp2970
        + 2.0 * tmp2971
        - 2_f64 / 3.0 * tmp4079
        + tmp4331
        + tmp4333;
    let tmp4335 = 2.0 * tmp2996
        + 2.0 * tmp2998
        + 2.0 * tmp3001
        + 2.0 * tmp3003
        + 2.0 * tmp3004
        + 2.0 * tmp3005
        + 2.0 * tmp3006
        + 2.0 * tmp3007
        + 2.0 * tmp3008
        + 2.0 * tmp3009
        + (2_f64 / 3.0) * tmp4076
        - 2_f64 / 3.0 * tmp4078
        + tmp4330
        + tmp4333;
    let tmp4336 = 0.666666666666667 * tmp4168;
    let tmp4337 = tmp1819 * tmp2106;
    let tmp4338 = tmp281 * tmp4337;
    let tmp4339 = tmp1819 * tmp2337;
    let tmp4340 = tmp4336 + tmp4338 - 1.0 * tmp4339;
    let tmp4341 = 0.666666666666667 * tmp4162;
    let tmp4342 = tmp1819 * tmp2339;
    let tmp4343 = tmp283 * tmp4337;
    let tmp4344 = tmp4341 + tmp4342 - 1.0 * tmp4343;
    let tmp4345 = tmp1829 + tmp1984 + tmp2036;
    let tmp4346 = 0.666666666666667 * tmp4174;
    let tmp4347 = tmp1819 * tmp4232;
    let tmp4348 = tmp217 * tmp4337;
    let tmp4349 = -1.0 * tmp4346 - 1.0 * tmp4347 - 1.0 * tmp4348;
    let tmp4350 = tmp1845 + tmp2028;
    let tmp4351 = 0.666666666666667 * tmp4172;
    let tmp4352 = tmp219 * tmp4337;
    let tmp4353 = tmp1819 * tmp4227;
    let tmp4354 = tmp4351 - 1.0 * tmp4352 - 1.0 * tmp4353;
    let tmp4355 = 0.666666666666667 * tmp4179;
    let tmp4356 = q_old[3] * tmp1819;
    let tmp4357 = tmp2107 * tmp4356;
    let tmp4358 = tmp1819 * tmp4176;
    let tmp4359 = tmp4355 + tmp4357 - 1.0 * tmp4358;
    let tmp4360 = 0.666666666666667 * tmp4182;
    let tmp4361 = tmp227 * tmp4337;
    let tmp4362 = tmp231 * tmp4337;
    let tmp4363 = -1.0 * tmp4360 + tmp4361 - 1.0 * tmp4362;
    let tmp4364 = -1.0 * tmp4355 - 1.0 * tmp4357 + tmp4358;
    let tmp4365 = (4_f64 / 3.0) * tmp2175;
    let tmp4366 = 0.333333333333333 * tmp4162;
    let tmp4367 = 0.333333333333333 * tmp4168;
    let tmp4368 = tmp1883 + tmp1999;
    let tmp4369 = 0.333333333333333 * tmp4179;
    let tmp4370 = 0.333333333333333 * tmp4182;
    let tmp4371 = -1.0 * tmp4370;
    let tmp4372 = 0.333333333333333 * tmp4172;
    let tmp4373 = -1.0 * tmp4372;
    let tmp4374 = 0.333333333333333 * tmp4174;
    let tmp4375 = tmp1504 * tmp4232;
    let tmp4376 = -1.0 * tmp4375;
    let tmp4377 = tmp217 * tmp3016;
    let tmp4378 = -1.0 * tmp4377;
    let tmp4379 = -1.0 * tmp4374 + tmp4376 + tmp4378;
    let tmp4380 = tmp219 * tmp3016;
    let tmp4381 = tmp1504 * tmp4227;
    let tmp4382 = tmp4380 + tmp4381;
    let tmp4383 = -1.0 * tmp4366;
    let tmp4384 = tmp3057 + tmp3061 + tmp3062 + tmp3066;
    let tmp4385 = (2_f64 / 3.0) * tmp2175;
    let tmp4386 = tmp1635 * tmp4156
        + tmp1640 * tmp4110
        + tmp1643 * tmp4150
        + tmp1649 * tmp4098
        + tmp1944 * tmp4216
        + tmp1945 * tmp4211
        + tmp1947 * tmp4104
        + tmp2100 * tmp4385
        + tmp2899 * tmp4129
        + tmp2904 * tmp4117
        + tmp2917 * tmp4114
        + tmp2921 * tmp4126
        - 1_f64 / 3.0 * tmp4120
        - 2_f64 / 3.0 * tmp4133
        + tmp646 * (tmp1925 + tmp1957 + tmp4373 + tmp4379 + tmp4382)
        + tmp700 * (tmp1914 + tmp2010 + tmp3019 + tmp3024 + tmp4369 + tmp4371)
        + tmp722 * (tmp1943 + tmp4367 + tmp4383 + tmp4384)
        + tmp740 * (tmp3067 + tmp3161 + tmp4366 + tmp4367 + tmp4368);
    let tmp4387 = tmp1951 + tmp2004;
    let tmp4388 = -1.0 * tmp4380;
    let tmp4389 = -1.0 * tmp4381;
    let tmp4390 = tmp4388 + tmp4389;
    let tmp4391 = tmp4372 + tmp4390;
    let tmp4392 = tmp4375 + tmp4377;
    let tmp4393 = -1.0 * tmp4369;
    let tmp4394 = -1.0 * tmp4367;
    let tmp4395 = tmp3055 + tmp3056 + tmp3063 + tmp3160;
    let tmp4396 = tmp1753 * tmp4157
        + tmp1755 * tmp4112
        + tmp1757 * tmp4154
        + tmp1759 * tmp4105
        + tmp1973 * tmp4217
        + tmp1974 * tmp4216
        - 1.0 * tmp1975 * tmp4109
        + tmp2977 * tmp4127
        + tmp2980 * tmp4130
        + tmp2982 * tmp4115
        + tmp2985 * tmp4118
        + (1_f64 / 3.0) * tmp4124
        - 2_f64 / 3.0 * tmp4134
        + (2_f64 / 3.0) * tmp4135
        + tmp671 * (tmp1966 + tmp3014 + tmp3018 + tmp3022 + tmp3155 + tmp4371 + tmp4393)
        + tmp706 * (tmp1972 + tmp4366 + tmp4394 + tmp4395)
        + tmp729 * (tmp4379 + tmp4387 + tmp4391)
        + tmp744 * (tmp1928 + tmp1958 + tmp4374 + tmp4391 + tmp4392);
    let tmp4397 = -1.0 * tmp4341 - 1.0 * tmp4342 + tmp4343;
    let tmp4398 = -1.0 * tmp4336 - 1.0 * tmp4338 + tmp4339;
    let tmp4399 = tmp1831 + tmp1983 + tmp2037;
    let tmp4400 = tmp4360 - 1.0 * tmp4361 + tmp4362;
    let tmp4401 = tmp4346 + tmp4347 + tmp4348;
    let tmp4402 = tmp1880 + tmp1998;
    let tmp4403 = tmp4377 + tmp4380;
    let tmp4404 = tmp4375 + tmp4381;
    let tmp4405 = tmp1950 + tmp2003;
    let tmp4406 = tmp1679 * tmp4108
        + tmp1681 * tmp4155
        + tmp1684 * tmp4147
        + tmp1687 * tmp4091
        + tmp2020 * tmp4211
        + tmp2021 * tmp4217
        - 1.0 * tmp2022 * tmp4111
        + tmp2023 * tmp4097
        - 1.0 * tmp2098 * tmp4385
        + tmp2932 * tmp4113
        + tmp2936 * tmp4116
        + tmp2939 * tmp4125
        + tmp2941 * tmp4128
        + (2_f64 / 3.0) * tmp4132
        + tmp601 * (tmp3058 + tmp3064 + tmp4383 + tmp4394 + tmp4402)
        + tmp693 * (tmp4373 + tmp4374 + tmp4403 + tmp4404 + tmp4405)
        + tmp716 * (tmp1910 + tmp2008 + tmp3156 + tmp3158 + tmp4370 + tmp4393)
        + tmp736 * (tmp2018 + tmp3015 + tmp3017 + tmp3023 + tmp3154 + tmp4369 + tmp4370);
    let tmp4407 = tmp1841 + tmp2027;
    let tmp4408 = -1.0 * tmp4351 + tmp4352 + tmp4353;
    let tmp4409 = tmp3318 * tmp4083;
    let tmp4410 = 2.0 * tmp4409;
    let tmp4411 = tmp21 * tmp4410;
    let tmp4412 = tmp25 * tmp4410;
    let tmp4413 = -1.0 * tmp2225 + tmp23 * tmp4412 + tmp3374;
    let tmp4414 = tmp2222 + tmp2233 + tmp3377 + tmp3400 - 1.0 * tmp4411 + tmp4413;
    let tmp4415 = tmp4414 * tmp476;
    let tmp4416 = tmp23 * tmp4410;
    let tmp4417 = tmp20 * tmp4416;
    let tmp4418 = tmp19 * tmp4412 + tmp2324 + tmp3327;
    let tmp4419 = tmp2288 + tmp2317 - 1.0 * tmp3323 + tmp3347 + tmp4417 + tmp4418;
    let tmp4420 = tmp4419 * tmp476;
    let tmp4421 = tmp2290 + tmp2323 + tmp3325 + tmp3346 - 1.0 * tmp4417 + tmp4418;
    let tmp4422 = tmp4421 * tmp513;
    let tmp4423 = tmp2218 + tmp2224 + tmp3380 + tmp3405 + tmp4411 + tmp4413;
    let tmp4424 = tmp4423 * tmp495;
    let tmp4425 = tmp4419 * tmp537;
    let tmp4426 = tmp4414 * tmp545;
    let tmp4427 = tmp20 * tmp4412;
    let tmp4428 = tmp19 * tmp4416 + tmp2268;
    let tmp4429 = tmp2264 + tmp2265 + tmp2271 + tmp2272 + tmp3362 + tmp3367 - 1.0 * tmp4427 + tmp4428 + tmp593;
    let tmp4430 = tmp4429 * tmp495;
    let tmp4431 = tmp2266 + tmp2274 + tmp3360 - 1.0 * tmp3362 + tmp3366 + tmp4427 + tmp4428;
    let tmp4432 = tmp4431 * tmp513;
    let tmp4433 = tmp26 * tmp4409;
    let tmp4434 = tmp27 * tmp4409;
    let tmp4435 = tmp24 * tmp4409;
    let tmp4436 = -1.0 * tmp3337 + tmp4434 - 1.0 * tmp4435;
    let tmp4437 = tmp22 * tmp4409;
    let tmp4438 = -1.0 * tmp2210 - 1.0 * tmp4437;
    let tmp4439 = -1.0 * tmp2209 + tmp2242 + tmp3334 + tmp3342 + tmp4433 + tmp4436 + tmp4438;
    let tmp4440 = tmp4439 * tmp513;
    let tmp4441 = tmp4439 * tmp495;
    let tmp4442 = tmp4423 * tmp539;
    let tmp4443 = tmp4421 * tmp543;
    let tmp4444 = tmp2175 * tmp3241;
    let tmp4445 = tmp2098 * tmp3302;
    let tmp4446 = tmp2177 * tmp3242;
    let tmp4447 = tmp2099 * tmp3304;
    let tmp4448 = tmp2102 * tmp3299;
    let tmp4449 = tmp2172 * tmp3245;
    let tmp4450 = tmp2172 * tmp3240;
    let tmp4451 = tmp2097 * tmp3299;
    let tmp4452 = tmp2175 * tmp3243;
    let tmp4453 = tmp2100 * tmp3302;
    let tmp4454 = tmp2101 * tmp3304;
    let tmp4455 = tmp2177 * tmp3244;
    let tmp4456 = -1.0 * tmp4433;
    let tmp4457 = tmp2195 + tmp2209 + tmp2241 + tmp3337 + tmp3343 + tmp4434 + tmp4435 + tmp4438 + tmp4456;
    let tmp4458 = tmp4457 * tmp476;
    let tmp4459 = tmp2211 + tmp2243 + tmp3336 + tmp3341 + tmp4436 + tmp4437 + tmp4456;
    let tmp4460 = tmp4459 * tmp495;
    let tmp4461 = tmp4459 * tmp476;
    let tmp4462 = tmp4457 * tmp513;
    let tmp4463 = tmp4431 * tmp535;
    let tmp4464 = tmp4429 * tmp541;
    let tmp4465 = tmp3414 * tmp4159;
    let tmp4466 = tmp4465 * tmp5;
    let tmp4467 = tmp4 * tmp4466;
    let tmp4468 = tmp4465 * tmp7;
    let tmp4469 = tmp10 * tmp4468;
    let tmp4470 = tmp10 * tmp4466;
    let tmp4471 = tmp4 * tmp4468;
    let tmp4472 = -1.0 * tmp4471;
    let tmp4473 = tmp219 * tmp3422;
    let tmp4474 = -1.0 * tmp4473;
    let tmp4475 = tmp217 * tmp3422;
    let tmp4476 = -1.0 * tmp4475;
    let tmp4477 = tmp4474 + tmp4476;
    let tmp4478 = q_old[1] * tmp3252;
    let tmp4479 = tmp228 * tmp4478;
    let tmp4480 = -1.0 * tmp4479;
    let tmp4481 = q_old[2] * tmp3250;
    let tmp4482 = tmp228 * tmp4481;
    let tmp4483 = -1.0 * tmp4482;
    let tmp4484 = tmp4480 + tmp4483;
    let tmp4485 = -1.0 * tmp4470;
    let tmp4486 = tmp4479 + tmp4482;
    let tmp4487 = tmp4473 + tmp4475;
    let tmp4488 = -1.0 * tmp4467;
    let tmp4489 = -1.0 * tmp4469;
    let tmp4490 = tmp4474 + tmp4475 + tmp4480 + tmp4482;
    let tmp4491 = tmp3436 + tmp3445;
    let tmp4492 = tmp3441 + tmp3447;
    let tmp4493 = tmp4473 + tmp4476 + tmp4479 + tmp4483;
    let tmp4494 = tmp4466 * tmp7;
    let tmp4495 = tmp2384 + tmp2431 + tmp4494;
    let tmp4496 = tmp4465 * tmp87;
    let tmp4497 = tmp2388 + tmp2401 + tmp4496;
    let tmp4498 = -1.0 * tmp4494;
    let tmp4499 = tmp3250 * tmp634;
    let tmp4500 = -1.0 * tmp4499;
    let tmp4501 = tmp3249 * tmp651;
    let tmp4502 = -1.0 * tmp4501;
    let tmp4503 = tmp4465 * tmp9;
    let tmp4504 = q_old[0] * tmp3252;
    let tmp4505 = tmp214 * tmp4504;
    let tmp4506 = tmp3249 * tmp666;
    let tmp4507 = tmp4505 + tmp4506;
    let tmp4508 = tmp11 * tmp4465;
    let tmp4509 = tmp4465 * tmp8;
    let tmp4510 = tmp4193 + tmp4508 - 1.0 * tmp4509;
    let tmp4511 = tmp4465 * tmp6;
    let tmp4512 = tmp4189 - 1.0 * tmp4511;
    let tmp4513 = tmp4185 + tmp4500 + tmp4502 + tmp4503 + tmp4507 + tmp4510 + tmp4512;
    let tmp4514 = -1.0 * tmp4496;
    let tmp4515 = q_old[0] * tmp2441;
    let tmp4516 = tmp4500 - 1.0 * tmp4505;
    let tmp4517 = tmp4502 - 1.0 * tmp4506;
    let tmp4518 = tmp2440 * tmp683 + tmp4184 - 1.0 * tmp4503 + tmp4517;
    let tmp4519 = tmp4180 + tmp4510 + tmp4511 - 1.0 * tmp4515 + tmp4516 + tmp4518;
    let tmp4520 = tmp4177 + tmp4192 + tmp4499 + tmp4505 + tmp4508 + tmp4509 + tmp4512 + tmp4515 + tmp4518;
    let tmp4521 = 0.5 * tmp4470;
    let tmp4522 = tmp219 * tmp3672;
    let tmp4523 = q_old[1] * tmp3674;
    let tmp4524 = tmp4521 - 1.0 * tmp4522 - 1.0 * tmp4523;
    let tmp4525 = 0.5 * tmp4471;
    let tmp4526 = q_old[2] * tmp3679;
    let tmp4527 = tmp217 * tmp3672;
    let tmp4528 = -1.0 * tmp4525 - 1.0 * tmp4526 - 1.0 * tmp4527;
    let tmp4529 = tmp2773 + tmp4524 + tmp4528;
    let tmp4530 = tmp4525 + tmp4526 + tmp4527;
    let tmp4531 = tmp1370 + tmp2734 + tmp2771 + tmp4524 + tmp4530;
    let tmp4532 = 0.5 * tmp4494;
    let tmp4533 = tmp2677 - 1.0 * tmp4532;
    let tmp4534 = 0.5 * tmp4496;
    let tmp4535 = tmp2680 - 1.0 * tmp4534;
    let tmp4536 = tmp2731 + tmp3705 + tmp3725 + tmp3751 + tmp4533 + tmp4535;
    let tmp4537 = 0.5 * tmp4469;
    let tmp4538 = q_old[2] * tmp2648;
    let tmp4539 = tmp219 * tmp2642;
    let tmp4540 = -1.0 * tmp4537 + tmp4538 - 1.0 * tmp4539;
    let tmp4541 = 0.5 * tmp4467;
    let tmp4542 = tmp217 * tmp2642;
    let tmp4543 = q_old[1] * tmp2644;
    let tmp4544 = tmp4541 + tmp4542 - 1.0 * tmp4543;
    let tmp4545 = tmp1366 + tmp3729 + tmp3759 + tmp4540 + tmp4544;
    let tmp4546 = tmp186 * tmp4414;
    let tmp4547 = grad_phi_node_i[1] * tmp4546;
    let tmp4548 = tmp2721 * tmp3504;
    let tmp4549 = tmp2452 * tmp3713;
    let tmp4550 = tmp4457 * tmp495;
    let tmp4551 = tmp4431 * tmp495;
    let tmp4552 = tmp4421 * tmp839;
    let tmp4553 = tmp1410 * tmp4423;
    let tmp4554 = tmp2754 * tmp3213;
    let tmp4555 = tmp2068 * tmp3722;
    let tmp4556 = -1.0 * tmp4541 - 1.0 * tmp4542 + tmp4543;
    let tmp4557 = tmp3732 + tmp4540 + tmp4556;
    let tmp4558 = tmp2686 + tmp4534;
    let tmp4559 = tmp2728 + tmp2764 + tmp2785 + tmp3754 + tmp4533 + tmp4558;
    let tmp4560 = tmp186 * tmp4419;
    let tmp4561 = grad_phi_node_i[2] * tmp4560;
    let tmp4562 = tmp2721 * tmp3550;
    let tmp4563 = tmp2509 * tmp3713;
    let tmp4564 = tmp1041 * tmp4429;
    let tmp4565 = tmp4423 * tmp968;
    let tmp4566 = tmp4439 * tmp964;
    let tmp4567 = tmp1311 * tmp4421;
    let tmp4568 = tmp2061 * tmp3737;
    let tmp4569 = tmp2704 * tmp3207;
    let tmp4570 = tmp1305 * tmp4545
        + tmp1415 * tmp4513
        + tmp1416 * tmp4559
        + tmp1417 * tmp4531
        + tmp1421 * tmp4519
        + tmp2102 * tmp3792
        + tmp2251 * tmp3567
        + tmp2515 * tmp3717
        + tmp2519 * tmp3791
        + tmp2780 * tmp3245
        + tmp2783 * tmp3232
        + tmp2798 * tmp3560
        + tmp2800 * tmp3556
        + tmp2832 * tmp4459
        + tmp4557 * tmp729
        + tmp4561
        + tmp4562
        + tmp4563
        + tmp4564 * tmp513
        + tmp4565 * tmp513
        + tmp4566 * tmp513
        - 1.0 * tmp4567
        - 1.0 * tmp4568
        - 1.0 * tmp4569;
    let tmp4571 = tmp4537 - 1.0 * tmp4538 + tmp4539;
    let tmp4572 = tmp3761 + tmp4544 + tmp4571;
    let tmp4573 = tmp2690 + tmp4532;
    let tmp4574 = tmp2730 + tmp2765 + tmp2786 + tmp3707 + tmp3726 + tmp4535 + tmp4573;
    let tmp4575 = -1.0 * tmp4521 + tmp4522 + tmp4523;
    let tmp4576 = tmp1258 + tmp2735 + tmp2772 + tmp4528 + tmp4575;
    let tmp4577 = tmp1219 + tmp3731 + tmp3760 + tmp4556 + tmp4571;
    let tmp4578 = tmp4423 * tmp823;
    let tmp4579 = tmp4429 * tmp818;
    let tmp4580 = tmp2736 + tmp4530 + tmp4575;
    let tmp4581 = tmp2766 + tmp3706 + tmp3724 + tmp3753 + tmp4558 + tmp4573;
    let tmp4582 = tmp4439 * tmp960;
    let tmp4583 = tmp4421 * tmp958;
    let tmp4584 = tmp1397 * tmp4576
        + tmp1398 * tmp4581
        + tmp1399 * tmp4577
        + tmp1442 * tmp4520
        + tmp1443 * tmp4513
        + tmp2100 * tmp3768
        + tmp2254 * tmp3569
        + tmp2464 * tmp3774
        + tmp2715 * tmp3516
        + tmp2717 * tmp3510
        + tmp2748 * tmp3527
        + tmp2777 * tmp3236
        + tmp2779 * tmp3243
        - 1.0 * tmp4547
        - 1.0 * tmp4548
        - 1.0 * tmp4549
        + tmp4550 * tmp818
        + tmp4551 * tmp834
        + tmp4553
        + tmp4554
        + tmp4555
        + tmp4580 * tmp740
        + tmp4582 * tmp495
        + tmp4583 * tmp495;
    let tmp4585 = grad_phi_node_i[0] * tmp4546;
    let tmp4586 = tmp2721 * tmp3196;
    let tmp4587 = tmp2050 * tmp3713;
    let tmp4588 = tmp4439 * tmp476;
    let tmp4589 = tmp1281 * tmp4423;
    let tmp4590 = tmp2754 * tmp3501;
    let tmp4591 = tmp2449 * tmp3722;
    let tmp4592 = tmp1311 * tmp4431;
    let tmp4593 = tmp2704 * tmp3216;
    let tmp4594 = tmp2072 * tmp3737;
    let tmp4595 = tmp4414 * tmp830;
    let tmp4596 = tmp4419 * tmp826;
    let tmp4597 = tmp4459 * tmp839;
    let tmp4598 = tmp1318 * tmp4429;
    let tmp4599 = tmp2754 * tmp3560;
    let tmp4600 = tmp2519 * tmp3722;
    let tmp4601 = tmp1286 * tmp4519
        + tmp1287 * tmp4572
        + tmp1295 * tmp4529
        + tmp1305 * tmp4536
        + tmp1309 * tmp4520
        + tmp2101 * tmp3742
        + tmp2257 * tmp3521
        + tmp2509 * tmp3743
        + tmp2515 * tmp3719
        + tmp2711 * tmp3244
        + tmp2713 * tmp3238
        + tmp2749 * tmp3550
        + tmp2751 * tmp3556
        + tmp4462 * tmp823
        + tmp4574 * tmp744
        + tmp4592
        + tmp4593
        + tmp4594
        + tmp4595 * tmp513
        + tmp4596 * tmp513
        + tmp4597 * tmp513
        - 1.0 * tmp4598
        - 1.0 * tmp4599
        - 1.0 * tmp4600;
    let tmp4602 = tmp4419 * tmp960;
    let tmp4603 = tmp4414 * tmp964;
    let tmp4604 = tmp1188 * tmp4580
        + tmp1222 * tmp4581
        + tmp1246 * tmp4577
        + tmp1363 * tmp4520
        + tmp1374 * tmp4513
        + tmp1377 * tmp4431
        + tmp2098 * tmp3774
        + tmp2247 * tmp3527
        + tmp2455 * tmp3767
        + tmp2461 * tmp3768
        + tmp2715 * tmp3241
        + tmp2717 * tmp3225
        + tmp2777 * tmp3507
        + tmp2779 * tmp3513
        + tmp4458 * tmp818
        + tmp4576 * tmp716
        + tmp4582 * tmp476
        + tmp4583 * tmp476
        - 1.0 * tmp4585
        - 1.0 * tmp4586
        - 1.0 * tmp4587
        + tmp4589
        + tmp4590
        + tmp4591;
    let tmp4605 = tmp2702 * tmp4421;
    let tmp4606 = tmp2461 * tmp3737;
    let tmp4607 = tmp2704 * tmp3513;
    let tmp4608 = grad_phi_node_i[0] * tmp4560;
    let tmp4609 = tmp2721 * tmp3210;
    let tmp4610 = tmp2065 * tmp3713;
    let tmp4611 = tmp2739 * tmp4431;
    let tmp4612 = tmp2704 * tmp3516;
    let tmp4613 = tmp2464 * tmp3737;
    let tmp4614 = tmp1410 * tmp4429;
    let tmp4615 = tmp2754 * tmp3202;
    let tmp4616 = tmp2056 * tmp3722;
    let tmp4617 = tmp1397 * tmp4574
        + tmp1476 * tmp4519
        + tmp1477 * tmp4572
        + tmp1478 * tmp4529
        + tmp1482 * tmp4520
        + tmp2099 * tmp3743
        + tmp2248 * tmp3541
        + tmp2452 * tmp3742
        + tmp2711 * tmp3504
        + tmp2713 * tmp3510
        + tmp2748 * tmp3521
        + tmp2749 * tmp3242
        + tmp2751 * tmp3229
        + tmp4536 * tmp722
        + tmp4550 * tmp823
        + tmp4595 * tmp495
        + tmp4596 * tmp495
        + tmp4597 * tmp495
        + tmp4611
        + tmp4612
        + tmp4613
        - 1.0 * tmp4614
        - 1.0 * tmp4615
        - 1.0 * tmp4616;
    let tmp4618 = tmp1188 * tmp4557
        + tmp1274 * tmp4459
        + tmp1463 * tmp4559
        + tmp1464 * tmp4531
        + tmp1491 * tmp4513
        + tmp1492 * tmp4519
        + tmp2097 * tmp3791
        + tmp2253 * tmp3590
        + tmp2449 * tmp3792
        + tmp2455 * tmp3769
        + tmp2780 * tmp3501
        + tmp2783 * tmp3507
        + tmp2798 * tmp3240
        + tmp2800 * tmp3234
        + tmp4545 * tmp736
        + tmp4564 * tmp476
        + tmp4565 * tmp476
        + tmp4566 * tmp476
        - 1.0 * tmp4605
        - 1.0 * tmp4606
        - 1.0 * tmp4607
        + tmp4608
        + tmp4609
        + tmp4610;
    let tmp4619 = 0.666666666666667 * tmp4467;
    let tmp4620 = 0.666666666666667 * tmp4469;
    let tmp4621 = tmp4376 + tmp4620;
    let tmp4622 = 0.666666666666667 * tmp4470;
    let tmp4623 = 0.666666666666667 * tmp4471;
    let tmp4624 = -1.0 * tmp4623;
    let tmp4625 = tmp219 * tmp3963;
    let tmp4626 = -1.0 * tmp4625;
    let tmp4627 = tmp217 * tmp3963;
    let tmp4628 = -1.0 * tmp4627;
    let tmp4629 = tmp4626 + tmp4628;
    let tmp4630 = tmp1504 * tmp4478;
    let tmp4631 = -1.0 * tmp4630;
    let tmp4632 = tmp1504 * tmp4481;
    let tmp4633 = -1.0 * tmp4632;
    let tmp4634 = tmp4631 + tmp4633;
    let tmp4635 = 0.666666666666667 * tmp4494;
    let tmp4636 = tmp3030 + tmp4635;
    let tmp4637 = 0.666666666666667 * tmp4496;
    let tmp4638 = tmp3027 + tmp3042 + tmp3964 - 1.0 * tmp4637;
    let tmp4639 = tmp3031 + tmp3045 + tmp3966 - 1.0 * tmp4635;
    let tmp4640 = 0.333333333333333 * tmp4469;
    let tmp4641 = tmp219 * tmp3079;
    let tmp4642 = tmp1583 * tmp4232;
    let tmp4643 = tmp4640 + tmp4641 - 1.0 * tmp4642;
    let tmp4644 = 0.333333333333333 * tmp4467;
    let tmp4645 = tmp217 * tmp3079;
    let tmp4646 = tmp1583 * tmp4227;
    let tmp4647 = tmp4644 + tmp4645 - 1.0 * tmp4646;
    let tmp4648 = 0.333333333333333 * tmp4494;
    let tmp4649 = 0.333333333333333 * tmp4496;
    let tmp4650 = tmp3111 + tmp3133 - 1.0 * tmp4649;
    let tmp4651 = 0.333333333333333 * tmp4470;
    let tmp4652 = tmp219 * tmp4006;
    let tmp4653 = tmp1583 * tmp4478;
    let tmp4654 = -1.0 * tmp4651 + tmp4652 + tmp4653;
    let tmp4655 = 0.333333333333333 * tmp4471;
    let tmp4656 = tmp1583 * tmp4481;
    let tmp4657 = tmp217 * tmp4006;
    let tmp4658 = -1.0 * tmp4655 - 1.0 * tmp4656 - 1.0 * tmp4657;
    let tmp4659 = -1.0 * tmp4644 - 1.0 * tmp4645 + tmp4646;
    let tmp4660 = tmp1635 * tmp4441
        + tmp1640 * tmp4424
        + tmp1643 * tmp4460
        + tmp1649 * tmp4430
        + tmp1944 * tmp4513
        + tmp1945 * tmp4519
        + tmp1947 * tmp4421
        + tmp2185 * tmp3878
        + tmp2191 * tmp3872
        + tmp2248 * tmp3881
        + tmp2254 * tmp3867
        + tmp2899 * tmp3354
        + tmp2904 * tmp3314
        + tmp2917 * tmp3311
        + tmp2921 * tmp3351
        - 1_f64 / 3.0 * tmp4425
        - 1_f64 / 3.0 * tmp4446
        - 1_f64 / 3.0 * tmp4447
        + (1_f64 / 3.0) * tmp4452
        + (1_f64 / 3.0) * tmp4453
        + tmp646 * (tmp1969 + tmp3143 + tmp3177 + tmp4654 + tmp4658)
        + tmp700 * (tmp3103 + tmp3107 + tmp3114 + tmp3123 + tmp3173 + tmp4011 + tmp4648 + tmp4650)
        + tmp722 * (tmp1957 + tmp4015 + tmp4062 + tmp4643 + tmp4659)
        + tmp740 * (tmp4065 + tmp4643 + tmp4647);
    let tmp4661 = tmp4651 - 1.0 * tmp4652 - 1.0 * tmp4653;
    let tmp4662 = tmp4655 + tmp4656 + tmp4657;
    let tmp4663 = tmp3108 + tmp3122 - 1.0 * tmp4648;
    let tmp4664 = -1.0 * tmp4640 - 1.0 * tmp4641 + tmp4642;
    let tmp4665 = tmp1753 * tmp4462
        + tmp1755 * tmp4432
        + tmp1757 * tmp4440
        + tmp1759 * tmp4422
        + tmp1973 * tmp4520
        + tmp1974 * tmp4513
        - 1.0 * tmp1975 * tmp4423
        + tmp2188 * tmp3931
        + tmp2193 * tmp3933
        + tmp2251 * tmp3927
        + tmp2257 * tmp3929
        + tmp2977 * tmp3352
        + tmp2980 * tmp3355
        + tmp2982 * tmp3312
        + tmp2985 * tmp3315
        + (1_f64 / 3.0) * tmp4426
        - 1_f64 / 3.0 * tmp4448
        - 1_f64 / 3.0 * tmp4449
        + (1_f64 / 3.0) * tmp4454
        + (1_f64 / 3.0) * tmp4455
        + tmp671 * (tmp3117 + tmp3137 + tmp4010 + tmp4060 + tmp4650 + tmp4663)
        + tmp706 * (tmp1928 + tmp4019 + tmp4064 + tmp4647 + tmp4664)
        + tmp729 * (tmp3179 + tmp4658 + tmp4661)
        + tmp744 * (tmp1942 + tmp3145 + tmp3178 + tmp4661 + tmp4662);
    let tmp4666 = -1.0 * tmp4620;
    let tmp4667 = tmp4378 - 1.0 * tmp4619;
    let tmp4668 = tmp3165 + tmp4637;
    let tmp4669 = tmp4626 + tmp4627 + tmp4631 + tmp4632;
    let tmp4670 = tmp1817 + tmp1981;
    let tmp4671 = tmp3985 + tmp3986 + tmp3993 + tmp4053;
    let tmp4672 = tmp1679 * tmp4420
        + tmp1681 * tmp4461
        + tmp1684 * tmp4458
        + tmp1687 * tmp4415
        + tmp2020 * tmp4519
        + tmp2021 * tmp4520
        - 1.0 * tmp2022 * tmp4431
        + tmp2023 * tmp4429
        + tmp2183 * tmp3892
        + tmp2190 * tmp3894
        + tmp2247 * tmp3897
        + tmp2253 * tmp3899
        + tmp2932 * tmp3310
        + tmp2936 * tmp3313
        + tmp2939 * tmp3350
        + tmp2941 * tmp3353
        - 1_f64 / 3.0 * tmp4444
        - 1_f64 / 3.0 * tmp4445
        + (1_f64 / 3.0) * tmp4450
        + (1_f64 / 3.0) * tmp4451
        + tmp601 * (tmp4020 + tmp4659 + tmp4664)
        + tmp693 * (tmp3146 + tmp4654 + tmp4662)
        + tmp716 * (tmp3116 + tmp3121 + tmp3136 + tmp3140 + tmp3172 + tmp4061 + tmp4649 + tmp4663)
        + tmp736 * (tmp3124 + tmp3141 + tmp3174 + tmp4005 + tmp4059 + tmp4648 + tmp4649);
    let tmp4673 = -1.0 * tmp4622;
    let tmp4674 = tmp4630 + tmp4632;
    let tmp4675 = tmp4625 + tmp4627;
    let tmp4676 = tmp4625 + tmp4628 + tmp4630 + tmp4633;
    let tmp4677 = tmp1831 + tmp1984;
    let tmp4678 = tmp3987 + tmp3991 + tmp3992 + tmp3996;
    let tmp4679 = tmp3215 * tmp535;
    let tmp4680 = tmp3212 * tmp539;
    let tmp4681 = tmp3201 * tmp541;
    let tmp4682 = tmp3206 * tmp543;
    let tmp4683 = tmp3190 * tmp438;
    let tmp4684 = tmp4683 * tmp73;
    let tmp4685 = tmp4683 * tmp79;
    let tmp4686 = tmp3189 * tmp3318 * tmp447;
    let tmp4687 = 2.0 * tmp4686;
    let tmp4688 = tmp21 * tmp4687;
    let tmp4689 = tmp25 * tmp4687;
    let tmp4690 = tmp23 * tmp4689 - 1.0 * tmp4683 * tmp70 + tmp4683 * tmp77;
    let tmp4691 = tmp4089 + tmp4684 - 1.0 * tmp4685 - 1.0 * tmp4688 + tmp4690;
    let tmp4692 = tmp4691 * tmp476;
    let tmp4693 = tmp3190 * tmp441;
    let tmp4694 = tmp3190 * tmp454;
    let tmp4695 = tmp20 * tmp4689;
    let tmp4696 = tmp23 * tmp4687;
    let tmp4697 = tmp19 * tmp4696 + tmp3190 * tmp439 + tmp3190 * tmp452;
    let tmp4698 = tmp4693 + tmp4694 - 1.0 * tmp4695 + tmp4697 + tmp493;
    let tmp4699 = tmp4698 * tmp495;
    let tmp4700 = tmp3190 * tmp444;
    let tmp4701 = tmp3190 * tmp446;
    let tmp4702 = tmp20 * tmp4696;
    let tmp4703 = tmp19 * tmp4689 - 1.0 * tmp3190 * tmp443 + tmp3190 * tmp445;
    let tmp4704 = tmp3404 + tmp4700 - 1.0 * tmp4701 - 1.0 * tmp4702 + tmp4703 + tmp510;
    let tmp4705 = tmp4704 * tmp513;
    let tmp4706 = tmp4106 - 1.0 * tmp4700 + tmp4701 + tmp4702 + tmp4703;
    let tmp4707 = tmp4706 * tmp476;
    let tmp4708 = -1.0 * tmp4684 + tmp4685 + tmp4688 + tmp4690 + tmp521;
    let tmp4709 = tmp4708 * tmp495;
    let tmp4710 = -1.0 * tmp4693 - 1.0 * tmp4694 + tmp4695 + tmp4697 + tmp485 + tmp487 + tmp492;
    let tmp4711 = tmp4710 * tmp513;
    let tmp4712 = 4.0 * tmp3196;
    let tmp4713 = 4.0 * tmp3202;
    let tmp4714 = 4.0 * tmp3207;
    let tmp4715 = 4.0 * tmp3210;
    let tmp4716 = 4.0 * tmp3213;
    let tmp4717 = 4.0 * tmp3216;
    let tmp4718 = tmp4710 * tmp535;
    let tmp4719 = tmp4706 * tmp537;
    let tmp4720 = tmp4708 * tmp539;
    let tmp4721 = tmp4698 * tmp541;
    let tmp4722 = tmp4704 * tmp543;
    let tmp4723 = tmp4691 * tmp545;
    let tmp4724 = 4.0 * tmp3225;
    let tmp4725 = 4.0 * tmp3229;
    let tmp4726 = 4.0 * tmp3232;
    let tmp4727 = 4.0 * tmp3234;
    let tmp4728 = 4.0 * tmp3236;
    let tmp4729 = 4.0 * tmp3238;
    let tmp4730 = 2.0 * tmp3302;
    let tmp4731 = tmp3240 * tmp3299;
    let tmp4732 = tmp3242 * tmp3304;
    let tmp4733 = tmp3245 * tmp3299;
    let tmp4734 = tmp3244 * tmp3304;
    let tmp4735 = tmp27 * tmp4686;
    let tmp4736 = -1.0 * tmp3190 * tmp499;
    let tmp4737 = tmp24 * tmp4686;
    let tmp4738 = tmp3190 * tmp497;
    let tmp4739 = tmp26 * tmp4686;
    let tmp4740 = tmp3190 * tmp502;
    let tmp4741 = -1.0 * tmp4739 + tmp4740;
    let tmp4742 = tmp22 * tmp4686;
    let tmp4743 = tmp3190 * tmp503;
    let tmp4744 = -1.0 * tmp4742 - 1.0 * tmp4743 + tmp573;
    let tmp4745 = tmp4735 + tmp4736 + tmp4737 + tmp4738 + tmp4741 + tmp4744 + tmp570 + tmp582;
    let tmp4746 = tmp4745 * tmp476;
    let tmp4747 = tmp4735 + tmp4736 - 1.0 * tmp4737 - 1.0 * tmp4738;
    let tmp4748 = tmp4151 + tmp4741 + tmp4742 + tmp4743 + tmp4747 + tmp569 + tmp572 + tmp594;
    let tmp4749 = tmp4748 * tmp495;
    let tmp4750 = tmp4152 + tmp4739 - 1.0 * tmp4740 + tmp4744 + tmp4747;
    let tmp4751 = tmp4750 * tmp513;
    let tmp4752 = tmp4748 * tmp476;
    let tmp4753 = tmp4750 * tmp495;
    let tmp4754 = tmp4745 * tmp513;
    let tmp4755 = tmp3248 * tmp3414 * tmp605;
    let tmp4756 = tmp4755 * tmp5;
    let tmp4757 = tmp4 * tmp4756;
    let tmp4758 = tmp3486 + tmp4757;
    let tmp4759 = tmp4755 * tmp7;
    let tmp4760 = tmp10 * tmp4759;
    let tmp4761 = tmp3484 + tmp4760;
    let tmp4762 = tmp10 * tmp4756;
    let tmp4763 = tmp214 * tmp3249;
    let tmp4764 = tmp281 * tmp4763;
    let tmp4765 = tmp2418 * tmp3252;
    let tmp4766 = tmp4762 + tmp4764 - 1.0 * tmp4765;
    let tmp4767 = tmp4 * tmp4759;
    let tmp4768 = tmp283 * tmp4763;
    let tmp4769 = tmp2415 * tmp3250;
    let tmp4770 = tmp4767 + tmp4768 - 1.0 * tmp4769;
    let tmp4771 = tmp4756 * tmp7;
    let tmp4772 = tmp4499 + tmp4771;
    let tmp4773 = tmp4755 * tmp87;
    let tmp4774 = tmp4517 + tmp4773;
    let tmp4775 = -1.0 * tmp4767 - 1.0 * tmp4768 + tmp4769;
    let tmp4776 = tmp4516 - 1.0 * tmp4771;
    let tmp4777 = tmp3475 - 1.0 * tmp4757;
    let tmp4778 = tmp4501 - 1.0 * tmp4773;
    let tmp4779 = tmp3477 - 1.0 * tmp4760;
    let tmp4780 = -1.0 * tmp4762 - 1.0 * tmp4764 + tmp4765;
    let tmp4781 = tmp4755 * tmp6;
    let tmp4782 = tmp3252 * tmp4205;
    let tmp4783 = tmp11 * tmp4755;
    let tmp4784 = tmp4755 * tmp8;
    let tmp4785 = tmp3250 * tmp4200;
    let tmp4786 = -1.0 * tmp3249 * tmp3413;
    let tmp4787 = tmp4783 - 1.0 * tmp4784 - 1.0 * tmp4785 + tmp4786;
    let tmp4788 = tmp4755 * tmp9;
    let tmp4789 = tmp3249 * tmp3412;
    let tmp4790 = -1.0 * tmp4788 + tmp4789;
    let tmp4791 = tmp4213 + tmp4781 + tmp4782 + tmp4787 + tmp4790 + tmp775 + tmp782;
    let tmp4792 = -1.0 * tmp4781 - 1.0 * tmp4782;
    let tmp4793 = tmp4214 + tmp4787 + tmp4788 - 1.0 * tmp4789 + tmp4792;
    let tmp4794 = tmp4783 + tmp4784 + tmp4785 + tmp4786 + tmp4790 + tmp4792 + tmp787;
    let tmp4795 = 0.5 * tmp4762;
    let tmp4796 = 0.5 * tmp4767;
    let tmp4797 = -1.0 * tmp4796;
    let tmp4798 = tmp3455 + tmp4230 + tmp4795 + tmp4797;
    let tmp4799 = tmp1220 + tmp4491 + tmp4795 + tmp4796;
    let tmp4800 = 0.5 * tmp4771;
    let tmp4801 = -1.0 * tmp4800;
    let tmp4802 = 0.5 * tmp4773;
    let tmp4803 = -1.0 * tmp4802;
    let tmp4804 = tmp1231 + tmp1350 + tmp3456 + tmp4801 + tmp4803;
    let tmp4805 = 0.5 * tmp4760;
    let tmp4806 = -1.0 * tmp4805;
    let tmp4807 = 0.5 * tmp4757;
    let tmp4808 = tmp1254 + tmp1370 + tmp4490 + tmp4806 + tmp4807;
    let tmp4809 = tmp3196 * tmp3304;
    let tmp4810 = tmp186 * tmp4691;
    let tmp4811 = grad_phi_node_i[0] * tmp4810;
    let tmp4812 = tmp3299 * tmp3501;
    let tmp4813 = tmp4750 * tmp476;
    let tmp4814 = tmp4704 * tmp839;
    let tmp4815 = 4.0 * tmp3595;
    let tmp4816 = 4.0 * tmp3590;
    let tmp4817 = 4.0 * tmp3533;
    let tmp4818 = 4.0 * tmp3541;
    let tmp4819 = tmp1281 * tmp4708;
    let tmp4820 = tmp4265 + tmp4484 + tmp4487 + tmp4805 + tmp4807;
    let tmp4821 = tmp1299 + tmp3428 + tmp3454 + tmp4800 + tmp4803;
    let tmp4822 = tmp3216 * tmp3302;
    let tmp4823 = tmp1311 * tmp4710;
    let tmp4824 = tmp3299 * tmp3560;
    let tmp4825 = tmp4691 * tmp830;
    let tmp4826 = tmp4706 * tmp826;
    let tmp4827 = tmp4748 * tmp839;
    let tmp4828 = 4.0 * tmp3531;
    let tmp4829 = tmp1318 * tmp4698;
    let tmp4830 = tmp1286 * tmp4791
        + tmp1287 * tmp4820
        + tmp1295 * tmp4798
        + tmp1305 * tmp4804
        + tmp1309 * tmp4794
        + tmp3244 * tmp4828
        + tmp3521 * tmp4729
        + tmp3550 * tmp4817
        + tmp3556 * tmp4818
        + tmp4754 * tmp823
        + tmp4821 * tmp744
        + tmp4822
        + tmp4823
        - 1.0 * tmp4824
        + tmp4825 * tmp513
        + tmp4826 * tmp513
        + tmp4827 * tmp513
        - 1.0 * tmp4829;
    let tmp4831 = -1.0 * tmp4807;
    let tmp4832 = tmp4278 + tmp4477 + tmp4486 + tmp4806 + tmp4831;
    let tmp4833 = -1.0 * tmp4795;
    let tmp4834 = tmp3458 + tmp4280 + tmp4796 + tmp4833;
    let tmp4835 = tmp1344 + tmp3425 + tmp3453 + tmp4801 + tmp4802;
    let tmp4836 = tmp1235 + tmp1351 + tmp3457 + tmp4800 + tmp4802;
    let tmp4837 = tmp4706 * tmp960;
    let tmp4838 = tmp4691 * tmp964;
    let tmp4839 = 4.0 * tmp3569;
    let tmp4840 = 4.0 * tmp3563;
    let tmp4841 = 4.0 * tmp3573;
    let tmp4842 = 4.0 * tmp3567;
    let tmp4843 = tmp1368 + tmp4492 + tmp4797 + tmp4833;
    let tmp4844 = tmp1258 + tmp1371 + tmp4493 + tmp4805 + tmp4831;
    let tmp4845 = tmp4750 * tmp960;
    let tmp4846 = tmp4704 * tmp958;
    let tmp4847 = 4.0 * tmp3537;
    let tmp4848 = tmp1188 * tmp4834
        + tmp1222 * tmp4836
        + tmp1246 * tmp4844
        + tmp1363 * tmp4794
        + tmp1374 * tmp4793
        + tmp1377 * tmp4710
        + tmp3241 * tmp4847
        + tmp3507 * tmp4839
        + tmp3513 * tmp4840
        + tmp3527 * tmp4724
        + tmp4746 * tmp818
        + tmp476 * tmp4845
        + tmp476 * tmp4846
        - 1.0 * tmp4809
        - 1.0 * tmp4811
        + tmp4812
        + tmp4819
        + tmp4843 * tmp716;
    let tmp4849 = tmp3304 * tmp3504;
    let tmp4850 = grad_phi_node_i[1] * tmp4810;
    let tmp4851 = tmp3213 * tmp3299;
    let tmp4852 = tmp4745 * tmp495;
    let tmp4853 = tmp4710 * tmp495;
    let tmp4854 = tmp1410 * tmp4708;
    let tmp4855 = tmp3304 * tmp3550;
    let tmp4856 = tmp186 * tmp4706;
    let tmp4857 = grad_phi_node_i[2] * tmp4856;
    let tmp4858 = tmp3207 * tmp3302;
    let tmp4859 = tmp1041 * tmp4698;
    let tmp4860 = tmp4708 * tmp968;
    let tmp4861 = tmp4750 * tmp964;
    let tmp4862 = tmp1311 * tmp4704;
    let tmp4863 = tmp1305 * tmp4808
        + tmp1415 * tmp4793
        + tmp1416 * tmp4835
        + tmp1417 * tmp4799
        + tmp1421 * tmp4791
        + tmp2832 * tmp4748
        + tmp3232 * tmp4842
        + tmp3245 * tmp4841
        + tmp3556 * tmp4816
        + tmp3560 * tmp4815
        + tmp4832 * tmp729
        + tmp4855
        + tmp4857
        - 1.0 * tmp4858
        + tmp4859 * tmp513
        + tmp4860 * tmp513
        + tmp4861 * tmp513
        - 1.0 * tmp4862;
    let tmp4864 = tmp4708 * tmp823;
    let tmp4865 = tmp4698 * tmp818;
    let tmp4866 = 4.0 * tmp3521;
    let tmp4867 = 4.0 * tmp3527;
    let tmp4868 = 4.0 * tmp3510;
    let tmp4869 = tmp1397 * tmp4843
        + tmp1398 * tmp4836
        + tmp1399 * tmp4844
        + tmp1442 * tmp4794
        + tmp1443 * tmp4793
        + tmp3236 * tmp4839
        + tmp3243 * tmp4840
        + tmp3516 * tmp4847
        + tmp3527 * tmp4868
        + tmp4834 * tmp740
        + tmp4845 * tmp495
        + tmp4846 * tmp495
        - 1.0 * tmp4849
        - 1.0 * tmp4850
        + tmp4851
        + tmp4852 * tmp818
        + tmp4853 * tmp834
        + tmp4854;
    let tmp4870 = tmp3302 * tmp3513;
    let tmp4871 = tmp2702 * tmp4704;
    let tmp4872 = tmp3210 * tmp3304;
    let tmp4873 = grad_phi_node_i[0] * tmp4856;
    let tmp4874 = tmp3302 * tmp3516;
    let tmp4875 = tmp2739 * tmp4710;
    let tmp4876 = tmp3202 * tmp3299;
    let tmp4877 = tmp1410 * tmp4698;
    let tmp4878 = tmp1397 * tmp4821
        + tmp1476 * tmp4791
        + tmp1477 * tmp4820
        + tmp1478 * tmp4798
        + tmp1482 * tmp4794
        + tmp3229 * tmp4818
        + tmp3242 * tmp4817
        + tmp3504 * tmp4828
        + tmp3521 * tmp4868
        + tmp4804 * tmp722
        + tmp4825 * tmp495
        + tmp4826 * tmp495
        + tmp4827 * tmp495
        + tmp4852 * tmp823
        + tmp4874
        + tmp4875
        - 1.0 * tmp4876
        - 1.0 * tmp4877;
    let tmp4879 = tmp1188 * tmp4832
        + tmp1274 * tmp4748
        + tmp1463 * tmp4835
        + tmp1464 * tmp4799
        + tmp1491 * tmp4793
        + tmp1492 * tmp4791
        + tmp3234 * tmp4816
        + tmp3240 * tmp4815
        + tmp3501 * tmp4841
        + tmp3507 * tmp4842
        + tmp476 * tmp4859
        + tmp476 * tmp4860
        + tmp476 * tmp4861
        + tmp4808 * tmp736
        - 1.0 * tmp4870
        - 1.0 * tmp4871
        + tmp4872
        + tmp4873;
    let tmp4880 = 2.0 * tmp3868
        + 2.0 * tmp3873
        + 2.0 * tmp3879
        + 2.0 * tmp3882
        + 2.0 * tmp3883
        + 2.0 * tmp3884
        + 2.0 * tmp3885
        + 2.0 * tmp3886
        + 2.0 * tmp3887
        + 2.0 * tmp3888
        + tmp3925
        - 1.0 * tmp3926;
    let tmp4881 = 2.0 * tmp3893
        + 2.0 * tmp3895
        + 2.0 * tmp3898
        + 2.0 * tmp3900
        + 2.0 * tmp3901
        + 2.0 * tmp3902
        + 2.0 * tmp3903
        + 2.0 * tmp3904
        + 2.0 * tmp3905
        + 2.0 * tmp3906
        + tmp3958
        - 1.0 * tmp3959;
    let tmp4882 = tmp3861 - 1.0 * tmp3862
        + 2.0 * tmp3928
        + 2.0 * tmp3930
        + 2.0 * tmp3932
        + 2.0 * tmp3934
        + 2.0 * tmp3935
        + 2.0 * tmp3936
        + 2.0 * tmp3937
        + 2.0 * tmp3938
        + 2.0 * tmp3939
        + 2.0 * tmp3940;
    let tmp4883 = 0.666666666666667 * tmp4760;
    let tmp4884 = tmp1819 * tmp3249;
    let tmp4885 = tmp219 * tmp4884;
    let tmp4886 = tmp1819 * tmp4481;
    let tmp4887 = tmp4883 + tmp4885 - 1.0 * tmp4886;
    let tmp4888 = 0.666666666666667 * tmp4757;
    let tmp4889 = tmp217 * tmp4884;
    let tmp4890 = tmp1819 * tmp4478;
    let tmp4891 = tmp4888 + tmp4889 - 1.0 * tmp4890;
    let tmp4892 = 0.666666666666667 * tmp4767;
    let tmp4893 = tmp1819 * tmp3433;
    let tmp4894 = tmp283 * tmp4884;
    let tmp4895 = -1.0 * tmp4892 + tmp4893 - 1.0 * tmp4894;
    let tmp4896 = 0.666666666666667 * tmp4762;
    let tmp4897 = tmp281 * tmp4884;
    let tmp4898 = tmp1819 * tmp3437;
    let tmp4899 = tmp4896 + tmp4897 - 1.0 * tmp4898;
    let tmp4900 = 0.666666666666667 * tmp4771;
    let tmp4901 = tmp1819 * tmp4504;
    let tmp4902 = tmp3250 * tmp4356;
    let tmp4903 = tmp4900 + tmp4901 + tmp4902;
    let tmp4904 = 0.666666666666667 * tmp4773;
    let tmp4905 = tmp256 * tmp4884;
    let tmp4906 = tmp263 * tmp4884;
    let tmp4907 = -1.0 * tmp4904 + tmp4905 + tmp4906;
    let tmp4908 = -1.0 * tmp4900 - 1.0 * tmp4901 - 1.0 * tmp4902;
    let tmp4909 = (4_f64 / 3.0) * tmp3302;
    let tmp4910 = 0.333333333333333 * tmp4757;
    let tmp4911 = 0.333333333333333 * tmp4760;
    let tmp4912 = 0.333333333333333 * tmp4771;
    let tmp4913 = 0.333333333333333 * tmp4773;
    let tmp4914 = -1.0 * tmp4913;
    let tmp4915 = 0.333333333333333 * tmp4762;
    let tmp4916 = -1.0 * tmp4915;
    let tmp4917 = 0.333333333333333 * tmp4767;
    let tmp4918 = -1.0 * tmp4917;
    let tmp4919 = -1.0 * tmp4910;
    let tmp4920 = (2_f64 / 3.0) * tmp3302;
    let tmp4921 = tmp1635 * tmp4753
        + tmp1640 * tmp4709
        + tmp1643 * tmp4749
        + tmp1649 * tmp4699
        + tmp1944 * tmp4793
        + tmp1945 * tmp4791
        + tmp1947 * tmp4704
        + tmp3243 * tmp4920
        + tmp3867 * tmp4728
        + tmp3872 * tmp4716
        + tmp3878 * tmp4713
        + tmp3881 * tmp4725
        - 1_f64 / 3.0 * tmp4719
        - 2_f64 / 3.0 * tmp4732
        + tmp646 * (tmp1929 + tmp4671 + tmp4916 + tmp4918)
        + tmp700 * (tmp1915 + tmp4058 + tmp4073 + tmp4912 + tmp4914)
        + tmp722 * (tmp1940 + tmp1969 + tmp4676 + tmp4911 + tmp4919)
        + tmp740 * (tmp4368 + tmp4634 + tmp4675 + tmp4910 + tmp4911);
    let tmp4922 = -1.0 * tmp4912;
    let tmp4923 = -1.0 * tmp4911;
    let tmp4924 = tmp1753 * tmp4754
        + tmp1755 * tmp4711
        + tmp1757 * tmp4751
        + tmp1759 * tmp4705
        + tmp1973 * tmp4794
        + tmp1974 * tmp4793
        - 1.0 * tmp1975 * tmp4708
        + tmp3927 * tmp4726
        + tmp3929 * tmp4729
        + tmp3931 * tmp4714
        + tmp3933 * tmp4717
        + (1_f64 / 3.0) * tmp4723
        - 2_f64 / 3.0 * tmp4733
        + (2_f64 / 3.0) * tmp4734
        + tmp671 * (tmp1965 + tmp2017 + tmp3983 + tmp4914 + tmp4922)
        + tmp706 * (tmp1942 + tmp1971 + tmp4669 + tmp4910 + tmp4923)
        + tmp729 * (tmp3988 + tmp4054 + tmp4387 + tmp4915 + tmp4918)
        + tmp744 * (tmp1959 + tmp4678 + tmp4915 + tmp4917);
    let tmp4925 = -1.0 * tmp4888 - 1.0 * tmp4889 + tmp4890;
    let tmp4926 = -1.0 * tmp4883 - 1.0 * tmp4885 + tmp4886;
    let tmp4927 = tmp4904 - 1.0 * tmp4905 - 1.0 * tmp4906;
    let tmp4928 = tmp4892 - 1.0 * tmp4893 + tmp4894;
    let tmp4929 = tmp1679 * tmp4707
        + tmp1681 * tmp4752
        + tmp1684 * tmp4746
        + tmp1687 * tmp4692
        + tmp2020 * tmp4791
        + tmp2021 * tmp4794
        - 1.0 * tmp2022 * tmp4710
        + tmp2023 * tmp4698
        - 1.0 * tmp3241 * tmp4920
        + tmp3892 * tmp4712
        + tmp3894 * tmp4715
        + tmp3897 * tmp4724
        + tmp3899 * tmp4727
        + (2_f64 / 3.0) * tmp4731
        + tmp601 * (tmp4402 + tmp4629 + tmp4674 + tmp4919 + tmp4923)
        + tmp693 * (tmp3994 + tmp3997 + tmp4405 + tmp4916 + tmp4917)
        + tmp716 * (tmp2011 + tmp4057 + tmp4072 + tmp4913 + tmp4922)
        + tmp736 * (tmp1963 + tmp2016 + tmp4070 + tmp4912 + tmp4913);
    let tmp4930 = -1.0 * tmp4896 - 1.0 * tmp4897 + tmp4898;
}
