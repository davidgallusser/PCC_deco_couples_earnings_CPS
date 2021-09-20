###########################################
Counterfactual::InferenceTestingEval()
function (method, qte_cov_booti, reps, taus, robust, qte_cov_obs, 
          last, first, alpha, cons_test, treatment, decomposition) 
{
  nqs = length(taus)
  reps = nrow(qte_cov_booti)
  sel0 = ((sum(taus < first)) + 1):(sum(taus <= last))
  if (is.element(0.5, taus)) {
    median_sel = which(taus == 0.5)
  }
  else {
    print("Q_0.5 is not evaluated. This is the median of the specified rather than the Q_0.5")
    median_sel = ceiling(nqs/2)
  }
  sel_median = c((sum(taus < first) + 1):(median_sel - 1), 
                 (median_sel + 1):sum(taus <= last))
  selall = 1:length(sel_median)
  quantile_effect_bootstrap = qte_cov_booti[, 1:nqs]
  marginal_obs_bootstrap = qte_cov_booti[, (nqs + 1):(2 * nqs)]
  marginal_fitted_bootstrap = qte_cov_booti[, (2 * nqs + 1):(3 * 
                                                               nqs)]
  marginal_counterfactual_bootstrap = qte_cov_booti[, (3 * 
                                                         nqs + 1):(4 * nqs)]
  quantile_effect = qte_cov_obs[1:nqs]
  marginal_obs = qte_cov_obs[(nqs + 1):(2 * nqs)]
  marginal_fitted = qte_cov_obs[(2 * nqs + 1):(3 * nqs)]
  marginal_counterfactual = qte_cov_obs[(3 * nqs + 1):(4 * 
                                                         nqs)]
  if (treatment & decomposition) {
    composition_effect_bootstrap = qte_cov_booti[, (4 * nqs + 
                                                      1):(5 * nqs)]
    total_effect_bootstrap = qte_cov_booti[, (5 * nqs + 1):(6 * 
                                                              nqs)]
    marginal_obs_1_bootstrap = qte_cov_booti[, (6 * nqs + 
                                                  1):(7 * nqs)]
    marginal_fitted_1_bootstrap = qte_cov_booti[, (7 * nqs + 
                                                     1):(8 * nqs)]
    F_ms_bootstrap = qte_cov_booti[, -(1:8 * nqs)]
    composition_effect = qte_cov_obs[(4 * nqs + 1):(5 * nqs)]
    total_effect = qte_cov_obs[(5 * nqs + 1):(6 * nqs)]
    marginal_obs_1 = qte_cov_obs[(6 * nqs + 1):(7 * nqs)]
    marginal_fitted_1 = qte_cov_obs[(7 * nqs + 1):(8 * nqs)]
    F_ms_obs = qte_cov_obs[-(1:8 * nqs)]
    sel_ms = 1:length(F_ms_obs)
  }
  else if (treatment & !decomposition) {
    marginal_obs_1_bootstrap = qte_cov_booti[, (4 * nqs + 
                                                  1):(5 * nqs)]
    marginal_fitted_1_bootstrap = qte_cov_booti[, (5 * nqs + 
                                                     1):(6 * nqs)]
    F_ms_bootstrap = qte_cov_booti[, -(1:6 * nqs)]
    marginal_obs_1 = qte_cov_obs[(4 * nqs + 1):(5 * nqs)]
    marginal_fitted_1 = qte_cov_obs[(5 * nqs + 1):(6 * nqs)]
    F_ms_obs = qte_cov_obs[-(1:6 * nqs)]
    sel_ms = 1:length(F_ms_obs)
  }
  pe_boot = VarianceEval(quantile_effect_bootstrap, quantile_effect, 
                         reps, nqs, alpha, robust, sel0)
  pe_lb_point = pe_boot$qte_cov - pe_boot$se * qnorm(1 - alpha/2)
  pe_ub_point = pe_boot$qte_cov + pe_boot$se * qnorm(1 - alpha/2)
  if (!treatment & !decomposition) {
    resCE = cbind(pe_boot$qte_cov, pe_boot$se, pe_lb_point, 
                  pe_ub_point, pe_boot$lb, pe_boot$ub)
    colnames(resCE) = c("composition.effect", "se.ce", 
                        "lb.ce.point", "ub.ce.point", "lb.ce", 
                        "ub.ce")
  }
  else if (treatment & !decomposition) {
    resSE = cbind(pe_boot$qte_cov, pe_boot$se, pe_lb_point, 
                  pe_ub_point, pe_boot$lb, pe_boot$ub)
    colnames(resSE) = c("structure.effect", "se.se", 
                        "lb.se.point", "ub.se.point", "lb.se", 
                        "ub.se")
  }
  else if (treatment & decomposition) {
    resSE = cbind(pe_boot$qte_cov, pe_boot$se, pe_lb_point, 
                  pe_ub_point, pe_boot$lb, pe_boot$ub)
    colnames(resSE) = c("structure.effect", "se.se", 
                        "lb.se.point", "ub.se.point", "lb.se", 
                        "ub.se")
    ce_boot = VarianceEval(composition_effect_bootstrap, 
                           composition_effect, reps, nqs, alpha, robust, sel0)
    ce_lb_point = ce_boot$qte_cov - ce_boot$se * qnorm(1 - 
                                                         alpha/2)
    ce_ub_point = ce_boot$qte_cov + ce_boot$se * qnorm(1 - 
                                                         alpha/2)
    resCE = cbind(ce_boot$qte_cov, ce_boot$se, ce_lb_point, 
                  ce_ub_point, ce_boot$lb, ce_boot$ub)
    colnames(resCE) = c("composition.effect", "se.ce", 
                        "lb.ce.point", "ub.ce.point", "lb.ce", 
                        "ub.ce")
    te_boot = VarianceEval(total_effect_bootstrap, total_effect, 
                           reps, nqs, alpha, robust, sel0)
    te_lb_point = te_boot$qte_cov - te_boot$se * qnorm(1 - 
                                                         alpha/2)
    te_ub_point = te_boot$qte_cov + te_boot$se * qnorm(1 - 
                                                         alpha/2)
    resTE = cbind(te_boot$qte_cov, te_boot$se, te_lb_point, 
                  te_ub_point, te_boot$lb, te_boot$ub)
    colnames(resTE) = c("total.effect", "se.te", 
                        "lb.te.point", "ub.te.point", "lb.te", 
                        "ub.te")
  }
  sample_quantile_ref0 = VarianceEval(marginal_obs_bootstrap, 
                                      marginal_obs, reps, nqs, alpha, robust, sel0)
  colnames(sample_quantile_ref0) = c("ME.ref0", "se.ref0", 
                                     "lb.ref0", "ub.ref0")
  model_quantile_ref0 = VarianceEval(marginal_fitted_bootstrap, 
                                     marginal_fitted, reps, nqs, alpha, robust, sel0)
  colnames(model_quantile_ref0) = c("ME.fitted0", "se.fitted0", 
                                    "lb.fitted0", "ub.fitted0")
  model_quantile_counter = VarianceEval(marginal_counterfactual_bootstrap, 
                                        marginal_counterfactual, reps, nqs, alpha, robust, sel0)
  colnames(model_quantile_counter) = c("ME.counter", 
                                       "se.counter", "lb.counter", "ub.counter")
  if (treatment) {
    sample_quantile_ref1 = VarianceEval(marginal_obs_1_bootstrap, 
                                        marginal_obs_1, reps, nqs, alpha, robust, sel0)
    model_quantile_ref1 = VarianceEval(marginal_fitted_1_bootstrap, 
                                       marginal_fitted_1, reps, nqs, alpha, robust, sel0)
    colnames(sample_quantile_ref1) = c("ME.ref1", "se.ref1", 
                                       "lb.ref1", "ub.ref1")
    colnames(model_quantile_ref1) = c("ME.fitted1", 
                                      "se.fitted1", "lb.fitted1", "ub.fitted1")
  }
  if (method == "logit" | method == "lpm") {
    test_MS = c(NA, NA)
  }
  else {
    qte_cov_boot_ms = marginal_obs_bootstrap - marginal_fitted_bootstrap
    qte_cov_ms = marginal_obs - marginal_fitted
    test_MS = TestingEval((qte_cov_boot_ms - kronecker(matrix(qte_cov_ms, 
                                                              1, nqs), matrix(1, reps, 1)))[, sel0], qte_cov_ms[sel0], 
                          qte_cov_boot_ms, sel0, reps, robust)
  }
  test_0 = TestingEval((quantile_effect_bootstrap - kronecker(matrix(quantile_effect, 
                                                                     1, nqs), matrix(1, reps, 1)))[, sel0], quantile_effect[sel0], 
                       quantile_effect_bootstrap, sel0, reps, robust)
  constest1 = sort(unique(cons_test))
  constestNo0 = setdiff(constest1, 0)
  nc = length(constestNo0)
  if (nc > 0) {
    test_const = matrix(0, nc, 2)
    for (i in 1:nc) {
      test_const[i, ] = TestingEval((quantile_effect_bootstrap - 
                                       kronecker(matrix(quantile_effect, 1, nqs), matrix(1, 
                                                                                         reps, 1)))[, sel0], quantile_effect[sel0] - 
                                      constestNo0[i], quantile_effect_bootstrap, sel0, 
                                    reps, robust)
    }
  }
  else test_const = NULL
  qte_cov_boot_median = quantile_effect_bootstrap[, sel_median] - 
    quantile_effect_bootstrap[, median_sel]
  qte_cov_def_median = quantile_effect[sel_median] - quantile_effect[median_sel]
  test_median = TestingEval(qte_cov_boot_median - kronecker(matrix(qte_cov_def_median, 
                                                                   1, ncol(qte_cov_boot_median)), matrix(1, reps, 1)), qte_cov_def_median, 
                            qte_cov_boot_median, selall, reps, robust)
  qte_cov_boot_SD_ex = (quantile_effect_bootstrap - kronecker(matrix(quantile_effect, 
                                                                     1, nqs), matrix(1, reps, 1)))[, sel0]
  test_SD = TestingEval((qte_cov_boot_SD_ex * (qte_cov_boot_SD_ex <= 
                                                 0)), (quantile_effect * (quantile_effect <= 0))[sel0], 
                        quantile_effect_bootstrap, sel0, reps, robust)
  test_SDD = TestingEval((qte_cov_boot_SD_ex * (qte_cov_boot_SD_ex >= 
                                                  0)), (quantile_effect * (quantile_effect >= 0))[sel0], 
                         quantile_effect_bootstrap, sel0, reps, robust)
  if (!treatment & !decomposition) {
    testCE = rbind(test_MS, test_0, test_const, test_median, 
                   test_SD, test_SDD)
  }
  else if (treatment & !decomposition) {
    testSE = rbind(test_MS, test_0, test_const, test_median, 
                   test_SD, test_SDD)
  }
  else if (treatment & decomposition) {
    testSE = rbind(test_MS, test_0, test_const, test_median, 
                   test_SD, test_SDD)
    if (method == "logit" | method == "lpm") {
      test_MS_CE = c(NA, NA)
    }
    else {
      qte_cov_boot_ms_CE = marginal_obs_1_bootstrap - marginal_fitted_1_bootstrap
      qte_cov_ms_CE = marginal_obs_1 - marginal_fitted_1
      test_MS_CE = TestingEval((qte_cov_boot_ms_CE - kronecker(matrix(qte_cov_ms_CE, 
                                                                      1, length(qte_cov_ms_CE)), matrix(1, reps, 1)))[, 
                                                                                                                      sel0], qte_cov_ms_CE[sel0], qte_cov_boot_ms_CE, 
                               sel0, reps, robust)
    }
    test_0_CE = TestingEval((composition_effect_bootstrap - 
                               kronecker(matrix(composition_effect, 1, nqs), matrix(1, 
                                                                                    reps, 1)))[, sel0], composition_effect[sel0], 
                            composition_effect_bootstrap, sel0, reps, robust)
    if (nc > 0) {
      test_const_CE = matrix(0, nc, 2)
      for (i in 1:nc) {
        test_const_CE[i, ] = TestingEval((composition_effect_bootstrap - 
                                            kronecker(matrix(composition_effect, 1, nqs), 
                                                      matrix(1, reps, 1)))[, sel0], composition_effect[sel0] - 
                                           constestNo0[i], composition_effect_bootstrap, 
                                         sel0, reps, robust)
      }
    }
    else test_const_CE = NULL
    ce_cov_boot_median = composition_effect_bootstrap[, sel_median] - 
      composition_effect_bootstrap[, median_sel]
    ce_cov_def_median = composition_effect[sel_median] - 
      composition_effect[median_sel]
    test_median_CE = TestingEval(ce_cov_boot_median - kronecker(matrix(ce_cov_def_median, 
                                                                       1, ncol(ce_cov_boot_median)), matrix(1, reps, 1)), 
                                 ce_cov_def_median, ce_cov_boot_median, selall, reps, 
                                 robust)
    ce_cov_boot_SD_ex = (composition_effect_bootstrap - kronecker(matrix(composition_effect, 
                                                                         1, nqs), matrix(1, reps, 1)))[, sel0]
    test_SD_CE = TestingEval((ce_cov_boot_SD_ex * (ce_cov_boot_SD_ex <= 
                                                     0)), (composition_effect * (composition_effect <= 
                                                                                   0))[sel0], composition_effect_bootstrap, sel0, reps, 
                             robust)
    test_SDD_CE = TestingEval((ce_cov_boot_SD_ex * (ce_cov_boot_SD_ex >= 
                                                      0)), (composition_effect * (composition_effect >= 
                                                                                    0))[sel0], composition_effect_bootstrap, sel0, reps, 
                              robust)
    testCE = rbind(test_MS_CE, test_0_CE, test_const_CE, 
                   test_median_CE, test_SD_CE, test_SDD_CE)
    if (method == "logit" | method == "lpm") {
      test_MS_TE = c(NA, NA)
    }
    else {
      qte_cov_boot_ms_TE = F_ms_bootstrap
      qte_cov_ms_TE = F_ms_obs
      test_MS_TE = TestingEval((qte_cov_boot_ms_TE - kronecker(matrix(qte_cov_ms_TE, 
                                                                      1, length(qte_cov_ms_TE)), matrix(1, reps, 1))), 
                               qte_cov_ms_TE, qte_cov_boot_ms_TE, sel_ms, reps, 
                               robust)
    }
    test_0_TE = TestingEval((total_effect_bootstrap - kronecker(matrix(total_effect, 
                                                                       1, nqs), matrix(1, reps, 1)))[, sel0], total_effect[sel0], 
                            total_effect_bootstrap, sel0, reps, robust)
    if (nc > 0) {
      test_const_TE = matrix(0, nc, 2)
      for (i in 1:nc) {
        test_const_TE[i, ] = TestingEval((total_effect_bootstrap - 
                                            kronecker(matrix(total_effect, 1, nqs), matrix(1, 
                                                                                           reps, 1)))[, sel0], total_effect[sel0] - 
                                           constestNo0[i], total_effect_bootstrap, sel0, 
                                         reps, robust)
      }
    }
    else test_const_TE = NULL
    te_cov_boot_median = total_effect_bootstrap[, sel_median] - 
      total_effect_bootstrap[, median_sel]
    te_cov_def_median = total_effect[sel_median] - total_effect[median_sel]
    test_median_TE = TestingEval(te_cov_boot_median - kronecker(matrix(te_cov_def_median, 
                                                                       1, ncol(te_cov_boot_median)), matrix(1, reps, 1)), 
                                 te_cov_def_median, te_cov_boot_median, selall, reps, 
                                 robust)
    te_cov_boot_SD_ex = (total_effect_bootstrap - kronecker(matrix(total_effect, 
                                                                   1, nqs), matrix(1, reps, 1)))[, sel0]
    test_SD_TE = TestingEval((te_cov_boot_SD_ex * (te_cov_boot_SD_ex <= 
                                                     0)), (total_effect * (total_effect <= 0))[sel0], 
                             total_effect_bootstrap, sel0, reps, robust)
    test_SDD_TE = TestingEval((te_cov_boot_SD_ex * (te_cov_boot_SD_ex >= 
                                                      0)), (total_effect * (total_effect >= 0))[sel0], 
                              total_effect_bootstrap, sel0, reps, robust)
    testTE = rbind(test_MS_TE, test_0_TE, test_const_TE, 
                   test_median_TE, test_SD_TE, test_SDD_TE)
  }
  res_boot = NULL
  res_boot$sample_quantile_ref0 = sample_quantile_ref0
  res_boot$model_quantile_ref0 = model_quantile_ref0
  res_boot$model_quantile_counter = model_quantile_counter
  if (!treatment & !decomposition) {
    res_boot$resCE = resCE
    res_boot$testCE = testCE
  }
  else if (treatment & !decomposition) {
    res_boot$sample_quantile_ref1 = sample_quantile_ref1
    res_boot$model_quantile_ref1 = model_quantile_ref1
    res_boot$resSE = resSE
    res_boot$testSE = testSE
  }
  else if (treatment & decomposition) {
    res_boot$sample_quantile_ref1 = sample_quantile_ref1
    res_boot$model_quantile_ref1 = model_quantile_ref1
    res_boot$resSE = resSE
    res_boot$testSE = testSE
    res_boot$resCE = resCE
    res_boot$testCE = testCE
    res_boot$resTE = resTE
    res_boot$testTE = testTE
  }
  return(res_boot)
}

###########################################
Counterfactual::VarianceEval
function (qte_cov_boot, qte_cov, reps, nqs, alpha, robust, sel) 
{
  if (robust) {
    seuqf = apply(qte_cov_boot, 2, function(x) {
      (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, 
                                                  na.rm = TRUE))/1.34
    })
  }
  else {
    seuqf = sqrt(diag(var(qte_cov_boot)))
  }
  Vuqf = seuqf^2
  Kuqf = sqrt((qte_cov_boot - kronecker(matrix(qte_cov, 1, 
                                               nqs), matrix(1, reps, 1)))^2/kronecker(matrix(Vuqf, 1, 
                                                                                             nqs), matrix(1, reps, 1)))[, sel]
  Kuqfsel = Kuqf[, apply(Kuqf, 2, function(x) all(is.finite(x)))]
  Kmaxuqf = apply(Kuqfsel, 1, max)
  Kalpha = quantile(Kmaxuqf, 1 - alpha, na.rm = TRUE, names = FALSE)
  lb = qte_cov - seuqf * Kalpha
  ub = qte_cov + seuqf * Kalpha
  simbootres = data.frame(qte_cov = qte_cov, se = seuqf, lb = lb, 
                          ub = ub)
  return(simbootres)
}
<bytecode: 0x0000021524a78798>
  <environment: namespace:Counterfactual>

###########################################
Counterfactual::TestingEval
function (boot_test_numerator, obs_test_numerator, variable_for_variance, 
          sel, reps, robust) 
{
  if (robust) {
    seboot = apply(variable_for_variance, 2, function(x) {
      (quantile(x, 0.75, na.rm = TRUE) - quantile(x, 0.25, 
                                                  na.rm = TRUE))/1.34
    })
    Vuqf = (seboot)^2
  }
  else {
    Vuqf = diag(var(variable_for_variance))
  }
  Vuqf = Vuqf[sel] + 1e-09
  nqs = ncol(boot_test_numerator)
  Kuqf = sqrt(boot_test_numerator^2/kronecker(matrix(Vuqf, 
                                                     1, nqs), matrix(1, reps, 1)))
  Kmaxuqf = apply(Kuqf, 1, max)
  KSstat = max(sqrt(obs_test_numerator^2/Vuqf))
  Kuqf2 = Kuqf^2
  Kmeanuqf = rowMeans(Kuqf2)
  CMSstat = mean(obs_test_numerator^2/Vuqf)
  testboot = c(mean(Kmaxuqf > KSstat), mean(Kmeanuqf > CMSstat))
  return(testboot)
}
<bytecode: 0x0000021516a96930>
  <environment: namespace:Counterfactual>