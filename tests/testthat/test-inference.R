#-----------------------------------------------------------------------------
#
# These are tests for the inference procedure for computing ATT(g,t)'s and
# different aggregation methods.  These tests are meant to run fast ---
# they only confirm that standard errors are not changing across different
# iterations of the code.  More extensive/direct tests are available at
# did/tests/att_gt_inference_tests.Rmd.
#
#-----------------------------------------------------------------------------

library(DRDID)
library(BMisc)
library(ggplot2)
library(ggpubr)
library(withr)

with_did_version_2 <- function(fn) {
  # with_did_version_2 sets the current 'did' version to 2.0.0 and executes the
  # passed function. If necessary, it installs did version 2.0.0 in a temp
  # directory. with_did_version_2 returns the return value of the passed
  # function.
  tryCatch(detach("package:did"), error=function(e) "")

  old_did_path <- paste(tempdir(), "old_packages", sep="/")

  if (!dir.exists(old_did_path)) {
    dir.create(old_did_path)
    with_libpaths(new = old_did_path, devtools::install_version("did", version="2.0.0"))
  }

  library(did, lib.loc=old_did_path)

  package_version <- packageVersion("did")
  if (package_version != "2.0.0") {
    expect_true(
      package_version == "2.0.0",
      paste("wrong version of package. Expected: 2.0.0, Actual:", package_version)
    )
  }
  ret <- fn()
  detach("package:did")
  ret
}

with_local_did <- function(fn) {
  # with_local_did sets the current 'did' version to the 'did' package in the
  # current working directory (i.e. getwd()) and executes the passed function.
  # with_local_did returns the return value of the passed function.
  devtools::load_all(getwd())
  expect_true(
    packageVersion("did") != "2.0.0",
    "Expected 'did' package to not be version 2.0.0, but it was"
  )
  ret <- fn()
  detach("package:did")
  ret
}

expect_agg_gts_equal <- function(att_gt_new, att_gt_2.0) {
  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(att_gt_new$inffunc == att_gt_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(att_gt_new$se[1], att_gt_2.0$se[1], tol=.01)
}

expect_aggregations_equal <- expect_agg_gts_equal
  
test_that("inference with balanced panel data and aggregations", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)

  set.seed(1234)
  dr_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr")
  })
  reg_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg")
  })
  ipw_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw")
  })

  # aggregations
  dyn_2.0 <- with_did_version_2(function() {
    aggte(ipw_2.0, type="dynamic")
  })
  group_2.0 <- with_did_version_2(function() {
    aggte(reg_2.0, type="group")
  })
  cal_2.0 <- with_did_version_2(function() {
    aggte(dr_2.0, type="calendar")
  })
  
  set.seed(1234)
  dr_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr")
  })
  reg_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg")
  })
  ipw_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw")
  })

  # aggregations
  dyn_new <- with_local_did(function() {aggte(ipw_new, type="dynamic")})
  group_new <- with_local_did(function() {aggte(reg_new, type="group")})
  cal_new <- with_local_did(function() {aggte(dr_new, type="calendar")})

  expect_agg_gts_equal(dr_new, dr_2.0)
  expect_agg_gts_equal(reg_new, reg_2.0)
  expect_agg_gts_equal(ipw_new, ipw_2.0)

  expect_aggregations_equal(dyn_new, dyn_2.0)
  expect_aggregations_equal(group_new, group_2.0)
  expect_aggregations_equal(cal_new, cal_2.0)
})  

test_that("inference with clustering", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)

  set.seed(1234)
  dr_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr", clustervars="cluster")
  })
  reg_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg", clustervars="cluster")
  })
  ipw_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw", clustervars="cluster")
  })

  # aggregations
  dyn_2.0 <- with_did_version_2(function() {
    aggte(ipw_2.0, type="dynamic")
  })
  group_2.0 <- with_did_version_2(function() {
    aggte(reg_2.0, type="group")
  })
  cal_2.0 <- with_did_version_2(function() {
    aggte(dr_2.0, type="calendar")
  })

  set.seed(1234)
  dr_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr", clustervars="cluster")
  })
  reg_new <- with_local_did(
    function() {
      att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
             gname="G", est_method="reg", clustervars="cluster")
  })
  ipw_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw", clustervars="cluster")
  })
  
  # aggregations
  dyn_new <- with_local_did(function() {
    aggte(ipw_new, type="dynamic")
  })
  group_new <- with_local_did(function() {
    aggte(reg_new, type="group")
  })
  cal_new <- with_local_did(function() {
    aggte(dr_new, type="calendar")
  })

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
})  

test_that("same inference with unbalanced panel and panel data", {
  skip(message="While fixing the inference unit tests, I noticed this test is failing. Should revisit.")
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  
  res_factor <- att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id", 
                       gname="G", est_method="dr", clustervars="cluster")

  #-----------------------------------------------------------------------------
  # clustered standard errors with unbalanced panel
  res_ub <- att_gt(yname="Y",
                   tname="period",
                   idname="id",
                   gname="G",
                   xformla=~X,
                   data=data,
                   panel=TRUE,
                   allow_unbalanced_panel=TRUE,
                   clustervars="cluster")

  # check that influence function is the same if we use unbalanced panel approach
  # vs. balanced panel approach
  # Note: differences are showing up in never-treated units (that are towards end of sample)
  expect_true(all(res_factor$inffunc == res_ub$inffunc))
})

test_that("inference with repeated cross sections", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp, panel=FALSE)

  set.seed(1234)
  dr_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr", panel=FALSE)
  })
  reg_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg", panel=FALSE)
  })
  ipw_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw", panel=FALSE)
  })

  # aggregations
  dyn_2.0 <- with_did_version_2(function() {aggte(ipw_2.0, type="dynamic")})
  group_2.0 <- with_did_version_2(function() {aggte(reg_2.0, type="group")})
  cal_2.0 <- with_did_version_2(function() {aggte(dr_2.0, type="calendar")})

  set.seed(1234)
  #dr 
  dr_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr", panel=FALSE)
  })
  # reg
  reg_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg", panel=FALSE)
  })
  # reg
  ipw_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw", panel=FALSE)
  })

  # aggregations
  dyn_new <- with_local_did(function() {aggte(ipw_new, type="dynamic")})
  group_new <- with_local_did(function() {aggte(reg_new, type="group")})
  cal_new <- with_local_did(function() {aggte(dr_new, type="calendar")})

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
})

test_that("inference with repeated cross sections and clustering", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp, panel=FALSE)

  set.seed(1234)
  dr_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr", clustervars="cluster", panel=FALSE)
  })
  reg_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg", clustervars="cluster", panel=FALSE)
  })
  ipw_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw", clustervars="cluster", panel=FALSE)
  })

  # aggregations
  dyn_2.0 <- with_did_version_2(function() {aggte(ipw_2.0, type="dynamic")})
  group_2.0 <- with_did_version_2(function() {aggte(reg_2.0, type="group")})
  cal_2.0 <- with_did_version_2(function() {aggte(dr_2.0, type="calendar")})

  set.seed(1234)
  #dr 
  dr_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr", clustervars="cluster", panel=FALSE)
  })
  # reg
  reg_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg", clustervars="cluster", panel=FALSE)
  })
  # reg
  ipw_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw", clustervars="cluster", panel=FALSE)
  })

  # aggregations
  dyn_new <- with_local_did(function() {aggte(ipw_new, type="dynamic")})
  group_new <- with_local_did(function() {aggte(reg_new, type="group")})
  cal_new <- with_local_did(function() {aggte(dr_new, type="calendar")})

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
})  

test_that("inference with unbalanced panel", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  data <- data[-3,]

  set.seed(1234)
  dr_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr", panel=TRUE, allow_unbalanced_panel=TRUE)
  })
  reg_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg", panel=TRUE, allow_unbalanced_panel=TRUE)
  })
  ipw_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw", panel=TRUE, allow_unbalanced_panel=TRUE)
  })

  # aggregations
  dyn_2.0 <- with_did_version_2(function() {
    aggte(ipw_2.0, type="dynamic")
  })
  group_2.0 <- with_did_version_2(function() {
    aggte(reg_2.0, type="group")
  })
  cal_2.0 <- with_did_version_2(function() {
    aggte(dr_2.0, type="calendar")
  })
  
  set.seed(1234)
  dr_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr", panel=TRUE, allow_unbalanced_panel=TRUE)
  })
  reg_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg", panel=TRUE, allow_unbalanced_panel=TRUE)
  })
  ipw_new <- with_local_did(function() {
  att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
         gname="G", est_method="ipw", panel=TRUE, allow_unbalanced_panel=TRUE)
  })

  # aggregations
  dyn_new <- with_local_did(function() {aggte(ipw_new, type="dynamic")})
  group_new <- with_local_did(function() {aggte(reg_new, type="group")})
  cal_new <- with_local_did(function() {aggte(dr_new, type="calendar")})

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
})

test_that("inference with unbalanced panel and clustering", {
  sp <- reset.sim()
  data <- build_sim_dataset(sp)
  data <- data[-3,]

  set.seed(1234)
  dr_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr", clustervars="cluster",
           allow_unbalanced_panel=TRUE)
  })
  reg_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg", clustervars="cluster",
           allow_unbalanced_panel=TRUE)
  })
  ipw_2.0 <- with_did_version_2(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw", clustervars="cluster",
           allow_unbalanced_panel=TRUE)
  })

  # aggregations
  dyn_2.0 <- with_did_version_2(function() {
    aggte(ipw_2.0, type="dynamic")
  })
  group_2.0 <- with_did_version_2(function() {
    aggte(reg_2.0, type="group")
  })
  cal_2.0 <- with_did_version_2(function() {
    aggte(dr_2.0, type="calendar")
  })

  set.seed(1234)
  dr_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="dr", clustervars="cluster", allow_unbalanced_panel=TRUE)
  })
  reg_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="reg", clustervars="cluster", allow_unbalanced_panel=TRUE)
  })
  ipw_new <- with_local_did(function() {
    att_gt(yname="Y", xformla=~X, data=data, tname="period", idname="id",
           gname="G", est_method="ipw", clustervars="cluster", allow_unbalanced_panel=TRUE)
  })

  # aggregations
  dyn_new <- with_local_did(function() {aggte(ipw_new, type="dynamic")})
  group_new <- with_local_did(function() {aggte(reg_new, type="group")})
  cal_new <- with_local_did(function() {aggte(dr_new, type="calendar")})

  # checks for ATT(g,t)'s
  # check that the influence function is the same
  expect_true(all(dr_new$inffunc == dr_2.0$inffunc))
  expect_true(all(reg_new$inffunc == reg_2.0$inffunc))
  expect_true(all(ipw_new$inffunc == ipw_2.0$inffunc))

  # standard errors should be close
  # not totally sure, but I think slight differences are expected
  # perhaps from implementing the multiplier on the C++ side
  # in newer versions of the code
  expect_equal(dr_2.0$se[1], dr_new$se[1], tol=.01)
  expect_equal(reg_2.0$se[1], reg_new$se[1], tol=.01)
  expect_equal(ipw_2.0$se[1], ipw_new$se[1], tol=.01)

  # checks for aggregations
  expect_true(all(dyn_2.0$inffunc == dyn_new$inffunc))
  expect_true(all(group_2.0$inffunc == group_new$inffunc))
  expect_true(all(cal_2.0$inffunc == cal_new$inffunc))

  # standard errors for aggregations
  expect_equal(dyn_2.0$se[1], dyn_new$se[1], tol=.01)
  expect_equal(group_2.0$se[1], group_new$se[1], tol=.01)
  expect_equal(cal_2.0$se[1], cal_new$se[1], tol=.01)
})  
