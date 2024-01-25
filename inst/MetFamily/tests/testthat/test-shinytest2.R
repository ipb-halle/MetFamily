library(shinytest2)

test_that("{shinytest2} recording: MetFamilyLoadExample", {
  app <- AppDriver$new(variant = platform_variant(), name = "MetFamilyLoadExample",
      height = 1113, width = 1549)
  #app$expect_screenshot()
  app$set_inputs(fileInputSelection = "Example data")
  app$click("loadExampleData")
  app$expect_values(screenshot_args = FALSE)
})
