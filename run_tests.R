library('RUnit')

#source('sample.R')

test.suite <- defineTestSuite("Main",
                              dirs = file.path("tests"),
                              testFileRegexp = '^\\d+\\.R')

test.result <- runTestSuite(test.suite)

printTextProtocol(test.result)