# Test .checkTechReplicate ------------------------------------------------

# Test 1: No technical replicates
input = data.table::data.table(GROUP = c("A", "A", "B", "B"),
                               RUN = c(1, 2, 3, 4),
                               SUBJECT = c(1, 2, 3, 4))
expect_false(MSstats:::.checkTechReplicate(input))

# Test 2: Repeated Measures - No technical replicates
input = data.table::data.table(GROUP = c("A", "A", "B", "B"),
                               RUN = c(1, 2, 3, 4),
                               SUBJECT = c(1, 2, 1, 2))
expect_false(MSstats:::.checkTechReplicate(input))

# Test 3: Technical replicates for one subject
input = data.table::data.table(GROUP = c("A", "A", "B", "B"),
                               RUN = c(1, 2, 3, 4),
                               SUBJECT = c(1, 1, 2, 3))
expect_true(MSstats:::.checkTechReplicate(input))

# Test 4: Technical replicates for all subjects
input = data.table::data.table(GROUP = c("A", "A", "B", "B"),
                               RUN = c(1, 2, 3, 4),
                               SUBJECT = c(1, 1, 2, 2))
expect_true(MSstats:::.checkTechReplicate(input))

# Test 5: Repeated Measures - Technical replicates
input = data.table::data.table(GROUP = c("A", "A", "B", "B", "B"),
                               RUN = c(1, 2, 3, 4, 5),
                               SUBJECT = c(1, 2, 1, 2, 2))
expect_true(MSstats:::.checkTechReplicate(input))

