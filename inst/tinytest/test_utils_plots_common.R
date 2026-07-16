# Test theme_msstats -------------------------------------------------------

theme_condition = theme_msstats("CONDITIONPLOT")
expect_true(inherits(theme_condition, "theme"))

theme_comparison = theme_msstats("COMPARISONPLOT")
expect_true(inherits(theme_comparison, "theme"))

theme_other = theme_msstats("PROFILEPLOT")
expect_true(inherits(theme_other, "theme"))

# Test getSelectedProteins --------------------------------------------------

all_proteins = c("ProtA", "ProtB", "ProtC")

# Test 1: character selection, all found
expect_equal(getSelectedProteins(c("ProtA", "ProtC"), all_proteins),
             c("ProtA", "ProtC"))

# Test 2: character selection with a missing protein errors
expect_error(getSelectedProteins(c("ProtA", "NotAProtein"), all_proteins))

# Test 3: numeric selection within range
expect_equal(getSelectedProteins(c(1, 3), all_proteins),
             c("ProtA", "ProtC"))

# Test 4: numeric selection out of range errors
expect_error(getSelectedProteins(c(1, 10), all_proteins))

# Test savePlot / getFileName ------------------------------------------------

tmp_dir = tempfile("msstats_saveplot_")
dir.create(tmp_dir)
old_wd = getwd()
setwd(tmp_dir)

# Test 5: savePlot with name_base = FALSE is a no-op, returns NULL, no file created
result_noop = savePlot(FALSE, "ProfilePlot", width = 800, height = 600)
expect_true(is.null(result_noop))
expect_equal(length(list.files(tmp_dir)), 0)

# Test 6: savePlot with a real name_base creates a pdf file
result_plot = savePlot("", "ProfilePlot", width = 800, height = 600)
dev.off()
expect_true(is.null(result_plot))
expect_true(file.exists(file.path(tmp_dir, "ProfilePlot.pdf")))

# Test 7: calling savePlot again with the same file_name appends a numeric suffix
result_plot2 = savePlot("", "QCPlot", width = 800, height = 600)
dev.off()
result_plot3 = savePlot("", "QCPlot", width = 800, height = 600)
dev.off()
expect_true(file.exists(file.path(tmp_dir, "QCPlot_2.pdf")))

setwd(old_wd)
unlink(tmp_dir, recursive = TRUE)

# Test .saveTable ------------------------------------------------------------

tmp_dir2 = tempfile("msstats_savetable_")
dir.create(tmp_dir2)
tbl = data.table::data.table(Protein = c("A", "B"), Value = c(1, 2))

# Test 8: .saveTable with name_base = FALSE is a no-op
res_noop_table = MSstats:::.saveTable(tbl, FALSE, "MyTable")
expect_true(is.null(res_noop_table))
expect_equal(length(list.files(tmp_dir2)), 0)

# Test 9: .saveTable with a real name_base writes a file
res_table = MSstats:::.saveTable(tbl, tmp_dir2, "MyTable")
expect_true(is.null(res_table))
expect_true(file.exists(file.path(tmp_dir2, "MyTable.pdf")))

unlink(tmp_dir2, recursive = TRUE)
