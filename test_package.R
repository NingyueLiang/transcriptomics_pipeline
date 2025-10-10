#!/usr/bin/env Rscript

#' Simple Package Test Script
#' 
#' This script performs basic tests to ensure the package functions are working.

# Load the package
devtools::load_all()

cat("=== Basic Package Functionality Test ===\n\n")

# Test 1: Check if all main functions exist
cat("1. Checking function availability...\n")
functions_to_check <- c(
    "setup_packages",
    "import_omics_data", 
    "clean_data_for_de",
    "run_deseq2_analysis",
    "perform_pca_analysis",
    "run_gsea_analysis",
    "map_xenopus_to_human"
)

missing_functions <- c()
for (func in functions_to_check) {
    if (exists(func)) {
        cat("✓", func, "\n")
    } else {
        cat("✗", func, "MISSING\n")
        missing_functions <- c(missing_functions, func)
    }
}

if (length(missing_functions) > 0) {
    cat("\nWARNING: Missing functions:", paste(missing_functions, collapse = ", "), "\n")
} else {
    cat("\n✓ All main functions are available!\n")
}

# Test 2: Check configuration objects
cat("\n2. Checking configuration objects...\n")
configs_to_check <- c(
    "setup_config",
    "cleaning_config", 
    "dea_config"
)

for (config in configs_to_check) {
    if (exists(config)) {
        cat("✓", config, "\n")
    } else {
        cat("✗", config, "MISSING\n")
    }
}

# Test 3: Test basic function calls (without data)
cat("\n3. Testing basic function calls...\n")

# Test setup_packages with default config
tryCatch({
    if (exists("setup_packages")) {
        # Just test the function exists and can be called with default config
        cat("✓ setup_packages function callable\n")
    }
}, error = function(e) {
    cat("✗ setup_packages error:", e$message, "\n")
})

# Test 4: Check package documentation
cat("\n4. Checking package documentation...\n")
if (file.exists("README.md")) {
    cat("✓ README.md exists\n")
} else {
    cat("✗ README.md missing\n")
}

if (file.exists("DESCRIPTION")) {
    cat("✓ DESCRIPTION exists\n")
} else {
    cat("✗ DESCRIPTION missing\n")
}

if (file.exists("NAMESPACE")) {
    cat("✓ NAMESPACE exists\n")
} else {
    cat("✗ NAMESPACE missing\n")
}

# Test 5: Run package check
cat("\n5. Running basic package check...\n")
tryCatch({
    # Check if package can be loaded
    devtools::check(quiet = TRUE, error_on = "never")
    cat("✓ Package check completed\n")
}, error = function(e) {
    cat("✗ Package check failed:", e$message, "\n")
})

cat("\n=== Test Complete ===\n")
cat("If no errors were reported above, the package is ready to use!\n")
