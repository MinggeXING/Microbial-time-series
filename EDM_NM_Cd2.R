# EDM_NM_Cd2.R - EDM calculation script

# 0. Initialization ----
workspace_file <- "NM_Cd2_workspace.RData"
progress_file <- "NM_Cd2_progress.rds"
log_file <- "edm_nm_cd2.log" 

# Ensure log directory exists
if(!dir.exists(dirname(log_file))) dir.create(dirname(log_file), recursive = TRUE)

# Initialize log
cat(paste(Sys.time(), "===== Script started =====\n"), file = log_file, append = TRUE)
cat(paste(Sys.time(), "Current working directory:", getwd(), "\n"), file = log_file, append = TRUE)
cat(paste(Sys.time(), "R version:", R.version.string, "\n"))

# 1. Load workspace ----
load("NM_Cd2_workspace.RData")


# Custom log function
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste(timestamp, msg, sep = " - ")
  cat(full_msg, file = log_file, append = TRUE)
  cat(full_msg)  # Also output to console
}

# 2. Restore workspace ----
if(file.exists(workspace_file)) {
  load(workspace_file, .GlobalEnv)
  log_message("Workspace restored\n")
} else {
  log_message("Initializing new workspace\n")
}

# 3. Load required packages ----
required_packages <- c("dplyr", "ggplot2", "agricolae", "ggpubr", "vegan", "rEDM", 
                       "igraph", "quantreg", "doParallel", "foreach", 
                       "Kendall", "MASS", "glmnet", "tidyr", "dplyr")

for(pkg in required_packages){
  if(!require(pkg, character.only = TRUE, quietly = TRUE)){
    install.packages(pkg)
    library(pkg, character.only = TRUE)
    log_message(paste("Loaded package:", pkg, "\n"))
  }
}


# 4. Set control flags ----
# Define whether to execute compilation and MDR analysis
if(!exists("CompileCV")) CompileCV <- TRUE
if(!exists("do.MDR")) do.MDR <- TRUE
if(!exists("ptype")) ptype <- 'aenet'  # Set default parameter type

# 5. Check and restore progress ----
if(file.exists(progress_file)) {
  completed_alfs <- readRDS(progress_file)
  log_message(paste("Current progress:", paste(completed_alfs, collapse=","), "\n"))
} else {
  completed_alfs <- integer(0)
  log_message("Initializing progress tracking\n")
}

# Ensure key variables exist
if(!exists("completed_alfs")) completed_alfs <- integer(0)

# 6. Custom function current_alf: safely run task ----
safe_run_alf <- function(current_alf) {
  tryCatch({
    # Get current parameter range
    alpha.s <- alpha.so[afsp[current_alf,1]:afsp[current_alf,2]]
    
    log_message(paste("Parameter range: alpha =", 
                      round(min(alpha.s), 3), "-", 
                      round(max(alpha.s), 3), "\n"))
    
    # Execute calculation
    cv.ind <- cv.MDR.demo(ds, ds_tp1, dmatrix.list = dmatrix.train.mvx, 
                          parall = TRUE, ncore = 4, keep_intra = TRUE, 
                          alpha.seq = alpha.s)
    
    # Save result
    out_file <- file.path('out', 
                          paste0(da.name, '_nin', nin, '_cvunit', cv.unit,
                                 '_alph', formatC(alpha.s[1]*100, width=3, flag="0"),
                                 '_cvout_Nmvx_Rallx.csv'))
    write.csv(cv.ind, out_file, row.names = FALSE)
    log_message(paste("Result saved to:", out_file, "\n"))
    
    # Mark task complete
    return(TRUE)
  }, error = function(e) {
    log_message(paste("Error processing alf =", current_alf, ":", e$message, "\n"))
    return(FALSE)
  })
}

# 7. Progress management function alf ----
update_progress <- function(alf) {
  if(!exists("completed_alfs")) {
    completed_alfs <<- integer(0)
  }
  completed_alfs <<- unique(c(completed_alfs, alf))
  saveRDS(completed_alfs, progress_file)
  
  # Compact save: only save key objects
  save_list <- c("completed_alfs", "dmatrix.train.mvx", "ds", "ds_tp1", 
                 "cv.unit", "alpha.so", "sub.da", "afsp", "da.name", "nin")
  
  # Add optional objects (if they exist)
  optional_objs <- c("esele", "ccm.sig", "ccm.rho", "Ed")
  for(obj in optional_objs) {
    if(exists(obj)) save_list <- c(save_list, obj)
  }
  
  save(list = save_list, file = workspace_file)
  
  log_message(paste("Progress updated: alf", alf, "completed\n"))
}

# Set default values (if not exist)
if(!exists("cv.unit")) cv.unit <- 0.1
if(!exists("alpha.so")) alpha.so <- seq(0, 1, cv.unit)
if(!exists("sub.da")) sub.da <- 6
if(!exists("afsp")) afsp <- eqsplit(1:length(alpha.so), sub.da)

# 8. Main calculation loop ----
log_message(paste("Starting CV tasks, total segments:", sub.da, "current progress:", paste(completed_alfs, collapse=","), "\n"))

for(current_alf in 1:sub.da) {
  # Check if current task is already completed
  if(current_alf %in% completed_alfs) {
    log_message(paste("Skipping already completed alf =", current_alf, "\n"))
    next
  }
  
  log_message(paste("\n===== Starting processing alf =", current_alf, "=====\n"))
  
  # Execute and save safely
  success <- safe_run_alf(current_alf)
  
  if(success) {
    # Update progress
    update_progress(current_alf)
    log_message(paste("===== alf =", current_alf, "processing completed =====\n\n"))
  } else {
    log_message(paste("===== alf =", current_alf, "processing failed =====\n\n"))
  }
}

save.image(file = "NM_Cd2_workspace.RData")

# 9. Compile CV results ----
if (CompileCV && length(completed_alfs) == sub.da) {
  log_message("Starting compilation of CV results\n")
  
  # Initialize result data frame
  cv.ind <- NULL
  
  # Collect all segment results
  for(alf in 1:sub.da) {
    # Build file name
    alpha_start_val <- alpha.so[afsp[alf, 1]] * 100
    alpha_start_val_formatted <- formatC(alpha_start_val, width=3, flag="0")
    
    file_name <- paste0(da.name, '_nin', nin, '_cvunit', cv.unit,
                        '_alph', alpha_start_val_formatted,
                        '_cvout_Nmvx_Rallx.csv')
    file_path <- file.path('out', file_name)
    
    if (file.exists(file_path)) {
      cv_part <- read.csv(file_path, header = TRUE)
      cv.ind <- rbind(cv.ind, cv_part)
      log_message(paste("Read file:", file_name, "\n"))
    } else {
      log_message(paste("Warning: file does not exist -", file_name, "\n"))
    }
  }
  
  if (!is.null(cv.ind)) {
    # Select optimal parameters
    paracv.demo <- secv.demo(cv.ind)
    optimal_file <- paste0(da.name, '_nin', nin, '_cvunit', cv.unit,
                           '_OptimalCV_Nmvx_Rallx_NEW.csv')
    write.csv(paracv.demo, file.path('out', optimal_file), row.names = FALSE)
    log_message(paste("CV result compilation complete, optimal parameters saved to:", optimal_file, "\n"))
    
    # Save paracv.demo to environment
    assign("paracv.demo", paracv.demo, envir = .GlobalEnv)
  } else {
    log_message("Error: No CV result files read, compilation failed\n")
  }
}

save.image(file = "NM_Cd2_workspace.RData")

# 10. Execute MDR analysis ----
if (do.MDR && exists("paracv.demo")) {
  log_message("Starting MDR analysis\n")
  
  # Check if necessary objects exist
  required_objs <- c("ds", "ds_tp1", "dmatrix.train.mvx", "ds.test", "dst_tp1", "dmatrix.test.mvx")
  missing_objs <- required_objs[!sapply(required_objs, exists)]
  
  if(length(missing_objs) == 0) {
    # Call MDRsmap.demo function
    smap.demo <- MDRsmap.demo(paracv = paracv.demo, 
                              ptype = ptype, 
                              keep_intra = TRUE, 
                              out.sample = out.sample,
                              ds = ds, 
                              ds_tp1 = ds_tp1, 
                              ds.test = ds.test, 
                              dst_tp1 = dst_tp1,
                              dmatrix.list = dmatrix.train.mvx,
                              dmatrix.test.list = dmatrix.test.mvx)
    
    # Save results
    if(save) {
      # Save prediction skill
      if(!is.null(smap.demo[['nr.out']])) {
        nrout_file <- paste0(da.name, '_nin', nin, '_cvunit', cv.unit, '_', ptype, '_nrout_Nmvx_Rallx_NEW.csv')
        write.csv(smap.demo[['nr.out']], file.path('out', nrout_file), row.names = FALSE)
        log_message(paste("Prediction skill saved to:", nrout_file, "\n"))
      }
      
      # Save Jacobian matrix
      if(!is.null(smap.demo[['jcof']])) {
        jcof_file <- paste0(da.name, '_nin', nin, '_cvunit', cv.unit, '_', ptype, '_jcof_Nmvx_Rallx_demo_NEW.csv')
        write.csv(smap.demo[['jcof']], file.path('out', jcof_file), row.names = FALSE)
        log_message(paste("Jacobian matrix saved to:", jcof_file, "\n"))
      }
    }
    log_message("MDR analysis completed\n")
  } else {
    log_message(paste("Error: Missing required objects -", paste(missing_objs, collapse=", "), "\n"))
  }
} else if(do.MDR && !exists("paracv.demo")) {
  log_message("Warning: Cannot execute MDR analysis, missing paracv.demo object\n")
}

# 11. Save and exit ----
save.image(file = "NM_Cd2_workspace.RData")
log_message("===== Script execution finished =====\n")