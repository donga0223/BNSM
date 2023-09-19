
myschool <- c("ARTS02", "ARTS03", "PREP02", "TECH03")
print(getwd())

num_cores <- "8"
memory <- "8000"

for(i in myschool){
  filename <- paste("BNSM_application_", i, sep = "") ## add T,L i
  
  requestCmds <- "#!/bin/bash\n"
  requestCmds <- paste0(requestCmds, "#SBATCH --job-name=", filename, "# Name of the job\n",
                        "#SBATCH --ntasks=1 # Number of tasks \n",
                        "#SBATCH -c 1 # Number of Cores per Task\n",
                        "#SBATCH --nodes=1 # Requested number of nodes\n",
                        "#SBATCH --mem=16000 # Requested Memory; this is 16 GB\n",
                        "#SBATCH --output=", paste("application/BNSM_fit_logs/", filename, ".out", sep = ""), "## add T,L i
 # File for output \n",
                        "#SBATCH --error=", paste("application/BNSM_fit_logs/", filename, ".err", sep = ""), "## add T,L i
 # File for error log \n",
                        "#SBATCH --partition cpu-long # Partition\n",
                        "#SBATCH --time 100:00:00 # Job time limit\n")
  
  cat(requestCmds, file = paste("application/shellfiles/", filename, ".sh", sep = ""))
  cat("module load r-rocker-ml-verse/4.2.3+apptainer\n",  file = paste("application/shellfiles/", filename, ".sh", sep = ""), append = TRUE)  
  cat("shopt -s expand_aliases\n",  file = paste("application/shellfiles/", filename, ".sh", sep = ""), append = TRUE)  
  cat(paste(
    "Rscript BNSM_application.R", i, "10 Mannini ankle y_intensity parametricBoostCRF FALSE"),
    file = paste("application/shellfiles/", filename, ".sh", sep = ""), append = TRUE)
  
  bsubCmd <- paste("sbatch ", paste("application/shellfiles/", filename, ".sh", sep = ""))
  system(bsubCmd)
}










