
# Benchmarking Setup Documentation

## 1. benchmark.R

This R script contains the benchmarking logic for the MSstats-based proteomics data analysis. 
It executes specific workflows to analyze proteomics datasets, compute metrics (e.g., False Discovery Rate (FDR)), 
and ensure the validity of the MSstats library's updates. The script reads input data, processes it, and outputs benchmarking metrics.

## 2. config.slurm

This SLURM configuration file automates the execution of the R benchmarking script on an HPC system. 
It includes directives for resource allocation, job naming, and runtime limits. 
It ensures efficient utilization of HPC resources for running computationally intensive workflows.

## 3. benchmark.yml

This YAML configuration file is part of a GitHub Actions pipeline. 
It defines workflows for automating benchmarking tasks. 
The file contains instructions for setting up the R environment, pulling the required repositories, and executing the benchmarks.


## Setup Instructions for New Users

### 1. Prerequisites

Ensure you have access to the following:
- An HPC account with SLURM job scheduler.
- Required R dependencies installed (check `benchmark.R` for library imports).
- A GitHub account with access to the repository containing these files.

### 2. Setup HPC Environment

1. Transfer the `benchmark.R` and `config.slurm` files to your HPC environment.
2. Modify the `config.slurm` file to include your job-specific parameters (e.g., email, account name, partitions).
3. Submit the job using `sbatch config.slurm`.

### 3. Setup GitHub Actions

1. Place the `benchmark.yml` file in the `.github/workflows/` directory of your repository.
2. Configure the `benchmark.yml` file with appropriate paths and repository settings.
3. Push the changes to your repository to trigger the pipeline.

### 4. Verify Execution

1. Check the SLURM job output logs for successful execution of the `benchmark.R` script.
2. Validate that the benchmarking metrics are generated correctly.
3. Monitor the GitHub Actions logs to ensure the workflows execute without errors.

---

# SSH Access Setup for a New User

## Steps to Set Up SSH Access for a New User

### 1. Generate SSH Key Pair
On the new user's local machine, generate an SSH key pair (if not already created):
```bash
ssh-keygen -t rsa -b 4096 -C "new_user_email@example.com"

Example : current user email configured is : raina.ans@login-00.discovery.neu.edu
You can check this by navigating to shell through Discovery Cluster Dashboard > Clusters > Discovery Shell Access

```
- When prompted, specify a file to save the key or press `Enter` to use the default (`~/.ssh/id_rsa`).
- Set a passphrase for additional security.

### 2. Copy the Public Key to the Remote Server

#### Manually Copy the Key

1. SSH into the remote server using an existing account with sufficient privileges:
   ```bash
   ssh existing_user@remote_server (raina.ans@login-00.discovery.neu.edu)
   ```

2. Switch to the new user account (if created):
   ```bash
   sudo su - new_user
   ```

3. Append the public key to the `authorized_keys` file:
   ```bash
   mkdir -p ~/.ssh
   echo "paste_the_public_key_here" >> ~/.ssh/authorized_keys
   chmod 600 ~/.ssh/authorized_keys
   chmod 700 ~/.ssh
   ```
---

### 3. Verify the New User's SSH Access
From the new user's local machine, attempt to log in to the remote server:
```bash
ssh new_user@remote_server
```

If successful, the new user should be logged into the remote server.


### Notes
- Ensure proper permissions are set on the `.ssh` directory and `authorized_keys` file.
- Regularly review and update SSH configurations for security.
---
