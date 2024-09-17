# Nextflow Documentation: `./resource_handling.config`

This Nextflow configuration file implements a memory escalation strategy for processes labeled `mem_k2r_escalate`, alongside a customizable retry strategy for handling errors. Below is an explanation of the key sections of the configuration:

`params` Block:

- `max_attempts`: Defines the maximum number of retries for a task before it is terminated.

- `default_error_strategy`: Specifies the default error handling strategy, either "retry" or "terminate".

- **Memory Escalation Parameters**:
  - `mem_k2r_b0_offset`, `mem_k2r_b0`, and `mem_k2r_b1` define the coefficients for the linear regression model used to predict memory requirements.
  - `mem_k2r_f1` is a multiplication factor used during the second retry attempt.
  - `mem_k2r_a2` is an additional memory allocation factor used for retries beyond the second attempt.

These parameters allow fine control over how memory requirements are scaled depending on task retry attempts.

`functions`:

- `bytesToSize()`: Converts a given number of bytes into the requested memory unit (GB, MB, etc.) for easier memory management.
- `linear_regression_fit()`: Implements a linear regression formula to predict memory requirements based on input size and coefficients (`b0` and `b1`).

- `retry_strategy()`: A function designed to handle retry logic. It analyzes task exit codes and retries based on the nature of the error. Non-scalable errors (e.g., `SIGKILL`, `SIGTERM`) terminate the workflow, while scalable errors (e.g., `SIGUSR2`, `SIGINT`) trigger a retry with more memory. If no exit code is available or if it's the first attempt, it also triggers a retry.
- `k2r_escalate_memory_strategy()`: Implements the memory escalation strategy based on task retry attempts and input data size. The initial memory allocation is predicted using a linear regression model, and further attempts apply increasing memory based on factors (`f1`, `a2`) and task attempt number.

`process` Block:

- **Error Strategy**: The default error strategy is defined by the retry_strategy() function. If a task fails, it retries up to params.max_attempts based on the task exit status and error code.

- **Memory Escalation (withLabel)**: For tasks labeled `mem_k2r_escalate`, memory is dynamically calculated using the `k2r_escalate_memory_strategy()` function. The function determines the required memory based on the size of the input data and adjusts it based on retry attempts, allowing for controlled memory escalation to handle larger tasks or errors.
  - This configuration ensures efficient resource management and robust error handling, minimizing pipeline interruptions due to memory or execution limits while making dynamic adjustments for failed tasks.

> **DEV NOTE**: this function was initially implemented as a response for `kraken2ref` requesting too much memory for big fqs files, which stops the pipeline until the user figures out how much memory it needs. This function gives a good prediction based on the memory usage vs fastq input size profile of one of our runs and it was supposed to be a temporary fix for it. Since then, improvements on `kraken2ref` and the spliting fastq files implemented on the pipeline made this unnecessary. It is currently only used on `run_k2r_sort_reads` step of  or, at least, it should be updated for that step process.