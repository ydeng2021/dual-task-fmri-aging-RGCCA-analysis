## Behavioral Summary Measures

This folder includes processed behavioral summary files containing subject-level performance measures derived from the motor and cognitive task.

The behavioral summary data are provided in tabular format (one row per participant). All measures were computed from the original task log files and do not contain absolute time stamps or personal identifiers.

## Task-Specific Behavioral Summary Files

Behavioral summary measures are provided separately for the cognitive task and the motor task:

-   `cognitive_behavior_summary.csv`
-   `motor_behavior_summary.csv`

Both files share the same column structure and variable names to ensure consistency across task domains. Each file contains subject-level behavioral measures derived from the corresponding task.

The distinction between cognitive and motor performance is encoded at the file level rather than the column level.

### Variables (applies to both files)

-   **subjectID**: Anonymized participant identifier
-   **S_RT_mean**: Mean reaction time at single-task
-   **RT_median**: Median reaction time at single-task
-   **RT_Std**: Standard deviation of reaction time at single-task
-   **RT_Error**: Proportion of incorrect responses at single-task
-   **Dual_RT_mean**: Mean reaction time at dual-task
-   **Dual_RT_median**: Median reaction time at dual-task
-   **Dual_RT_Std**: Standard deviation of reaction time at dual-task
-   **Dual_RT_Error**: Proportion of incorrect responses at dual-task
-   **DTC_RT**: Dual-task cost of reaction time
-   **DTC_Std**: Dual-task cost of reaction time variability

## Dual-Task Cost (DTC)

Dual-task cost (DTC) was computed as:

DTC = ((Dual − Single) / Single) × 100

where *Dual* and *Single* refer to the corresponding reaction time measure (mean RT or RT standard deviation).

DTC was calculated for reaction time measures only. Error rates are reported descriptively and were not used to compute dual-task cost.

For reaction time measures, positive DTC values indicate performance slowing (i.e., increased reaction time or variability) under dual-task relative to single-task conditions.

All reaction time measures are reported in seconds. Dual-task cost measures are reported in percentage.
