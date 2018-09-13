import os
import pandas as pd
import numpy as np
import itertools

from numba import njit, guvectorize


class OutputTable(object):
    def __init__(self, out_dir, sampling_method, imputation_method, snap, quantiles, stats_intervals):
        # Initialize path
        self.dir = out_dir
        self.base = "bf_results_{}_{}_{}.csv".format(sampling_method, imputation_method, snap)
        self.path = os.path.join(self.dir, self.base)
        self.quantiles = quantiles
        self.intervals = stats_intervals

        # Initialize header
        self.header = ["Site", "Year", "Sample Interval", "Realization"] + \
                      ["{}-day Bias Factor".format(interval) for interval in stats_intervals] + \
                      ["{}P {}".format(quantile, var) for var in ("Slope", "RMSE") for quantile in quantiles]

        # Initialize file
        self.initialize_table()

    def initialize_table(self):
        # Create output directory if it does not exist
        if not os.path.exists(self.dir):
            os.mkdir(self.dir)

        # Create file and write header
        with open(self.path, 'w') as f:
            f.write(",".join(self.header) + "\n")

    def update_table(self, site, year, sample_interval, realization, bias_factors, slopes, residuals):
        """ Insert run results into table """
        # Arrange data into row
        out_data = [site, year, sample_interval, realization] + list(bias_factors) + list(slopes) + list(residuals)

        # Update file
        with open(self.path, 'a') as f:
            f.write(",".join(map(str, out_data)) + "\n")


def check_inputs(chemograph_dir, sampling_intervals, start_date):
    # Ensure quarterly sampling
    assert all(1 <= sampling_interval <= 90 for sampling_interval in sampling_intervals), \
        "Sampling interval must be between 1 and 90"
    assert start_date[0] <= 3, "Start date must be prior to April 1"

    # Check to see if input files exist
    assert os.path.exists(chemograph_dir), "Directory {} does not exist".format(chemograph_dir)


def compile_strategies(sampling, imputation, snap):
    """ This function compiles valid combinations of all run parameters """
    valid_fields = [("sample selection", sampling, ["random", "stratified"]),
                    ("imputation", imputation, ["linear", "log_linear", "stair"]),
                    ("snap", snap, ["yes", "no"])]

    output_sets = []
    for label, strategy, strategies in valid_fields:
        assert strategy in strategies + ["all"], \
            "{} strategy must be '{}', or 'all'".format(label, "', '".join(strategies))
        output_sets.append((strategy,) if strategy != 'all' else strategies)

    combinations = list(itertools.product(*output_sets))
    for i, invalid_set in enumerate(filter(lambda x: x[0] == "random" and x[2] == "yes", combinations)):
        if not i:
            print("Cannot run random sampling with date snap")
        combinations.remove(invalid_set)

    assert len(combinations) > 0, "No valid combinations of sampling method, imputation method, and snapping provided"

    return combinations


@guvectorize(['void(float64[:], float64[:, :], float64[:], float64[:])'], '(n),(o, n),(p)->(p)')
def compute_bias_factor(chemograph, simulated_chemographs, stats_intervals, results):
    # Loop through intervals
    for i in range(stats_intervals.size):
        stats_interval = stats_intervals[i]

        # Get the max moving-window mean value from the chemograph
        true_value = max_mean(chemograph, stats_interval)

        # Get the max moving-window mean for each of the 10,000 simulations
        sampling_values = np.zeros(simulated_chemographs.shape[0])
        for j in range(simulated_chemographs.shape[0]):
            sampling_values[j] = max_mean(simulated_chemographs[j], stats_interval)

        # The bias factor is chemograph max-mean divided by the 5th percentile simulation max-mean
        p5 = np.percentile(sampling_values, 5).round(10)
        results[i] = true_value / p5


@njit
def interpolate(chemograph, n_simulations, n_periods, sampling_interval, tail, sample_method, imputation_method, snap):
    """ Generate a set of 10,000 simulated chemographs based on random sampling """

    # Initialize the output array
    results = np.zeros((n_simulations, chemograph.size)) - 1

    # If doing log-linear imputation, log-transform the chemograph
    if imputation_method == 3:
        chemograph = np.log10(chemograph)

    # Generate sample period boundaries for use if using stratified sampling
    if sample_method > 0:
        sample_periods = np.arange(0, chemograph.size, sampling_interval)

    # Iterate simulations
    for i in range(n_simulations):

        # Set sampling dates
        if sample_method == 1:  # random
            sample_dates = np.sort(np.random.choice(chemograph.size, n_periods + 1, replace=False))
        else:  # stratified
            sample_dates = sample_periods + np.random.randint(0, sampling_interval, sample_periods.shape)

        # bound upper sample to be less than 365. only affects stair-step (JCH - when does this happen?)
        if sample_dates[-1] > 364:
            sample_dates[-1] = 364

        # Iterate sampling dates
        for j in range(n_periods - 1):
            if snap:  # assign the sampled value to the first day in the period
                period_start, period_end = j * sampling_interval, ((j * sampling_interval) + sampling_interval)
            else:  # don't do that
                period_start, period_end = sample_dates[j], sample_dates[j + 1]

            start_val = chemograph[sample_dates[j]]
            end_val = chemograph[sample_dates[j + 1]]

            # Iterate days in sampling period
            for k in range(period_start, period_end + 1):
                if imputation_method == 1:  # stair step
                    results[i, k] = start_val
                    if k == period_end:
                        break
                else:
                    day = k - period_start
                    if imputation_method == 2:  # linear
                        results[i, k] = start_val + ((end_val - start_val) / (period_end - period_start)) * day
                    if imputation_method == 3:  # log-linear
                        results[i, k] = 10 ** (start_val + ((end_val - start_val) / (period_end - period_start)) * day)

        # Extend values across last period if not interpolating
        if imputation_method == 1 and sample_method == 2 and snap == 1:
            if chemograph.size - period_end > sampling_interval:
                for k in range(sampling_interval):
                    results[i, period_end + k] = chemograph[sample_dates[-2]]

            for k in range(tail):
                results[i, period_end + sampling_interval + k] = chemograph[sample_dates[-1]]

    return results


@njit
def max_mean(chemograph, stats_interval):
    """ Use a moving window to get a series of sample means and return the maximum """
    maximum = 0
    for j in range(chemograph.size - stats_interval + 1):
        window_total = 0
        window_count = 0
        for k in range(stats_interval):
            value = chemograph[j + k]
            if value >= 0:
                window_total += value
                window_count += 1
        if window_count > 0:
            window_mean = window_total / window_count
            if window_mean > maximum:
                maximum = window_mean
    return maximum


def percentiles(values, quantiles):
    return [np.percentile(values, quantile) for quantile in quantiles]


@njit
def process_residuals(simulated_chemographs, chemograph):
    """ Calculate slopes in a batch """
    n_repetitions = simulated_chemographs.shape[0]
    slope, rmse = np.zeros(n_repetitions), np.zeros(n_repetitions)
    for i in range(simulated_chemographs.shape[0]):

        # Initialize the x and y arrays to feed to the linear regression function
        x = np.ones((chemograph.size, 2))
        residual = np.zeros(chemograph.size)

        # Loop through each day in chemograph, get residual and add to x, y arrays if the value > 0 (not no-data)
        index = 0
        for j in range(chemograph.size):
            simulated_value = simulated_chemographs[i, j]
            if simulated_value >= 0:
                x[index, 0] = j
                residual[index] = chemograph[j] - simulated_value
                index += 1

        slope[i] = np.linalg.lstsq(x[:index], residual[:index])[0][0]  # linear regression function
        rmse[i] = np.sqrt((residual[:index] ** 2).mean())  # RMSE function

    return slope, rmse


def print_chemographs(chemograph, simulated_chemographs, site_id, n=10):
    """ Write sampled chemographs to file for QC purposes """

    out_dir = os.path.join("QC")  # CAP - os.path.join("..", "usgs files", "atrazine")
    outfile = os.path.join(out_dir, site_id + "chemograph.csv")
    sample = pd.DataFrame(data=np.vstack([chemograph, simulated_chemographs[:n]]).T)
    sample.to_csv(outfile)


def process_chemograph(chemograph, n_simulations, sampling_interval, quantiles,
                       sampling_method, imputation_method, snap, stats_intervals,
                       site=None, write_chemographs=False, qc_n=10):
    """ The main routine: generate 10,000 simulations, and calculate bias factors and such """

    # numba compiled functions (interpolate) don't take strings: convert methods to numbers
    sample_method = {'random': 1, 'stratified': 2}[sampling_method]
    imputation_method = {'stair': 1, 'linear': 2, 'log_linear': 3}[imputation_method]
    snap = {'yes': 1, 'no': 0}[snap]
    stats_intervals = np.array(stats_intervals)

    # The number of sampling periods in a year, and the length of the last incomplete period
    n_periods, tail = divmod(chemograph.size, sampling_interval)

    # Generate randomly-sampled chemograph simulations
    simulated_chemographs = interpolate(chemograph, n_simulations, n_periods, sampling_interval, tail,
                                        sample_method, imputation_method, snap)

    # Print out chemograph and simulated chemographs to QC
    if write_chemographs:
        print_chemographs(chemograph, simulated_chemographs, site, qc_n)

    # Get slopes of input and simulated chemographs
    simulated_slopes, simulation_residuals = process_residuals(simulated_chemographs, chemograph)

    # Calculate bias factors
    bias_factors = compute_bias_factor(chemograph, simulated_chemographs, stats_intervals)

    return percentiles(simulated_slopes, quantiles), percentiles(simulation_residuals, quantiles), bias_factors


def read_chemographs(chemograph_dir, start_date=(1, 1), end_date=(12, 31)):
    # Loop through chemograph files
    for f in filter(lambda x: x.endswith(".txt"), os.listdir(chemograph_dir)):

        # Load tables
        site_name = f.rstrip(".csv")
        table_path = os.path.join(chemograph_dir, f)
        table = pd.read_csv(table_path, sep="\t")

        # Initialize date index and filter out dates outside the specified range
        dates = pd.DatetimeIndex(table.date)  # JCH - originally coerces dtype to np.datetime64 - not necessary?

        # Adjust dates if manual start or end dates are used
        if start_date != (1, 1) or end_date != (12, 31):
            (start_month, start_day), (end_month, end_day) = start_date, end_date
            date_filter = ((dates.month > start_month) | ((dates.month == start_month) & (dates.day >= start_day))) & \
                          ((dates.month < end_month) | ((dates.month == end_month) & (dates.day <= end_day)))
            table = table[date_filter]
            dates = dates[date_filter]

        # Iterate through years to get site-year tables
        chemograph_fields = [f for f in table.columns if f.startswith("csim")]
        for year in sorted(np.unique(dates.year)):
            yield site_name.rstrip(".txt"), year, table[dates.year == year][chemograph_fields].T.as_matrix()


def main():
    # Data paths
    chemograph_dir = os.path.join("AtzRockV2100")  # CAP - os.path.join("..", "usgs_files", "atrazine")
    output_dir = os.path.join("AtzRockV2100","Output")  # CAP - os.path.join("..", "usgs_files", "atrazine", "Output")

    # Parameters
    stats_intervals = [1, 4, 7, 21, 60, 90]  # any positive number between 1 and 90
    start_date = (1, 1)  # date to pick first sample from, must be before Mar 31 (ensures we have at least quarterly)
    end_date = (12, 31)
    sampling_strategy = "random"  # can be 'random', 'stratified', or 'all'
    imputation = "log_linear"  # can be 'linear', 'log_linear', or 'stair'
    snap = "no"  # can be "yes", "no", or "both"
    n_simulations = 10000  # number of random samples
    sampling_intervals = [7, 14, 21, 28]
    quantiles = (0.1, 10, 50, 90, 99.9)

    # Diagnostic
    qc_mode = False  # If True, will write sampled chemographs to file
    qc_sample = 10  # Number of sample chemographs to be written if QC mode is on

    """ All path and parameter variables are set before this line. The rest is program functionality """

    # Inspect inputs and get all combinations of run parameters
    check_inputs(chemograph_dir, sampling_intervals, start_date)
    strategies = compile_strategies(sampling_strategy, imputation, snap)

    # Iterate through all combinations of sampling, imputation, and snap methods
    for sampling_method, imputation_method, snap in strategies:

        # Iterate through all 50 or so chemographs for the site
        output = OutputTable(output_dir, sampling_method, imputation_method, snap, quantiles, stats_intervals)

        # Iterate through sampling intervals:
        for sampling_interval in sampling_intervals:

            # Loop through each site-year and process
            for site, year, chemographs in read_chemographs(chemograph_dir, start_date, end_date):

                # Loop through all chemographs in file
                for i, chemograph in enumerate(chemographs):
                    print("Processing {}, {}, {}-day {} sampling, {} imputation, {} snap, simulation #{} ...".format(
                        site, year, sampling_interval, sampling_method, imputation_method, snap, i + 1))

                    # Analyze the chemograph and perform Monte Carlo sampling
                    slopes, residuals, bias_factors = \
                        process_chemograph(chemograph, n_simulations, sampling_interval, quantiles,
                                           sampling_method, imputation_method, snap, stats_intervals,
                                           site="{}_{}".format(site, i), write_chemographs=qc_mode, qc_n=qc_sample)

                    # Write results to output table
                    output.update_table(site, year, sampling_interval, i+1, bias_factors, slopes, residuals)


main()
