import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat

# Load data
base = "/home/runnalja/git/bgb-jelle/data/CCISTEP1_timeSeriesAndNullInfo/Lake_297_327_GLWD00000327/"
index = 586

# Load time series info
info = loadmat(base + "CCI_NullInfo_chla_mean_731335_738886.mat")
time = pd.to_datetime(np.array(info["timeSeries_dateNR"][0]) - 719529, unit='D', origin='unix')

# Load time series data
data = loadmat(base + "CCI_Data_chla_mean_731335_738886_pixNRs_1_982.mat")
value = np.array(data["timeSeries_data"][index])
qa = np.array(data["timeSeries_QA"][index]) if "timeSeries_QA" in data else None

# Load phenology results
result = loadmat(
    "/home/runnalja/git/bgb-jelle/data/CCISTEP3_phenoMETs/Lake_297_327_GLWD00000327/pixFileDataList_chla_mean_731335_738886/cubicSpline_1_14_0_0_31_0_365_5_pixNRs_1_982.mat")

# Print available keys and shapes
print("Available keys in result file:")
for key in result.keys():
    if not key.startswith('__'):
        print(
            f"  {key}: {type(result[key])}, shape: {np.array(result[key]).shape if hasattr(result[key], 'shape') else 'N/A'}")
print()

# The pixel range is pixNRs_1_982, so pixel 586 is at index 586-1=585 (0-based)
# Or it could be that index 586 in the original data corresponds to index 586 in this file
# Let's check the startPixCoorNR and endPixCoorNR from the original data file
start_pix = data['startPixCoorNR'][0, 0] if 'startPixCoorNR' in data else 1
end_pix = data['endPixCoorNR'][0, 0] if 'endPixCoorNR' in data else 982

print(f"Pixel range in data file: {start_pix} to {end_pix}")
print(f"Looking for pixel index: {index}")
print(f"Adjusted index (0-based): {index - start_pix}")
print()

# Adjust the pixel index
pix_idx = index - start_pix


# Helper function to safely extract data
def safe_extract(result, key, pix_idx):
    """Safely extract data from result dictionary"""
    if key not in result:
        return None

    data = result[key]

    # Handle different data structures
    if isinstance(data, np.ndarray):
        # Check if it's a cell array (object dtype)
        if data.dtype == object:
            # Flatten if needed
            if data.ndim > 1:
                data = data.flatten()
            # Check bounds
            if pix_idx >= len(data):
                print(f"Warning: {key} - index {pix_idx} out of bounds (size: {len(data)})")
                return None
            # Extract the cell content
            cell_data = data[pix_idx]
            # Convert to array if it's a matrix
            if isinstance(cell_data, np.ndarray) and cell_data.size > 0:
                return cell_data.flatten()
            else:
                return None
        else:
            # Regular array
            if data.ndim == 1:
                if pix_idx < len(data):
                    return data[pix_idx]
            elif data.ndim == 2:
                if pix_idx < data.shape[1]:
                    return data[:, pix_idx]
                elif pix_idx < data.shape[0]:
                    return data[pix_idx, :]

    return None


# Get smoothed spline data
smooth_x_axis = result['cubicSpline_smoothXaxis'][0] if 'cubicSpline_smoothXaxis' in result else None

# Get peaks
pks_x = safe_extract(result, 'phenoMETs_pks_X', pix_idx)
pks_y = safe_extract(result, 'phenoMETs_pks_Y', pix_idx)
pks_qa = safe_extract(result, 'phenoMETs_pks_QA', pix_idx)

# Get troughs
trgs_x = safe_extract(result, 'phenoMETs_trgs_X', pix_idx)
trgs_y = safe_extract(result, 'phenoMETs_trgs_Y', pix_idx)
trgs_qa = safe_extract(result, 'phenoMETs_trgs_QA', pix_idx)

# Get green-up metrics
greenup_onset_x = safe_extract(result, 'phenoMETs_GreenUp_Onset_X', pix_idx)
greenup_onset_y = safe_extract(result, 'phenoMETs_GreenUp_Onset_Y', pix_idx)
greenup_mid_x = safe_extract(result, 'phenoMETs_GreenUp_Mid_X', pix_idx)
greenup_mid_y = safe_extract(result, 'phenoMETs_GreenUp_Mid_Y', pix_idx)
greenup_adv_x = safe_extract(result, 'phenoMETs_GreenUp_Advanced_X', pix_idx)
greenup_adv_y = safe_extract(result, 'phenoMETs_GreenUp_Advanced_Y', pix_idx)

# Get green-down metrics
greendown_onset_x = safe_extract(result, 'phenoMETs_GreenDown_Onset_X', pix_idx)
greendown_onset_y = safe_extract(result, 'phenoMETs_GreenDown_Onset_Y', pix_idx)
greendown_mid_x = safe_extract(result, 'phenoMETs_GreenDown_Mid_X', pix_idx)
greendown_mid_y = safe_extract(result, 'phenoMETs_GreenDown_Mid_Y', pix_idx)
greendown_adv_x = safe_extract(result, 'phenoMETs_GreenDown_Advanced_X', pix_idx)
greendown_adv_y = safe_extract(result, 'phenoMETs_GreenDown_Advanced_Y', pix_idx)

# Get data gaps
datagap_start = safe_extract(result, 'cubicSpline_dataGap_startDAYs', pix_idx)
datagap_end = safe_extract(result, 'cubicSpline_dataGap_endDAYs', pix_idx)

# Get quality metrics
rsquare = safe_extract(result, 'cubicSpline_Rsquare', pix_idx)
rmse = safe_extract(result, 'cubicSpline_RMSE', pix_idx)
smpar = safe_extract(result, 'cubicSpline_smPar', pix_idx)

print(f"Extracted data summary:")
print(f"  Peaks: {len(pks_x) if pks_x is not None else 0}")
print(f"  Troughs: {len(trgs_x) if trgs_x is not None else 0}")
print(f"  Green-up phases: {len(greenup_onset_x) if greenup_onset_x is not None else 0}")
print(f"  Green-down phases: {len(greendown_onset_x) if greendown_onset_x is not None else 0}")
print()


# Convert date numbers to datetime for plotting
def datenum_to_datetime(datenum):
    """Convert MATLAB datenum to pandas datetime"""
    if datenum is None or len(datenum) == 0:
        return None
    return pd.to_datetime(datenum - 719529, unit='D', origin='unix')


# Create figure with subplots
fig, axes = plt.subplots(2, 1, figsize=(16, 10))

# ========== PLOT 1: Full time series with peaks and troughs ==========
ax1 = axes[0]

# Filter out NaN values for plotting
valid_mask = ~np.isnan(value)
time_valid = time[valid_mask]
value_valid = value[valid_mask]

# Color code by QA if available
if qa is not None:
    qa_valid = qa[valid_mask]
    colors = ['green' if q == 0 else 'orange' if q == 1 else 'red' for q in qa_valid]
    ax1.scatter(time_valid, value_valid, c=colors, s=30, alpha=0.6, zorder=2)
    # Create custom legend for QA
    from matplotlib.patches import Patch

    legend_elements = [Patch(facecolor='green', alpha=0.6, label='QA=0 (Good)'),
                       Patch(facecolor='orange', alpha=0.6, label='QA=1 (Fair)'),
                       Patch(facecolor='red', alpha=0.6, label='QA=2 (Poor)')]
else:
    ax1.scatter(time_valid, value_valid, c='blue', s=30, alpha=0.6, label='Original data', zorder=2)
    legend_elements = []

# Plot peaks
if pks_x is not None and len(pks_x) > 0:
    pks_time = datenum_to_datetime(pks_x)
    if pks_qa is not None and len(pks_qa) > 0:
        # Color peaks by QA
        qa_colors = {0: 'darkgreen', 1: 'gold', 2: 'darkred'}
        qa_markers = {0: '^', 1: '^', 2: 'x'}
        for q in np.unique(pks_qa):
            mask = pks_qa == q
            ax1.scatter(pks_time[mask], pks_y[mask], c=qa_colors.get(q, 'gray'),
                        marker=qa_markers.get(q, '^'), s=200,
                        edgecolors='black', linewidth=1.5, zorder=5,
                        label=f'Peak (QA={int(q)})')
    else:
        ax1.scatter(pks_time, pks_y, c='green', marker='^', s=200,
                    edgecolors='black', linewidth=1.5, label='Peaks', zorder=5)

# Plot troughs
if trgs_x is not None and len(trgs_x) > 0:
    trgs_time = datenum_to_datetime(trgs_x)
    if trgs_qa is not None and len(trgs_qa) > 0:
        # Color troughs by QA
        qa_colors = {0: 'darkorange', 1: 'gold', 2: 'darkred'}
        qa_markers = {0: 'v', 1: 'v', 2: 'x'}
        for q in np.unique(trgs_qa):
            mask = trgs_qa == q
            ax1.scatter(trgs_time[mask], trgs_y[mask], c=qa_colors.get(q, 'gray'),
                        marker=qa_markers.get(q, 'v'), s=200,
                        edgecolors='black', linewidth=1.5, zorder=5,
                        label=f'Trough (QA={int(q)})')
    else:
        ax1.scatter(trgs_time, trgs_y, c='orange', marker='v', s=200,
                    edgecolors='black', linewidth=1.5, label='Troughs', zorder=5)

# Plot data gaps as shaded regions
if datagap_start is not None and datagap_end is not None and len(datagap_start) > 0:
    for i, (start, end) in enumerate(zip(datagap_start, datagap_end)):
        start_time = datenum_to_datetime(np.array([start]))[0]
        end_time = datenum_to_datetime(np.array([end]))[0]
        ax1.axvspan(start_time, end_time, alpha=0.2, color='red',
                    label='Data gap' if i == 0 else '')

ax1.set_xlabel('Date', fontsize=12, fontweight='bold')
ax1.set_ylabel('Chlorophyll-a (mean)', fontsize=12, fontweight='bold')
ax1.set_title(f'Time Series with Detected Peaks and Troughs (Pixel {index})', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(loc='best', fontsize=9)

# ========== PLOT 2: Phenology metrics detail ==========
ax2 = axes[1]

# Plot the same base data
ax2.scatter(time_valid, value_valid, c='blue', s=20, alpha=0.4, label='Original data', zorder=2)

# Plot peaks and troughs again
if pks_x is not None and len(pks_x) > 0:
    pks_time = datenum_to_datetime(pks_x)
    ax2.scatter(pks_time, pks_y, c='green', marker='^', s=150,
                edgecolors='black', linewidth=1.5, label='Peaks', zorder=5)

if trgs_x is not None and len(trgs_x) > 0:
    trgs_time = datenum_to_datetime(trgs_x)
    ax2.scatter(trgs_time, trgs_y, c='orange', marker='v', s=150,
                edgecolors='black', linewidth=1.5, label='Troughs', zorder=5)

# Plot green-up metrics
if greenup_onset_x is not None and len(greenup_onset_x) > 0:
    onset_time = datenum_to_datetime(greenup_onset_x)
    ax2.scatter(onset_time, greenup_onset_y, c='lightgreen', marker='o', s=100,
                edgecolors='darkgreen', linewidth=2, label='Green-up Onset (25%)', zorder=4)

if greenup_mid_x is not None and len(greenup_mid_x) > 0:
    mid_time = datenum_to_datetime(greenup_mid_x)
    ax2.scatter(mid_time, greenup_mid_y, c='lime', marker='s', s=100,
                edgecolors='darkgreen', linewidth=2, label='Green-up Mid (50%)', zorder=4)

if greenup_adv_x is not None and len(greenup_adv_x) > 0:
    adv_time = datenum_to_datetime(greenup_adv_x)
    ax2.scatter(adv_time, greenup_adv_y, c='forestgreen', marker='D', s=100,
                edgecolors='darkgreen', linewidth=2, label='Green-up Advanced (75%)', zorder=4)

# Plot green-down metrics
if greendown_onset_x is not None and len(greendown_onset_x) > 0:
    onset_time = datenum_to_datetime(greendown_onset_x)
    ax2.scatter(onset_time, greendown_onset_y, c='lightsalmon', marker='o', s=100,
                edgecolors='darkred', linewidth=2, label='Green-down Onset (75%)', zorder=4)

if greendown_mid_x is not None and len(greendown_mid_x) > 0:
    mid_time = datenum_to_datetime(greendown_mid_x)
    ax2.scatter(mid_time, greendown_mid_y, c='coral', marker='s', s=100,
                edgecolors='darkred', linewidth=2, label='Green-down Mid (50%)', zorder=4)

if greendown_adv_x is not None and len(greendown_adv_x) > 0:
    adv_time = datenum_to_datetime(greendown_adv_x)
    ax2.scatter(adv_time, greendown_adv_y, c='darkred', marker='D', s=100,
                edgecolors='darkred', linewidth=2, label='Green-down Advanced (25%)', zorder=4)

ax2.set_xlabel('Date', fontsize=12, fontweight='bold')
ax2.set_ylabel('Chlorophyll-a (mean)', fontsize=12, fontweight='bold')
ax2.set_title(f'Phenology Metrics Detail (Pixel {index})', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(loc='best', fontsize=9, ncol=2)

plt.tight_layout()
plt.savefig(f'phenology_metrics_pixel_{index}.png', dpi=150, bbox_inches='tight')
plt.show()

# Print summary statistics
print(f"\n{'=' * 60}")
print(f"Phenology Metrics Summary for Pixel {index}")
print(f"{'=' * 60}")
print(f"Number of peaks detected: {len(pks_x) if pks_x is not None else 0}")
print(f"Number of troughs detected: {len(trgs_x) if trgs_x is not None else 0}")
print(f"Number of green-up phases: {len(greenup_onset_x) if greenup_onset_x is not None else 0}")
print(f"Number of green-down phases: {len(greendown_onset_x) if greendown_onset_x is not None else 0}")
print(f"Number of data gaps: {len(datagap_start) if datagap_start is not None else 0}")

if rsquare is not None and len(rsquare) > 0:
    print(f"R-squared of fit: {rsquare[0]:.4f}")

if rmse is not None and len(rmse) > 0:
    print(f"RMSE of fit: {rmse[0]:.4f}")

if smpar is not None and len(smpar) > 0:
    print(f"Smoothing parameter: {smpar[0]:.6f}")