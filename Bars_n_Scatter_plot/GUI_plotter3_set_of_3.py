%reset -f

import tkinter as tk
from tkinter import messagebox
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random

def convert_to_numpy():
    try:
        text = text_box.get("1.0", tk.END).strip()
        if not text:
            messagebox.showwarning("Warning", "No data pasted!")
            return

        # Convert tab-separated text into a DataFrame
        data = [row.split('\t') for row in text.split('\n') if row]
        df = pd.DataFrame(data)

        # Convert all values to numeric, coercing errors to NaN and dropping them
        df = df.apply(pd.to_numeric, errors='coerce').dropna()

        # Convert DataFrame to NumPy array
        global data2plot
        data2plot = df.to_numpy()

        messagebox.showinfo("Success", "Data stored in variable 'data2plot'")

        # Close the GUI window
        root.destroy()

        # Plot the data
        plot_data(data2plot)
    except Exception as e:
        messagebox.showerror("Error", f"Failed to convert to numpy array: {e}")

def plot_data(data):
    # Convert to DataFrame for easier handling
    df = pd.DataFrame(data)
    total_cols = df.shape[1]
    num_full_groups = total_cols // 3  # number of complete groups of 3 columns
    remainder = total_cols % 3         # number of extra columns (0, 1, or 2)

    # Pre-calculate means and standard errors for each column
    means = df.mean()
    std_errors = df.std() / np.sqrt(df.count())

    # Parameters for plotting
    dot_spread = 0.05

    # Set up the figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Excel-like colors (we use the first three for full groups)
    excel_colors = ['#4472C4', '#ED7D31', '#A5A5A5', '#FFC000', '#70AD47', '#FF0000']

    # --- Plot full groups of 3 ---
    group_width = 0.8          # total width for a full group
    bar_width_full = group_width / 3  # each bar in a triplet
    # x positions for groups: 0, 1, 2, ..., num_full_groups - 1
    x_full = np.arange(num_full_groups)
    
    # For each full group, plot the three bars with offsets:
    # We use offsets: [-bar_width_full, 0, +bar_width_full] relative to the group center.
    for i in range(num_full_groups):
        x_center = x_full[i]
        offsets = [-bar_width_full, 0, bar_width_full]
        for j in range(3):
            col_index = 3 * i + j
            x_pos = x_center + offsets[j]
            bar = ax.bar(x_pos, means.iloc[col_index], width=bar_width_full, color=excel_colors[j % len(excel_colors)])
            # Plot scatter points for this column
            for measurement in df.iloc[:, col_index].dropna():
                ax.plot(x_pos + random.uniform(-dot_spread, dot_spread),
                        measurement,
                        marker='o', markersize=4, color='black', linestyle='None')
            # Add error bar for this column
            ax.errorbar(x_pos, means.iloc[col_index], yerr=std_errors.iloc[col_index], fmt='none',
                        ecolor='red', elinewidth=5, capsize=8, capthick=3, zorder=3)

    # --- Plot extra group if there are leftover columns ---
    if remainder > 0:
        x_extra = np.array([num_full_groups])  # place extra group at the next integer x position
        if remainder == 2:
            # For a group of 2, use the same bar width as in full groups.
            bar_width_extra = group_width / 3
            # Offsets: center the two bars; use offsets: -bar_width_extra/2 and +bar_width_extra/2.
            offsets = [-bar_width_extra/2, bar_width_extra/2]
            for j in range(2):
                col_index = 3 * num_full_groups + j
                x_pos = x_extra[0] + offsets[j]
                bar = ax.bar(x_pos, means.iloc[col_index], width=bar_width_extra, color=excel_colors[j % len(excel_colors)])
                for measurement in df.iloc[:, col_index].dropna():
                    ax.plot(x_pos + random.uniform(-dot_spread, dot_spread),
                            measurement,
                            marker='o', markersize=4, color='black', linestyle='None')
                ax.errorbar(x_pos, means.iloc[col_index], yerr=std_errors.iloc[col_index], fmt='none',
                            ecolor='red', elinewidth=5, capsize=8, capthick=3, zorder=3)
        elif remainder == 1:
            # For a single extra column, center it.
            col_index = 3 * num_full_groups
            bar_width_extra = 0.4  # an arbitrary width for a single bar
            x_pos = x_extra[0]
            bar = ax.bar(x_pos, means.iloc[col_index], width=bar_width_extra, color=excel_colors[0])
            for measurement in df.iloc[:, col_index].dropna():
                ax.plot(x_pos + random.uniform(-dot_spread, dot_spread),
                        measurement,
                        marker='o', markersize=4, color='black', linestyle='None')
            ax.errorbar(x_pos, means.iloc[col_index], yerr=std_errors.iloc[col_index], fmt='none',
                        ecolor='red', elinewidth=5, capsize=8, capthick=3, zorder=3)
    
    # Remove all text elements except enlarged Y-axis tick labels
    ax.set_xticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    
    # Increase Y-axis tick label size
    ax.tick_params(axis='y', labelsize=20)
    
    plt.show()

# Create main window
root = tk.Tk()
root.title("Excel to NumPy Converter")

# Input text box
text_box = tk.Text(root, height=10, width=50)
text_box.pack(pady=10)

# Convert button
convert_button = tk.Button(root, text="Convert to NumPy", command=convert_to_numpy)
convert_button.pack()

# Run the application
root.mainloop()
